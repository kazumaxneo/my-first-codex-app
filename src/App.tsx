import { useEffect, useMemo, useRef, useState } from 'react';
import type { PointerEvent } from 'react';
import {
  Activity,
  Dna,
  Gauge,
  Pause,
  Play,
  RotateCcw,
  SlidersHorizontal,
  Sparkles,
  UploadCloud,
} from 'lucide-react';
import { parseBamFile, parseReferenceFile, type BamAnalysis, type ReferenceSequence } from './bam';

type Settings = {
  readCount: number;
  readLength: number;
  errorRate: number;
  gcBias: number;
  duplicateRate: number;
  variantAlleleFraction: number;
  variantPosition: number;
  seed: number;
  playing: boolean;
};

type Read = {
  id: number;
  start: number;
  end: number;
  row: number;
  strand: 1 | -1;
  duplicate: boolean;
  hasAlt: boolean;
  quality: number;
  mismatchPositions: number[];
  indelPositions: number[];
};

type Simulation = {
  source: 'synthetic' | 'bam';
  fileName?: string;
  referenceName?: string;
  referenceSequence?: string;
  genomeLength: number;
  variantPosition: number;
  reads: Read[];
  coverage: number[];
  altCoverage: number[];
  gc: number[];
  maxCoverage: number;
  observedAltFraction: number;
  variantDepth: number;
  confidence: 'Low' | 'Watch' | 'Strong';
  mappedCount?: number;
};

const initialSettings: Settings = {
  readCount: 180,
  readLength: 76,
  errorRate: 1.2,
  gcBias: 34,
  duplicateRate: 8,
  variantAlleleFraction: 50,
  variantPosition: 690,
  seed: 18,
  playing: false,
};

const genomeLength = 1200;

function mulberry32(seed: number) {
  return () => {
    let t = (seed += 0x6d2b79f5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function clamp(value: number, min: number, max: number) {
  return Math.min(max, Math.max(min, value));
}

function makeGcProfile(length: number) {
  return Array.from({ length }, (_, i) => {
    const wave = Math.sin(i / 62) * 0.17 + Math.sin(i / 23 + 1.4) * 0.07;
    const island = Math.exp(-Math.pow((i - 420) / 95, 2)) * 0.28;
    const desert = Math.exp(-Math.pow((i - 910) / 130, 2)) * -0.23;
    return clamp(0.46 + wave + island + desert, 0.16, 0.82);
  });
}

function simulate(settings: Settings): Simulation {
  const random = mulberry32(settings.seed);
  const gc = makeGcProfile(genomeLength);
  const reads: Read[] = [];
  const coverage = new Array(genomeLength).fill(0);
  const altCoverage = new Array(genomeLength).fill(0);
  const rowEnds: number[] = [];
  const duplicates = Math.round(settings.readCount * (settings.duplicateRate / 100));
  const uniqueReads = settings.readCount - duplicates;

  const placeRead = (id: number, forcedStart?: number, duplicate = false) => {
    let start = forcedStart ?? Math.floor(random() * (genomeLength - settings.readLength));
    if (forcedStart === undefined && settings.gcBias > 0) {
      for (let tries = 0; tries < 10; tries += 1) {
        const candidate = Math.floor(random() * (genomeLength - settings.readLength));
        const localGc = gc[candidate + Math.floor(settings.readLength / 2)];
        const distanceFromIdeal = Math.abs(localGc - 0.52);
        const rejection = distanceFromIdeal * (settings.gcBias / 55);
        if (random() > rejection) {
          start = candidate;
          break;
        }
      }
    }

    const end = start + settings.readLength;
    const row = rowEnds.findIndex((rowEnd) => start > rowEnd + 5);
    const finalRow = row === -1 ? rowEnds.length : row;
    rowEnds[finalRow] = end;
    const hasAlt = start <= settings.variantPosition && end >= settings.variantPosition && random() * 100 < settings.variantAlleleFraction;
    const mismatchPositions: number[] = [];
    const indelPositions: number[] = [];
    const errorEvents = Math.floor((settings.errorRate / 100) * settings.readLength + random() * 2);
    for (let i = 0; i < errorEvents; i += 1) {
      const eventPosition = start + Math.floor(random() * settings.readLength);
      if (random() < 0.24) indelPositions.push(eventPosition);
      else mismatchPositions.push(eventPosition);
    }

    const read: Read = {
      id,
      start,
      end,
      row: finalRow,
      strand: random() > 0.5 ? 1 : -1,
      duplicate,
      hasAlt,
      quality: Math.round(20 + random() * 18 - settings.errorRate * 0.7),
      mismatchPositions,
      indelPositions,
    };
    reads.push(read);

    for (let pos = start; pos < end && pos < genomeLength; pos += 1) {
      coverage[pos] += 1;
      if (hasAlt && Math.abs(pos - settings.variantPosition) < 2) altCoverage[pos] += 1;
    }
  };

  for (let i = 0; i < uniqueReads; i += 1) placeRead(i);
  for (let i = 0; i < duplicates; i += 1) {
    const source = reads[Math.floor(random() * Math.max(reads.length, 1))];
    placeRead(uniqueReads + i, source?.start, true);
  }

  const maxCoverage = Math.max(...coverage, 1);
  const variantDepth = coverage[settings.variantPosition];
  const observedAltFraction = variantDepth ? altCoverage[settings.variantPosition] / variantDepth : 0;
  const confidence =
    variantDepth >= 18 && observedAltFraction > 0.22
      ? 'Strong'
      : variantDepth >= 8 && observedAltFraction > 0.12
        ? 'Watch'
        : 'Low';

  return {
    source: 'synthetic',
    genomeLength,
    variantPosition: settings.variantPosition,
    reads,
    coverage,
    altCoverage,
    gc,
    maxCoverage,
    observedAltFraction,
    variantDepth,
    confidence,
  };
}

function simulationFromBam(analysis: BamAnalysis, variantPosition: number): Simulation {
  const position = clamp(variantPosition, 0, analysis.genomeLength - 1);
  const reads: Read[] = analysis.reads.map((read) => {
    const supportsEvent =
      read.mismatchPositions.some((pos) => Math.abs(pos - position) <= 1) ||
      read.indelPositions.some((pos) => Math.abs(pos - position) <= 1);
    return {
      id: read.id,
      start: read.start,
      end: read.end,
      row: read.row,
      strand: read.strand,
      duplicate: read.duplicate,
      hasAlt: supportsEvent,
      quality: read.mapq,
      mismatchPositions: read.mismatchPositions,
      indelPositions: read.indelPositions,
    };
  });
  const eventReads = reads.filter(
    (read) =>
      read.start <= position &&
      read.end >= position &&
      (read.mismatchPositions.some((pos) => Math.abs(pos - position) <= 1) ||
        read.indelPositions.some((pos) => Math.abs(pos - position) <= 1)),
  ).length;
  const variantDepth = analysis.coverage[position] ?? 0;
  const observedAltFraction = variantDepth ? eventReads / variantDepth : 0;
  const confidence =
    variantDepth >= 18 && observedAltFraction > 0.22
      ? 'Strong'
      : variantDepth >= 8 && observedAltFraction > 0.12
        ? 'Watch'
        : 'Low';

  return {
    source: 'bam',
    fileName: analysis.fileName,
    referenceName: analysis.referenceName,
    referenceSequence: analysis.referenceSequence,
    genomeLength: analysis.genomeLength,
    variantPosition: position,
    reads,
    coverage: analysis.coverage,
    altCoverage: new Array(analysis.genomeLength).fill(0),
    gc: new Array(analysis.genomeLength).fill(0.5),
    maxCoverage: analysis.maxCoverage,
    observedAltFraction,
    variantDepth,
    confidence,
    mappedCount: analysis.mappedCount,
  };
}

function drawRoundedRect(
  context: CanvasRenderingContext2D,
  x: number,
  y: number,
  width: number,
  height: number,
  radius: number,
) {
  const r = Math.min(radius, width / 2, height / 2);
  context.beginPath();
  context.moveTo(x + r, y);
  context.arcTo(x + width, y, x + width, y + height, r);
  context.arcTo(x + width, y + height, x, y + height, r);
  context.arcTo(x, y + height, x, y, r);
  context.arcTo(x, y, x + width, y, r);
  context.closePath();
}

function CoverageCanvas({
  simulation,
  onVariantPositionChange,
}: {
  simulation: Simulation;
  onVariantPositionChange: (position: number) => void;
}) {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const [isDragging, setIsDragging] = useState(false);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const context = canvas.getContext('2d');
    if (!context) return;

    const draw = () => {
      const rect = canvas.getBoundingClientRect();
      const pixelRatio = window.devicePixelRatio || 1;
      canvas.width = Math.floor(rect.width * pixelRatio);
      canvas.height = Math.floor(rect.height * pixelRatio);
      context.setTransform(pixelRatio, 0, 0, pixelRatio, 0, 0);

      const width = rect.width;
      const height = rect.height;
      const left = 58;
      const right = 28;
      const trackWidth = width - left - right;
      const xFor = (pos: number) => left + (pos / simulation.genomeLength) * trackWidth;

      context.clearRect(0, 0, width, height);
      context.fillStyle = '#ffffff';
      context.fillRect(0, 0, width, height);

      for (let tick = 0; tick <= 12; tick += 1) {
        const x = left + (tick / 12) * trackWidth;
        context.strokeStyle = '#eef1f3';
        context.lineWidth = 1;
        context.beginPath();
        context.moveTo(x, 22);
        context.lineTo(x, height - 30);
        context.stroke();
      }

      context.font = '600 11px Inter, system-ui, sans-serif';
      context.fillStyle = '#66727a';
      context.fillText('GC%', 16, 45);
      context.fillText('DEPTH', 16, 122);
      context.fillText('READS', 16, 193);

      const gcY = 28;
      const gcHeight = 62;
      context.strokeStyle = '#dde4e8';
      context.lineWidth = 1;
      context.strokeRect(left, gcY, trackWidth, gcHeight);
      context.beginPath();
      if (simulation.source === 'synthetic') {
        simulation.gc.forEach((gc, pos) => {
          const x = xFor(pos);
          const y = gcY + gcHeight - gc * gcHeight;
          if (pos === 0) context.moveTo(x, y);
          else context.lineTo(x, y);
        });
        context.strokeStyle = '#45a36f';
        context.lineWidth = 2;
        context.stroke();
      } else {
        context.fillStyle = '#f7f9fa';
        context.fillRect(left, gcY, trackWidth, gcHeight);
        context.fillStyle = '#69757c';
        context.font = '600 12px Inter, system-ui, sans-serif';
        const label = simulation.referenceSequence
          ? `${simulation.referenceName ?? 'reference'} (${simulation.genomeLength.toLocaleString()} bp, FASTA loaded)`
          : `${simulation.referenceName ?? 'reference'} (${simulation.genomeLength.toLocaleString()} bp)`;
        context.fillText(label, left + 12, gcY + 36);
      }

      const coverageY = 110;
      const coverageHeight = 56;
      context.fillStyle = '#f7f9fa';
      context.fillRect(left, coverageY, trackWidth, coverageHeight);
      context.strokeStyle = '#dfe5e8';
      context.strokeRect(left, coverageY, trackWidth, coverageHeight);
      for (let i = 0; i < simulation.coverage.length; i += 2) {
        const depth = simulation.coverage[i];
        const normalized = depth / simulation.maxCoverage;
        const barHeight = normalized * (coverageHeight - 8);
        context.fillStyle = '#cbd2d6';
        context.fillRect(xFor(i), coverageY + coverageHeight - 4 - barHeight, Math.max(1, trackWidth / simulation.genomeLength * 2), barHeight);
      }

      const readY = 185;
      const rowHeight = 12;
      const maxRows = Math.floor((height - 230) / rowHeight);
      simulation.reads.forEach((read) => {
        if (read.row >= maxRows) return;
        const x = xFor(read.start);
        const w = Math.max(2, xFor(read.end) - x);
        const y = readY + read.row * rowHeight;
        const alpha = read.duplicate ? 0.38 : 0.86;
        context.fillStyle = read.hasAlt ? `rgba(244, 177, 66, ${alpha})` : `rgba(152, 160, 165, ${alpha})`;
        drawRoundedRect(context, x, y, w, 6, 2);
        context.fill();
        if (read.strand === -1) {
          context.fillStyle = 'rgba(82, 90, 96, 0.32)';
          context.fillRect(x + 2, y + 3, Math.max(1, w - 4), 1);
        }
        context.fillStyle = '#d94a4a';
        read.mismatchPositions.slice(0, 5).forEach((pos) => {
          const mx = xFor(pos);
          if (mx >= x && mx <= x + w) context.fillRect(mx, y - 2, 2, 10);
        });
        read.indelPositions.slice(0, 3).forEach((pos) => {
          const mx = xFor(pos);
          if (mx < x || mx > x + w) return;
          context.strokeStyle = '#7b61d8';
          context.lineWidth = 2;
          context.beginPath();
          context.moveTo(mx - 4, y - 2);
          context.lineTo(mx, y + 7);
          context.lineTo(mx + 4, y - 2);
          context.stroke();
        });
      });

      const vx = xFor(simulation.variantPosition);
      context.strokeStyle = isDragging ? '#d08a00' : '#f2b928';
      context.lineWidth = isDragging ? 3 : 2;
      context.beginPath();
      context.moveTo(vx, 16);
      context.lineTo(vx, height - 32);
      context.stroke();
      context.fillStyle = isDragging ? '#d08a00' : '#f2b928';
      context.beginPath();
      context.arc(vx, 17, 8, 0, Math.PI * 2);
      context.fill();
      context.fillStyle = '#ffffff';
      context.font = '800 10px Inter, system-ui, sans-serif';
      context.fillText('S', vx - 3, 21);
      context.fillStyle = '#8d6512';
      context.fillText(`${simulation.variantPosition} bp`, Math.min(vx + 10, width - 82), 23);

      context.fillStyle = '#69757c';
      context.font = '600 11px Inter, system-ui, sans-serif';
      context.fillText('0 bp', left, height - 20);
      context.fillText(`${simulation.genomeLength} bp`, width - 82, height - 20);
    };

    draw();
    const resizeObserver = new ResizeObserver(draw);
    resizeObserver.observe(canvas);
    return () => resizeObserver.disconnect();
  }, [simulation, isDragging]);

  const positionFromPointer = (event: PointerEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return simulation.variantPosition;
    const rect = canvas.getBoundingClientRect();
    const left = 58;
    const right = 28;
    const trackWidth = rect.width - left - right;
    const x = clamp(event.clientX - rect.left, left, rect.width - right);
    return Math.round(((x - left) / trackWidth) * simulation.genomeLength);
  };

  return (
    <canvas
      className={`coverage-canvas ${isDragging ? 'dragging' : ''}`}
      ref={canvasRef}
      aria-label="Genome coverage visualization"
      onPointerDown={(event) => {
        setIsDragging(true);
        event.currentTarget.setPointerCapture(event.pointerId);
        onVariantPositionChange(positionFromPointer(event));
      }}
      onPointerMove={(event) => {
        if (isDragging) onVariantPositionChange(positionFromPointer(event));
      }}
      onPointerUp={(event) => {
        setIsDragging(false);
        if (event.currentTarget.hasPointerCapture(event.pointerId)) {
          event.currentTarget.releasePointerCapture(event.pointerId);
        }
      }}
      onPointerLeave={(event) => {
        if (isDragging) {
          setIsDragging(false);
          if (event.currentTarget.hasPointerCapture(event.pointerId)) {
            event.currentTarget.releasePointerCapture(event.pointerId);
          }
        }
      }}
    />
  );
}

function Metric({ label, value, tone }: { label: string; value: string; tone?: 'green' | 'amber' | 'rose' }) {
  return (
    <div className={`metric ${tone ?? ''}`}>
      <span>{label}</span>
      <strong>{value}</strong>
    </div>
  );
}

function Slider({
  label,
  value,
  min,
  max,
  step = 1,
  unit,
  onChange,
}: {
  label: string;
  value: number;
  min: number;
  max: number;
  step?: number;
  unit?: string;
  onChange: (value: number) => void;
}) {
  return (
    <label className="control">
      <span>
        {label}
        <strong>
          {value}
          {unit}
        </strong>
      </span>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(event) => onChange(Number(event.target.value))}
      />
    </label>
  );
}

export default function App() {
  const [settings, setSettings] = useState<Settings>(initialSettings);
  const [bamAnalysis, setBamAnalysis] = useState<BamAnalysis | null>(null);
  const [reference, setReference] = useState<ReferenceSequence | null>(null);
  const [uploadMessage, setUploadMessage] = useState('BAM: not loaded');
  const [referenceMessage, setReferenceMessage] = useState('Reference: not loaded');

  useEffect(() => {
    if (!settings.playing) return;
    const interval = window.setInterval(() => {
      setSettings((current) => ({ ...current, seed: current.seed + 1 }));
    }, 900);
    return () => window.clearInterval(interval);
  }, [settings.playing]);

  const simulation = useMemo(
    () => (bamAnalysis ? simulationFromBam(bamAnalysis, settings.variantPosition) : simulate(settings)),
    [bamAnalysis, settings],
  );
  const duplicateReads = simulation.reads.filter((read) => read.duplicate).length;
  const confidenceTone = simulation.confidence === 'Strong' ? 'green' : simulation.confidence === 'Watch' ? 'amber' : 'rose';

  const update = (patch: Partial<Settings>) => setSettings((current) => ({ ...current, ...patch }));

  const handleBamUpload = async (file: File | undefined) => {
    if (!file) return;
    setUploadMessage('BAMを解析中...');
    try {
      const analysis = await parseBamFile(file, reference);
      setBamAnalysis(analysis);
      update({ variantPosition: Math.min(settings.variantPosition, analysis.genomeLength - 1), playing: false });
      setUploadMessage(`${analysis.fileName} / ${analysis.referenceName} / ${analysis.mappedCount} reads`);
    } catch (error) {
      setBamAnalysis(null);
      setUploadMessage(error instanceof Error ? error.message : 'BAMの解析に失敗しました。');
    }
  };

  const handleReferenceUpload = async (file: File | undefined) => {
    if (!file) return;
    setReferenceMessage('FASTAを読み込み中...');
    try {
      const parsed = await parseReferenceFile(file);
      setReference(parsed);
      setBamAnalysis(null);
      setReferenceMessage(`${parsed.fileName} / ${parsed.name} / ${parsed.sequence.length.toLocaleString()} bp`);
      setUploadMessage('BAM: not loaded');
      update({ variantPosition: Math.min(settings.variantPosition, parsed.sequence.length - 1), playing: false });
    } catch (error) {
      setReference(null);
      setBamAnalysis(null);
      setReferenceMessage(error instanceof Error ? error.message : 'FASTAの読み込みに失敗しました。');
    }
  };

  return (
    <main className="app-shell">
      <section className="workspace">
        <header className="topbar">
          <div className="brand">
            <div className="brand-mark">
              <Dna size={24} />
            </div>
            <div>
              <h1>Genome Coverage Visualizer</h1>
              <p>short-read sequencing depth and variant evidence</p>
            </div>
          </div>
          <div className="toolbar">
            <button
              className="icon-button"
              type="button"
              title={settings.playing ? 'Pause simulation' : 'Animate simulation'}
              onClick={() => update({ playing: !settings.playing })}
            >
              {settings.playing ? <Pause size={18} /> : <Play size={18} />}
            </button>
            <button
              className="icon-button"
              type="button"
              title="New random library"
              onClick={() => {
                setBamAnalysis(null);
                setReference(null);
                setReferenceMessage('Reference: not loaded');
                setUploadMessage('BAM: not loaded');
                update({ seed: settings.seed + 11 });
              }}
            >
              <Sparkles size={18} />
            </button>
            <button
              className="icon-button"
              type="button"
              title="Reset controls"
              onClick={() => {
                setBamAnalysis(null);
                setReference(null);
                setReferenceMessage('Reference: not loaded');
                setUploadMessage('BAM: not loaded');
                setSettings(initialSettings);
              }}
            >
              <RotateCcw size={18} />
            </button>
          </div>
        </header>

        <div className="visual-stage">
          <CoverageCanvas simulation={simulation} onVariantPositionChange={(variantPosition) => update({ variantPosition })} />
        </div>

        <section className="metrics-row" aria-label="Simulation metrics">
          <Metric
            label="Mean depth"
            value={`${(simulation.coverage.reduce((sum, depth) => sum + depth, 0) / simulation.genomeLength).toFixed(1)}x`}
            tone="green"
          />
          <Metric label="Max depth" value={`${simulation.maxCoverage}x`} />
          <Metric label="SNV depth" value={`${simulation.variantDepth}x`} tone={confidenceTone} />
          <Metric label="Alt fraction" value={`${Math.round(simulation.observedAltFraction * 100)}%`} tone={confidenceTone} />
          <Metric label="Duplicates" value={`${duplicateReads}`} tone="amber" />
        </section>
      </section>

      <aside className="side-panel">
        <div className="upload-card">
          <div className="panel-heading compact">
            <UploadCloud size={18} />
            <h2>Input Files</h2>
          </div>
          <label className="file-picker">
            <input type="file" accept=".fa,.fasta,.fna,text/plain" onChange={(event) => handleReferenceUpload(event.target.files?.[0])} />
            <span>Choose FASTA reference</span>
          </label>
          <p>{referenceMessage}</p>
          <label className="file-picker secondary">
            <input type="file" accept=".bam,application/octet-stream" onChange={(event) => handleBamUpload(event.target.files?.[0])} />
            <span>Choose BAM</span>
          </label>
          <p>{uploadMessage}</p>
          <small>FASTAは10kb以下、BAMは5MB以下に限定しています。BAIは不要です。FASTAを先に選ぶとSNV検出がより正確になります。</small>
        </div>

        <div className="panel-heading">
          <SlidersHorizontal size={18} />
          <h2>{bamAnalysis ? 'Synthetic Model' : 'Library Model'}</h2>
        </div>
        <div className="controls">
          <Slider label="Reads" value={settings.readCount} min={40} max={460} step={10} onChange={(readCount) => update({ readCount })} />
          <Slider label="Read length" value={settings.readLength} min={38} max={160} step={2} unit=" bp" onChange={(readLength) => update({ readLength })} />
          <Slider label="Sequencing error" value={settings.errorRate} min={0} max={8} step={0.2} unit="%" onChange={(errorRate) => update({ errorRate })} />
          <Slider label="GC bias" value={settings.gcBias} min={0} max={85} step={1} unit="%" onChange={(gcBias) => update({ gcBias })} />
          <Slider label="PCR duplicates" value={settings.duplicateRate} min={0} max={50} step={1} unit="%" onChange={(duplicateRate) => update({ duplicateRate })} />
          <Slider
            label="Alt allele"
            value={settings.variantAlleleFraction}
            min={0}
            max={100}
            step={5}
            unit="%"
            onChange={(variantAlleleFraction) => update({ variantAlleleFraction })}
          />
        </div>

        <div className="insight-block">
          <div className="panel-heading compact">
            <Activity size={18} />
            <h2>Evidence</h2>
          </div>
          <div className={`callout ${confidenceTone}`}>
            <strong>{simulation.confidence}</strong>
            <span>
              depth {simulation.variantDepth}x / alt {Math.round(simulation.observedAltFraction * 100)}%
            </span>
          </div>
          <p>
            黄色の縦線はドラッグできる候補位置です。赤い縦マークはSNV、紫のV字マークはindelを表します。周囲のpileupとalt allele比率から、低深度や高エラーが判定に与える影響を比較できます。
          </p>
        </div>

        <div className="legend">
          <div>
            <span className="swatch read" />
            Reference read
          </div>
          <div>
            <span className="swatch alt" />
            Alt-supporting read
          </div>
          <div>
            <span className="swatch error" />
            SNV / mismatch
          </div>
          <div>
            <span className="swatch indel" />
            Indel
          </div>
          <div>
            <span className="swatch depth" />
            Coverage depth
          </div>
        </div>

        <div className="footer-note">
          <Gauge size={16} />
          <span>Simulated 1.2 kb locus · synthetic reads · teaching mode</span>
        </div>
      </aside>
    </main>
  );
}
