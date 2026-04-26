import { useEffect, useMemo, useRef, useState } from 'react';
import {
  Activity,
  Dna,
  Gauge,
  Pause,
  Play,
  RotateCcw,
  SlidersHorizontal,
  Sparkles,
} from 'lucide-react';

type Settings = {
  readCount: number;
  readLength: number;
  errorRate: number;
  gcBias: number;
  duplicateRate: number;
  variantAlleleFraction: number;
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
};

type Simulation = {
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
};

const initialSettings: Settings = {
  readCount: 180,
  readLength: 76,
  errorRate: 1.2,
  gcBias: 34,
  duplicateRate: 8,
  variantAlleleFraction: 50,
  seed: 18,
  playing: false,
};

const genomeLength = 1200;
const variantPosition = 690;

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
    const hasAlt = start <= variantPosition && end >= variantPosition && random() * 100 < settings.variantAlleleFraction;
    const mismatchPositions: number[] = [];
    const errorEvents = Math.floor((settings.errorRate / 100) * settings.readLength + random() * 2);
    for (let i = 0; i < errorEvents; i += 1) {
      mismatchPositions.push(start + Math.floor(random() * settings.readLength));
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
    };
    reads.push(read);

    for (let pos = start; pos < end && pos < genomeLength; pos += 1) {
      coverage[pos] += 1;
      if (hasAlt && Math.abs(pos - variantPosition) < 2) altCoverage[pos] += 1;
    }
  };

  for (let i = 0; i < uniqueReads; i += 1) placeRead(i);
  for (let i = 0; i < duplicates; i += 1) {
    const source = reads[Math.floor(random() * Math.max(reads.length, 1))];
    placeRead(uniqueReads + i, source?.start, true);
  }

  const maxCoverage = Math.max(...coverage, 1);
  const variantDepth = coverage[variantPosition];
  const observedAltFraction = variantDepth ? altCoverage[variantPosition] / variantDepth : 0;
  const confidence =
    variantDepth >= 18 && observedAltFraction > 0.22
      ? 'Strong'
      : variantDepth >= 8 && observedAltFraction > 0.12
        ? 'Watch'
        : 'Low';

  return {
    genomeLength,
    variantPosition,
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

function CoverageCanvas({ simulation }: { simulation: Simulation }) {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);

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
      const left = 34;
      const right = 22;
      const trackWidth = width - left - right;
      const xFor = (pos: number) => left + (pos / simulation.genomeLength) * trackWidth;

      context.clearRect(0, 0, width, height);
      const bg = context.createLinearGradient(0, 0, width, height);
      bg.addColorStop(0, '#07110f');
      bg.addColorStop(0.5, '#11130f');
      bg.addColorStop(1, '#150d17');
      context.fillStyle = bg;
      context.fillRect(0, 0, width, height);

      for (let i = 0; i < 90; i += 1) {
        const x = (i * 79) % width;
        const y = (i * 47) % height;
        context.fillStyle = `rgba(230, 255, 238, ${0.04 + (i % 5) * 0.015})`;
        context.fillRect(x, y, 1, 1);
      }

      context.font = '600 11px Inter, system-ui, sans-serif';
      context.fillStyle = '#8ca096';
      context.fillText('GC%', 8, 41);
      context.fillText('READS', 8, 136);
      context.fillText('DEPTH', 8, 344);

      const gcY = 28;
      const gcHeight = 62;
      context.strokeStyle = 'rgba(180, 202, 188, 0.18)';
      context.lineWidth = 1;
      context.strokeRect(left, gcY, trackWidth, gcHeight);
      context.beginPath();
      simulation.gc.forEach((gc, pos) => {
        const x = xFor(pos);
        const y = gcY + gcHeight - gc * gcHeight;
        if (pos === 0) context.moveTo(x, y);
        else context.lineTo(x, y);
      });
      context.strokeStyle = '#7df0b2';
      context.lineWidth = 2;
      context.stroke();

      const readY = 108;
      const rowHeight = 9;
      const maxRows = Math.floor((height - 245) / rowHeight);
      simulation.reads.forEach((read) => {
        if (read.row >= maxRows) return;
        const x = xFor(read.start);
        const w = Math.max(2, xFor(read.end) - x);
        const y = readY + read.row * rowHeight;
        const alpha = read.duplicate ? 0.34 : 0.78;
        context.fillStyle = read.hasAlt ? `rgba(255, 202, 87, ${alpha})` : `rgba(88, 199, 178, ${alpha})`;
        drawRoundedRect(context, x, y, w, 5, 2);
        context.fill();
        if (read.strand === -1) {
          context.fillStyle = 'rgba(10, 12, 11, 0.5)';
          context.fillRect(x + 2, y + 2, Math.max(1, w - 4), 1);
        }
        context.fillStyle = '#ff637d';
        read.mismatchPositions.slice(0, 5).forEach((pos) => {
          const mx = xFor(pos);
          if (mx >= x && mx <= x + w) context.fillRect(mx, y - 1, 2, 7);
        });
      });

      const coverageBase = height - 76;
      const coverageHeight = 118;
      context.fillStyle = 'rgba(255, 255, 255, 0.04)';
      context.fillRect(left, coverageBase - coverageHeight, trackWidth, coverageHeight);
      for (let i = 0; i < simulation.coverage.length; i += 3) {
        const depth = simulation.coverage[i];
        const normalized = depth / simulation.maxCoverage;
        const barHeight = normalized * coverageHeight;
        const hue = 165 - normalized * 110;
        context.fillStyle = `hsl(${hue}, 78%, ${38 + normalized * 20}%)`;
        context.fillRect(xFor(i), coverageBase - barHeight, Math.max(1, trackWidth / simulation.genomeLength * 3), barHeight);
      }

      const vx = xFor(simulation.variantPosition);
      context.strokeStyle = '#ffcf5a';
      context.lineWidth = 2;
      context.beginPath();
      context.moveTo(vx, 18);
      context.lineTo(vx, height - 32);
      context.stroke();
      context.fillStyle = '#ffcf5a';
      context.beginPath();
      context.arc(vx, coverageBase - 126, 7, 0, Math.PI * 2);
      context.fill();
      context.fillStyle = '#160d06';
      context.font = '800 10px Inter, system-ui, sans-serif';
      context.fillText('SNP', vx - 9, coverageBase - 122);

      context.fillStyle = 'rgba(230, 255, 238, 0.72)';
      context.font = '600 11px Inter, system-ui, sans-serif';
      context.fillText('0 bp', left, height - 20);
      context.fillText(`${simulation.genomeLength} bp`, width - 82, height - 20);
    };

    draw();
    const resizeObserver = new ResizeObserver(draw);
    resizeObserver.observe(canvas);
    return () => resizeObserver.disconnect();
  }, [simulation]);

  return <canvas className="coverage-canvas" ref={canvasRef} aria-label="Genome coverage visualization" />;
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

  useEffect(() => {
    if (!settings.playing) return;
    const interval = window.setInterval(() => {
      setSettings((current) => ({ ...current, seed: current.seed + 1 }));
    }, 900);
    return () => window.clearInterval(interval);
  }, [settings.playing]);

  const simulation = useMemo(() => simulate(settings), [settings]);
  const duplicateReads = simulation.reads.filter((read) => read.duplicate).length;
  const confidenceTone = simulation.confidence === 'Strong' ? 'green' : simulation.confidence === 'Watch' ? 'amber' : 'rose';

  const update = (patch: Partial<Settings>) => setSettings((current) => ({ ...current, ...patch }));

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
              onClick={() => update({ seed: settings.seed + 11 })}
            >
              <Sparkles size={18} />
            </button>
            <button className="icon-button" type="button" title="Reset controls" onClick={() => setSettings(initialSettings)}>
              <RotateCcw size={18} />
            </button>
          </div>
        </header>

        <div className="visual-stage">
          <CoverageCanvas simulation={simulation} />
        </div>

        <section className="metrics-row" aria-label="Simulation metrics">
          <Metric label="Mean depth" value={`${(settings.readCount * settings.readLength / genomeLength).toFixed(1)}x`} tone="green" />
          <Metric label="Max depth" value={`${simulation.maxCoverage}x`} />
          <Metric label="SNP depth" value={`${simulation.variantDepth}x`} tone={confidenceTone} />
          <Metric label="Alt fraction" value={`${Math.round(simulation.observedAltFraction * 100)}%`} tone={confidenceTone} />
          <Metric label="Duplicates" value={`${duplicateReads}`} tone="amber" />
        </section>
      </section>

      <aside className="side-panel">
        <div className="panel-heading">
          <SlidersHorizontal size={18} />
          <h2>Library Model</h2>
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
            黄色の縦線は候補SNPです。周囲のリード数、alt alleleを持つリードの比率、赤いミスマッチの増え方で、低深度や高エラーが判定に与える影響を比較できます。
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
            Sequencing error
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
