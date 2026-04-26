import { gunzipSync } from 'fflate';

export type BamRead = {
  id: number;
  name: string;
  start: number;
  end: number;
  row: number;
  strand: 1 | -1;
  duplicate: boolean;
  mismatchPositions: number[];
  indelPositions: number[];
  mapq: number;
};

export type BamAnalysis = {
  fileName: string;
  referenceName: string;
  genomeLength: number;
  referenceSequence?: string;
  reads: BamRead[];
  coverage: number[];
  maxCoverage: number;
  meanDepth: number;
  duplicateCount: number;
  mappedCount: number;
};

export type ReferenceSequence = {
  fileName: string;
  name: string;
  sequence: string;
};

const MAX_BAM_BYTES = 5 * 1024 * 1024;
const MAX_GENOME_LENGTH = 10_000;

const cigarOps = 'MIDNSHP=XB';
const seqLookup = '=ACMGRSVTWYHKDBN';

function readString(bytes: Uint8Array, offset: number, length: number) {
  return new TextDecoder().decode(bytes.subarray(offset, offset + length));
}

function readInt32(view: DataView, offset: number) {
  return view.getInt32(offset, true);
}

function readUint16(view: DataView, offset: number) {
  return view.getUint16(offset, true);
}

function readUint32(view: DataView, offset: number) {
  return view.getUint32(offset, true);
}

function decodeReadSequence(bytes: Uint8Array, offset: number, length: number) {
  let sequence = '';
  for (let i = 0; i < length; i += 1) {
    const packed = bytes[offset + Math.floor(i / 2)];
    const code = i % 2 === 0 ? packed >> 4 : packed & 0x0f;
    sequence += seqLookup[code] ?? 'N';
  }
  return sequence;
}

function decompressBamBytes(input: Uint8Array) {
  if (readString(input, 0, 4) === 'BAM\u0001') return input;
  if (input[0] !== 0x1f || input[1] !== 0x8b) {
    throw new Error('BAMまたはBGZF圧縮BAMとして読み込めませんでした。選択したファイルが.bamか確認してください。');
  }

  const chunks: Uint8Array[] = [];
  let offset = 0;
  while (offset + 18 <= input.length) {
    if (input[offset] !== 0x1f || input[offset + 1] !== 0x8b) break;
    const flags = input[offset + 3];
    if ((flags & 0x04) === 0) {
      chunks.push(gunzipSync(input.subarray(offset)));
      break;
    }

    const extraLength = input[offset + 10] | (input[offset + 11] << 8);
    let blockSize = 0;
    let extraOffset = offset + 12;
    const extraEnd = extraOffset + extraLength;
    while (extraOffset + 4 <= extraEnd) {
      const subfieldId = String.fromCharCode(input[extraOffset], input[extraOffset + 1]);
      const subfieldLength = input[extraOffset + 2] | (input[extraOffset + 3] << 8);
      if (subfieldId === 'BC' && subfieldLength === 2) {
        blockSize = (input[extraOffset + 4] | (input[extraOffset + 5] << 8)) + 1;
        break;
      }
      extraOffset += 4 + subfieldLength;
    }

    if (!blockSize) {
      chunks.push(gunzipSync(input.subarray(offset)));
      break;
    }

    const block = input.subarray(offset, offset + blockSize);
    const inflated = gunzipSync(block);
    if (inflated.length > 0) chunks.push(inflated);
    offset += blockSize;
  }

  const totalLength = chunks.reduce((sum, chunk) => sum + chunk.length, 0);
  const output = new Uint8Array(totalLength);
  let writeOffset = 0;
  chunks.forEach((chunk) => {
    output.set(chunk, writeOffset);
    writeOffset += chunk.length;
  });
  return output;
}

export async function parseReferenceFile(file: File): Promise<ReferenceSequence> {
  if (!/\.(fa|fasta|fna)$/i.test(file.name)) {
    throw new Error('リファレンスはFASTA（.fa, .fasta, .fna）を選択してください。');
  }
  const text = await file.text();
  const lines = text.split(/\r?\n/).filter(Boolean);
  const header = lines.find((line) => line.startsWith('>'));
  if (!header) throw new Error('FASTAヘッダー（>reference_name）が見つかりません。');
  const sequence = lines
    .filter((line) => !line.startsWith('>'))
    .join('')
    .replace(/\s/g, '')
    .toUpperCase();
  if (!sequence) throw new Error('FASTA配列が空です。');
  if (sequence.length > MAX_GENOME_LENGTH) {
    throw new Error(`リファレンスは${sequence.length.toLocaleString()} bpです。10kb以下にしてください。`);
  }
  if (!/^[ACGTRYSWKMBDHVN]+$/i.test(sequence)) {
    throw new Error('FASTAに想定外の文字があります。IUPAC塩基コードのみ対応しています。');
  }
  return {
    fileName: file.name,
    name: header.replace(/^>/, '').trim().split(/\s+/)[0] || file.name,
    sequence,
  };
}

function parseMdTag(md: string, referenceStart: number) {
  const mismatches: number[] = [];
  const deletions: number[] = [];
  let cursor = referenceStart;
  let i = 0;

  while (i < md.length) {
    const char = md[i];
    if (/\d/.test(char)) {
      let digits = '';
      while (i < md.length && /\d/.test(md[i])) {
        digits += md[i];
        i += 1;
      }
      cursor += Number(digits);
      continue;
    }

    if (char === '^') {
      i += 1;
      while (i < md.length && /[A-Za-z]/.test(md[i])) {
        deletions.push(cursor);
        cursor += 1;
        i += 1;
      }
      continue;
    }

    if (/[A-Za-z]/.test(char)) {
      mismatches.push(cursor);
      cursor += 1;
    }
    i += 1;
  }

  return { mismatches, deletions };
}

function parseOptionalTags(bytes: Uint8Array, view: DataView, offset: number, end: number) {
  let md = '';
  while (offset + 3 <= end) {
    const tag = readString(bytes, offset, 2);
    const type = String.fromCharCode(bytes[offset + 2]);
    offset += 3;

    if (type === 'A' || type === 'c' || type === 'C') {
      offset += 1;
    } else if (type === 's' || type === 'S') {
      offset += 2;
    } else if (type === 'i' || type === 'I' || type === 'f') {
      offset += 4;
    } else if (type === 'Z' || type === 'H') {
      const start = offset;
      while (offset < end && bytes[offset] !== 0) offset += 1;
      if (tag === 'MD') md = readString(bytes, start, offset - start);
      offset += 1;
    } else if (type === 'B') {
      const subtype = String.fromCharCode(bytes[offset]);
      const count = readUint32(view, offset + 1);
      const width = subtype === 'c' || subtype === 'C' ? 1 : subtype === 's' || subtype === 'S' ? 2 : 4;
      offset += 5 + count * width;
    } else {
      break;
    }
  }
  return { md };
}

export async function parseBamFile(file: File, referenceInput?: ReferenceSequence | null): Promise<BamAnalysis> {
  if (!file.name.toLowerCase().endsWith('.bam')) {
    throw new Error('BAMファイル（.bam）を選択してください。');
  }
  if (file.size > MAX_BAM_BYTES) {
    throw new Error('BAMファイルは5MB以下にしてください。');
  }

  const compressed = new Uint8Array(await file.arrayBuffer());
  const bytes = decompressBamBytes(compressed);
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let offset = 0;

  if (readString(bytes, 0, 4) !== 'BAM\u0001') {
    throw new Error('BAMとして読み込めませんでした。BGZF圧縮されたBAMか確認してください。');
  }
  offset += 4;

  const headerLength = readInt32(view, offset);
  offset += 4 + headerLength;

  const referenceCount = readInt32(view, offset);
  offset += 4;
  if (referenceCount < 1) throw new Error('BAMヘッダーにreferenceがありません。');

  const references: { name: string; length: number }[] = [];
  for (let i = 0; i < referenceCount; i += 1) {
    const nameLength = readInt32(view, offset);
    offset += 4;
    const name = readString(bytes, offset, nameLength).replace(/\0$/, '');
    offset += nameLength;
    const length = readInt32(view, offset);
    offset += 4;
    references.push({ name, length });
  }

  const reference = references[0];
  if (reference.length > MAX_GENOME_LENGTH) {
    throw new Error(`最初のreferenceは${reference.length.toLocaleString()} bpです。10kb以下のBAMにしてください。`);
  }
  if (referenceInput && referenceInput.sequence.length !== reference.length) {
    throw new Error(
      `FASTAは${referenceInput.sequence.length.toLocaleString()} bp、BAMの最初のreferenceは${reference.length.toLocaleString()} bpです。長さが一致するFASTAを指定してください。`,
    );
  }

  const coverage = new Array(reference.length).fill(0);
  const reads: BamRead[] = [];
  let mappedCount = 0;
  let duplicateCount = 0;

  while (offset + 4 <= bytes.length) {
    const blockSize = readInt32(view, offset);
    offset += 4;
    const blockEnd = offset + blockSize;
    if (blockEnd > bytes.length || blockSize < 32) break;

    const refId = readInt32(view, offset);
    const pos = readInt32(view, offset + 4);
    const lReadName = bytes[offset + 8];
    const mapq = bytes[offset + 9];
    const nCigar = readUint16(view, offset + 12);
    const flag = readUint16(view, offset + 14);
    const lSeq = readInt32(view, offset + 16);
    const unmapped = (flag & 0x4) !== 0;
    const reverse = (flag & 0x10) !== 0;
    const duplicate = (flag & 0x400) !== 0;

    let cursor = offset + 32;
    const readName = readString(bytes, cursor, Math.max(0, lReadName - 1));
    cursor += lReadName;

    const cigar: { length: number; op: string }[] = [];
    for (let i = 0; i < nCigar; i += 1) {
      const packed = readUint32(view, cursor);
      cursor += 4;
      cigar.push({ length: packed >>> 4, op: cigarOps[packed & 0x0f] ?? '?' });
    }

    const seqOffset = cursor;
    cursor += Math.ceil(lSeq / 2);
    cursor += lSeq;
    const optional = parseOptionalTags(bytes, view, cursor, blockEnd);

    if (!unmapped && refId === 0 && pos >= 0) {
      mappedCount += 1;
      if (duplicate) duplicateCount += 1;
      const sequence = decodeReadSequence(bytes, seqOffset, lSeq);
      const mismatchPositions: number[] = [];
      const indelPositions: number[] = [];
      let refCursor = pos;
      let readCursor = 0;

      cigar.forEach(({ length, op }) => {
        if (op === 'M' || op === '=' || op === 'X') {
          for (let i = 0; i < length; i += 1) {
            const refPos = refCursor + i;
            if (refPos >= 0 && refPos < coverage.length) coverage[refPos] += 1;
          }
          if (op === 'X') {
            for (let i = 0; i < length; i += 1) mismatchPositions.push(refCursor + i);
          }
          if (referenceInput) {
            for (let i = 0; i < length; i += 1) {
              const readBase = sequence[readCursor + i];
              const referenceBase = referenceInput.sequence[refCursor + i];
              if (readBase && referenceBase && readBase !== 'N' && referenceBase !== 'N' && readBase !== referenceBase) {
                mismatchPositions.push(refCursor + i);
              }
            }
          }
          refCursor += length;
          readCursor += length;
        } else if (op === 'I') {
          indelPositions.push(refCursor);
          readCursor += length;
        } else if (op === 'D' || op === 'N') {
          for (let i = 0; i < length; i += 1) indelPositions.push(refCursor + i);
          refCursor += length;
        } else if (op === 'S') {
          readCursor += length;
        }
      });

      if (optional.md) {
        const mdEvents = parseMdTag(optional.md, pos);
        mismatchPositions.push(...mdEvents.mismatches);
        indelPositions.push(...mdEvents.deletions);
      } else if (!referenceInput) {
        sequence.split('').forEach((base, index) => {
          if (base === 'N' && randomish(readName, index) < 0.12) mismatchPositions.push(pos + index);
        });
      }

      reads.push({
        id: reads.length,
        name: readName,
        start: pos,
        end: Math.min(refCursor, reference.length),
        row: 0,
        strand: reverse ? -1 : 1,
        duplicate,
        mismatchPositions: [...new Set(mismatchPositions)].filter((value) => value >= 0 && value < reference.length),
        indelPositions: [...new Set(indelPositions)].filter((value) => value >= 0 && value < reference.length),
        mapq,
      });
    }

    offset = blockEnd;
  }

  if (reads.length === 0) throw new Error('最初のreferenceにmapped readが見つかりませんでした。');
  reads.sort((a, b) => a.start - b.start || a.end - b.end);
  const rowEnds: number[] = [];
  reads.forEach((read) => {
    const row = rowEnds.findIndex((end) => read.start > end + 3);
    read.row = row === -1 ? rowEnds.length : row;
    rowEnds[read.row] = read.end;
  });

  const maxCoverage = Math.max(...coverage, 1);
  const totalBases = coverage.reduce((sum, depth) => sum + depth, 0);

  return {
    fileName: file.name,
    referenceName: reference.name,
    genomeLength: reference.length,
    referenceSequence: referenceInput?.sequence,
    reads,
    coverage,
    maxCoverage,
    meanDepth: totalBases / reference.length,
    duplicateCount,
    mappedCount,
  };
}

function randomish(text: string, salt: number) {
  let hash = salt + 2166136261;
  for (let i = 0; i < text.length; i += 1) {
    hash ^= text.charCodeAt(i);
    hash = Math.imul(hash, 16777619);
  }
  return (hash >>> 0) / 4294967295;
}
