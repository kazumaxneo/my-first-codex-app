# Genome Coverage Visualizer

An interactive teaching app for understanding short-read NGS coverage.

The first version simulates a 1.2 kb locus, places synthetic reads, and shows how read count, read length, sequencing error, GC bias, duplicate rate, and allele fraction change the visual evidence around a candidate SNP.

## Run locally

```bash
npm install
npm run dev
```

Open the local URL shown by Vite.

## Build

```bash
npm run build
```

## Deploy on Vercel

Import this GitHub repository into Vercel and keep the default Vite settings:

- Build command: `npm run build`
- Output directory: `dist`
