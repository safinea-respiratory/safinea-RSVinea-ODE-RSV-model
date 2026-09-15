/* ========================================================== *
 * REPORT TABLES -> WORD
 *
 * Builds a .docx of the values behind figures 7.4 to 7.13, one
 * table per figure, ready to copy into the report.
 *
 * Input : <report_figures>/figure_tables_wide.csv
 *         (written by report_figures.R)
 * Output: <report_figures>/figure_tables.docx
 *
 * Usage : node report_tables_docx.js [path/to/report_figures]
 * ========================================================== */

const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, Table, TableRow, TableCell,
  TextRun, HeadingLevel, WidthType, AlignmentType, ShadingType, BorderStyle,
} = require('docx');

const FIGDIR = process.argv[2] ||
  'output/IRL-RespiCompass-2026/3_results/IE/report_figures';
const SRC = path.join(FIGDIR, 'figure_tables_wide.csv');
const OUT = path.join(FIGDIR, 'figure_tables.docx');

// ---- minimal CSV reader (quoted fields, no embedded newlines) -------------
function parseCsv(text) {
  const lines = text.replace(/\r\n/g, '\n').trim().split('\n');
  const split = (line) => {
    const out = []; let cur = ''; let q = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (c === '"') { if (q && line[i + 1] === '"') { cur += '"'; i++; } else q = !q; }
      else if (c === ',' && !q) { out.push(cur); cur = ''; }
      else cur += c;
    }
    out.push(cur); return out;
  };
  const head = split(lines[0]);
  return lines.slice(1).map((l) => Object.fromEntries(split(l).map((v, i) => [head[i], v])));
}

const rows = parseCsv(fs.readFileSync(SRC, 'utf8'));

// Figure titles, and the order they appear in the report
const TITLES = {
  '7.4':  'Seasonal burden by model and age group, compared with the observed data',
  '7.5':  'Seasonal burden under the no vaccination scenario',
  '7.6':  'Relative change under the no vaccination scenario, compared with baseline',
  '7.7':  'Relative risk under the no vaccination scenario, compared with baseline',
  '7.8':  'Seasonal burden under the high vaccination uptake scenario',
  '7.9':  'Relative change under the high vaccination uptake scenario, compared with baseline',
  '7.10': 'Relative risk under the high vaccination uptake scenario, compared with baseline',
  '7.11': 'Seasonal burden under the catch-up scenarios',
  '7.12': 'Relative change under the catch-up scenarios, compared with baseline',
  '7.13': 'Relative risk under the catch-up scenarios, compared with baseline',
};
const ORDER = ['7.4', '7.5', '7.6', '7.7', '7.8', '7.9', '7.10',
               '7.11', '7.12', '7.13'];

// ---- table construction ---------------------------------------------------
const TABLE_W = 9360;                       // DXA, fits A4 with 1" margins
const HEAD_BG = 'D9E2F3';
const thin = { style: BorderStyle.SINGLE, size: 4, color: '9CA3AF' };
const BORDERS = { top: thin, bottom: thin, left: thin, right: thin,
                  insideHorizontal: thin, insideVertical: thin };

// Range separator. An en dash is fine for positive bounds (693–694) but
// unreadable once a bound is negative (-0.3--0.1), which is common in the
// relative-change tables, so use an explicit "to" in that case.
const range = (s) => {
  const m = s.match(/^(.*) \((-?[\d.]+)-(-?[\d.]+)\)$/);
  if (!m) return s;
  const [, med, lo, hi] = m;
  const sep = (lo.startsWith('-') || hi.startsWith('-')) ? ' to ' : '–';
  return `${med} (${lo}${sep}${hi})`;
};

function cell(text, { bold = false, header = false, width, align } = {}) {
  return new TableCell({
    width: { size: width, type: WidthType.DXA },
    shading: header ? { type: ShadingType.CLEAR, fill: HEAD_BG } : undefined,
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children: [new Paragraph({
      alignment: align || (header ? AlignmentType.CENTER : AlignmentType.LEFT),
      children: [new TextRun({ text, bold: bold || header, size: 18 })],
    })],
  });
}

function buildTable(figRows) {
  // Only include model columns that this figure actually has
  const models = ['Dynamic', 'Static', 'Observed']
    .filter((m) => figRows.some((r) => r[m] && r[m].length));

  // Key columns: show Season and/or Scenario only where they vary within the
  // figure. The season figures (7.4-7.10) vary by season at one scenario; the
  // catch-up figures (7.11-7.13) vary by scenario within one season.
  const varies = (k) => new Set(figRows.map((r) => r[k])).size > 1;
  const keys = ['age_band', ...(varies('season') ? ['season'] : []),
                             ...(varies('scenario') ? ['scenario'] : [])];
  const KEYLAB = { age_band: 'Age group', season: 'Season', scenario: 'Scenario' };
  const cols = [...keys.map((k) => KEYLAB[k]), ...models];

  const keyW = keys.map((k) => (k === 'age_band' ? 1900 : (k === 'season' ? 1500 : 2400)));
  const used = keyW.reduce((a, b) => a + b, 0);
  const rest = Math.floor((TABLE_W - used) / models.length);
  const widths = [...keyW, ...models.map(() => rest)];

  const header = new TableRow({
    tableHeader: true,
    children: cols.map((c, i) => cell(c, { header: true, width: widths[i] })),
  });

  const body = figRows.map((r, idx) => {
    // repeat the age group only on its first row, as in a published table
    const showBand = idx === 0 || figRows[idx - 1]['age_band'] !== r['age_band'];
    const keyCells = keys.map((k, i) =>
      cell(k === 'age_band' ? (showBand ? r[k] : '') : (r[k] || ''),
           { bold: k === 'age_band', width: widths[i] }));
    return new TableRow({
      children: [
        ...keyCells,
        ...models.map((m, i) => cell(range(r[m] || ''), {
          width: widths[keys.length + i], align: AlignmentType.RIGHT })),
      ],
    });
  });

  return new Table({
    columnWidths: widths,
    width: { size: TABLE_W, type: WidthType.DXA },
    borders: BORDERS,
    rows: [header, ...body],
  });
}

// ---- document -------------------------------------------------------------
const children = [
  new Paragraph({ text: 'RSV model comparison: figure values', heading: HeadingLevel.HEADING_1 }),
  new Paragraph({ children: [new TextRun({
    text: 'Values behind figures 7.4 to 7.13. Each cell is the median with the '
        + '5th–95th percentile in brackets. Observed data has no uncertainty '
        + 'interval. Generated from figure_tables_wide.csv by report_figures.R.',
    italics: true, size: 18 })] }),
  new Paragraph({ text: '' }),
];

for (const fig of ORDER) {
  const figRows = rows.filter((r) => r['figure'] === fig);
  if (!figRows.length) { console.warn(`no rows for figure ${fig}`); continue; }
  const quantity = figRows[0]['quantity'];
  const scenarios = [...new Set(figRows.map((r) => r['scenario']))];
  const seasons = [...new Set(figRows.map((r) => r['season']))];
  // Name the scenario only when the whole table is one scenario; otherwise it
  // is a column. Same for the season.
  const note = [quantity,
                scenarios.length === 1 ? scenarios[0] : null,
                seasons.length === 1 ? seasons[0] : null]
               .filter(Boolean).join(' — ');

  children.push(new Paragraph({
    text: `Figure ${fig}. ${TITLES[fig] || ''}`, heading: HeadingLevel.HEADING_2 }));
  children.push(new Paragraph({ children: [new TextRun({
    text: note, italics: true, size: 18 })] }));
  children.push(buildTable(figRows));
  children.push(new Paragraph({ text: '' }));
}

const doc = new Document({
  styles: { default: { document: { run: { font: 'Calibri', size: 20 } } } },
  sections: [{
    properties: { page: { margin: { top: 1000, bottom: 1000, left: 1000, right: 1000 } } },
    children,
  }],
});

Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync(OUT, buf);
  console.log(`written: ${OUT} (${rows.length} rows across ${ORDER.length} tables)`);
});
