import { useMemo, useRef, type ChangeEvent } from 'react';
import {
  flexRender,
  getCoreRowModel,
  useReactTable,
  type ColumnDef,
} from '@tanstack/react-table';
import { FileUp, Plus, Trash2 } from 'lucide-react';
import { Button } from '../../components/Button';
import type { SampleColumn, WorkbenchSchema } from './schemaContract';
import { parseSampleTsv, type DraftSampleRow } from './sampleTsv';
import type { InputDraftController } from './useInputDraft';

interface SampleEditorProps {
  schema: WorkbenchSchema;
  draft: InputDraftController;
}

interface SampleCellProps {
  column: SampleColumn;
  row: DraftSampleRow;
  rowNumber: number;
  onChange: (value: string) => void;
}

function SampleCell({ column, row, rowNumber, onChange }: SampleCellProps) {
  const label = `Sample ${rowNumber} ${column.key}`;
  const value = row.values[column.key] ?? '';
  const className =
    'w-full min-w-32 rounded border border-[var(--color-border)] bg-white px-2 py-1.5 text-sm focus:border-[var(--color-accent)] focus:outline-none focus:ring-1 focus:ring-[var(--color-accent)]';
  if (column.enumValues) {
    const hasEmpty = column.enumValues.includes('');
    return (
      <select
        aria-label={label}
        className={className}
        value={value}
        onChange={(event) => onChange(event.target.value)}
      >
        {!hasEmpty && <option value="">Select…</option>}
        {column.enumValues.map((option) => (
          <option key={option || '__empty__'} value={option}>
            {option || 'Empty'}
          </option>
        ))}
      </select>
    );
  }
  return (
    <input
      type="text"
      aria-label={label}
      className={className}
      value={value}
      onChange={(event) => onChange(event.target.value)}
    />
  );
}

export function SampleEditor({ schema, draft }: SampleEditorProps) {
  const fileSelectionRevision = useRef(0);
  const columns = useMemo<ColumnDef<DraftSampleRow>[]>(
    () => [
      {
        id: 'row-number',
        header: '#',
        cell: ({ row }) => (
          <span className="font-mono text-xs text-[var(--color-text-muted)]">
            {row.index + 1}
          </span>
        ),
      },
      ...schema.sampleColumns.map<ColumnDef<DraftSampleRow>>((column) => ({
        id: column.key,
        header: () => (
          <span title={column.description ?? undefined}>
            {column.label}
            {column.required && <span className="text-[var(--color-error)]"> *</span>}
          </span>
        ),
        cell: ({ row }) => (
          <SampleCell
            column={column}
            row={row.original}
            rowNumber={row.index + 1}
            onChange={(value) =>
              draft.updateSample(row.original.id, column.key, value)
            }
          />
        ),
      })),
      {
        id: 'actions',
        header: 'Actions',
        cell: ({ row }) => (
          <Button
            type="button"
            className="h-8 w-8 p-0"
            aria-label={`Remove sample ${row.index + 1}`}
            title={`Remove sample ${row.index + 1}`}
            onClick={() => draft.removeSample(row.original.id)}
          >
            <Trash2 aria-hidden="true" size={15} />
          </Button>
        ),
      },
    ],
    [draft.removeSample, draft.updateSample, schema.sampleColumns],
  );
  const table = useReactTable({
    data: draft.state.rows,
    columns,
    getCoreRowModel: getCoreRowModel(),
    getRowId: (row) => row.id,
  });

  async function handleFile(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    event.target.value = '';
    if (!file) return;
    const selectionRevision = ++fileSelectionRevision.current;
    if (file.size > schema.limits.max_request_bytes) {
      draft.failSampleImport({
        code: 'SAMPLE_TSV_LIMIT_EXCEEDED',
        message: 'The sample file is larger than the workflow authoring limit.',
        row: null,
        column: null,
      });
      return;
    }
    let text: string;
    try {
      text = await file.text();
    } catch {
      if (selectionRevision !== fileSelectionRevision.current) return;
      draft.failSampleImport({
        code: 'SAMPLE_TSV_INVALID',
        message: 'The browser could not read this sample file.',
        row: null,
        column: null,
      });
      return;
    }
    if (selectionRevision !== fileSelectionRevision.current) return;
    const result = parseSampleTsv(text, schema);
    if (!result.ok) {
      draft.failSampleImport(result.issue);
      return;
    }
    draft.replaceSamples(result.rows, file.name);
  }

  return (
    <section className="min-w-0 space-y-3" aria-labelledby="sample-editor-title">
      <div className="flex flex-col gap-2 sm:flex-row sm:items-start sm:justify-between">
        <div>
          <h3 id="sample-editor-title" className="text-sm font-semibold">
            Inline sample rows
          </h3>
          <p className="mt-1 text-xs text-[var(--color-text-muted)]">
            Import tab-separated strings or edit the adapter-declared columns.
            Values are not scientifically inferred.
          </p>
        </div>
        <div className="flex flex-wrap gap-2">
          <label className="inline-flex cursor-pointer items-center gap-1.5 rounded border border-[var(--color-border)] bg-white px-3 py-1.5 text-sm font-medium hover:bg-[var(--color-bg)] focus-within:ring-2 focus-within:ring-[var(--color-accent)]">
            <FileUp aria-hidden="true" size={16} />
            Import TSV
            <input
              type="file"
              className="sr-only"
              accept=".tsv,text/tab-separated-values,text/plain"
              aria-label="Import samples TSV"
              onChange={handleFile}
            />
          </label>
          <Button
            type="button"
            className="gap-1.5"
            aria-label="Add sample row"
            disabled={draft.state.rows.length >= schema.limits.max_sample_rows}
            onClick={draft.addSample}
          >
            <Plus aria-hidden="true" size={16} />
            Add row
          </Button>
        </div>
      </div>

      <div className="flex min-h-8 flex-wrap items-center gap-x-3 gap-y-1 text-xs text-[var(--color-text-muted)]" aria-live="polite">
        <span>{draft.state.rows.length} rows</span>
        <span>{schema.sampleColumns.length} declared columns</span>
        {draft.state.importedFileName && (
          <span className="break-all">Imported {draft.state.importedFileName}</span>
        )}
      </div>
      {draft.state.sampleImportIssue && (
        <p
          className="rounded border border-red-200 bg-[var(--color-error-bg)] px-3 py-2 text-sm text-[var(--color-error)]"
          role="alert"
        >
          {draft.state.sampleImportIssue.message}
          {draft.state.sampleImportIssue.row !== null &&
            ` Row ${draft.state.sampleImportIssue.row}.`}
          {draft.state.sampleImportIssue.column !== null && (
            <>
              {' Column '}
              <code className="break-all font-mono">
                {draft.state.sampleImportIssue.column}
              </code>
              .
            </>
          )}
        </p>
      )}

      {draft.state.rows.length === 0 ? (
        <div className="flex min-h-48 items-center justify-center rounded border border-dashed border-[var(--color-border)] px-4 text-center text-sm text-[var(--color-text-muted)]">
          Import a TSV or add a row. At least one sample is required before this
          draft is structurally ready.
        </div>
      ) : (
        <>
          <div
            data-testid="sample-desktop-table"
            className="hidden min-w-0 max-w-full overflow-x-auto rounded border border-[var(--color-border)] md:block"
          >
            <table className="w-max min-w-full border-collapse text-left text-xs">
              <thead className="bg-[var(--color-bg)]">
                {table.getHeaderGroups().map((headerGroup) => (
                  <tr key={headerGroup.id}>
                    {headerGroup.headers.map((header) => (
                      <th
                        key={header.id}
                        className="whitespace-nowrap border-b border-r border-[var(--color-border)] px-2 py-2 font-semibold last:border-r-0"
                      >
                        {flexRender(header.column.columnDef.header, header.getContext())}
                      </th>
                    ))}
                  </tr>
                ))}
              </thead>
              <tbody>
                {table.getRowModel().rows.map((row) => (
                  <tr key={row.id} className="align-top">
                    {row.getVisibleCells().map((cell) => (
                      <td
                        key={cell.id}
                        className="border-b border-r border-[var(--color-border)] p-1.5 last:border-r-0"
                      >
                        {flexRender(cell.column.columnDef.cell, cell.getContext())}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          <div data-testid="sample-mobile-list" className="space-y-3 md:hidden">
            {draft.state.rows.map((row, rowIndex) => (
              <section
                key={row.id}
                className="space-y-3 border-b border-[var(--color-border)] pb-4 last:border-b-0"
                aria-labelledby={`sample-row-${row.id}`}
              >
                <div className="flex items-center justify-between">
                  <h4 id={`sample-row-${row.id}`} className="text-sm font-semibold">
                    Sample {rowIndex + 1}
                  </h4>
                  <Button
                    type="button"
                    className="h-8 w-8 p-0"
                    aria-label={`Remove sample ${rowIndex + 1}`}
                    title={`Remove sample ${rowIndex + 1}`}
                    onClick={() => draft.removeSample(row.id)}
                  >
                    <Trash2 aria-hidden="true" size={15} />
                  </Button>
                </div>
                <div className="grid min-w-0 gap-3 sm:grid-cols-2">
                  {schema.sampleColumns.map((column) => (
                    <label key={column.key} className="min-w-0 space-y-1 text-xs font-medium">
                      <span className="block break-words">
                        {column.label}
                        {column.required && (
                          <span className="text-[var(--color-error)]"> *</span>
                        )}
                      </span>
                      <SampleCell
                        column={column}
                        row={row}
                        rowNumber={rowIndex + 1}
                        onChange={(value) =>
                          draft.updateSample(row.id, column.key, value)
                        }
                      />
                    </label>
                  ))}
                </div>
              </section>
            ))}
          </div>
        </>
      )}
    </section>
  );
}
