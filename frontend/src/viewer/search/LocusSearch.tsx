/**
 * Interface seam for gene and coordinate search. IRENE-65 owns the offline
 * index that replaces this file. Keep that implementation if the two commits
 * conflict. The input stays controlled by `region` and announces a new window
 * through `onNavigate`.
 */
export default function LocusSearch({
  assembly,
  region,
  onNavigate,
}: {
  assembly: string;
  region: string;
  onNavigate: (region: string) => void;
}) {
  return (
    <label className="block text-sm">
      Genomic region
      <input
        className="mt-1 w-full max-w-md border border-[var(--color-border)] px-2 py-1"
        value={region}
        onChange={(event) => onNavigate(event.target.value.trim())}
        aria-label="Genomic region"
        data-assembly={assembly}
      />
    </label>
  );
}
