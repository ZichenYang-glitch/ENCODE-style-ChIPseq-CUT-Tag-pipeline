"""Shared blacklist interval helpers for workflow scripts.

Peak files and BAM template filters both use 0-based half-open intervals after
conversion; blacklist BED files are already 0-based half-open.
"""

import bisect


class Blacklist:
    """Merged, sorted blacklist intervals with fast overlap lookup."""

    def __init__(self, path):
        raw = {}
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                line_str = line.strip()
                if not line_str or line_str.startswith(("#", "track", "browser")):
                    continue
                fields = line_str.split("\t")
                if len(fields) < 3:
                    fields = line_str.split()
                if len(fields) < 3:
                    continue
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError:
                    continue
                if start < 0 or end <= start:
                    continue
                raw.setdefault(fields[0], []).append((start, end))

        self.intervals = {}
        self.starts = {}
        for chrom, intervals in raw.items():
            intervals.sort()
            merged = []
            for start, end in intervals:
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self.intervals[chrom] = merged
            self.starts[chrom] = [item[0] for item in merged]
        self.total_intervals = sum(len(ivs) for ivs in self.intervals.values())

    def overlaps(self, chrom, start, end):
        """Return True when [start, end) overlaps the blacklist."""
        intervals = self.intervals.get(chrom)
        if not intervals or end <= start:
            return False
        starts = self.starts[chrom]
        index = bisect.bisect_left(starts, end) - 1
        if index < 0:
            return False
        return intervals[index][1] > start

    def peak_line_overlaps(self, line):
        """Return True when a BED-like peak line overlaps the blacklist."""
        line_str = line.strip()
        if not line_str:
            return False
        fields = line_str.split("\t")
        if len(fields) < 3:
            fields = line_str.split()
        if len(fields) < 3:
            return False
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError:
            return False
        return self.overlaps(fields[0], start, end)

    def filter_peak_file(self, input_path, output_path):
        """Write non-overlapping peak lines; return (kept, dropped)."""
        kept = 0
        dropped = 0
        with (
            open(input_path, encoding="utf-8") as source,
            open(output_path, "w", encoding="utf-8") as target,
        ):
            for line in source:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    target.write(line)
                    continue
                if self.peak_line_overlaps(line):
                    dropped += 1
                else:
                    target.write(line)
                    kept += 1
        return kept, dropped
