#!/usr/bin/env python3
"""broad_idr_summary.py — Broad-peak IDR reproducibility QC.

Counts peaks from true-replicate, pooled-pseudorep, and self-pseudorep IDR
thresholded broadPeak files, computes rescue/self-consistency ratios, and
writes a 15-column reproducibility summary TSV. The true-replicate thresholded
peak set is copied as the final validated peak output.

IDR broadPeak thresholded output has 17 columns (for 2 replicates):
  1-9: standard broadPeak (chrom, start, end, name, score, strand,
       signalValue, pValue, qValue)
  10-17: localIDR, globalIDR, rep1_start, rep1_end, rep1_signalValue,
         rep2_start, rep2_end, rep2_signalValue

Usage:
  python3 scripts/broad_idr_summary.py \
      --true-peaks <path> \
      --pooled-peaks <path> \
      --self1-peaks <path> \
      --self2-peaks <path> \
      --experiment <id> \
      --assay <chipseq|cuttag> \
      --caller macs3 \
      --peak-mode broad \
      --bio-rep-a <n> \
      --bio-rep-b <n> \
      --final-method idr \
      --final-output <path> \
      --output-tsv <path> \
      --output-peak <path>
"""

import argparse
import shutil


def count_peaks(path):
    """Count non-header, non-track lines in a broadPeak file."""
    count = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            count += 1
    return count


def compute_ratio(numerator, denominator):
    """Compute ratio with safe handling of zero denominators.

    Returns "NA" for 0/0, "inf" for N/0 (0 < N), or "{:.3f}" formatted float.
    """
    if denominator == 0:
        if numerator == 0:
            return "NA"
        return "inf"
    return "{:.3f}".format(numerator / denominator)


def main():
    parser = argparse.ArgumentParser(
        description="Broad-peak IDR reproducibility QC"
    )
    parser.add_argument("--true-peaks", required=True,
                        help="True-replicate IDR thresholded broadPeak")
    parser.add_argument("--pooled-peaks", required=True,
                        help="Pooled pseudorep IDR thresholded broadPeak")
    parser.add_argument("--self1-peaks", required=True,
                        help="Self-pseudorep A IDR thresholded broadPeak")
    parser.add_argument("--self2-peaks", required=True,
                        help="Self-pseudorep B IDR thresholded broadPeak")
    parser.add_argument("--experiment", required=True,
                        help="Experiment ID")
    parser.add_argument("--assay", required=True,
                        help="Assay type (chipseq or cuttag)")
    parser.add_argument("--caller", required=True,
                        help="Peak caller name")
    parser.add_argument("--peak-mode", required=True,
                        help="Peak mode (broad)")
    parser.add_argument("--bio-rep-a", required=True,
                        help="Biological replicate A label")
    parser.add_argument("--bio-rep-b", required=True,
                        help="Biological replicate B label")
    parser.add_argument("--final-method", default="idr",
                        help="Final method (default: idr)")
    parser.add_argument("--final-output", default="",
                        help="Path for final validated output")
    parser.add_argument("--output-tsv", required=True,
                        help="Path to write reproducibility summary TSV")
    parser.add_argument("--output-peak", required=True,
                        help="Path for final validated peak")
    args = parser.parse_args()

    # Count peaks
    Nt = count_peaks(args.true_peaks)
    Np = count_peaks(args.pooled_peaks)
    N1 = count_peaks(args.self1_peaks)
    N2 = count_peaks(args.self2_peaks)

    # Compute ratios
    rescue_numer = max(Np, Nt)
    rescue_denom = min(Np, Nt)
    rescue_ratio = compute_ratio(rescue_numer, rescue_denom)

    self_numer = max(N1, N2)
    self_denom = min(N1, N2)
    self_ratio = compute_ratio(self_numer, self_denom)

    # Reproducibility status: pass if both ratios < 2 and finite
    def _is_ok(r):
        return r not in ("NA", "inf") and float(r) < 2.0
    status = "pass" if _is_ok(rescue_ratio) and _is_ok(self_ratio) else "fail"

    # Copy true-replicate thresholded peaks to final output
    shutil.copyfile(args.true_peaks, args.output_peak)

    # Write 15-column reproducibility summary TSV
    header = [
        "experiment", "assay", "peak_mode", "caller",
        "bio_rep_a", "bio_rep_b",
        "true_peaks_Nt", "pooled_peaks_Np",
        "self1_peaks_N1", "self2_peaks_N2",
        "rescue_ratio", "self_consistency_ratio",
        "reproducibility_status", "final_method", "final_output",
    ]
    row = [
        args.experiment, args.assay, args.peak_mode, args.caller,
        args.bio_rep_a, args.bio_rep_b,
        str(Nt), str(Np), str(N1), str(N2),
        rescue_ratio, self_ratio, status,
        args.final_method, args.final_output,
    ]

    with open(args.output_tsv, "w") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(row) + "\n")

    print(
        "Broad-peak IDR reproducibility summary for %s (assay=%s): "
        "Nt=%d Np=%d N1=%d N2=%d "
        "rescue=%s self=%s status=%s"
        % (args.experiment, args.assay,
           Nt, Np, N1, N2,
           rescue_ratio, self_ratio, status)
    )


if __name__ == "__main__":
    main()
