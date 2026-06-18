#!/usr/bin/env python3
"""Stage 55 ATAC IDR reproducibility QC summary and final peak copy.

Reads ATAC IDR thresholded peak files, computes rescue ratio and
self-consistency ratio, writes a 15-column reproducibility summary TSV,
and copies the true-replicate IDR thresholded peak set as the final
replicate-validated output.

Usage:
    python3 scripts/atac_idr_summary.py \\
        --true-peaks <idr/thresholded.narrowPeak> \\
        --pooled-peaks <idr/pooled.thresholded.narrowPeak> \\
        --self1-peaks <idr/self1.thresholded.narrowPeak> \\
        --self2-peaks <idr/self2.thresholded.narrowPeak> \\
        --experiment <exp> \\
        --assay atac \\
        --caller macs3 \\
        --peak-mode narrow \\
        --bio-rep-a <label> --bio-rep-b <label> \\
        --final-method idr \\
        --final-output <final.narrowPeak> \\
        --output-tsv <summary.tsv> \\
        --output-peak <final.narrowPeak>
"""

import argparse
import shutil
import sys


def count_peaks(path):
    """Count lines in a peak file, excluding headers."""
    count = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#") and not line.startswith("track"):
                count += 1
    return count


def compute_ratio(numerator, denominator):
    """Return ratio as a formatted string: float, 'inf', or 'NA'."""
    if denominator == 0:
        if numerator == 0:
            return "NA"
        return "inf"
    return "{:.3f}".format(numerator / denominator)


def main():
    parser = argparse.ArgumentParser(
        description="Stage 55 ATAC IDR reproducibility summary"
    )
    parser.add_argument("--true-peaks", required=True)
    parser.add_argument("--pooled-peaks", required=True)
    parser.add_argument("--self1-peaks", required=True)
    parser.add_argument("--self2-peaks", required=True)
    parser.add_argument("--experiment", required=True)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--caller", required=True)
    parser.add_argument("--peak-mode", required=True)
    parser.add_argument("--bio-rep-a", required=True)
    parser.add_argument("--bio-rep-b", required=True)
    parser.add_argument("--final-method", default="idr")
    parser.add_argument("--final-output", default="")
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-peak", required=True)
    args = parser.parse_args()

    Nt = count_peaks(args.true_peaks)
    Np = count_peaks(args.pooled_peaks)
    N1 = count_peaks(args.self1_peaks)
    N2 = count_peaks(args.self2_peaks)

    # Rescue ratio: max(Np,Nt) / min(Np,Nt)
    rescue_num = Np if Np > Nt else Nt
    rescue_den = Nt if Np > Nt else Np
    rescue_ratio = compute_ratio(rescue_num, rescue_den)

    # Self-consistency ratio: max(N1,N2) / min(N1,N2)
    self_num = N1 if N1 > N2 else N2
    self_den = N2 if N1 > N2 else N1
    self_ratio = compute_ratio(self_num, self_den)

    # Reproducibility status: pass if both ratios < 2 and finite
    status = "fail"
    try:
        r = float(rescue_ratio)
        s = float(self_ratio)
        if r < 2 and s < 2:
            status = "pass"
    except ValueError:
        pass

    # Copy true-replicate IDR thresholded as final validated peak
    shutil.copyfile(args.true_peaks, args.output_peak)

    # Write 15-column summary TSV
    with open(args.output_tsv, "w") as fh:
        fh.write("\t".join([
            "experiment",
            "assay",
            "peak_mode",
            "caller",
            "bio_rep_a",
            "bio_rep_b",
            "true_peaks_Nt",
            "pooled_peaks_Np",
            "self1_peaks_N1",
            "self2_peaks_N2",
            "rescue_ratio",
            "self_consistency_ratio",
            "reproducibility_status",
            "final_method",
            "final_output",
        ]) + "\n")
        fh.write("\t".join([
            args.experiment,
            args.assay,
            args.peak_mode,
            args.caller,
            args.bio_rep_a,
            args.bio_rep_b,
            str(Nt),
            str(Np),
            str(N1),
            str(N2),
            rescue_ratio,
            self_ratio,
            status,
            args.final_method,
            args.final_output,
        ]) + "\n")

    print(
        f"ATAC IDR summary: exp={args.experiment} "
        f"Nt={Nt} Np={Np} N1={N1} N2={N2} "
        f"rescue={rescue_ratio} self={self_ratio} status={status}"
    )


if __name__ == "__main__":
    main()
