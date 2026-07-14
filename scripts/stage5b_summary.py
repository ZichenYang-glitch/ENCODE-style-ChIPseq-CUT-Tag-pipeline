#!/usr/bin/env python3
"""Stage 5b reproducibility QC summary and final peak set assembly.

Reads IDR thresholded peak files, computes rescue ratio and self-consistency
ratio, handles zero denominators, writes reproducibility_summary.tsv, and
copies conservative/optimal peak sets.

Usage:
    python3 scripts/stage5b_summary.py \
        --true-peaks true_replicates/idr.thresholded.narrowPeak \
        --pooled-peaks pooled_pseudoreps/idr.thresholded.narrowPeak \
        --self1-peaks self_pseudoreps/biorepA.idr.thresholded.narrowPeak \
        --self2-peaks self_pseudoreps/biorepB.idr.thresholded.narrowPeak \
        --experiment exp1 \
        --bio-rep-a 1 --bio-rep-b 2 \
        --output-tsv reproducibility_summary.tsv \
        --output-cons conservative.narrowPeak \
        --output-opt optimal.narrowPeak
"""

import argparse
import shutil


def count_peaks(path):
    """Count lines in a peak file, excluding header lines."""
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
        description="Stage 5b reproducibility summary"
    )
    parser.add_argument("--true-peaks", required=True)
    parser.add_argument("--pooled-peaks", required=True)
    parser.add_argument("--self1-peaks", required=True)
    parser.add_argument("--self2-peaks", required=True)
    parser.add_argument("--experiment", required=True)
    parser.add_argument("--bio-rep-a", required=True)
    parser.add_argument("--bio-rep-b", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-cons", required=True)
    parser.add_argument("--output-opt", required=True)
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

    # Reproducibility status: pass if both ratios < 2 and both finite
    status = "fail"
    try:
        r = float(rescue_ratio)
        s = float(self_ratio)
        if r < 2 and s < 2:
            status = "pass"
    except ValueError:
        pass  # NA or inf → fail

    # Copy final peak sets
    shutil.copyfile(args.true_peaks, args.output_cons)
    shutil.copyfile(args.pooled_peaks, args.output_opt)

    # Write summary TSV
    with open(args.output_tsv, "w") as fh:
        fh.write("\t".join([
            "experiment",
            "true_peaks(Nt)",
            "pooled_peaks(Np)",
            "self1_peaks(N1)",
            "self2_peaks(N2)",
            "rescue_ratio",
            "self_consistency_ratio",
            "reproducibility_status",
        ]) + "\n")
        fh.write("\t".join([
            args.experiment,
            str(Nt),
            str(Np),
            str(N1),
            str(N2),
            rescue_ratio,
            self_ratio,
            status,
        ]) + "\n")

    print(f"Summary: exp={args.experiment} Nt={Nt} Np={Np} "
          f"N1={N1} N2={N2} "
          f"rescue={rescue_ratio} self={self_ratio} status={status}")


if __name__ == "__main__":
    main()
