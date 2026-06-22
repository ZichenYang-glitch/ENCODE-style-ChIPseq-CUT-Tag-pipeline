"""Stage 55 ATAC IDR dry-run target tests.

These checks exercise Snakemake target expansion only. They do not execute
MACS3, IDR, or pseudoreplicate splitting.
"""

import os
import shutil
import subprocess
import sys
import tempfile

from _tool_resolver import resolve_tool


SNAKEFILE = "workflow/Snakefile"
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")


BASE_CONFIG = """\
samples: "{samples}"
outdir: "{outdir}"
use_control: false
threads: 1
stage4b: true
stage5: false
multiqc: false
"""


ATAC_2BIOREP = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
    "\tbowtie2_index\texperiment\tbiological_replicate\n"
    "A1\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\tatac_exp\t1\n"
    "A2\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\tatac_exp\t2\n"
)


MIXED_ATAC_CHIPSEQ = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
    "\tbowtie2_index\texperiment\tbiological_replicate\n"
    "A1\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\tatac_exp\t1\n"
    "A2\t{r1}\t{r2}\tPE\tatac\tT\tnarrow\ths\tidx\tatac_exp\t2\n"
    "C1\t{r1}\t{r2}\tPE\tchipseq\tTF\tnarrow\ths\tidx\tchip_exp\t1\n"
    "C2\t{r1}\t{r2}\tPE\tchipseq\tTF\tnarrow\ths\tidx\tchip_exp\t2\n"
)


def _write_case(tmp, config_extra, samples_tsv):
    r1 = os.path.join(tmp, "R1.fq")
    r2 = os.path.join(tmp, "R2.fq")
    open(r1, "w").close()
    open(r2, "w").close()

    samples_path = os.path.join(tmp, "samples.tsv")
    with open(samples_path, "w") as fh:
        fh.write(samples_tsv.format(r1=r1, r2=r2))

    outdir = os.path.join(tmp, "results")
    config_path = os.path.join(tmp, "config.yaml")
    with open(config_path, "w") as fh:
        fh.write(BASE_CONFIG.format(samples=samples_path, outdir=outdir))
        fh.write(config_extra)

    return config_path, outdir


def _dryrun(config_path):
    result = subprocess.run(
        [
            SNAKEMAKE,
            "-s",
            SNAKEFILE,
            "--configfile",
            config_path,
            "-n",
            "-p",
        ],
        capture_output=True,
        text=True,
    )
    return result.returncode, result.stdout + result.stderr


def _run_case(config_extra, samples_tsv):
    with tempfile.TemporaryDirectory(prefix="stage55_dryrun_") as tmp:
        config_path, outdir = _write_case(tmp, config_extra, samples_tsv)
        rc, output = _dryrun(config_path)
        output = output.replace(outdir, "{outdir}")
        return rc, output


def main():
    passed = 0
    total = 0

    def check(name, fn):
        nonlocal passed, total
        total += 1
        try:
            if fn():
                print("PASS: %s" % name)
                passed += 1
        except Exception as exc:
            print("FAIL: %s\n   Unexpected error: %s" % (name, exc))

    enabled = (
        "reproducibility:\n"
        "  enabled: true\n"
        "  idr:\n"
        "    atac_narrow: true\n"
    )
    disabled_flag = (
        "reproducibility:\n"
        "  enabled: true\n"
        "  idr:\n"
        "    atac_narrow: false\n"
    )
    disabled_block = (
        "reproducibility:\n"
        "  enabled: false\n"
        "  idr:\n"
        "    atac_narrow: true\n"
    )

    def d1():
        rc, output = _run_case(enabled, ATAC_2BIOREP)
        if rc != 0:
            print("   dry-run failed:\n%s" % output[-1000:])
            return False
        required = [
            "atac_macs3_idr_biorep",
            "atac_idr_true_replicates",
            "atac_split_pseudoreps",
            "atac_macs3_idr_pseudorep",
            "atac_idr_self_pseudoreps",
            "atac_idr_pooled_pseudoreps",
            "atac_idr_summary",
            "06_reproducibility/idr/",
        ]
        missing = [text for text in required if text not in output]
        if missing:
            print("   missing: %s" % missing)
            return False
        return True
    check("D1: ATAC narrow IDR targets appear", d1)

    def d2():
        rc, output = _run_case(enabled, MIXED_ATAC_CHIPSEQ)
        if rc != 0:
            print("   dry-run failed:\n%s" % output[-1000:])
            return False
        if "atac_exp/06_reproducibility/idr/" not in output:
            print("   ATAC IDR path not found")
            return False
        forbidden = [
            "chip_exp/06_idr/",
            "chip_exp/06_reproducibility/idr/",
        ]
        found = [text for text in forbidden if text in output]
        if found:
            print("   forbidden ChIP-seq IDR paths found: %s" % found)
            return False
        return True
    check("D2: mixed run creates ATAC IDR only", d2)

    def d3():
        rc, output = _run_case(disabled_flag, ATAC_2BIOREP)
        if rc != 0:
            print("   dry-run failed:\n%s" % output[-1000:])
            return False
        return "atac_idr_" not in output and "06_reproducibility/idr/" not in output
    check("D3: atac_narrow=false disables ATAC IDR targets", d3)

    def d4():
        rc, output = _run_case(disabled_block, ATAC_2BIOREP)
        if rc != 0:
            print("   dry-run failed:\n%s" % output[-1000:])
            return False
        return "atac_idr_" not in output and "06_reproducibility/idr/" not in output
    check("D4: reproducibility.enabled=false gates ATAC IDR", d4)

    def d5():
        rc, output = _run_case(enabled, ATAC_2BIOREP)
        if rc != 0:
            print("   dry-run failed:\n%s" % output[-1000:])
            return False
        required = [
            "06_reproducibility/final/atac_exp.atac.macs3.narrow.replicate_validated.idr.narrowPeak",
            "06_reproducibility/final/reproducibility_summary.tsv",
        ]
        missing = [text for text in required if text not in output]
        if missing:
            print("   missing final paths: %s" % missing)
            return False
        return True
    check("D5: ATAC IDR final paths", d5)

    print("\n%d/%d tests passed" % (passed, total))

    for path in ("__pycache__", "scripts/__pycache__", "test/__pycache__"):
        if os.path.isdir(path):
            shutil.rmtree(path)

    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
