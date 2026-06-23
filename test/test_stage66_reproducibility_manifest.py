"""Stage 66 reproducibility manifest tests.

Verifies mode-specific output_type rows, gating, and that
no chipseq narrow IDR final entries appear under 06_reproducibility/.
"""

import os
import subprocess
import sys
import tempfile

MANIFEST_SCRIPT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts", "make_manifest.py",
)


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
        except Exception as e:
            print("FAIL: %s\n   Unexpected error: %s" % (name, e))

    HEADER = (
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
        "\tbowtie2_index\texperiment\tbiological_replicate\n"
    )

    def ms(*rows):
        return HEADER + "\n".join(
            "{id}\t{r1}\t{r2}\t{layout}\t{assay}\tT\t{pm}\ths\tidx\t{exp}\t{br}".format(
                id=id, r1="{r1}", r2="{r2}", layout=layout,
                assay=assay, pm=pm, exp=exp, br=br,
            )
            for id, assay, pm, exp, br, layout in rows
        )

    CT_PE_N2 = ms(
        ("T1", "cuttag", "narrow", "exp_ct", "1", "PE"),
        ("T2", "cuttag", "narrow", "exp_ct", "2", "PE"),
    )
    CS_B2 = ms(
        ("CB1", "chipseq", "broad", "exp_csb", "1", "PE"),
        ("CB2", "chipseq", "broad", "exp_csb", "2", "PE"),
    )
    AT_N2 = ms(
        ("A1", "atac", "narrow", "exp_at", "1", "PE"),
        ("A2", "atac", "narrow", "exp_at", "2", "PE"),
    )
    CS_N2 = ms(
        ("C1", "chipseq", "narrow", "exp_cs", "1", "PE"),
        ("C2", "chipseq", "narrow", "exp_cs", "2", "PE"),
    )
    CT_PE_B2 = ms(
        ("TB1", "cuttag", "broad", "exp_ctb", "1", "PE"),
        ("TB2", "cuttag", "broad", "exp_ctb", "2", "PE"),
    )

    BASE = (
        "samples: \"{samples}\"\nuse_control: false\nthreads: 1\n"
        "stage4b: true\n"
    )

    def make_config(samples_tsv, extra=""):
        cfg = BASE + extra
        return cfg.format(samples=samples_tsv)

    def run_manifest(name, config_yaml, samples_tsv_str):
        with tempfile.TemporaryDirectory(prefix="stage66_") as tmp:
            r1 = os.path.join(tmp, "R1.fq")
            r2 = os.path.join(tmp, "R2.fq")
            open(r1, "w").close()
            open(r2, "w").close()

            samples_fmt = samples_tsv_str.format(r1=r1, r2=r2)
            samples_path = os.path.join(tmp, "samples.tsv")
            with open(samples_path, "w") as f:
                f.write(samples_fmt)

            config_fmt = config_yaml.format(samples=samples_path)
            config_path = os.path.join(tmp, "config.yaml")
            with open(config_path, "w") as f:
                f.write(config_fmt)

            out_path = os.path.join(tmp, "manifest.tsv")
            result = subprocess.run(
                [sys.executable, MANIFEST_SCRIPT,
                 "--config", config_path, "--output", out_path],
                capture_output=True, text=True, timeout=60,
            )
            if result.returncode != 0:
                print("FAIL: %s\n   make_manifest.py failed: %s"
                      % (name, result.stderr.strip()[:500]))
                return None

            lines = []
            with open(out_path) as f:
                for line in f:
                    line = line.rstrip("\n")
                    if line == "":
                        continue
                    if line.startswith("#"):
                        continue
                    lines.append(line)
            return lines

    def output_types(lines):
        header = lines[0].split("\t")
        ot_idx = header.index("output_type") if "output_type" in header else -1
        if ot_idx < 0:
            return set()
        return {row.split("\t")[ot_idx] for row in lines[1:]}

    # T1: consensus enabled → consensus output types appear
    def t1():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: true\n  consensus:\n    enabled: true\n"
        ))
        lines = run_manifest("T1", cfg, CS_B2)
        if lines is None:
            return False
        ots = output_types(lines)
        return "chipseq_macs3_broad_consensus_peak" in ots
    check("T1: consensus rows appear", t1)

    # T2: cuttag_narrow IDR → cuttag IDR final row
    def t2():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: true\n  idr:\n    cuttag_narrow: true\n"
        ))
        lines = run_manifest("T2", cfg, CT_PE_N2)
        if lines is None:
            return False
        ots = output_types(lines)
        return "cuttag_macs3_narrow_idr_final_peak" in ots
    check("T2: cuttag IDR final row", t2)

    # T3: atac_narrow IDR → atac IDR final row
    def t3():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: true\n  idr:\n    atac_narrow: true\n"
        ))
        lines = run_manifest("T3", cfg, AT_N2)
        if lines is None:
            return False
        ots = output_types(lines)
        return "atac_macs3_narrow_idr_final_peak" in ots
    check("T3: atac IDR final row", t3)

    # T4: broad experimental true → broad IDR final row
    def t4():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: true\n  idr:\n"
            "    chipseq_broad_experimental: true\n"
        ))
        lines = run_manifest("T4", cfg, CS_B2)
        if lines is None:
            return False
        ots = output_types(lines)
        return "chipseq_macs3_broad_idr_final_peak" in ots
    check("T4: broad IDR final row", t4)

    # T5: broad experimental false → no broad IDR final row
    def t5():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: true\n  idr:\n"
            "    chipseq_broad_experimental: false\n"
        ))
        lines = run_manifest("T5", cfg, CS_B2)
        if lines is None:
            return False
        ots = output_types(lines)
        if "chipseq_macs3_broad_idr_final_peak" in ots:
            print("   broad IDR final row present but flag is false")
            return False
        # But consensus should still be there
        if "chipseq_macs3_broad_consensus_final_peak" not in ots:
            print("   consensus final row missing when IDR is off")
            return False
        return True
    check("T5: broad IDR off → no IDR final, consensus final present", t5)

    # T6: no chipseq_macs3_narrow_idr_final_peak ever
    def t6():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: true\n  consensus:\n    enabled: true\n"
        ))
        lines = run_manifest("T6", cfg, CS_N2)
        if lines is None:
            return False
        ots = output_types(lines)
        if "chipseq_macs3_narrow_idr_final_peak" in ots:
            print("   chipseq narrow IDR final peak should never appear in 06_reproducibility/")
            return False
        return True
    check("T6: no chipseq_macs3_narrow_idr_final_peak", t6)

    # T7: no chipseq narrow replicate_validated.idr under 06_reproducibility/
    def t7():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: true\n  consensus:\n    enabled: true\n"
        ))
        lines = run_manifest("T7", cfg, CS_N2)
        if lines is None:
            return False
        for line in lines[1:]:
            cols = line.split("\t")
            path = cols[7] if len(cols) > 7 else ""
            ot = cols[6] if len(cols) > 6 else ""
            if ("replicate_validated.idr" in path
                and "06_reproducibility" in path
                and "chipseq" in ot and "narrow" in ot):
                print("   Found chipseq narrow replicate_validated.idr under"
                      " 06_reproducibility/: %s" % path)
                return False
            if ("chipseq_macs3_narrow_idr_final_peak" in ot):
                print("   Found chipseq_macs3_narrow_idr_final_peak output_type")
                return False
        return True
    check("T7: no chipseq narrow 06_reproducibility/ idr final path", t7)

    # T8: SEACR enabled → SEACR consensus rows
    def t8():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: true\n  consensus:\n    enabled: true\n"
            "cuttag:\n  seacr:\n    enabled: true\n    mode: stringent\n"
        ))
        lines = run_manifest("T8", cfg, CT_PE_N2)
        if lines is None:
            return False
        ots = output_types(lines)
        has_seacr_consensus = "cuttag_seacr_consensus_peak" in ots
        has_seacr_final = "cuttag_seacr_consensus_final_peak" in ots
        return has_seacr_consensus and has_seacr_final
    check("T8: SEACR consensus + final rows", t8)

    # T9: reproducibility disabled → zero 06_reproducibility rows
    def t9():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: false\n"
        ))
        lines = run_manifest("T9", cfg, CT_PE_N2)
        if lines is None:
            return False
        for line in lines[1:]:
            if "06_reproducibility" in line:
                print("   Found 06_reproducibility/ path when repro disabled: %s"
                      % line.split("\t")[7])
                return False
        return True
    check("T9: repro disabled → zero 06_reproducibility rows", t9)

    # T10: stage4b false → zero 06_reproducibility rows
    def t10():
        cfg = BASE.replace("stage4b: true", "stage4b: false") + (
            "reproducibility:\n  enabled: true\n  consensus:\n    enabled: true\n"
        )
        cfg = cfg.format(samples="{samples}")
        lines = run_manifest("T10", cfg, CT_PE_N2)
        if lines is None:
            return False
        for line in lines[1:]:
            if "06_reproducibility" in line:
                print("   Found 06_reproducibility/ path when stage4b false: %s"
                      % line.split("\t")[7])
                return False
        return True
    check("T10: stage4b false → zero 06_reproducibility rows", t10)

    # T11: consensus final row present for modes where IDR not enabled
    def t11():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: true\n  consensus:\n    enabled: true\n"
        ))
        lines = run_manifest("T11", cfg, CT_PE_B2)
        if lines is None:
            return False
        ots = output_types(lines)
        return "cuttag_macs3_broad_consensus_final_peak" in ots
    check("T11: consensus final for cuttag broad (no IDR)", t11)

    # T12: disabled mode omitted (not not_applicable)
    def t12():
        cfg = make_config("{samples}", (
            "reproducibility:\n  enabled: true\n"
        ))
        lines = run_manifest("T12", cfg, CT_PE_N2)
        if lines is None:
            return False
        for line in lines[1:]:
            cols = line.split("\t")
            status = cols[8] if len(cols) > 8 else ""
            ot = cols[6] if len(cols) > 6 else ""
            if "idr_final" in ot and status == "not_applicable":
                print("   IDR final row has not_applicable status (should be omitted)")
                return False
        return True
    check("T12: disabled modes omitted, not not_applicable", t12)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
