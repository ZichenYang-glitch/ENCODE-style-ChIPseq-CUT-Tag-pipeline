import subprocess
import os
import sys

# Path to the validator script
VALIDATOR = "scripts/validate_samples.py"

def run_test(name, config_yaml, samples_tsv, expect_fail=True, expected_error=""):
    """
    Helper to run a single validation test case.
    Writes temp config and samples files, runs the validator, and checks results.
    """
    with open("test_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_samples.tsv", "w") as f:
        f.write(samples_tsv)
    
    # Run the validator script as a subprocess
    result = subprocess.run(
        [sys.executable, VALIDATOR, "--config", "test_config.yaml"],
        capture_output=True, text=True
    )
    
    passed = False
    if expect_fail:
        if result.returncode == 0:
            print(f"❌ FAIL: {name}\n   Expected failure, but it passed.")
        elif expected_error not in result.stderr and expected_error not in result.stdout:
            print(f"❌ FAIL: {name}\n   Expected error '{expected_error}' not found.\n   Actual output: {result.stderr.strip()}")
        else:
            print(f"✅ PASS: {name}")
            passed = True
    else:
        if result.returncode != 0:
            print(f"❌ FAIL: {name}\n   Expected success, but it failed.\n   Output: {result.stderr.strip()}")
        else:
            print(f"✅ PASS: {name}")
            passed = True
            
    return passed

def main():
    print("Starting Stress Tests for Stage 2 (validate_samples.py)...\n")
    tests = 0
    passed = 0

    # Common valid snippets
    valid_samples = "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n" \
                    "S1\tR1.fq\tR2.fq\tPE\tchipseq\tT\tnarrow\ths\tidx\n"

    valid_config = "samples: 'test_samples.tsv'\nuse_control: false\n"

    def t(name, cfg, smp, expect_fail=True, err=""):
        nonlocal tests, passed
        tests += 1
        if run_test(name, cfg, smp, expect_fail, err):
            passed += 1

    # --- 1. Baseline ---
    t("Baseline valid", valid_config, valid_samples, expect_fail=False)

    # --- 2. Configuration Validation ---
    t("Config: threads not int", valid_config + "threads: abc\n", valid_samples, err="must be an integer")
    t("Config: threads zero", valid_config + "threads: 0\n", valid_samples, err="must be positive")
    t("Config: mapq negative", valid_config + "mapq: -5\n", valid_samples, err="must be non-negative")
    t("Config: invalid remove_dup", valid_config + "remove_dup: whatever\n", valid_samples, err="auto, yes, or no")
    
    # --- 3. Genome Resources ---
    t("Config: missing genome_resources egs", 
      valid_config + "genome_resources:\n  hs:\n    chrom_sizes: ''\n", 
      valid_samples, err="effective_genome_size is required")
      
    t("Config: invalid genome_resources egs (string instead of shortcut/numeric)", 
      valid_config + "genome_resources:\n  hs:\n    effective_genome_size: 'hg38'\n", 
      valid_samples, err="'hs', 'mm', or a positive integer")

    # --- 4. Sample Sheet Structure ---
    bad_samples_col = "sample\tfastq_2\tlayout\nS1\tR2\tSE\n"
    t("Sample: missing fastq_1 column", valid_config, bad_samples_col, err="missing required column: 'fastq_1'")
    
    # --- 5. Sample Validation ---
    bad_samples_id = "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n" \
                     "S1 Bad\tR1\tR2\tPE\tchipseq\tT\tnarrow\ths\tidx\n"
    t("Sample: invalid sample ID (space)", valid_config, bad_samples_id, err="invalid characters")

    bad_samples_pe = "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n" \
                     "S1\tR1\t\tPE\tchipseq\tT\tnarrow\ths\tidx\n"
    t("Sample: PE missing fastq_2", valid_config, bad_samples_pe, err="PE layout requires 'fastq_2'")

    bad_samples_assay = "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n" \
                        "S1\tR1\t\tSE\trnaseq\tT\tnarrow\ths\tidx\n"
    t("Sample: invalid assay", valid_config, bad_samples_assay, err="assay must be chipseq, cuttag, or atac")

    # --- 6. Control Logic Validation ---
    cfg_ctrl = "samples: 'test_samples.tsv'\nuse_control: true\n"
    
    # Missing control sample row
    smps_miss_ctrl = "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\trole\tcontrol_sample\n" \
                     "S1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tCtrl1\n"
    t("Control: referencing missing sample", cfg_ctrl, smps_miss_ctrl, err="not found in sample sheet")

    # Control sample with role=treatment
    smps_bad_role = "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\trole\tcontrol_sample\n" \
                    "S1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tCtrl1\n" \
                    "Ctrl1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\t\n"
    t("Control: reference has treatment role", cfg_ctrl, smps_bad_role, err="expected role=control")

    # Mutually exclusive control_sample and control_bam
    smps_mutex = "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\trole\tcontrol_sample\tcontrol_bam\n" \
                 "S1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tCtrl1\t/tmp/a.bam\n" \
                 "Ctrl1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\tcontrol\t\t\n"
    t("Control: mutually exclusive controls", cfg_ctrl, smps_mutex, err="cannot set both control_sample and control_bam")
    
    # Missing control bam file
    smps_bad_bam = "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\trole\tcontrol_bam\n" \
                   "S1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\t/does/not/exist.bam\n"
    t("Control: missing control_bam file", cfg_ctrl, smps_bad_bam, err="control_bam file not found")

    print(f"\nSummary: {passed}/{tests} tests passed.")
    
    # Cleanup temp files
    for f in ["test_config.yaml", "test_samples.tsv"]:
        if os.path.exists(f): 
            os.remove(f)
            
    if passed < tests:
        sys.exit(1)

if __name__ == "__main__":
    main()
