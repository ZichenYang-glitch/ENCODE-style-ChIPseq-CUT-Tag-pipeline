# ENCODE-style ChIP-seq & CUT&Tag Pipeline

This is an automated, high-performance bioinformatic pipeline designed for both **ChIP-seq** and **CUT&Tag** data analysis. It follows ENCODE best practices and features adaptive logic for different sequencing protocols.

<br>

## 📁 Repository Structure

- `chipseq.yml`: Conda environment recipe with all required tools (FastQC, Bowtie2, MACS3, etc.).
- `scripts/chipseq.sh`: The core master script with integrated QC, alignment, and peak calling.
- `.gitignore`: Configured to exclude large genomic files (FASTQ, BAM, BW).

<br>

## 🚀 Getting Started

<br>

### 1. Installation

```bash
# Clone the repository
git clone https://github.com/YangZiChen-glitch/ENCODE-style-ChIP-CUT-Tag-pipeline.git
cd ENCODE-style-ChIP-CUT-Tag-pipeline

# Create and activate environment
conda env create -f chipseq.yml
conda activate chipseq

# Optional: Install to your bin for global access
cp scripts/chipseq.sh $CONDA_PREFIX/bin/chipseq
chmod +x $CONDA_PREFIX/bin/chipseq
```

<br>

### 2. Usage Examples
#### Paired-End (PE) - CUT&Tag Mode
For CUT&Tag, the script automatically applies high-resolution centering (--nomodel --shift -100 --extsize 200).

```bash
chipseq -s MySample -x /path/to/bt2_index -g mm --pe -r1 R1.fq.gz -r2 R2.fq.gz --tech cuttag
```
#### Single-End (SE) - Standard ChIP-seq

```bash
chipseq -s MySample -x /path/to/bt2_index -g hg38 --se -r1 sample.fq.gz --tech chip
```

<br>

## 🛠 Script Parameters
| Option | Description | Default |
| :--- | :--- | :--- |
| `-s \| --sample` | Sample ID (used for naming outputs) | **Required** |
| `-x \| --index` | Path to Bowtie2 index prefix | **Required** |
| `-g \| --genome` | Genome build (e.g., `mm`, `hg`) | `mm` |
| `--pe / --se` | Paired-end or Single-end mode | **Required** |
| `--tech` | Protocol type: `chip` or `cuttag` | `chip` |
| `--peakMode` | Peak type: `narrow` or `broad` | `narrow` |
| `-p \| --threads` | Number of CPU threads | `8` |
| `-o \| --outdir` | Output directory | `./chip_out` |

<br>

## 🔄 Pipeline Workflow

1. **Quality Control**: `FastQC` generates raw data metrics.

2. **Trimming**: `Trim Galore` removes adapters and low-quality bases.

3. **Alignment**: `Bowtie2` maps reads with Read Group (RG) info for Picard compatibility.

4. **Filtering**: `Samtools` filters for MAPQ >= 30 and removes unmapped/secondary reads.

5. **Deduplication**: `Picard` or `Samtools` markdup removes PCR duplicates (auto-enabled for narrow peaks).

6. **Peak Calling**: `MACS3` calls peaks based on the chosen --tech and --peakMode.

7. **Visualization**: `deepTools` generates CPM-normalized BigWig files for IGV.

<br>

## 📊 Expected Outputs
```text
chip_out/
├── 00_raw/        # Trimmed Fastq files
├── 01_qc/         # FastQC, flagstat, and duplication metrics
├── 02_align/      # Sorted and filtered BAM files
├── 03_bigwig/     # CPM-normalized BigWig tracks
└── 04_peaks/      # MACS3 peak files (.narrowPeak or .broadPeak)
```

<br>

## 📝 License
This project is open-source under the MIT License.
