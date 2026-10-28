# SCI-7 source log

Query date: 2026-10-28. All queries are read-only ENCODE portal API calls.

| Query | Result | What it established |
|---|---|---|
| `https://www.encodeproject.org/search/?type=Experiment&assay_title=TF%20ChIP-seq&status=released&format=json&limit=3` | 3 experiments, total 5002 | TF ChIP-seq candidates |
| `.../experiments/ENCSR054ZMK/?format=json` | 1 bio replicate, 2 PE FASTQ, controls ENCSR142JFU / ENCSR441NSM | TF treatment + matched control |
| `.../experiments/ENCSR142JFU/?format=json` | 2 PE FASTQ, Control ChIP-seq | TF control files |
| `.../search/?type=Experiment&assay_title=Histone%20ChIP-seq&status=released&format=json&limit=3` | 3 experiments, total 3832 | histone candidates |
| `.../experiments/ENCSR817LUF/?format=json` | 1 bio replicate, 5 SE FASTQ, control ENCSR603VTC | histone treatment + matched control |
| `.../experiments/ENCSR603VTC/?format=json` | 2 SE FASTQ, Control ChIP-seq | histone control files |
| `.../search/?type=Experiment&assay_title=ATAC-seq&status=released&format=json&limit=3` | 3 experiments, total 559 | ATAC candidates |
| `.../experiments/ENCSR844TIU/?format=json` | mouse Treg, 2 bio replicates, 4 PE FASTQ, no control | ATAC PE candidate |
| `.../search/?type=Experiment&assay_title=MNase-seq&status=released&format=json&limit=3` | 2 experiments | ENCODE MNase candidates |
| `.../experiments/ENCSR000CXP/?format=json` | 9 bio replicates, SOLiD csfasta/csqual/BAM hg19, no standard FASTQ | unsuitable for this pipeline; gap |
| `.../search/?type=Experiment&searchTerm=CUT%26Tag&status=released&format=json&limit=3` | returned DNase-seq, not CUT&Tag | ENCODE CUT&Tag not exposed by this API path; gap |
| `.../search/?type=Experiment&assay_title=CUT%26Tag...` | HTTP 404 | no such assay_title facet |

Raw responses saved under `/tmp/sci7-work/metadata/` (`search_*.json`, `exp_*.json`).
No FASTQ/BAM was downloaded. File sizes/md5 are the official portal metadata values.
