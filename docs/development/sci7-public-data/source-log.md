# SCI-7 source log

Query date: 2026-10-28. The initial queries below were read-only ENCODE portal
and NCBI GEO API calls. Raw responses are stored in this directory's
[`metadata/`](metadata/) subdirectory; no FASTQ/BAM was downloaded.

## Discovery metadata

| Complete query | Saved response | What it established |
| --- | --- | --- |
| `https://www.encodeproject.org/search/?type=Experiment&assay_title=TF%20ChIP-seq&status=released&format=json&limit=3` | [`search_TF_ChIP-seq.json`](metadata/search_TF_ChIP-seq.json) | TF ChIP-seq search results (3 returned, total 5002) |
| `https://www.encodeproject.org/experiments/ENCSR054ZMK/?format=json` | [`exp_ENCSR054ZMK.json`](metadata/exp_ENCSR054ZMK.json) | One biological replicate, two PE FASTQ files, possible controls including ENCSR142JFU |
| `https://www.encodeproject.org/experiments/ENCSR142JFU/?format=json` | [`exp_ENCSR142JFU.json`](metadata/exp_ENCSR142JFU.json) | Two PE FASTQ files for the matched Control ChIP-seq experiment |
| `https://www.encodeproject.org/search/?type=Experiment&assay_title=Histone%20ChIP-seq&status=released&format=json&limit=3` | [`search_Histone_ChIP-seq.json`](metadata/search_Histone_ChIP-seq.json) | Histone ChIP-seq search results (3 returned, total 3832) |
| `https://www.encodeproject.org/experiments/ENCSR817LUF/?format=json` | [`exp_ENCSR817LUF.json`](metadata/exp_ENCSR817LUF.json) | One biological replicate and five SE FASTQ files belonging to the same library/technical replicate; matched control ENCSR603VTC |
| `https://www.encodeproject.org/experiments/ENCSR603VTC/?format=json` | [`exp_ENCSR603VTC.json`](metadata/exp_ENCSR603VTC.json) | Two SE FASTQ files for the matched control |
| `https://www.encodeproject.org/search/?type=Experiment&assay_title=ATAC-seq&status=released&format=json&limit=3` | [`search_ATAC-seq.json`](metadata/search_ATAC-seq.json) | ATAC-seq search results (3 returned, total 559) |
| `https://www.encodeproject.org/experiments/ENCSR844TIU/?format=json` | [`exp_ENCSR844TIU.json`](metadata/exp_ENCSR844TIU.json) | Mouse Treg, two biological replicates and four 26 bp PE FASTQ files |
| `https://www.encodeproject.org/search/?type=Experiment&assay_title=MNase-seq&status=released&format=json&limit=3` | [`search_MNase-seq.json`](metadata/search_MNase-seq.json) | ENCODE MNase-seq candidate search (two experiments) |
| `https://www.encodeproject.org/experiments/ENCSR000CXP/?format=json` | [`exp_ENCSR000CXP.json`](metadata/exp_ENCSR000CXP.json) | SOLiD csfasta/csqual/BAM hg19 rather than compatible standard FASTQ; rejected for this workflow |

For every selected ENCODE row, file accession, official `file_size`, official
`md5sum`, read length and `paired_with` were rechecked against the saved
experiment response. The manifest's accession-only `paired_with` field is a
normalization of ENCODE's `/files/<ACCESSION>/` reference.

## Assay search limitations

The ENCODE query
`https://www.encodeproject.org/search/?type=Experiment&searchTerm=CUT%26Tag&status=released&format=json&limit=3`
returned DNase-seq rather than a CUT&Tag experiment, and the attempted
`assay_title=CUT%26Tag` facet path returned HTTP 404. These results establish
only that those API paths did not yield a compatible CUT&Tag record; they do
**not** establish that ENCODE/GEO/SRA lacks CUT&Tag data.

The initial NCBI GEO DataSets esearch responses were:

- [`geo_CUT_Tag.json`](metadata/geo_CUT_Tag.json): query translation
  `CUT[All Fields] AND Tag[All Fields]`, 19,003 hits. This is a broad text
  search, not 19,003 verified CUT&Tag experiments.
- [`geo_MNase-seq.json`](metadata/geo_MNase-seq.json): query translation
  `MNase-seq[All Fields]`, 3,885 hits. These hits were not file-level
  candidates by themselves.

Follow-up file-level CUT&Tag/MNase selection and the official queries used are
documented in [`cuttag-mnase-selection.md`](cuttag-mnase-selection.md).
Confirmation queries are documented in
[`confirmation-review.md`](confirmation-review.md), and reference sources in
[`reference-identity.md`](reference-identity.md).
