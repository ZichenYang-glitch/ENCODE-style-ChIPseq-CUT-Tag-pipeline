
################################################################################################
This directory contains data files produced by the 

           GENCODE project

which is headed by Fergal Martin at the EMBL-EBI, UK.

For questions, please contact gencode-help@ebi.ac.uk
and check the website http://www.gencodegenes.org

################################################################################################


##################################
GENCODE Release Files Description
##################################

(X is the release version, eg. 21 for human, M4 for mouse):

#################
#Annotation files
#################

1. gencode.vX.annotation.{gtf,gff3}.gz:
  Main file, gene annotation on reference chromosomes in GTF and GFF3 file formats.
  These are the main GENCODE gene annotation files. They contain annotation (genes, 
  transcripts, exons, start_codon, stop_codon, UTRs, CDS) on the reference chromosomes,
  which are chr1-22, X, Y, M in human and chr1-19, X, Y, M in mouse.

2. gencode.vX.chr_patch_hapl_scaff.annotation.{gtf,gff3}.gz:
  Gene annotation on reference-chromosomes/patches/scaffolds/haplotypes in GTF and GFF3 file formats.

3. gencode.vX.primary_assembly.annotation.{gtf,gff3}.gz:
  Gene annotation on reference chromosomes and scaffolds in GTF and GFF3 file formats.

4. gencode.vX.basic.annotation.{gtf,gff3}.gz:
  Basic gene annotation on reference chromosomes in GTF and GFF3 file formats.
  This is a subset of the corresponding comprehensive annotation including only 
  those transcripts tagged as 'basic' in every gene.

5. gencode.vX.chr_patch_hapl_scaff.basic.annotation.{gtf,gff3}.gz:
  Basic gene annotation on reference-chromosomes/patches/scaffolds/haplotypes in 
  GTF and GFF3 file formats.

6. gencode.vX.primary_assembly.basic.annotation.{gtf,gff3}.gz:
  Basic gene annotation on reference chromosomes and scaffolds in
  GTF and GFF3 file formats.

7. gencode.vX.long_noncoding_RNAs.{gtf,gff3}.gz: 
  Long non-coding RNAs on reference chromosomes in GTF and GFF3 file formats.
  These files are a sub-set of the main annotation files on the reference chromosomes. 
  They contain only the lncRNA genes, which are those with any of these biotypes: 
  "processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", 
  "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", 
  "macro_lncRNA", "lncRNA".

8. gencode.vX.polyAs.{gtf,gff3}.gz:
  PolyA features annotated by Havana on reference chromosomes in GTF and GFF3 file formats.
  These files contain polyA signals, polyA sites and pseudo polyAs manually 
  annotated by HAVANA. They include only the reference chromosomes.
  The value of the 'gene_id', 'transcript_id', 'gene_name' and 'transcript name' fields
  corresponds to a random identifier. The polyA features are not directly associated to 
  any gene or transcript when annotated by Havana. 

9. gencode.vX.2wayconspseudos.{gtf,gff3}.gz:
  (Retrotransposed) pseudogenes predicted by the Yale & UCSC pipelines, but not by Havana on 
  reference chromosomes in GTF and GFF3 file formats.

10. gencode.vX.tRNAs.{gtf,gff3}.gz:
  tRNA structures predicted by tRNA-Scan on reference chromosomes in GTF and GFF3 file formats.


###############
#Sequence files
###############

11. gencode.vX.transcripts.fa.gz:
  All transcript sequences on reference chromosomes in Fasta format.
  The sequence headers contain the following fields:
    transcript-id|
    gene-id|
    Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
    Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
    transcript-name|
    gene-name|
    sequence-length|
    transcript biotype

12. gencode.vX.pc_transcripts.fa.gz:
  Protein-coding transcript sequences on reference chromosomes in Fasta format.
  The sequence headers contain the following fields:
    transcript-id|
    gene-id|
    Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
    Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
    transcript-name|
    gene-name|
    sequence-length|
    5'-UTR (3'-UTR if reverse strand) location in the transcript|
    CDS location in the transcript|
    3'-UTR (5'-UTR if reverse strand) location in the transcript

13. gencode.vX.pc_translations.fa.gz:
  Translations of protein-coding transcripts on reference chromosomes Fasta file.
  The sequence headers contain the following fields:
    transcript-id|
    gene-id|
    Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
    Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
    transcript-name|gene-name|sequence-length

14. gencode.vX.lncRNA_transcripts.fa.gz:
  Long non-coding RNA transcript sequences on reference chromosomes Fasta file.
  The sequence headers contain the following fields:
    transcript-id|
    gene-id|
    Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
    Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
    transcript-name|
    gene-name|
    sequence-length
  
  Long non-coding RNA transcripts are those with any of these biotypes: 
  "processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", 
  "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", 
  "macro_lncRNA", "lncRNA".

15. A.genome.fa.gz (where A is the current assembly name eg. GRCh38 for human, GRCm38 for mouse):
  Genome sequence fasta file (sequence region names are the same as in the GTF/GFF3 files).
  It includes reference chromosomes, scaffolds, assembly patches and haplotypes. 

16. A.primary_assembly.genome.fa.gz (where A is the current assembly name eg. GRCh38 for human, GRCm38 for mouse):
  Primary assembly genome sequence fasta file (sequence region names are the same as in the GTF/GFF3 files).
  It includes reference chromosomes and scaffolds only.


###############
#Metadata files
###############

17. gencode.vX.metadata.Annotation_remark.gz:
  Remarks made during the manual annotation of the transcript.
  1 - transcript id
  2 - annotation remark

18. gencode.vX.metadata.EntrezGene.gz:
  Entrez Gene id associated to the transcript.
  1 - transcript id
  2 - Entrez Gene id

19. gencode.vX.metadata.Exon_supporting_feature.gz:
  Piece of evidence used in the annotation of an exon (usually peptides, mRNAs, ESTs).
  1 - transcript id
  2 - external id of the feature supporting the exon annotation
  3 - external source of the supporting feature
  4 - exon id
  5 - exon coordinates

20. gencode.vX.metadata.Gene_source.gz:
  Source of the gene annotation (Ensembl, Havana, Ensembl-Havana merged model or 
  imported in the case of small RNA and mitochondrial genes).
  1 - gene id
  2 - gene source

21. gencode.vX.metadata.HGNC.gz:
  HGNC approved gene symbol.
  1 - transcript id
  2 - HGNC gene symbol
  3 - HGNC unique id

22. gencode.vX.metadata.PDB.gz:
  PDB entry associated to the transcript.
  1 - transcript id
  2 - PDB id

23. gencode.vX.metadata.PolyA_feature.gz:
  Manually annotated polyA feature overlapping the transcript 3'-end.
  1 - transcript id
  2 - transcript-based start coordinate of the polyA feature
  3 - transcript-based end coordinate of the polyA feature
  4 - polyA feature chromosome
  5 - polyA feature start coordinate
  6 - polyA feature end coordinate
  7 - polyA feature strand
  8 - polyA feature type ("polyA_site", "polyA_signal", "pseudo_polyA")

24. gencode.vX.metadata.Pubmed_id.gz:
  PubMed id of a publication associated to the transcript.
  1 - transcript id
  2 - PubMed id

25. gencode.vX.metadata.RefSeq.gz:
  RefSeq RNA and/or protein associated to the transcript.
  1 - transcript id
  2 - RefSeq RNA id
  3 - RefSeq protein id (optional)

26. gencode.vX.metadata.Selenocysteine.gz:
  Amino acid position of a selenocysteine residue in the transcript.
  1 - transcript id
  2 - selenocysteine position

27. gencode.vX.metadata.SwissProt.gz:
  UniProtKB/SwissProt entry associated to the transcript.
  1 - transcript id
  2 - UniProtKB/SwissProt accession number
  3 - UniProtKB/SwissProt accession number

28. gencode.vX.metadata.Transcript_source.gz:
  Source of the transcript annotation (Ensembl, Havana, etc).
  1 - transcript id
  2 - transcript source

29. gencode.vX.metadata.Transcript_supporting_feature.gz:
  Piece of evidence used in the annotation of the transcript.
  1 - transcript id
  2 - external id of the feature supporting the transcript annotation
  3 - external source of the supporting feature
  
30. gencode.vX.metadata.TrEMBL.gz:
  UniProtKB/TrEMBL entry associated to the transcript.
  1 - transcript id
  2 - UniProtKB/TrEMBL accession number
  3 - UniProtKB/TrEMBL accession number




######################################
General format of the annotation files
######################################

We supply genome-wide features on three different confidence levels.
Level 1 + 2 should be used for high-quality local analysis.
1 + 2 + 3 should be used for genome-wide analysis.

* Level 1: validated 

Pseudogene loci, that were predicted by the analysis-pipelines from YALE, UCSC 
as well as by HAVANA manual annotation from WTSI. 
Other transcripts, that were verified experimentally by RT-PCR and sequencing
through the GENCODE experimental pipeline.

* Level 2: manual annotation 

HAVANA manual annotation from WTSI (and Ensembl annotation where it is identical to Havana).
The following regions are considered "fully annotated" although they will still be updated:
chromosomes 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, 
ENCODE pilot regions.

* Level 3: automated annotation 

ENSEMBL loci where they are different from the HAVANA annotation or where no annotation 
can be found. 

This data is supplied in GTF and GFF3 format as defined here:
 http://www.gencodegenes.org/data_format.html
with the following tags added to the attributes column where appropriate:

* level [1,2,3]: validation status as described.
* tag "3_nested_supported_extension": 3' end extended based on RNA-seq data.
* tag "3_standard_supported_extension": 3' end extended based on RNA-seq data.
* tag "454_RNA_Seq_supported": annotated based on RNA-seq data.
* tag "5_nested_supported_extension": 5' end extended based on RNA-seq data.
* tag "5_standard_supported_extension": 5' end extended based on RNA-seq data.
* tag "alternative_3_UTR": shares an identical CDS but has alternative 3' UTR with respect to a 
  reference variant.
* tag "alternative_5_UTR": shares an identical CDS but has alternative 5' UTR with respect to a 
  reference variant.
* tag "appris_principal": transcript expected to code for the main functional isoform based on a 
  range of protein features (APPRIS pipeline, Nucleic Acids Res. 2013 Jan;41(Database issue):D110-7). 
  (this tag is not found after Gencode 21)
* tag "appris_candidate": where there is no single 'appris_principal' variant the main functional 
  isoform will be translated from one of the 'appris_candidate' genes.  (this tag is not found after 
  Gencode 21)
* tag "appris_candidate_ccds": the "appris_candidate" transcript that has an unique CCDS.  
  (this tag is not found after Gencode 21)
* tag "appris_candidate_longest_ccds": the "appris_candidate" transcripts where there are several 
  CCDS, in this case APPRIS labels the longest CCDS.  (this tag is not found after Gencode 21)
* tag "appris_candidate_longest_seq": where there is no "appris_candidate_ccds" or 
  "appris_candidate_longest_ccds" variant, the longest protein of the "appris_candidate" variants is 
  selected as the primary variant.  (this tag is not found after Gencode 21)
* tag "appris_candidate_highest_score": where there is no 'appris_principal' variant, the candidate 
  with highest APPRIS score is selected as the primary variant. (this tag is not found after Gencode 
  20)
* tag "appris_candidate_longest": where there is no 'appris_principal' variant, the longest of the 
  'appris_candidate' variants is selected as the primary variant. (this tag is not found after 
  Gencode 20)
* tag "appris_principal_1": (This flag corresponds to the older flag "appris_principal") where the 
  transcript expected to code for the main functional isoform based solely on the core modules in the 
  APPRIS database. The APPRIS core modules map protein structural and functional information and 
  cross-species conservation to the annotated variants.
* tag "appris_principal_2": (This flag corresponds to the older flag "appris_candidate_ccds") Where 
  the APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human 
  protein coding genes), the database chooses two or more of the CDS variants as "candidates" to be 
  the principal variant.
  If one (but no more than one) of these candidates has a distinct CCDS identifier it is selected as 
  the principal variant for that gene. A CCDS identifier shows that there is consensus between RefSeq 
  and GENCODE/Ensembl for that variant, guaranteeing that the variant has cDNA support.
* tag "appris_principal_3": Where the APPRIS core modules are unable to choose a clear principal 
  variant and there more than one of the variants have distinct CCDS identifiers, APPRIS selects the 
  variant with lowest CCDS identifier as the principal variant. The lower the CCDS identifier, the 
  earlier it was annotated.
  Consensus CDS annotated earlier are likely to have more cDNA evidence. Consecutive CCDS identifiers 
  are not included in this flag, since they will have been annotated in the same release of CCDS. 
  These are distinguished with the next flag. 
* tag "appris_principal_4": (This flag corresponds to the Ensembl 78 flag 
  "appris_candidate_longest_ccds") Where the APPRIS core modules are unable to choose a clear 
  principal CDS and there is more than one variant with a distinct (but consecutive) CCDS 
  identifiers, APPRIS selects the longest CCDS isoform as the principal variant.
* tag "appris_principal_5": (This flag corresponds to the Ensembl 78 flag 
  "appris_candidate_longest_seq") Where the APPRIS core modules are unable to choose a clear 
  principal variant and none of the candidate variants are annotated by CCDS, APPRIS selects the 
  longest of the candidate isoforms as the principal variant.
* tag "appris_alternative_1": Candidate transcript(s) models that are conserved in at least three 
  tested non-primate species.
* tag "appris_alternative_2": Candidate transcript(s) models that appear to be conserved in fewer 
  than three tested non-primate species.
* tag "artifactual_duplication": annotated on an artifactual duplicate region of the genome assembly. 
* tag "basic": identifies a subset of representative transcripts for each gene; prioritises 
  full-length protein coding transcripts over partial or non-protein coding transcripts within the 
  same gene, and intends to highlight those transcripts that will be useful to the majority of users. 
* tag "bicistronic": transcript contains two confidently annotated CDSs. Support may come from eg 
  proteomic data, cross-species conservation or published experimental work.
* tag "CAGE_supported_TSS": transcript 5' end overlaps ENCODE or Fantom CAGE cluster.
* tag "CCDS": member of the consensus CDS gene set, confirming coding regions between ENSEMBL, 
  UCSC, NCBI and HAVANA.
* tag "cds_end_NF": the coding region end could not be confirmed.
* tag "cds_start_NF": the coding region start could not be confirmed.
* tag "dotter_confirmed": transcript QC checked using dotplot to identify features eg splice 
  junctions, end of homology.
* tag "downstream_ATG": downstream ATG assessed as less likely to initiate the translation of the 
  functional protein due to eg experimental evidence, poor cross-species conservation, weak Kozak 
  context, interference with signal peptides or targeting signals.
* tag "Ensembl_canonical": most representative transcript of the gene. This will be the MANE_Select 
  transcript if there is one, or a transcript chosen by an Ensembl algorithm otherwise.
* tag "exp_conf": transcript was tested and confirmed experimentally.
* tag "fragmented_locus": locus consists of non-overlapping transcript fragments either because of 
  genome assembly issues (i.e., gaps or mis-assemblies), or because supporting transcripts (e.g., 
  from another species) cannot be completely mapped, or because the supporting transcripts are 
  non-overlapping end pairs (i.e., 5' and 3' ESTs from a single cDNA).
* tag "GENCODE_Primary": belongs to a minimal set that contains MANE Select, MANE Plus Clinical and
  Ensembl Canonical transcripts and transcripts containing any conserved exons and common alternative
  splicing events (including exons skips) that are absent from the MANE and Ensembl Canonical
  transcripts for protein-coding genes. Other biotypes will have the GENCODE_Primary flag added to
  the Ensembl Canonical transcript and, for lncRNA genes only, this will be the transcript with the
  longest genomic span. 
* tag "inferred_exon_combination": transcript model contains all possible in-frame exons supported 
  by homology, experimental evidence or conservation, but the exon combination is not directly 
  supported by a single piece of evidence and may not be biological. Used for large genes with 
  repetitive exons (e.g. titin (TTN)) to represent all the exons individual transcript variants can 
  pool from.
* tag "inferred_transcript_model": transcript model is not supported by a single piece of 
  transcript evidence. May be supported by multiple fragments of transcript evidence or by combining 
  different evidence sources e.g. protein homology, RNA-seq data, published experimental data.
* tag "low_sequence_quality": transcript supported by transcript evidence that, while mapping 
  best-in-genome, shows regions of poor sequence quality.
* tag "MANE_Select": the transcript belongs to the MANE Select data set. The Matched Annotation 
  from NCBI and EMBL-EBI project (MANE) is a collaboration between Ensembl-GENCODE and RefSeq 
  to select a default transcript per human protein coding locus that is representative of biology, 
  well-supported, expressed and conserved. This transcript set matches GRCh38 and is 100% identical
  between RefSeq and Ensembl-GENCODE for 5' UTR, CDS, splicing and 3' UTR.
* tag "MANE_Plus_Clinical": the transcript belongs to the MANE Plus Clinical data set. Within the 
  MANE project, these are additional transcripts per locus necessary to support clinical variant 
  reporting, for example transcripts containing known pathogenic or likely pathogenic clinical 
  variants not reportable using the MANE Select data set. This transcript set matches GRCh38 
  and is 100% identical between RefSeq and Ensembl-GENCODE for 5' UTR, CDS, splicing and 3' UTR.
* tag "mRNA_end_NF": the mRNA end could not be confirmed.
* tag "mRNA_start_NF": the mRNA start could not be confirmed.
* tag "NAGNAG_splice_site": in-frame type of variation where, at the acceptor site, some variants 
  splice after the first AG and others after the second AG.
* tag "ncRNA_host": the locus is a host for small non-coding RNAs. 
* tag "nested_454_RNA_Seq_supported": annotated based on RNA-seq data.
* tag "NMD_exception": the transcript looks like it is subject to NMD but publications, experiments 
  or conservation support the translation of the CDS.
* tag "NMD_likely_if_extended": codon if the transcript were longer but cannot currently be 
  annotated as NMD as does not fulfil all criteria - most commonly lack of an intron downstream of 
  the stop codon.
* tag "non_ATG_start": the CDS has a non-ATG start and its validity is supported by publication or 
  conservation.
* tag "non_canonical_conserved": the transcript has a non-canonical splice site conserved in other 
  species.
* tag "non_canonical_genome_sequence_error": the transcript has a non-canonical splice site 
  explained by a genomic sequencing error.
* tag "non_canonical_other": the transcript has a non-canonical splice site explained by other 
  reasons.
* tag "non_canonical_polymorphism": the transcript has a non-canonical splice site explained by a 
  SNP.
* tag "non_canonical_TEC": the transcript has a non-canonical splice site that needs experimental 
  confirmation.
* tag "non_canonical_U12": the transcript has a non-canonical splice site explained by a U12 intron 
  (i.e. AT-AC splice site).
* tag "non_submitted_evidence": a splice variant for which supporting evidence has not been 
  submitted to databases, i.e. the model is based on literature or collaborator evidence.
* tag "not_best_in_genome_evidence": a transcript is supported by evidence from same species 
  paralogous loci.
* tag "not_organism_supported": evidence from other species was used to build model.
* tag "orphan": protein-coding locus with no paralogues or orthologs.
* tag "overlapping locus": exon(s) of the locus overlap exon(s) of a readthrough transcript or a 
  transcript belonging to another locus.
* tag "overlapping_uORF": a low confidence upstream ATG existing in other coding variant would lead 
  to NMD in this trancript, that uses the high confidence canonical downstream ATG.
* tag "PAR": annotation in the pseudo-autosomal region, which is duplicated between X & Y.
* tag "pseudo_consens": member of the pseudogene set predicted by YALE, UCSC and HAVANA.
* tag "readthrough_gene": protein-coding gene that has a readthrough transcript.
* tag "readthrough_transcript": a transcript that overlaps two or more independent loci but is 
  considered to belong to a 3rd, separate locus.
* tag "reference_genome_error": locus overlaps a sequence error or an assembly error in the 
  reference genome that affects its annotation (e.g., 1 or 2bp insertion/deletion, substitution 
  causing premature stop codon). The main effect is that affected transcripts that would have had a 
  CDS are currently annotated without one.
* tag "retained_intron_CDS": internal intron of CDS portion of transcript is retained.
* tag "retained_intron_final": final intron of CDS portion of transcript is retained.
* tag "retained_intron_first": first intron of CDS portion of transcript is retained.
* tag "retrogene": protein-coding locus created via retrotransposition.
* tag "RNA_Seq_supported_only": transcript supported by RNAseq data and not supported by mRNA or 
  EST evidence.
* tag "RNA_Seq_supported_partial": transcript annotated based on mixture of RNA-seq data and 
  EST/mRNA/protein evidence.
* tag "RP_supported_TIS": transcript that contains a CDS that has a translation initiation site 
  supported by Ribosomal Profiling data.
* tag "seleno": contains a selenocysteine.
* tag "semi_processed": a processed pseudogene with one or more introns still present. These are 
  likely formed through the retrotransposition of a retained intron transcript.
* tag "sequence_error": transcript contains ≥ 1 non-canonical splice junction that is associated 
  with a known or novel genome sequence error
* tag "stop_codon_readthrough": transcript whose coding sequence contains an internal stop codon 
  that does not cause the translation termination.
* tag "TAGENE": a transcript created or extended using assembled RNA-seq long reads.
* tag "upstream_ATG": an upstream ATG exists when a downstream ATG is better supported.
* tag "upstream_uORF": a low confidence upstream ATG existing in other coding variant would lead to 
  NMD in this trancript, that uses the high confidence canonical downstream ATG.


Please note: if start codons are split between two exons, two start-codon features will be listed.
Please note: pre-release 4, "cds_start_NF" was listed as "cds start not found", etc.
Please note: pre-release 6, "seleno" tags included the selenocystein position as the amino acid 
number within the protein, these are now given as genomic coordinates as separate GTF 
features.
Please note: the stable ids with the ENSGR0000XXXXXX, ENSTR0000XXXXXX format (until release 24) or
the ENSG00000XXXXXX.X_PAR_Y, ENST00000XXXXXX.X_PAR_Y format (from release 25 onwards) are genes and 
transcripts in the pseudoautosomal regions (PAR regions) of human chromosome Y. These genes/transcripts 
are tagged with "PAR" and have a different stable_id than their counterpart in chromosome X to avoid 
redundancy. For example: ENST00000431238.7 in chrX and ENST00000431238.7_PAR_Y in chrY.





#################################################################################
Release 50 (June 2026)
#################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2025.
This release corresponds to Ensembl version 116.


#################################################################################################
Release 49 (September 2025)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to February 2025.
This release corresponds to Ensembl version 115.



#################################################################################################
Release 48 (May 2025)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2024.
This release corresponds to Ensembl version 114.

The following files, which were restricted to reference chromosomes, now contain annotation on
primary assembly and alternate regions:
-gencode.v48.long_noncoding_RNAs.gtf, gencode.v48.long_noncoding_RNAs.gff3
-gencode.v48.transcripts.fa
-gencode.v48.pc_transcripts.fa
-gencode.v48.pc_translations.fa
-gencode.v48.lncRNA_transcripts.fa

The gencode.v48.pc_translations.fa file now includes all translation sequences. Previously,
sequences with internal stop codons were normally excluded.

The first release of the GENCODE promoter window file (gencode.v48.promoter_windows.gff3)
is now available (see https://www.gencodegenes.org/pages/cls_lncrnas/promoter_windows/).


#################################################################################################
Release 47 (Oct 2024)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to April 2024.
This release corresponds to Ensembl version 113.

The following tags or attributes have been added:
-GENCODE_Primary



#################################################################################################
Release 46 (May 2024)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to September 2023.
This release corresponds to Ensembl version 112.



#################################################################################################
Release 45 (January 2024)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to March 2023.
This release corresponds to Ensembl version 111.



#################################################################################################
Release 44 (July 2023)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to December 2022.
This release corresponds to Ensembl version 110.

The assembly version has been updated to GRCh38.p14.
Annotation on the chrY PAR has separate Ensembl stable ids.



#################################################################################################
Release 43 (February 2023)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2022
This release corresponds to Ensembl version 109.



#################################################################################################
Release 42 (October 2022)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to April 2022
This release corresponds to Ensembl version 108.

The following biotype has been added:
-protein_coding_CDS_not_defined: Transcript that belongs to a protein_coding gene 
  and doesn't contain an ORF. Replaces the processed_transcript transcript biotype
  in protein_coding genes.

Annotation files in gtf and gff3 format having the basic set of transcripts in the 
primary assembly are now included in the release.



#################################################################################################
Release 41 (July 2022)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to January 2022
This release corresponds to Ensembl version 107.

The following biotypes have been added:
-protein_coding_LoF: Not translated in the reference genome owing to a SNP/DIP 
  but in other individuals/haplotypes/strains the transcript is translated. 
  This biotype replaces the polymorphic_pseudogene transcript biotype.
-artifact: Annotated on artifactual regions of the genome assembly.

The following tags or attributes have been added:
-artifactual_duplication / artif_dupl
-readthrough_gene



#################################################################################################
Release 40 (April 2022)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2021
This release corresponds to Ensembl version 106.



#################################################################################################
Release 39 (December 2021)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to May 2021
This release corresponds to Ensembl version 105.

Clone-based gene names have been retired from the release files. From this release onwards, 
genes without a name in HGNC, EntrezGene, RFAM or miRBase will get their gene id as default name.   



#################################################################################################
Release 38 (May 2021)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to December 2020
This release corresponds to Ensembl version 104.

New tag added: Ensembl_canonical.



#################################################################################################
Release 37 (February 2021)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2020
This release corresponds to Ensembl version 103.

New tag added: MANE_Plus_Clinical.



#################################################################################################
Release 36 (November 2020)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to May 2020
This release corresponds to Ensembl version 102.



#################################################################################################
Release 35 (August 2020)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to March 2020
This release corresponds to Ensembl version 101.



#################################################################################################
Release 34 (April 2020)
##################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to November 2019.
This release corresponds to Ensembl version 100.



#################################################################################################
Release 33 (January 2020)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2019.
This release corresponds to Ensembl version 99.

NOTE: the following genes had their biotype erroneously changed from protein_coding to polymorphic_pseudogene. The correct biotypes will be restored in our next release.

ENSG00000158887.18  MPZ
ENSG00000014641.20  MDH1
ENSG00000211456.13  SACM1L
ENSG00000109339.24  MAPK10
ENSG00000112715.23  VEGFA
ENSG00000123505.18  AMD1
ENSG00000082556.13  OPRK1
ENSG00000134575.13  ACP2
ENSG00000111716.14  LDHB
ENSG00000111424.12  VDR
ENSG00000184992.13  BRI3BP
ENSG00000171885.17  AQP4
ENSG00000125510.18  OPRL1



#################################################################################################
Release 32 (September 2019)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to May 2019.
This release corresponds to Ensembl version 98.

New tag added: stop_codon_readthrough.



#################################################################################################
Release 31 (July 2019)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to February 2019.
This release corresponds to Ensembl version 97.

The following biotypes have been replaced by the "lncRNA" biotype: 
-3prime_overlapping_ncRNA
-antisense
-bidirectional_promoter_lncRNA
-lincRNA
-macro_lncRNA
-non_coding
-processed_transcript
-sense_intronic
-sense_overlapping

The following tags or attributes have been added:
-TAGENE: Transcript created or extended using assembled RNA-seq long reads.
-HGNC_id: Unique stable id provided by the HGNC for each gene with an approved symbol.



#################################################################################################
Release 30 (April 2019)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to November 2018.
This release corresponds to Ensembl version 96.



#################################################################################################
Release 29 (October 2018)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to May 2018.
This release corresponds to Ensembl version 94.



#################################################################################################
Release 28 (April 2018)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to November 2017.
This release corresponds to Ensembl version 92.



#################################################################################################
Release 27 (August 2017)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to January 2017.
This release corresponds to Ensembl version 90.

The antisense biotype is now called "antisense_RNA".



#################################################################################################
Release 26 (March 2017)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to October 2016.
This release corresponds to Ensembl version 88.

The "gene_status" and "transcript_status" attributes have been removed from the GENCODE GTF and 
GFF3 files.

New tags added:
-3_nested_supported_extension
-3_standard_supported_extension
-454_RNA_Seq_supported
-5_nested_supported_extension
-5_standard_supported_extension
-nested_454_RNA_Seq_supported



#################################################################################################
Release 25 (July 2016)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to March 2016.
This release corresponds to Ensembl version 85.

Genes and transcripts on the chrY PAR regions now have "_PAR_Y" appended to their ids. Until 
release 24 these ids had the "ENSGR00..." and "ENSTR00..." formats.

The "UTR" features in the GFF3 files have been replaced with "five_prime_UTR" and "three_prime_UTR".



#################################################################################################
Release 24 (December 2015)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2015.
This release corresponds to Ensembl version 83.



#################################################################################################
Release 23 (July 2015)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to March 2015.
This release corresponds to Ensembl version 81.

New files added in this release:

a)Basic annotation on the reference chromosomes and on all sequence regions:
-gencode.v23.basic.annotation.{gtf,gff3}.gz
-gencode.v23.basic.annotation.gff3.gz
-gencode.v23.chr_patch_hapl_scaff.basic.annotation.{gtf,gff3}.gz
-gencode.v23.chr_patch_hapl_scaff.basic.annotation.gff3.gz

b)Nucleotide sequences of all annotated transcripts on the reference chromosomes:
-gencode.v23.transcripts.fa.gz



#################################################################################################
Release 22 (March 2015)
#################################################################################################
This is a merge between a full new Ensembl gene build and updates from HAVANA up to October 2014.
This release corresponds to Ensembl version 79.

Three new APPRIS tags have been added in this release:

-appris_principal_1: (This flag corresponds to the older flag "appris_principal") where the 
transcript expected to code for the main functional isoform based solely on the core modules in the 
APPRIS database. The APPRIS core modules map protein structural and functional information and 
cross-species conservation to the annotated variants.

-appris_principal_2: (This flag corresponds to the older flag "appris_candidate_ccds") Where the 
APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human 
protein coding genes), the database chooses two or more of the CDS variants as "candidates" to be 
the principal variant.
If one (but no more than one) of these candidates has a distinct CCDS identifier it is selected as 
the principal variant for that gene. A CCDS identifier shows that there is consensus between RefSeq 
and GENCODE/Ensembl for that variant, guaranteeing that the variant has cDNA support.

-appris_principal_3: Where the APPRIS core modules are unable to choose a clear principal variant 
and there more than one of the variants have distinct CCDS identifiers, APPRIS selects the variant 
with lowest CCDS identifier as the principal variant. The lower the CCDS identifier, the earlier it 
was annotated.
Consensus CDS annotated earlier are likely to have more cDNA evidence. Consecutive CCDS identifiers 
are not included in this flag, since they will have been annotated in the same release of CCDS. 
These are distinguished with the next flag. 

-appris_principal_4: (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_ccds") 
Where the APPRIS core modules are unable to choose a clear principal CDS and there is more than one 
variant with a distinct (but consecutive) CCDS identifiers, APPRIS selects the longest CCDS isoform 
as the principal variant.

-appris_principal_5: (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_seq") 
Where the APPRIS core modules are unable to choose a clear principal variant and none of the 
candidate variants are annotated by CCDS, APPRIS selects the longest of the candidate isoforms as 
the principal variant.

-appris_alternative_1: Candidate transcript(s) models that are conserved in at least three tested 
non-primate species.

-appris_alternative_2: Candidate transcript(s) models that appear to be conserved in fewer than 
three tested non-primate species.

In contrast, five APPRIS tags present in the previous release have been dropped:
appris_principal 
appris_candidate_ccds
appris_candidate_longest_ccds
appris_candidate_longest_seq
appris_candidate


Also, in this release the transcript attributes were removed from the gene lines in all GTF and 
GFF3 annotation files. ie. the transcript_id transcript_type, transcript_status, transcript_name, 
havana_transcript, etc.


NOTE: The GFF3 files listed below were replaced on 7 May 2015 to fix non-unique IDs of some UTR 
features as well as to extend the CDS until the last nucleotide of the stop codon in stop codons 
split between two exons.
 - gencode.v22.chr_patch_hapl_scaff.annotation.gff3.gz
 - gencode.v22.annotation.gff3.gz

NOTE: GTF and GFF3 files containing gene annotation on the primary assembly only (main chromosomes 
and unplaced/unlocalized scaffolds) were added on 7 May 2015:
 - gencode.v22.primary_assembly.annotation.{gtf,gff3}.gz
 - gencode.v22.primary_assembly.annotation.gff3.gz
These files are meant to be used mainly for NGS analyses in conjunction with the equivalent primary 
assembly genome sequence provided:
 - GRCh38.primary_assembly.genome.fa.gz



#################################################################################################
Release 21 (October 2014)
#################################################################################################
This is a merge between a full new Ensembl gene build and updates from HAVANA up to June 2014.
This release corresponds to Ensembl version 77.

The transcript support levels imported from UCSC have been introduced in the main annotation files. 
Transcripts are scored according to how well mRNA and EST alignments match over its full length:
- 1 (all splice junctions of the transcript are supported by at least one non-suspect mRNA)
- 2 (the best supporting mRNA is flagged as suspect or the support is from multiple ESTs)
- 3 (the only support is from a single EST)
- 4 (the best supporting EST is flagged as suspect)
- 5 (no single transcript supports the model structure)
- NA (the transcript was not analyzed) 

Three new APPRIS tags have been added in this release:
-appris_candidate_ccds: the "appris_candidate" transcript that has an unique CCDS.
-appris_candidate_longest_ccds: the "appris_candidate" transcripts where there are several CCDS, in 
this case APPRIS labels the longest CCDS.
-appris_candidate_longest_seq: where there is no "appris_candidate_ccds" or 
"appris_candidate_longest_ccds" variant, the longest protein of the "appris_candidate" variants is 
selected as the primary variant.

In contrast, two APPRIS tags present in the previous release have been dropped:
-appris_candidate_highest_score
-appris_candidate_longest


NOTE: The files listed below were replaced on 12-11-2014 to fix the haplotype chromosome names, 
which should correspond to the GRC accessions as in previous GENCODE releases.
 - gencode.v21.chr_patch_hapl_scaff.annotation.{gtf,gff3}.gz
 - gencode.v21.chr_patch_hapl_scaff.annotation.gff3.gz
 - gencode.v21.metadata.PolyA_feature.gz
 - GRCh38.genome.fa.gz



#################################################################################################
Release 20 (August 2014)
#################################################################################################
This is a merge between a full new Ensembl gene build and updates from HAVANA up to April 2014.
This release corresponds to Ensembl version 76 and is in the new human assembly GRCh38.

A new annotation file type is added in the Gencode releases which is GFF3.
Note:
In GFF3 the start and stop codon are included in the CDS.
In GTF the start codon is included in the CDS but the stop codon is included in the UTR.  

1 extra appris tag is added in this release in the transcript lines:

* tag "appris_candidate_highest_score": where there is no 'appris_principal' variant, the candidate 
with highest APPRIS score is selected as the primary variant.

1 extra attribute is added in the annotation files in all protein coding transcripts: the protein_id

For example:

chr1	HAVANA	transcript	69091	70008	.	+	.	gene_id 
"ENSG00000186092.4"; transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status 
"KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; 
transcript_name "OR4F5-001"; level 2; protein_id "ENSP00000334393.3"; tag "basic"; tag 
"appris_principal"; tag "CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.1"; 
havana_transcript "OTTHUMT00000003223.1";


1 extra metadata file is added:

gencode.v20.metadata.EntrezGene: EntrezGene gene id associated to the transcript

NOTE: The files listed below were updated on 27-08-2014 to correct the following minor issues:
- some CCDS tags had been wrongly assigned to transcripts with the same start and end coordinates 
as a CCDS model but with different amino acid sequence
- absence of semicolon after "havana_transcript" attributes in the GTF files.

gencode.v20.annotation.{gtf,gff3}.gz
gencode.v20.chr_patch_hapl_scaff.annotation.{gtf,gff3}.gz
gencode.v20.long_noncoding_RNAs.{gtf,gff3}.gz
gencode.v20.annotation.gff3.gz
gencode.v20.chr_patch_hapl_scaff.annotation.gff3.gz


#################################################################################################
Release 19 (December 2013)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to July 2013.
This release corresponds to Ensembl version 74.

3 extra tags are added in this release in the transcript lines:

* tag "appris_principal": transcript expected to code for the main functional isoform based on a 
range of protein features (APPRIS pipeline).
* tag "appris_candidate": where there is no single 'appris_principal' variant the main functional 
isoform will be translated from one of the 'appris_candidate' genes.
* tag "appris_candidate_longest": where there is no 'appris_principal' variant, the longest of the 
'appris_candidate' variants is selected as the primary variant.


NOTE: The tag "basic" was added to the gencode.v19.annotation.gff3 file on 2 Feb 2016.

NOTE: The protein_id attribute was added to the gencode.v19.annotation.gtf file on 1 March 2016.


#################################################################################################
Release 18 (September 2013)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to April 2013.
This release corresponds to Ensembl version 73.

NOTE: File GRCh37.p12.genome.fa.gz with the whole genome fasta sequence was replaced in the server 
on the 4th October 2013 because it was faulty.

There is an extra field in the exon lines in the GTFs: the exon_id, which is placed 
after the exon_number and before the level.


#################################################################################################
Release 17 (June 2013)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to February 2013.
This release corresponds to Ensembl version 72.

There is an extra file in this release:
GRCh37.p11.genome.fa.gz
which has the genome sequence in fasta format.
In the label there are 2 names of the sequence region:
the official GRCh37.p11 assembly names then tab and then the Ensembl name.
*The names of the the official GRCh37.p11 assembly used in the GRCh37.p11.genome.fa.gz file 
are the same names used the GTF files.

There are 2 extra transcript tags in the GTF:

RP_supported_TIS
inferred_exon_combination

Please see above in the tags description for more details.


#################################################################################################
Release 16 (April 2013)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to November 2012.
This release corresponds to Ensembl version 71.

(2 files were re-uploaded in the server on the 28 May 2013; they had in the scaffolds name an extra 
 '.1' extension, these are :
1. gencode.v16.chr_patch_hapl_scaff.annotation.{gtf,gff3}.gz 
2. gencode.v16.patch_hapl_scaff.annotation.{gtf,gff3}.gz )

(3 files were re-uploaded in the server on the 30 April 2013; they had as mitochondrial chromosome 
'chrMT' instead of 'chrM'; now they have 'chrM'  these are:
1. gencode.v16.annotation.gtf
2. gencode.v16.tRNAs.gtf
3. gencode.v16.chr_patch_hapl_scaff.annotation.gtf)

This release has 2 extra annotation GTFS:
1) a GTF that has not only the annotation in the main chromosomes but also the annotation in 
patches/scaffolds/haplotypes. It's called "gencode.v16.chr_patch_hapl_scaff.annotation.gtf".
2) a GTF that has only the annotation in patches/scaffolds/haplotypes. It's called 
"gencode.v16.patch_hapl_scaff.annotation.gtf".

The coordinates of the patches/scaffolds/haplotypes start from 1 , i.e. they are not genomic 
coordinates.
The coordinates of the annotation (genes,transcripts,exons etc) in patches/scaffolds/haplotypes are 
relative to the patches/scaffolds/haplotypes.
The naming of the patches/scaffolds/haplotypes are from assembly GRCh37.p10 and they are accession 
numbers which are stable. ie. for the patches can be found here:
ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p10/PATCHES/alt_scaffolds/alt_scaffold_placement.txt
More info about GRCh37.p10 can be found here: 
ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p10/


Files which contain metadata associated to transcripts and genes, usually displayed in the UCSC 
browser, are now part of the public release:
-gencode.v16.metadata.Annotation_remark: remarks made during the manual annotation of the transcript
-gencode.v16.metadata.Exon_supporting_feature: piece of evidence used in the annotation of an exon 
(usually peptides, mRNAs, ESTs)
-gencode.v16.metadata.Gene_source: source of the gene annotation (Ensembl, Havana, Ensembl-Havana 
merged model or imported in the case of small RNA and mitochondrial genes)
-gencode.v16.metadata.HGNC: HGNC approved gene symbol
-gencode.v16.metadata.PDB: PDB entry associated to the transcript
-gencode.v16.metadata.PolyA_feature: manually annotated polyA feature overlapping the transcript 
3'-end (polyA_signal, polyA_site, pseudo_polyA) - Fields are 1)transcript_id, 2-3)polyA feature 
coordinates relative to the transcript, 4-7)polyA feature genomic coordinates, 8)type of polyA 
feature 
-gencode.v16.metadata.Pubmed_id: Pubmed id of a publication associated to the transcript
-gencode.v16.metadata.RefSeq: RefSeq RNA and/or protein associated to the transcript
-gencode.v16.metadata.Selenocysteine: amino acid position of a selenocysteine residue in the 
transcript
-gencode.v16.metadata.SwissProt: UniProtKB/SwissProt entry associated to the transcript
-gencode.v16.metadata.Transcript_source: source of the transcript annotation (Ensembl, Havana, 
Ensembl-Havana merged model or imported in the case of small RNA and mitochondrial genes)
-gencode.v16.metadata.Transcript_supporting_feature: piece of evidence used in the annotation of the 
transcript (usually peptides, mRNAs, ESTs)
-gencode.v16.metadata.TrEMBL: UniProtKB/TrEMBL entry associated to the transcript

All the metadata files include annotation on the main chromosomes/patches/scaffolds/haplotypes.

Also, there are 9 extra transcript tags in the GTF files:

NMD_likely_if_extended
RNA_Seq_supported_only
bicistronic
downstream_ATG
low_sequence_quality
retained_intron_CDS
retained_intron_first
retained_intron_final
sequence_error

Please see above in the tag specification for details about these 9 tags.


#################################################################################################
Release 15 (Jan 2013)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2012.
This release corresponds to Ensembl version 70.

There is an extra attribute on the exon lines called exon_number, that indicates the order of the
exon inside the transcript,starting from the 5'end; i.e. the closest to the 5'end will have 
"exon_number 1", the second next to it will have "exon_number 2" etc.

Also, there is an extra tag introduced in the GTF files, the tag 'readthrough_transcript'.This tag 
is found only on the transcript lines and indicates whether the transcript is read-through:
(a transcript that overlaps two or more independent loci but is considered to belong to a 3rd, 
separate locus)


#################################################################################################
Release 14 (November 2012)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to June 2012.
This release corresponds to Ensembl version 69.

There is an extra FASTA file called 'gencode.v14.lncRNA_transcripts.fa' with the long non-coding 
RNA transcripts in FASTA format.
Also, an extra tag is introduced in the GTF files, the tag 'basic'. This tag corresponds to the 
transcripts belonging to the GENCODE basic annotation set in the UCSC Genome Browser.


#################################################################################################
Release 13 (Aug 2012)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to March 2012.
This release corresponds to Ensembl version 68.

There is an extra GTF file called 'gencode13_corrected_statuses.gtf'.
Ensembl 68 database had errors in the statuses and the GTF was manually corrected after the release.
The statuses were corrected using the Havana database which includes only the manual annotations.
Therefore, there are some statuses "NULL" which exist in the Havana db but not in Ensembl as they 
become "KNOWN" or "NOVEL" according to what sources of evidence are found when the merge happens 
and Ensembl db is created.
Ensembl 69 release will have the correct statuses.


#################################################################################################
Release 12 (May 2012)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to December 2011.
This release corresponds to Ensembl version 67.

As in the previous release the "ncrna_host" biotype was replaced with "processed_transcript" in 
21 genes and 80 transcripts in the Gencode GTF although it is still present in Ensembl 67.


#################################################################################################
Release 11 (February 2012)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to October 2011.
This release corresponds to Ensembl version 66.

The locus-level tag "ncRNA_host" has been introduced in this Gencode release. A total of 2041 genes
were automatically tagged as ncRNA_host if they encompassed a small ncRNA. This is expected to be 
substituted by Havana manual annotation in future releases. In addition the "ncrna_host" biotype 
has been replaced by "processed_transcript" in 19 genes and 80 transcripts. The "ncrna_host" biotype
has now been removed from the Havana annotation but this change will not be effective in Ensembl 
until version 68, so please be aware of this discrepancy between the Gencode GTF and the Ensembl 
database. 

The files pc_transcripts.fa and pc_translations.fa were updated with an extended range of transcript
biotypes (see above) on 2012-04-04.

The following GTF files were updated on 2012-06-01 to ammend 19958 gene names and 27058 transcript 
names that had a wrongly appended version number (mostly ".1"). Please note that the original names 
in the Ensembl core database (e66) remain unchanged.
gencode.v11.annotation.level.gtf
gencode.v11.long_noncoding_RNAs.gtf
gencode.v11.annotation.level_1_2.gtf (gencode11_GRCh37.tgz)
gencode.v11.annotation.level_3.gtf (gencode11_GRCh37.tgz)


#################################################################################################
Release 10 (December 2011)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to July 2011.
This release corresponds to Ensembl version 65.

The files pc_transcripts.fa and pc_translations.fa were updated with an extended range of transcript
biotypes (see above) on 2012-04-04.


#################################################################################################
Release 9 (September 2011)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to May 2011.
Ensembl annotation: 8 genes and 27 transcripts were removed from the previous version.
Havana annotation: updated and first pass of chromosome 14 completed, including the ABO gene on a 
GRC patch. In the merge set, the GAGE cluster is imported from Havana only.
This release corresponds to Ensembl version 64.


#################################################################################################
Release 8 (June 2011)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to March 2011.
Poorly supported Ensembl models were removed and Ensembl annotation on the haplotypes, which were 
missing in Gencode 7, was also added. Annotation by Havana of chromosome 12 has been completed 
and that of chromosome 14 is nearly finished.
This release corresponds to Ensembl version 63.


#################################################################################################
Release 7 (March 2011)
#################################################################################################
This is a merge between a full new Ensembl gene build and updates from HAVANA up to December 2010.
The gene and transcript version numbers have been appended to their stable ids in the release files.
This release corresponds to Ensembl version 62.
This is the release used for the ENCODE integration analysis.


#################################################################################################
Release 6 (January 2011)
#################################################################################################

This is a remerge between the Ensembl annotation and updates from HAVANA up to November 2010.
This includes mainly data on chromosomes 8, 9. The level 1 pseudogenes, 2-way pseudogenes and 
tRNA scans are the same as in version 5.
This release corresponds to Ensembl version 61.


#################################################################################################
Release 5 (December 2010)
#################################################################################################

This is a remerge between the Ensembl annotation and updates from HAVANA up to April 2010.
This includes mainly data on chromosomes 6,7,8. The level 1 pseudogenes and the 2-way pseudogenes
were also updated.
This release corresponds to Ensembl version 60.


#################################################################################################
Release 4 (May 2010)
#################################################################################################

This is a remerge between the Ensembl annotation and updates from HAVANA up to February 2010.
This includes mainly data on chromosomes 4 & 5. The level 1 pseudogenes and the 2-way pseudogenes
were not updated from version 3c.
This release corresponds to Ensembl version 58 (with minor fixes).
It was submitted for the June 2010 ENCODE data freeze.

!!Please note that this file was updated with new transcript status fields on 02.07.2010, 17:00 
BST!!


#################################################################################################
Release 3d (February 2010)
#################################################################################################

This is an update of freeze 3(c), in which ENSEMBL removed a number of dubious transcripts.
Tags cds_start_NF, cds_end_NF, mRNA_start_NF, mRNA_end_NF, exp_conf, non_org_supp are introduced.
This release corresponds to Ensembl version 57.
Only the gencode.v3d.annotation.GRCh37.gtf file was updated.


#################################################################################################
Release 3c (October 2009 freeze)
#################################################################################################

This is an update of freeze 3(b), mainly for chromosomes 3 & 4 for which the latest annotation 
was held back and QC'ed again to be used in the RNASeq Genome Annotation Assessment Project.
Only gencode.v3b.annotation.GRCh37.gtf and gencode.v3b.annotation.NCBI36.gtf were updated.
This release corresponds to Ensembl version 56.
This is the release used for the ENCODE integration analysis.


#################################################################################################
Release 3 (July 2009 freeze)
#################################################################################################

New full merge between HAVANA and ENSEMBL.
Native on genome assembly GRCh37, features are projected back to NCBI36 were possible.

1. all loci on both assemblies:
   gencode.v3b.annotation.GRCh37.gtf (updated on 3.10.09)
   gencode.v3b.annotation.NCBI36.gtf (updated on 3.10.09)

2. polyA features:
   gencode.v3.polyAs.GRCh37.gtf
   gencode.v3.polyAs.NCBI36.gtf

3. tRNA predictions:
   gencode.v3.tRNAs.GRCh37.gtf
   gencode.v3.tRNAs.NCBI36.gtf

4. pseudogenes predicted by UCSC and Yale:
   gencode.v3.2wayconspseudos.GRCh37.gtf
   gencode.v3.2wayconspseudos.NCBI36.gtf


################################################################################################
Updated dump (23.06.2009)
#################################################################################################

1. Updated version of the January 2009 freeze file produced for use in the 1000 Genomes projects 
and others.
gencode_data.rel2b.{gtf,gff3}.gz

The field "CCDSID" with a valid CCDS id has been added to tha 9th column.
The field "CCDSOL" has been added indicating that the transcript overlaps a CCDS transcript, but 
was not flagged as such directly.
A few genes have been removed and added respectively. Some formatting issues have been resolved.

2. protein-coding transcripts, their sequences and translations
gencode_data.rel2b.pc_transcripts.fa.gz
gencode_data.rel2b.pc_translations.fa.gz


################################################################################################
Updated dumps (21.04.2009)
#################################################################################################

<directory release_2>

1. protein-coding transcripts, their sequences and translations
gencode_data.rel2.pc_transcripts.{gtf,gff3}.gz
gencode_data.rel2.pc_transcripts_cdnas.fa.gz
gencode_data.rel2.pc_transcripts_translations.fa.gz

2. re-dump (as January freeze) without external annotation 
gencode_data.rel2a.{gtf,gff3}.gz

################################################################################################
For the analysis data freeze of January 2009 there are the following files in this directory:
################################################################################################

<directory release_2>

1. gencode_data.rel2.{gtf,gff3}.gz
  Data file in GTF format, compressed with gzip, containing annotation on three levels.
  Same format as described below, with the addition of one line for every gene and transcript.
  In case you don not want these, do something like 
  awk '{if($3 !~ "gene|transcript"){print $0}}' gencode_data.rel2.gtf > gencode_data.rel2_mod.gtf 

2. gencode_tRNAscans.rel2.{gtf,gff3}.gz
   tRNAscan predictions from the ENSEMBL simpleFeature table (level 3).

3. gencode_polyAs.rel2.{gtf,gff3}.gz
   poly signals from the loutre database (polyA_site, pseudo_polyA) (seperate level).

################################################################################################
For the initial data freeze of October 1st 2008 there are the following files in this directory:
################################################################################################

<directory release_1>

1.gencode_data.rel1.v2.{gtf,gff3}.gz
  Data file in GTF format, compressed with gzip, containing annotation on three levels.
  Here is what the first 11 lines look like (containing first transcript):

  ##description: evidence-based annotation of the human genome (NCBI36)
  ##provider: GENCODE
  ##contact: fsk@sanger.ac.uk
  ##format: gtf 2.2
  ##date: 2008-10-02
  chr1	HAVANA	exon	1873	1920	.	+	.	 gene_id "OTTHUMG00000000961"; transcript_id "OTTH  UMT00000002844"; transcript_type "unprocessed_pseudogene"; transcript_status "UNKNOWN"; gene_type "unprocessed_pseudogene"; gene_status "NOVEL"; gene_name "RP11-34P13.1"; transcript_name "RP11-34P13.1-001"; level 2;
  chr1	HAVANA	exon	2042	2090	.	+	.	 gene_id "OTTHUMG00000000961"; transcript_id "OTTHUMT00000002844"; transcript_type "unprocessed_pseudogene"; transcript_status "UNKNOWN"; gene_type "unprocessed_pseudogene"; gene_status "NOVEL"; gene_name "RP11-34P13.1"; transcript_name "RP11-34P13.1-001"; level 2;
  chr1	HAVANA	exon	2476	2560	.	+	.	 gene_id "OTTHUMG00000000961"; transcript_id "OTTHUMT00000002844"; transcript_type "unprocessed_pseudogene"; transcript_status "UNKNOWN"; gene_type "unprocessed_pseudogene"; gene_status "NOVEL"; gene_name "RP11-34P13.1"; transcript_name "RP11-34P13.1-001"; level 2;
  chr1	HAVANA	exon	2838	2915	.	+	.	 gene_id "OTTHUMG00000000961"; transcript_id "OTTHUMT00000002844"; transcript_type "unprocessed_pseudogene"; transcript_status "UNKNOWN"; gene_type "unprocessed_pseudogene"; gene_status "NOVEL"; gene_name "RP11-34P13.1"; transcript_name "RP11-34P13.1-001"; level 2;
  chr1	HAVANA	exon	3084	3237	.	+	.	 gene_id "OTTHUMG00000000961"; transcript_id "OTTHUMT00000002844"; transcript_type "unprocessed_pseudogene"; transcript_status "UNKNOWN"; gene_type "unprocessed_pseudogene"; gene_status "NOVEL"; gene_name "RP11-34P13.1"; transcript_name "RP11-34P13.1-001"; level 2;
  chr1	HAVANA	exon	3316	3533	.	+	.	 gene_id "OTTHUMG00000000961"; transcript_id "OTTHUMT00000002844"; transcript_type "unprocessed_pseudogene"; transcript_status "UNKNOWN"; gene_type "unprocessed_pseudogene"; gene_status "NOVEL"; gene_name "RP11-34P13.1"; transcript_name "RP11-34P13.1-001"; level 2;

2.gencode_data.rel1.v2.regions.txt
  list of regions (genomic coordinates) where only HAVANA annotation (level 1 & 2) can be found.

3.gencode_data.rel1.v2.regions_with_ids.txt
  list of regions (genomic coordinates) where only HAVANA annotation (level 1 & 2) can be found, 
with all OTT-ids from the region listed.


