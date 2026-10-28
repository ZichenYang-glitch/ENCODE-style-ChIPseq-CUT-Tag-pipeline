
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

21. gencode.vX.metadata.MGI.gz:
  MGI approved gene symbol.
  1 - transcript id
  2 - MGI gene symbol
  3 - MGI unique id

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
chromosomes 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, X, Y.

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
* tag "upstream_ATG": an upstream ATG exists when a downstream ATG is better supported.
* tag "upstream_uORF": a low confidence upstream ATG existing in other coding variant would lead to 
  NMD in this trancript, that uses the high confidence canonical downstream ATG.


Please note: if start codons are split between two exons, two start-codon features will be listed.
Please note: pre-release 4, "cds_start_NF" was listed as "cds start not found", etc.
Please note: pre-release 6, "seleno" tags included the selenocysteine position as the amino acid 
number within the protein, these are now given as genomic coordinates as separate GTF 
features.





#################################################################################
Release M39 (June 2026)
#################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2025.
This release corresponds to Ensembl version 116.


#################################################################################################
Release M38 (September 2025)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to February 2025.
This release corresponds to Ensembl version 115.



#################################################################################################
Release M37 (May 2025)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2024.
This release corresponds to Ensembl version 114.

The following files, which were restricted to reference chromosomes, now contain annotation on
primary assembly and alternate regions:
-gencode.vM37.long_noncoding_RNAs.gtf, gencode.vM37.long_noncoding_RNAs.gff3
-gencode.vM37.transcripts.fa
-gencode.vM37.pc_transcripts.fa
-gencode.vM37.pc_translations.fa
-gencode.vM37.lncRNA_transcripts.fa

The gencode.vM37.pc_translations.fa file now includes all translation sequences. Previously,
sequences with internal stop codons were normally excluded.



#################################################################################################
Release M36 (October 2024)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to May 2024.
This release corresponds to Ensembl version 113 in the GRCm39 assembly.



#################################################################################################
Release M35 (May 2024)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to September 2023.
This release corresponds to Ensembl version 112 in the GRCm39 assembly.



#################################################################################################
Release M34 (January 2024)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to March 2023.
This release corresponds to Ensembl version 111 in the GRCm39 assembly.



#################################################################################################
Release M33 (July 2023)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to December 2022.
This release corresponds to Ensembl version 110 in the GRCm39 assembly.



#################################################################################################
Release M32 (February 2023)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2022
This release corresponds to Ensembl version 109 in the GRCm39 assembly.



#################################################################################################
Release M31 (October 2022)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to April 2022
This release corresponds to Ensembl version 108 in the GRCm39 assembly.

The following biotype has been added:
-protein_coding_CDS_not_defined: Transcript that belongs to a protein_coding gene 
  and doesn't contain an ORF. Replaces the processed_transcript transcript biotype
  in protein_coding genes.

Annotation files in gtf and gff3 format having the basic set of transcripts in the 
primary assembly are now included in the release.



#################################################################################################
Release M30 (July 2022)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to January 2022
This release corresponds to Ensembl version 107 in the GRCm39 assembly.

The following biotypes have been added:
-protein_coding_LoF: Not translated in the reference genome owing to a SNP/DIP 
  but in other individuals/haplotypes/strains the transcript is translated. 
  This biotype replaces the polymorphic_pseudogene transcript biotype.
-artifact: Annotated on artifactual regions of the genome assembly.

The following tags or attributes have been added:
-artifactual_duplication / artif_dupl
-readthrough_gene



################################################################################################
Release M29 (April 2022)
################################################################################################
This release contains updates to the Ensembl-HAVANA merged annotation up to August 2021.
This release corresponds to Ensembl version 106 in the GRCm39 assembly.



################################################################################################
Release M28 (December 2021)
################################################################################################
This release contains updates to the Ensembl-HAVANA merged annotation up to May 2021.
This release corresponds to Ensembl version 105 in the GRCm39 assembly.

Clone-based gene names have been retired from the release files. From this release onwards, 
genes without a name in MGI, EntrezGene, RFAM or miRBase will get their gene id as default name.   



################################################################################################
Release M27 (May 2021)
################################################################################################
This release contains updates to the Ensembl-HAVANA merged annotation up to December 2020.
This release corresponds to Ensembl version 104 in the GRCm39 assembly.



#################################################################################################
Release M26 (February 2021)
#################################################################################################
This release contains updates to the Ensembl-HAVANA merged annotation up to August 2020.
This release corresponds to Ensembl version 103 in the GRCm39 assembly.



#################################################################################################
Release M25 (April 2020)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to November 2019.
This release corresponds to Ensembl version 100 in the GRCm38 assembly.



#################################################################################################
Release M24 (January 2020)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2019.
This release corresponds to Ensembl version 99 in the GRCm38 assembly.



#################################################################################################
Release M23 (September 2019)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to May 2019.
This release corresponds to Ensembl version 98 in the GRCm38 assembly.

New tag added: stop_codon_readthrough.



#################################################################################################
Release M22 (July 2019)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to February 2019.
This release corresponds to Ensembl version 97 in the GRCm38 assembly.

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

The following attribute has been added:
-MGI_id: Unique stable id provided by the MGI database for each gene with an approved symbol.



#################################################################################################
Release M21 (April 2019)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to November 2018.
This release corresponds to Ensembl version 96 in the GRCm38 assembly.



#################################################################################################
Release M20 (January 2019)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2018.
This release corresponds to Ensembl version 95 in the GRCm38 assembly.



#################################################################################################
Release M19 (October 2018)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to May 2018.
This release corresponds to Ensembl version 94 in the GRCm38 assembly.



#################################################################################################
Release M18 (July 2018)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to March 2018.
This release corresponds to Ensembl version 93 in the GRCm38 assembly.



#################################################################################################
Release M17 (April 2018)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to November 2017.
This release corresponds to Ensembl version 92 in the GRCm38 assembly.



#################################################################################################
Release M16 (December 2017)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2017.
This release corresponds to Ensembl version 91 in the GRCm38 assembly.



#################################################################################################
Release M15 (August 2017)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to May 2017.
This release corresponds to Ensembl version 90 in the GRCm38 assembly.

The antisense biotype is now called "antisense_RNA".



#################################################################################################
Release M14 (May 2017)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to January 2017.
This release corresponds to Ensembl version 89 in the GRCm38 assembly.



#################################################################################################
Release M13 (March 2017)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to October 2016.
This release corresponds to Ensembl version 88 in the GRCm38 assembly.

New tags added:
-3_nested_supported_extension
-3_standard_supported_extension
-454_RNA_Seq_supported
-5_nested_supported_extension
-5_standard_supported_extension



#################################################################################################
Release M12 (December 2016)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2016.
This release corresponds to Ensembl version 87 in the GRCm38 assembly.

The 'gene_status' and 'transcript_status' attributes have been removed from the GENCODE GTF and 
GFF3 files.



#################################################################################################
Release M11 (October 2016)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to March 2016.
This release corresponds to Ensembl version 86 in the GRCm38 assembly.



#################################################################################################
Release M10 (July 2016)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to January 2016.
This release corresponds to Ensembl version 85 in the GRCm38 assembly.



#################################################################################################
Release M9 (March 2016)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to October 2015.
This release corresponds to Ensembl version 84 in the GRCm38 assembly.

Two changes have been introduced in the GFF3 files:
- the 'UTR' feature type has been replaced with 'five_prime_UTR' and 'three_prime_UTR';
- CDS, start_codon and stop_codon features that are split across different exons of the same 
  transcript now have a shared ID to conform with the SO GFF3 specification.



#################################################################################################
Release M8 (December 2015)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2015.
This release corresponds to Ensembl version 83 in the GRCm38 assembly.



#################################################################################################
Release M7 (September 2015)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to June 2015.
This release corresponds to Ensembl version 82 in the GRCm38 assembly.

New gene and transcript biotype introduced in this release:
-bidirectional_promoter_lncrna
-transcribed_unitary_pseudogene

New tags introduced in this release:
-overlapping_locus
-retrogene
-inferred_transcript_model



#################################################################################################
Release M6 (July 2015)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to March 2015.
This release corresponds to Ensembl version 81 in the GRCm38 assembly.


New files added in this release:

a)Basic annotation on the reference chromosomes and on all sequence regions:
-gencode.vM6.basic.annotation.gtf.gz
-gencode.vM6.basic.annotation.gff3.gz
-gencode.vM6.chr_patch_hapl_scaff.basic.annotation.gtf.gz
-gencode.vM6.chr_patch_hapl_scaff.basic.annotation.gff3.gz

b)Nucleotide sequences of all annotated transcripts on the reference chromosomes:
-gencode.vM6.transcripts.fa.gz



#################################################################################################
Release M5 (May 2015)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to December 2014.
This release corresponds to Ensembl version 80 in the GRCm38 assembly.

New gene and transcript biotypes introduced in M5:
-macro_lncRNA
-ribozyme
-sRNA
-scaRNA 


Seven new APPRIS tags have been added in this release:

-appris_principal_1: (This flag corresponds to the older flag "appris_principal") where the transcript expected to code for the main functional isoform based solely on the core modules in the APPRIS database. The APPRIS core modules map protein structural and functional information and cross-species conservation to the annotated variants.

-appris_principal_2: (This flag corresponds to the older flag "appris_candidate_ccds") Where the APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human protein coding genes), the database chooses two or more of the CDS variants as "candidates" to be the principal variant.
If one (but no more than one) of these candidates has a distinct CCDS identifier it is selected as the principal variant for that gene. A CCDS identifier shows that there is consensus between RefSeq and GENCODE/Ensembl for that variant, guaranteeing that the variant has cDNA support.

-appris_principal_3: Where the APPRIS core modules are unable to choose a clear principal variant and there more than one of the variants have distinct CCDS identifiers, APPRIS selects the variant with lowest CCDS identifier as the principal variant. The lower the CCDS identifier, the earlier it was annotated.
Consensus CDS annotated earlier are likely to have more cDNA evidence. Consecutive CCDS identifiers are not included in this flag, since they will have been annotated in the same release of CCDS. These are distinguished with the next flag. 

-appris_principal_4: (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_ccds") Where the APPRIS core modules are unable to choose a clear principal CDS and there is more than one variant with a distinct (but consecutive) CCDS identifiers, APPRIS selects the longest CCDS isoform as the principal variant.

-appris_principal_5: (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_seq") Where the APPRIS core modules are unable to choose a clear principal variant and none of the candidate variants are annotated by CCDS, APPRIS selects the longest of the candidate isoforms as the principal variant.

-appris_alternative_1: Candidate transcript(s) models that are conserved in at least three tested non-primate species.

-appris_alternative_2: Candidate transcript(s) models that appear to be conserved in fewer than three tested non-primate species.


The UTR lines in the GTF/GFF3 files now include 'exon_number' and 'exon_id' attributes.

Also, in this release the transcript attributes were removed from the gene lines in all GTF and GFF3 annotation files. ie. the transcript_id transcript_type , transcript_status, transcript_name , havana_transcript etc.

NOTE: GTF and GFF3 files containing gene annotation on the primary assembly only (main chromosomes and unplaced/unlocalized scaffolds) were added:
 - gencode.vM5.primary_assembly.annotation.gtf.gz
 - gencode.vM5.primary_assembly.annotation.gff3.gz
These files are meant to be used mainly for NGS analyses in conjunction with the equivalent primary assembly genome sequence provided:
 - GRCm38.primary_assembly.genome.fa.gz



#################################################################################################
Release M4 (December 2014)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to August 2014.
This release corresponds to Ensembl version 78 in the GRCm38 assembly.

Some new gene and transcript biotypes are introduced in M4:

TEC
IG_C_pseudogene
IG_D_pseudogene
TR_C_gene
TR_D_gene
TR_J_gene
TR_J_pseudogene
processed_pseudogene
transcribed_processed_pseudogene
transcribed_unprocessed_pseudogene
unitary_pseudogene
unprocessed_pseudogene

There are some transcript tag changes in M4:

appris_candidate_longest_seq tag is introduced in M4
and appris_candidate_highest_score and appris_canditate_longest are removed from M4
exp_conf tag is introduced in M4

NOTE: The files listed below were updated on 5-12-2014 to correct a small number of gene and transcript names containing a semicolon:
-gencode.vM4.annotation.gff3
-gencode.vM4.annotation.gtf
-gencode.vM4.chr_patch_hapl_scaff.annotation.gff3
-gencode.vM4.chr_patch_hapl_scaff.annotation.gtf
-gencode.vM4.pc_transcripts.fa
-gencode.vM4.pc_translations.fa



#################################################################################################
Release M3 (August 2014)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to April 2014.
This release corresponds to Ensembl version 76 in the GRCm38 assembly.

A new annotation file type is added in the Gencode releases which is GFF3.
Note:
In GFF3 the start and stop codon are included in the CDS.
In GTF the start codon is included in the CDS but the stop codon is included in the UTR.  

1 extra appris tag is added in this release in the transcript lines:

* tag "appris_candidate_highest_score": where there is no 'appris_principal' variant, the candidate with highest APPRIS score is selected as the primary variant.

1 extra attribute is added in the annotation files in all protein coding transcripts: the protein_id

For example:

chr1	HAVANA	transcript	3214482	3671498	.	-	.	gene_id "ENSMUSG00000051951.5"; transcript_id "ENSMUST00000070533.4"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "Xkr4"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "Xkr4-001"; level 2; protein_id "ENSMUSP00000070648.4"; tag "basic"; tag "appris_principal"; tag "CCDS"; ccdsid "CCDS14803.1"; havana_gene "OTTMUSG00000026353.2"; havana_transcript "OTTMUST00000065166.1";


1 extra metadata file is added:

gencode.vM3.metadata.EntrezGene: EntrezGene gene id associated to the transcript

NOTE: The files listed below were updated on 28-08-2014 to correct the following minor issues:
- some CCDS tags had been wrongly assigned to transcripts with the same start and end coordinates as a CCDS model but with different amino acid sequence.
- absence of semicolon after "havana_transcript" attributes in the GTF files.

gencode.vM3.annotation.gff3.gz
gencode.vM3.annotation.gtf.gz
gencode.vM3.chr_patch_hapl_scaff.annotation.gff3.gz
gencode.vM3.chr_patch_hapl_scaff.annotation.gtf.gz
gencode.vM3.long_noncoding_RNAs.gtf.gz 


#################################################################################################
Release M2 (December 2013)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to July 2013.
This GENCODE release corresponds to Ensembl version 74 in the GRCm38 assembly.


3 extra tags are added in this release in the trasncript lines:

* tag "appris_principal": transcript expected to code for the main functional isoform based on a range of protein features (APPRIS pipeline).
* tag "appris_candidate": where there is no single 'appris_principal' variant the main functional isoform will be translated from one of the 'appris_candidate' genes.
* tag "appris_candidate_longest": where there is no 'appris_principal' variant, the longest of the 'appris_candidate' variants is selected as the primary variant.


#################################################################################################
Release M1 (October 2013)
#################################################################################################
This is a remerge between the Ensembl annotation and updates from HAVANA up to July 2011.
This GENCODE release corresponds to Ensembl version 65 which was released in Ensembl browser in December 2011 in the NCBIM37 assembly.

This Ensembl release has been chosen as the first Mouse GENCODE release, Gencode M1, because this dataset was used in mouse ENCODE analysis (http://www.ncbi.nlm.nih.gov/pubmed/22889292)

This is the last Ensembl release in assembly version NCBIM37.


Files which contain metadata associated to transcripts and genes, usually displayed in the UCSC browser, are part of the public release:

gencode.M1.metadata.Annotation_remark: remarks made during the manual annotation of the transcript
gencode.M1.metadata.Exon_supporting_feature: piece of evidence used in the annotation of an exon (usually peptides, mRNAs, ESTs)
gencode.M1.metadata.Gene_source: source of the gene annotation (Ensembl, Havana, Ensembl-Havana merged model or imported in the case of small RNA and mitochondrial genes)
gencode.M1.metadata.MGI: MGI approved gene symbol
gencode.M1.metadata.PDB: PDB entry associated to the transcript
gencode.M1.metadata.PolyA_feature: manually annotated polyA feature overlapping the transcript 3'-end (polyA_signal, polyA_site, pseudo_polyA) - Fields are 1)transcript_id, 2-3)polyA feature coordinates relative to the transcript, 4-7)polyA feature genomic coordinates, 8)type of polyA feature 
gencode.M1.metadata.Pubmed_id: Pubmed id of a publication associated to the transcript
gencode.M1.metadata.RefSeq: RefSeq RNA and/or protein associated to the transcript
gencode.M1.metadata.Selenocysteine: amino acid position of a selenocysteine residue in the transcript
gencode.M1.metadata.SwissProt: UniProtKB/SwissProt entry associated to the transcript
gencode.M1.metadata.Transcript_source: source of the transcript annotation (Ensembl, Havana, Ensembl-Havana merged model or imported in the case of small RNA and mitochondrial genes)
gencode.M1.metadata.Transcript_supporting_feature: piece of evidence used in the annotation of the transcript (usually peptides, mRNAs, ESTs)
gencode.M1.metadata.TrEMBL: UniProtKB/TrEMBL entry associated to the transcript






