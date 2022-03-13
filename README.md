# RNA-Seq-Pipeline - a pipeline used to process bulk-RNA-Seq data. 
## Current Implementation:
-  Quality control on raw sequence data (FastQC).
- Revert and merge BAM files (Picard + Samtools). 
- Alignment and quantification for both linear and circRNAs (STAR).
- Post-Alignment QC (Picard tools: MarkDuplicates, CollectRNASeqMetrics, CollectAlignmentSummaryMetrics)
- Transcript Integrity (RSeQC, TIN).

## Genome Reference & Gene Models (GRCh38)
The genome reference and gene models for the entire pipeline were selected similar to the TOPMed pipeline.
- Reference file:
  - The GRCh38 reference genome and GENCODE 33 annotation, including the addition of ERCC spike-in annotations were used.
  - All ALT, HLA, and Decoy contigs were excluded from the reference genome (cannot be properly handled).
- Reference flat file:
  - The file was generated based on GENCODE 33 + ERCC spike-in annotations for Picard tools (CollectRnaSeqMetrics).
- A ribosomal RNA interval:
  - A list of rRNAs was generated based on GENCODE 33 + ERCC spike-in annotations for Picard tools (CollectRnaSeqMetrics).
- SimpleRepeats & RepeatMasker:
  - repeats were downloaded from UCSC Genome Browser. This file was used with the DCC too for circRNA detection.
- TIN calculations:
  - All protein coding, tag basic, and level 1 transcripts were extracted from the GENCODE 33 + ERCC spike-in annotations.
- Salmon Fasta file:
  - The file was generated based on the reference genome as well as the GENCODE 33 + ERCC spike-in annotations.
 
