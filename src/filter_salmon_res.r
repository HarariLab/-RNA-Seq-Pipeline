# Filtering Salmon results

library(stringr)

usage = '
Post-filtering of Salmon results. Output will be written to Salmon result directory:
	quant.genes.level_1.sf
	quant.level_1.sf

Usage:  Rscript filter_salmon_res.r
        <Salmon_result_dir>
        <GTF file>

'

args = commandArgs(T)

if (length(args) !=2 ){stop(usage)}

in.res.dir = args[1]
gtf = args[2]

cat('\nExtracting level 1 genes/transcripts...')

## read GTF file
#cat('Reading annotation file...')
gtf <- read.table(gtf, header =F, sep="\t", stringsAsFactors = F)
#cat('done.\n')

## read quantification files 
quant <- read.table(paste0(in.res.dir, '/quant.sf'), header =T, sep="\t", stringsAsFactors = F)
colnames(quant)[colnames(quant) == "Name"] = "TranscriptID"
colnames(quant)[colnames(quant) == "Length"] = "TranscriptLength"

quant.genes <- read.table(paste0(in.res.dir, '/quant.genes.sf'), header =T, sep="\t", stringsAsFactors = F)
colnames(quant.genes)[colnames(quant.genes) == "Name"] = "GeneID"
colnames(quant.genes)[colnames(quant.genes) == "Length"] = "GeneLength"

## extract genes levels
#gtf.genes <- gtf[gtf$V3 == 'gene', ]
gtf.genes <- gtf[gtf$V3 == 'gene' & !grepl('ERCC_', gtf$V1), ]
gtf.genes.ercc <- gtf[gtf$V3 == 'gene' & grepl('ERCC_', gtf$V1), ]
gtf.gg <- as.data.frame(str_split_fixed(gtf.genes$V9, ";", 6)[,c(1,3,2,4)])
gtf.gg.ercc <- as.data.frame(str_split_fixed(gtf.genes.ercc$V9, ";", 9)[,c(1,5,3,9)])
gtf.gg <- rbind(gtf.gg, gtf.gg.ercc)
colnames(gtf.gg) <- c('GeneID', 'GeneName', 'GeneBiotype', 'GeneAnnotationtLevel')

gtf.gg$GeneID <- gsub("gene_id ", "", gtf.gg$GeneID, fixed = TRUE)
gtf.gg$GeneID <- gsub(" ", '', gtf.gg$GeneID, fixed = TRUE)

gtf.gg$GeneAnnotationtLevel <- gsub("level ", "", gtf.gg$GeneAnnotationtLevel, fixed = TRUE)
gtf.gg$GeneAnnotationtLevel <- gsub(" ", "", gtf.gg$GeneAnnotationtLevel, fixed = TRUE)
gtf.gg$GeneAnnotationtLevel <- gsub(";", "", gtf.gg$GeneAnnotationtLevel, fixed = TRUE)

gtf.gg$GeneName <- gsub("gene_name ", "", gtf.gg$GeneName, fixed = TRUE)
gtf.gg$GeneName <- gsub(" ", "", gtf.gg$GeneName, fixed = TRUE)

gtf.gg$GeneBiotype <- gsub("gene_type ", "", gtf.gg$GeneBiotype, fixed = TRUE)
gtf.gg$GeneBiotype <- gsub(" ", "", gtf.gg$GeneBiotype, fixed = TRUE)

## create gene full name 
gtf.gg$GeneFullName <- paste0(gtf.gg$GeneID, '_', gtf.gg$GeneName)
## merge all
quant.genes <- merge(quant.genes, gtf.gg, all.x = T, sort = F)
quant.genes <- quant.genes[, c('GeneFullName','GeneID','GeneName','GeneBiotype','GeneAnnotationtLevel','GeneLength','EffectiveLength','TPM', 'NumReads')] 

## write updated result
write.table(quant.genes, file=paste0(in.res.dir, '/quant.genes.sf'), sep='\t', quote = F, row.names = F)

## write level1 genes 
quant.genes.level1 <- quant.genes[!is.na(quant.genes$GeneAnnotationtLevel) & quant.genes$GeneAnnotationtLevel == 1, ]
write.table(quant.genes.level1, file=paste0(in.res.dir, '/quant.genes.level_1.sf'), sep='\t', quote = F, row.names = F)


## extract transcripts 
gtf.tran <- gtf[gtf$V3 == 'transcript' & !grepl('ERCC_', gtf$V1), ]
gtf.tran.ercc <- gtf[gtf$V3 == 'transcript' & grepl('ERCC_', gtf$V1), ]
gtf.tt <- as.data.frame(str_split_fixed(gtf.tran$V9, ";", 8)[,c(2,6,5,7)])
gtf.tt.ercc <- as.data.frame(str_split_fixed(gtf.tran.ercc$V9, ";", 9)[,c(2,8,6,9)])
gtf.tt <- rbind(gtf.tt, gtf.tt.ercc)
colnames(gtf.tt) <- c('TranscriptID', 'TranscriptName', 'TranscriptBiotype', 'TranscriptAnnotationtLevel')

gtf.tt$TranscriptID <- gsub("transcript_id ", "", gtf.tt$TranscriptID, fixed = TRUE)
gtf.tt$TranscriptID <- gsub(" ", '', gtf.tt$TranscriptID, fixed = TRUE)

gtf.tt$TranscriptAnnotationtLevel <- gsub("level ", "", gtf.tt$TranscriptAnnotationtLevel, fixed = TRUE)
gtf.tt$TranscriptAnnotationtLevel <- gsub(" ", "", gtf.tt$TranscriptAnnotationtLevel, fixed = TRUE)
gtf.tt$TranscriptAnnotationtLevel <- gsub(";", "", gtf.tt$TranscriptAnnotationtLevel, fixed = TRUE)

gtf.tt$TranscriptName <- gsub("transcript_name ", "", gtf.tt$TranscriptName, fixed = TRUE)
gtf.tt$TranscriptName <- gsub(" ", "", gtf.tt$TranscriptName, fixed = TRUE)

gtf.tt$TranscriptBiotype <- gsub("transcript_type ", "", gtf.tt$TranscriptBiotype, fixed = TRUE)
gtf.tt$TranscriptBiotype <- gsub(" ", "", gtf.tt$TranscriptBiotype, fixed = TRUE)

## create gene full name 
gtf.tt$TranscriptFullName <- paste0(gtf.tt$TranscriptID, '_', gtf.tt$TranscriptName)
## merge all
quant.tran <- merge(quant, gtf.tt, all.x = T, sort = F)
quant.tran <- quant.tran[, c('TranscriptFullName','TranscriptID','TranscriptName','TranscriptBiotype','TranscriptAnnotationtLevel','TranscriptLength','EffectiveLength','TPM', 'NumReads')] 

## write updated result
write.table(quant.tran, file=paste0(in.res.dir, '/quant.sf'), sep='\t', quote = F, row.names = F)

## write level1 tran 
quant.tran.level1 <- quant.tran[!is.na(quant.tran$TranscriptAnnotationtLevel) & quant.tran$TranscriptAnnotationtLevel == 1, ]
write.table(quant.tran.level1, file=paste0(in.res.dir, '/quant.level_1.sf'), sep='\t', quote = F, row.names = F)

cat('done.\n')

