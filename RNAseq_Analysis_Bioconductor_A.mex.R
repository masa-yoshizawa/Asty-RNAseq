### RNAseq analysis by following the protocol 
## http://www.bioconductor.org/help/workflows/rnaseqGene/

### Please note that this script is not only for R but also include shell script.####################
### These are labeled as "Shell" or "R" in followings.

#### RNAseq analysis using SRA data set

## Shell ################################################################
## Shell in followings ################################################################

## SRA tool
## download SRA tool kit at http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
## manual: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc
## fastq-dump or illumina-dump will be useful to transfer the SRA data 
## fastq-dump to download
## test

######## Folder structure
/yourLocalFolder/RNAseq/
/yourLocalFolder/RNAseq/fastq		### a folder for SRA fast-q files
/yourLocalFolder/RNAseq/aligned 	### a folder for aligned files
/yourLocalFolder/RNAseq/Astyanax_mexicanus.AstMex102.85.gtf		### GTF file from ensembl.org
/yourLocalFolder/AstMex102.dna_sm.toplevel/		### a folder for annotated genome sequences used by STAR
/yourLocalFolder/Astyanax_mexicanus.AstMex102.dna_sm.toplevel.fa	### FASTA file of genome sequence downloaded from ensembl.org
/yourLocalFolder/RNAseq/SraRunTable_emb.csv		### Run table for Astyanax trascriptome (BioProject Acc#: PRJNA258661)
/yourLocalFolder/RNAseq/SraRunTable_72hpf.csv	### Only 72hpf Run table for Astyanax trascriptome (BioProject Acc#: PRJNA258661)
/yourLocalFolder/RNAseq/files.txt		### a text file containing SRR numbers


############################################################################
############################################################################

### download SRA files in the list.txt
cd /yourLocalFolder/RNAseq/fastq
for f in `cat ../files.txt`; do /Users/masato/sratoolkit.2.5.7-mac64/bin/fastq-dump $f; done

### prepare genome sequence for STAR, aligning to the genomic sequence of "Astyanax_mexicanus.AstMex102.dna_sm.toplevel.fa" downloaded from ensembl.org
cd /yourLocalFolder/RNAseq
/Users/masato/STAR --runThreadN 5 \
--runMode genomeGenerate \
--sjdbGTFfile /yourLocalFolder/RNAseq/Astyanax_mexicanus.AstMex102.84.gtf \
--genomeDir /yourLocalFolder/AstMex102.dna_sm.toplevel/ \
--genomeFastaFiles /yourLocalFolder/Astyanax_mexicanus.AstMex102.dna_sm.toplevel.fa
## run the alignments
for f in `cat files.txt`; do /Users/masato/STAR --genomeDir /yourLocalFolder/AstMex102.dna_sm.toplevel/ \
--outSAMtype BAM Unsorted \
--readFilesIn fastq/$f.fastq \
--runThreadN 5 --outFileNamePrefix aligned/$f.; done

## R ################################################################
## R scripts in followings ################################################################
####### 
## in R  convert .sam into .bam
source("https://bioconductor.org/biocLite.R")
biocLite("Rsamtools")
library("Rsamtools")
asBam("/yourLocalFolder/RNAseq/aligned/SRR1556193.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556193"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556198.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556198"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556201.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556201"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556204.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556204"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556205.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556205"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556258.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556258"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556273.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556273"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556275.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556275"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556277.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556277"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556301.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556301"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556303.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556303"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556304.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556304"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556305.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556305"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556306.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556306"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556307.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556307"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556308.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556308"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556309.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556309"); asBam("/yourLocalFolder/RNAseq/aligned/SRR1556316.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR1556316"); asBam("/yourLocalFolder/RNAseq/aligned/SRR639083.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR639083"); asBam("/yourLocalFolder/RNAseq/aligned/SRR639085.Aligned.out.sam", "/yourLocalFolder/RNAseq/aligned/SRR639085")

library("GenomicFeatures")
mydir <- "/yourLocalFolder/RNAseq"
gtffile <- file.path(mydir,"Astyanax_mexicanus.AstMex102.84.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())

ebg <- exonsBy(txdb, by="gene")
biocLite("GenomicAlignments"); biocLite("BiocParallel")
library("GenomicAlignments")

mydir <- "/yourLocalFolder/RNAseq"
csvfile <- file.path(mydir,"SraRunTable_emb.csv")
mysampleTable <- read.csv(csvfile,header = T)
mydir <- "/yourLocalFolder/RNAseq/aligned"
library("DESeq2")
myfilenames <- file.path(mydir, paste0(mysampleTable$Run_s, ".Aligned.out.bam"))
file.exists(myfilenames)
bamfiles <- BamFileList(myfilenames)
seqinfo(bamfiles[1])
##### check the bamfiles info

se <- summarizeOverlaps(features=ebg, reads=bamfiles, mode="Union", singleEnd=TRUE, ignore.strand=TRUE)
colData(se) <- DataFrame(mysampleTable)
countdata <- assay(se)
coldata <- colData(se)

se
##class: RangedSummarizedExperiment 
##dim: 23772 24 
##metadata(0):
##assays(1): counts
##rownames(23772): ENSAMXG00000000001 ENSAMXG00000000002 ... ENSAMXG00000027264 ENSAMXG00000027265
##rowData names(0):
##colnames(24): SRR1555606.Aligned.out.bam SRR1556144.Aligned.out.bam ... SRR1556309.Aligned.out.bam SRR1556316.Aligned.out.bam
##colData names(0):

### change the order of column level
se$Sample_Name_s <- relevel(se$Sample_Name_s, "Surface_fish")


##################################################################################################################################################################
###### set "design= ~ Sample_Name_s + age_s + Sample_Name_s:age_s" for full comparizon in two-way, and drop "Sample_Name_s + age_s" to calculate the interaction.
##################################################################################################################################################################
ddsTC <- DESeqDataSet(se, design=~ Sample_Name_s + age_s + Sample_Name_s:age_s)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ Sample_Name_s + age_s)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
resultsNames(ddsTC)
## check the ddsTC info
res_SfCf_age <- results(ddsTC, name="Sample_Name_sPachon_Cavefish.age_s72.hpf", test="Wald")
res_SfCf_age[which.min(resTC$padj),]
### check the result of res_SfCf_age. Followings are output of it
#log2 fold change (MLE): Sample Name sPachon Cavefish.age s72.hpf 
#Wald test p-value: Sample Name sPachon Cavefish.age s72.hpf 
#DataFrame with 1 row and 6 columns
#                    baseMean log2FoldChange     lfcSE      stat        pvalue         padj
#                   <numeric>      <numeric> <numeric> <numeric>     <numeric>    <numeric>
#ENSAMXG00000021623   1241.15       -4.97287  0.210131 -23.66557 8.160012e-124 8.98499e-120
#...

write.csv(res_SfCf_age, file="/yourLocalFolder/RNAseq/results_SfCf_age.csv") ### save the result to analyze further in MS Excel

##################################################################################################################################################################
###### import only 72hpf files (written down in SraRunTable_72hpf.csv
##################################################################################################################################################################
csvfile <- file.path("/yourLocalFolder/RNAseq","SraRunTable_72hpf.csv")
mysampleTable <- read.csv(csvfile,header = T)
myfilenames <- file.path(mydir, paste0(mysampleTable$Run_s, ".Aligned.out.bam"))
file.exists(myfilenames)
bamfiles <- BamFileList(myfilenames)
seqinfo(bamfiles[1])

se <- summarizeOverlaps(features=ebg, reads=bamfiles, mode="Union", singleEnd=TRUE, ignore.strand=TRUE)
colData(se) <- DataFrame(mysampleTable)
countdata <- assay(se)
coldata <- colData(se)

se
### check 'se' info. Followings are output of it
# class: RangedSummarizedExperiment 
# dim: 23772 6 
# metadata(0):
# assays(1): counts
# rownames(23772): ENSAMXG00000000001 ENSAMXG00000000002 ... ENSAMXG00000027264 ENSAMXG00000027265
# rowData names(0):
# colnames: NULL
# colData names(31): BioSample_s Experiment_s ... source_s tissue_s

se$Sample_Name_s
### check 'se$Sample_Name_s' info. Followings are output of it
# [1] Surface_fish    Surface_fish    Surface_fish    Pachon_Cavefish Pachon_Cavefish Pachon_Cavefish
# Levels: Pachon_Cavefish Surface_fish
se$Sample_Name_s <- relevel(se$Sample_Name_s, "Surface_fish")
se$Sample_Name_s
### check 'se$Sample_Name_s' info after changing the level. Followings are output of it
# Levels: Surface_fish Pachon_Cavefish
se$age_s
### check 'se$age_s' info. Followings are output of it
# [1] 72 hpf 72 hpf 72 hpf 72 hpf 72 hpf 72 hpf
# Levels: 72 hpf   --> only 72hpf data is available
dds <- DESeqDataSet(se, design= ~ Sample_Name_s)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ Sample_Name_s)
ddsMat
### check 'ddsMat' info. Followings are output of it
# class: DESeqDataSet 
# dim: 23772 6 
# metadata(1): version
# assays(1): counts
# rownames(23772): ENSAMXG00000000001 ENSAMXG00000000002 ... ENSAMXG00000027264 ENSAMXG00000027265
# rowData names(0):
# colnames: NULL
# colData names(31): BioSample_s Experiment_s ... source_s tissue_s

dds <- dds[ rowSums(counts(dds)) > 1, ]  ## removing no-count genes

rld <- rlog(dds, blind=FALSE)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
resSfCf <- results(dds, contrast=c("Sample_Name_s", "Surface_fish","Pachon_Cavefish"))
write.csv(resSfCf, file="/yourLocalFolder/RNAseq/results_SfCf_72hpf.csv")  ### save the result to analyze further in MS Excel
