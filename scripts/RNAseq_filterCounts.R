args <- commandArgs()

annoFile = snakemake@params[['anno']]

biotypes <- snakemake@params[['biotypes']]

countsFile <- snakemake@input[['countsFile']]

mito <- snakemake@params[['mito']]

##----------load counts------------#
print("Loading counts table")
print(countsFile)

## check if an rda file or tab sep
if(grepl('rda|RData|Rdata',countsFile)){
    counts <- get(load(file=countsFile))
}
if(grepl('txt|tsv',countsFile)){
    counts <- read.delim(file=countsFile)
}

##----------load counts------------#
print("Loading counts table")
print(countsFile)

## must be a tsv or txt tab sep file

counts <- read.delim(file=countsFile)

##----------load anno------------#
print("Loading annotation table")
print(annoFile)

## check if an rda file or tab sep
if(grepl('rda|RData|Rdata',annoFile)){
    anno <- get(load(file=annoFile))
}
if(grepl('txt|tsv',annoFile)){
    anno <- read.delim(file=annoFile)
}

if(strsplit(biotypes, split='\\,')[[1]]!=""){
    anno.sub <- anno[paste(anno$gene_biotype) %in% strsplit(biotypes, split='\\,')[[1]] ,]
    counts.sub <- counts[paste(counts$Genes) %in% unique(paste(anno.sub$ensembl_gene_id)) , ]
}else{
    print("no biotypes provided")
    counts.sub <- counts
}

if(mito==1){
    library(biomaRt)
    assembly = snakemake@params[["assembly"]]
    ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
    if (assembly == "mm10") {
      ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
    } else if (assembly == "hg38") {
      ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl) 
    } else {
      stop("Must use either mm10 or hg38 assembly")
    }
    gene_id_name <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl) # cols = gene ID, gene name

    print("tossing MT- genes")
    mito_index <- grep("^MT-", gene_id_name$external_gene_name, ignore.case=TRUE)
    mito_genes <- gene_id_name[mito_index, ]
    counts.sub <- counts.sub[ !(counts.sub$Genes %in% mito_genes$ensembl_gene_id), ]
}

write.table(counts.sub, file=sub(".txt", ".filt.txt", countsFile), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
