degFile = snakemake@input[['degFile']]

assembly <- snakemake@params[['assembly']]

FC <- snakemake@params[['FC']]

adjp <- snakemake@params[['adjp']]

printTree <- snakemake@params[['printTree']]

library(GO.db)
library(topGO)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(GenomicFeatures)
library(Rgraphviz)

##----------load differentially expressed genes --------#
print("Loading differential expressed gene table")
print(degFile)

if(grepl('txt$|tsv$',degFile)){
    deg <- read.delim(file=degFile,header=TRUE,sep="\t")

    ### Since gene ID and gene name are in the diffexp file, push gene ID to row name and delete gene name col.
    rownames(deg) <- deg$Gene.stable.ID
    deg <- deg[,-c(1,2)]
    ###
}

##---------load correct Biomart------------------------#
print(getwd())
if (assembly == "hg19") {
    organismStr <- "hsapiens"
    geneID2GO <- get(load("./anno/biomaRt/hg19.Ens_75.biomaRt.GO.external.geneID2GO.RData"))
    xx <- get(load("./anno/biomaRt/GO.db.Term.list.rda"))
}
if (assembly == "hg38.89") {
    organismStr <- "hsapiens"
    ### to get to hg38 mappings ensembl 89!
    geneID2GO <- get(load("./anno/biomaRt/hg38.Ens_89.biomaRt.GO.external.geneID2GO.RData"))
    xx <- get(load("./anno/biomaRt/GO.db.Term.list.rda"))
}
if (assembly == "hg38.90") {
    organismStr <- "hsapiens"
    ### to get to hg38 mappings ensembl 90!
    geneID2GO <- get(load("./anno/biomaRt/hg38.Ens_90.biomaRt.GO.external.geneID2GO.RData"))
    xx <- get(load("./anno/biomaRt/GO.db.Term.list.rda"))
}
if (assembly == "mm10") {
    organismStr <- "mmusculus"
    geneID2GO <- get(load("./anno/biomaRt/mm10.Ens_78.biomaRt.geneAnno.Rdata.external.geneID2GO.RData"))
    xx <- get(load("./anno/biomaRt/GO.db.Term.list.rda"))
}

### Convert ENSMUSG ID into gene names with biomaRt ###
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")
if (assembly == "mm10") {
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
} else if (assembly == "hg38") {
  ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl) 
} else {
  stop("Must use either mm10 or hg38 assembly")
}

gene_id_name <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)

##-----------------------------------Functions--------------------------------------#
runGO <- function(geneList,xx=xx,otype,setName){
    setLength       <- sum(as.numeric(levels(geneList))[geneList]) 
    fname           <- paste(Dir, paste(setName, otype, "GO.txt", sep="_"), sep="/")
    GOData          <- new("topGOdata", ontology=otype, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher    <- runTest(GOData, algorithm = "classic", statistic = "fisher")## statistical test for topGO
    x               <- GenTable(GOData, classicFisher=resultFisher, topNodes=length(names(resultFisher@score)))## make go table for all terms
    x               <- data.frame(x)
    pVal            <- data.frame(pval=signif(resultFisher@score, 6)) ## get unrounded pvalue
    x$enrich        <- x$Significant/x$Expected ## calculate enrichment based on what you expect by chance
    x$p.unround     <- pVal[x$GO.ID,"pval"]## put unrounded pvalue in the table
    x$p.adj         <- signif(p.adjust(x$p.unround, method="BH"), 6)## calculate the adjusted pvalue with Benjamini & Hochberg correction
    x$log.p.adj     <- -log10(x$p.adj) ## convert adjusted p value to -log10 for plot magnitude
    #x$Term.full     <- sapply(x$GO.ID, FUN=function(n){Term(xx[[n]])}) ## get the full term name
    x <- x[order(x$GO.ID),]
    write.table(x, file=fname, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE) ## save the table
    ## you can print the tree if you want, but since I keep the list of all of them skip
    if(printTree>0){
        printGraph(GOData,## make the tree for the go data
                   resultFisher,
                   firstSigNodes = 5,
                   fn.prefix = sub("_GO.txt$", "", fname),
                   useInfo = "all",
                   pdfSW = TRUE
                   )
    }
    return(x)  
}

## function to make barplot of -log10 adjusted pvalues colored by enrichment

drawBarplot <- function(go, ontology, setName){
    go <- go[!go$p.adj > 0.01,]
    if(nrow(go)>1){
        #go$Term.full <- make.unique(paste(sapply(strsplit(as.character(substring(go$Term.full,1,50)), "\\,"), `[`, 1)))
        go$Term <- make.unique(paste(sapply(strsplit(as.character(substring(go$Term,1,50)), "\\,"), `[`, 1)))
        print(setName)
        go <- go[with(go, order(p.adj, -enrich)),]
        ## Currently there is a discrepency between xx and x, so we only use Term right now, not Term.full
        #go$Term.full <-factor(paste(go$Term.full), levels=rev(paste(go$Term.full))) ## sort table by adjusted p-value
        go$Term <-factor(paste(go$Term), levels=rev(paste(go$Term))) ## sort table by adjusted p-value
        ptitle <- paste(ontology, setName) ## plot title
        ptitle <- gsub("^.*/","",ptitle)
        pfname <- paste(setName,ontology,"pdf",sep=".")## name of png file
        if(nrow(go) < 20 ){
            toprange <- 1:nrow(go)
        }else{
            toprange <- 1:20
        }
        top <- go[toprange,]
        col <- colorRampPalette(c("white","navy"))(16)   
        pdf(file=paste(Dir, pfname, sep="/"),height=5,width=7)
        print({
           p <- ggplot(top, aes(y=log.p.adj, x=Term, fill=enrich)) + ## ggplot barplot function
               geom_bar(stat="identity",colour="black") +
               ggtitle(ptitle) +
               xlab("") + ylab("-log10(fdr)") +
               scale_fill_gradient(low=col[2], high=col[15], name="enrichment", limits=c(0,ceiling(max(top$enrich))))+
               coord_flip()+
               theme(panel.grid.major = element_line(colour = "grey"), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
               theme(text = element_text(size=8),
                     axis.text.x = element_text(vjust=1,color="black",size=8),
                     axis.text.y = element_text(color="black",size=8),
                     plot.title=element_text(size=10))   
       })
       dev.off()
    }
}

print("get up genes and make geneList")
up <- deg$padj < adjp & deg$log2FoldChange >= FC
up <- unique(rownames(deg[up,]))
up <- toupper(up)
up <- subset(gene_id_name, subset = gene_id_name$ensembl_gene_id %in% up)[,"external_gene_name"] ### take gene ID and convert to gene name

all <-unique(names(geneID2GO))
up.geneList <-  factor(as.integer(all %in% up))
names(up.geneList) <- all

up.setsize <- sum(as.numeric(levels(up.geneList))[up.geneList])
print("setsize for significant genes") 
up.setsize

adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.tsv$|\\.txt$|\\.rda$|\\.RData$","",degFile)

Dir <- sub("$", "/GOterms", dirname(comparison))
if(!(file.exists(Dir))) {
      dir.create(Dir,FALSE,TRUE)
}

### if upregulated gene name cannot be found in gene ID, print them out.
no_name_index = which(up %in% all == FALSE)
if (length(no_name_index) > 0) {
  print("Upregulated gene name not found in GO database!")
  print(gene_id_name[no_name_index, ])
}

if (up.setsize >= 2){

print("make GO table for the up genes")
go.UP.BP <- runGO(geneList=up.geneList,xx=xx,otype="BP",setName=paste(basename(comparison),"upFC",FC, "adjp", adjp, sep="."))
drawBarplot(go=go.UP.BP,ontology="BP",setName=paste(basename(comparison),"upFC",FC, "adjp", adjp, sep="."))

}else{
up_out = snakemake@output[[grep("diffexp.upFC", snakemake@output)]]
write.table('No Significant Genes', file=up_out)
}

print("get down genes and make geneList")
dn <- deg$padj < adjp & deg$log2FoldChange <= -FC
dn <- unique(rownames(deg[dn,]))
dn <- toupper(dn)
dn <- subset(gene_id_name, subset = gene_id_name$ensembl_gene_id %in% dn)[,"external_gene_name"] ### take gene ID and convert to gene name

all <-unique(names(geneID2GO))
dn.geneList <-  factor(as.integer(all %in% dn))
names(dn.geneList) <- all

dn.setsize <- sum(as.numeric(levels(dn.geneList))[dn.geneList])
print("setsize for significant genes") 
dn.setsize

### if downregulated gene name cannot be found in gene ID, print them out.
no_name_index = which(dn %in% all == FALSE)
if (length(no_name_index) > 0) {
  print("Upregulated gene name not found in GO database!")
  print(gene_id_name[no_name_index, ])
}

if(dn.setsize >= 2){

print("make GO table for down genes")
go.DN.BP <- runGO(geneList=dn.geneList,xx=xx,otype="BP",setName=paste(basename(comparison),"downFC",FC, "adjp", adjp, sep="."))
print("make barplot for down genes")
drawBarplot(go=go.DN.BP,ontology="BP",setName=paste(basename(comparison),"downFC",FC, "adjp", adjp, sep="."))

}else{
down_out = snakemake@output[[grep("diffexp.downFC", snakemake@output)]]
write.table('No Significant Gene', file=down_out)
}
