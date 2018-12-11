args <- commandArgs()

help <- function(){
    cat("runGOforDESeq2.R :
- This script generates the annotation files necessary for GO Term analysis.
- It takes no file input other than specifying the assembly you are using.
- Options for this include: hg19, hg38.89, hg38.90, mm9 and mm10. \n")
    cat("Usage: \n")
    cat("--assembly: genome assembly                             [ required ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args))){
    help()
} else {
    assembly <-sub('--assembly=', '', args[grep('--assembly=', args)])
}


library(biomaRt)
library(GO.db)

#GO.to.Term  <- as.data.frame(do.call(rbind, lapply(as.list(GOTERM), function(x){c(x@GOID, x@Term)})), stringsAsFactors = FALSE)
GO.to.Term <- get(load("/home/groups/CEDAR/anno/biomaRt2/GO.to.Term.df.rda"))
#save(GO.to.Term, file = "/home/groups/CEDAR/anno/biomaRt2/GO.to.Term.df.rda")

#xx          <- as.list(GOTERM)
xx <- get(load("/home/groups/CEDAR/anno/biomaRt2/GO.db.Term.list.rda"))
#save(xx, file = "/home/groups/CEDAR/anno/biomaRt2/GO.db.Term.list.rda")

if (assembly == "hg19") {
    annoFile <- "/home/groups/CEDAR/anno/biomaRt/hg19.Ens_75.biomaRt.GO"
    if(!(file.exists(annoFile))) {
        
        bm                    <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37) 
        EG2GO                 <- getBM(mart=bm, attributes=c('ensembl_gene_id','ensembl_transcript_id','external_gene_name','go_id'))
        EG2GO                 <- EG2GO[EG2GO$go_id != '',]
        allTermsWithIds        <- merge(EG2GO, GO.to.Term, by.x="go_id", by.y="V1")
        names(allTermsWithIds) <- sub("V2", "go_term_name", names(allTermsWithIds))
        save(allTermsWithIds, file=paste(annoFile, "allTermsWithIds.RData", sep="."))
        save(EG2GO, file=paste(annoFile, "goIds.RData", sep="."))
        geneID2GO              <- by(EG2GO$go_id,EG2GO$ensembl_gene_id,
                                     function(x) as.character(x)
                                     )
        save(geneID2GO,          file=paste(annoFile, "geneID2GO.RData", sep="."))
    
    }
}

if (assembly == "hg38.90") {
    annoFile <- "/home/groups/CEDAR/anno/biomaRt/hg38.Ens_90.biomaRt.GO"
    if(!(file.exists(annoFile))) {  
        bm                    <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="aug2017.archive.ensembl.org")
        EG2GO                 <- getBM(mart=bm, attributes=c('ensembl_gene_id','ensembl_transcript_id','external_gene_name','go_id'))
        EG2GO                 <- EG2GO[EG2GO$go_id != '',]
        allTermsWithIds        <- merge(EG2GO, GO.to.Term, by.x="go_id", by.y="V1")
        names(allTermsWithIds) <- sub("V2", "go_term_name", names(allTermsWithIds))
        save(allTermsWithIds, file=paste(annoFile, "allTermsWithIds.RData", sep="."))
        save(EG2GO, file=paste(annoFile, "goIds.RData", sep="."))
        geneID2GO              <- by(EG2GO$go_id,EG2GO$ensembl_gene_id,
                                     function(x) as.character(x)
                                     )
        save(geneID2GO,          file=paste(annoFile, "geneID2GO.RData", sep="."))
    }
}

if (assembly == "hg38.89") {
    annoFile <- "/home/groups/CEDAR/anno/biomaRt/hg38.Ens_89.biomaRt.GO"
    if(!(file.exists(annoFile))) {  
        bm                    <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="may2017.archive.ensembl.org")
        EG2GO                 <- getBM(mart=bm, attributes=c('ensembl_gene_id','ensembl_transcript_id','external_gene_name','go_id'))
        EG2GO                 <- EG2GO[EG2GO$go_id != '',]
        allTermsWithIds        <- merge(EG2GO, GO.to.Term, by.x="go_id", by.y="V1")
        names(allTermsWithIds) <- sub("V2", "go_term_name", names(allTermsWithIds))
        save(allTermsWithIds, file=paste(annoFile, "allTermsWithIds.RData", sep="."))
        save(EG2GO, file=paste(annoFile, "goIds.RData", sep="."))
        geneID2GO              <- by(EG2GO$go_id,EG2GO$ensembl_gene_id,
                                     function(x) as.character(x)
                                     )
        save(geneID2GO,          file=paste(annoFile, "geneID2GO.RData", sep="."))
    }
}


