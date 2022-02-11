#!/usr/bin/env Rscript
library(plyr)
library(XML)
library(fgsea)
library(msigdbr)
library(reshape2)
library(R.utils)
### OptParse ########
library(optparse)
option_list <- list(
  make_option(c("-s", "--species"), type="character", default="hsa",
              help="Species: hsa, mmu, rno or dre. (Default: hsa"),
  make_option( "--database", type="character", default="kegg",
               help="Database to get gmt from. Supported: kegg, GO:BP, GO:MF, GO:CC, wikipathways, reactome, panglaodb. (Default: kegg)"),
  make_option(c("-o", "--outtable"), type="character", default="pathways.gmt",
              help="Output file name. (default= pathways.gmt)")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

species <- opt$species
db      <- opt$database
out     <- opt$outtable
#####################
### Functions #######
get_gmt <- function(species = "hsa", db ="kegg"){
  ### 1 - Get Entrez ID Mapping file (kegg, wikipathways)
  if (db=="kegg"|| db=="wikipathways") {
    print("Downloading entrez mapping file from NCBI...")
    if (species=="hsa"){entrezIDmapping.url = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"; ncbi.filename="Homo_sapiens.gene_info.gz"} else
      if (species=="mmu"){entrezIDmapping.url = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz"; ncbi.filename="Mus_musculus.gene_info.gz"} else
        if (species=="rno"){entrezIDmapping.url = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz"; ncbi.filename="Rattus_norvegicus.gene_info.gz"} else
          if (species=="dre"){entrezIDmapping.url = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Non-mammalian_vertebrates/Danio_rerio.gene_info.gz"; ncbi.filename="Danio_rerio.gene_info.gz"}
    temp <- tempfile()
    download.file(entrezIDmapping.url, temp)
    ncbi <- read.csv(gzfile(temp, ncbi.filename),
                     header = T, sep = "\t",
                     comment.char = "", quote = "", stringsAsFactors = F)
    closeAllConnections()
    mapping.data <- ncbi[, c("GeneID", "Symbol", "Synonyms")] # only Entrez IDs, Gene symbols and Synonyms
    map = setNames(mapping.data$Symbol, mapping.data$GeneID)
    rm(temp, ncbi)
    print("Entrez mapping file successfully downloaded.")
  }
  ### 2 - Get Pathway File
  #### 2.1 - KEGG
  if (db=="kegg") {
    print("Downloading Kegg pathway annotations...")
    kegg.pathways.url = paste0("http://rest.kegg.jp/link/",species,"/pathway")
    kid2name.url = paste0("http://rest.kegg.jp/list/pathway/",species)
    
    kegg.pathways <- read.csv(file = kegg.pathways.url, header = F, sep = "\t")
    colnames(kegg.pathways) <- c("PathwayID", "Gene")
    kid2name <- read.csv(file = kid2name.url, header = F, sep = "\t", quote="")
    colnames(kid2name) <- c("PathwayID", "PathwayName")
    print("Kegg pathway annotations successfully downloaded.")
    # map KEGG pathway id to pathway name
    print("Transforming to .gmt file...")
    kegg.pathways.new <- merge(kegg.pathways, kid2name, by = "PathwayID")
    # remove hsa prefix
    kegg.pathways.new$EntrezID = as.character(lapply(strsplit(as.character(kegg.pathways.new$Gene), split = ":"), "[", 2))
    kegg.pathways.new$PathwayID <- NULL
    kegg.pathways.new$Gene <- NULL
    kegg.pathways.mapped <- merge(kegg.pathways.new, mapping.data, by.x="EntrezID", by.y="GeneID", sort = F)
    kegg.pathways.mapped$PathwayName <- gsub(" ","_",kegg.pathways.mapped$PathwayName)
    pathways.mapped.aslist  <- split(kegg.pathways.mapped$Symbol, kegg.pathways.mapped$PathwayName)  # generate list of lists with pathway names as list names
  }
  #### 2.2 - WikiPathways
  if (db=="wikipathways"){
    print("Downloading Wikipathways annotations...")
    if (species=="hsa"){ wikipaths.suffix.filename="*Homo_sapiens.gmt"} else
      if (species=="mmu"){ wikipaths.suffix.filename="*Mus_musculus.gmt"} else
        if (species=="rno"){ wikipaths.suffix.filename="*Rattus_norvegicus.gmt"} else
          if (species=="dre"){ wikipaths.suffix.filename="*Danio_rerio.gmt"}
    wikipathways.gmt.url <- "https://wikipathways-data.wmcloud.org/current/gmt/"
    doc <- readLines(wikipathways.gmt.url)
    links <- doc[grep("*Homo_sapiens.gmt", doc)]
    # doc <- htmlParse(wikipathways.gmt.url)
    # links <- xpathSApply(doc, "//a/@href")
    # free(doc)
    wanted <- strsplit(links[grepl(wikipaths.suffix.filename, links)][1],
                       split = "./")[[1]][2]
    wanted <- strsplit(wanted, split = ">")[[1]][2]
    GetMe <- paste(wikipathways.gmt.url, wanted, sep = "")
    temp <- tempfile()
    download.file(GetMe, temp)
    wikipaths.file <- gmtPathways(temp)
    closeAllConnections()
    print("Wikipathways annotations successfully downloaded.")
    rm(wikipathways.gmt.url, doc, links, wanted, GetMe, temp)
    # mapping entrez ids to synonyms
    print("Transforming to .gmt file...")
    pathways.mapped.aslist <- lapply(wikipaths.file, function(x) mapvalues(x, from=mapping.data$GeneID, to=mapping.data$Symbol, warn_missing = F))
    names(pathways.mapped.aslist) <- unlist(lapply(strsplit(names(pathways.mapped.aslist), "%"), function(x) gsub(" ","_",x[1]) ))
  }
  #### 2.3 - Reactome
  if (db=="reactome"){
    if (species=="hsa"){ reactome.url="http://www.reactome.org/download/current/ReactomePathways.gmt.zip"} else
      if (species=="mmu"){ stop(paste("Reactome pathways for",species,"are not available.",sep = " "))} else
        if (species=="rno"){ stop(paste("Reactome pathways for",species,"are not available.",sep = " "))} else
          if (species=="dre"){ stop(paste("Reactome pathways for",species,"are not available.",sep = " "))}
    temp <- tempfile()
    download.file("http://www.reactome.org/download/current/ReactomePathways.gmt.zip", temp)
    pathways.mapped.aslist <- gmtPathways(unz(temp, "ReactomePathways.gmt"))
    names(pathways.mapped.aslist) <- unlist(lapply(strsplit(names(pathways.mapped.aslist), "%"),
                                                   function(x) gsub(" ","_",x[1]) ))
    rm(temp)
    closeAllConnections()
    print("Reactome annotations successfully downloaded.")
    print("Transforming to .gmt file...")
  }
  #### 2.4 - PanglaoDB
  if (db=="panglaodb"){
    #temp <- tempfile()
    path_to_tsv="/data/gmt/"
    ct.anno <- read.table(paste0(path_to_tsv,"panglaodb20190311.tsv"), sep="\t", header = T, quote = '')
    #rm(temp)
    print("PanglaoDB Cell type gene expression markers annotations successfully downloaded.")
    print("Transforming to .gmt file...")
    if (species=="hsa"){ ct.anno <- ct.anno[ct.anno$species %in% c("Mm Hs", "Hs"), ]} else
      if (species == "mmu"){ct.anno <- ct.anno[ct.anno$species %in% c("Mm Hs", "Mm"), ]} else
        if (species=="rno"){ stop("PanglaoDB is only available for Human and Mouse.")} else
          if (species=="dre"){stop("PanglaoDB is only available for Human and Mouse.")}
    
    
    ct <- as.list(dcast(ct.anno[, c("official.gene.symbol", "cell.type")],
                        official.gene.symbol ~ cell.type,
                        value.var = "official.gene.symbol")[, -1])
    ct <- merge(ct.anno[!duplicated(ct.anno[, c("cell.type", "organ")]), c("cell.type", "organ")],
                ldply(lapply(ct, function(x) x[!is.na(x)]), rbind),
                by.x="cell.type", by.y=".id", all=T)
    closeAllConnections()
    ct$cell.type <- gsub(" ","_", ct$cell.type)
    rownames(ct) <- ct$cell.type
    ct2 <- as.matrix( ct[,c(-1,-2)])
    pathways.mapped.aslist <- Map(Filter, list(Negate(is.na)), split(ct2, row(ct2)))
    names(pathways.mapped.aslist)<- ct$cell.type
  }
  #### 2.5 - GO
  if (db %in% c("GO:BP", "GO:MF", "GO:CC")){
    if (species=="hsa"){ species="Homo sapiens" }
    if (species=="mmu"){ species="Mus musculus" }
    if (species=="rno"){ species="Rattus norvegicus" }
    if (species=="dre"){ species="Danio rerio" }
    # taxids <- c(9606, 10090, 10116, 7955)
    # names(taxids) <- c("hsa", "mmu", "rno", "dre")
    col <- msigdbr_collections()
    GO = msigdbr(species = species,
                 category = "C5",
                 subcategory = db)
    pathways.mapped.aslist <- split(GO$gene_symbol, GO$gs_name)
    print(paste0(db, " annotations successfully imported."))
    print("Transforming to .gmt file...")
  }
  ### 3 - Transform nested lists into a gmt-file:
  # gmt1 <- ldply(lapply(pathways.mapped.aslist, function(x) x[!is.na(x)]), rbind)
  # gmt1$Description <- paste0(species, " from ",db)
  # gmt2 <- gmt1[,c(1,length(gmt1),3:length(gmt1)-1)]
  return(pathways.mapped.aslist)
}
#write.gmt <- function(df, filename){write.table(df, file = filename, sep = "\t", na = "", row.names = F, col.names = F, quote = F)}
#write.gmt2 <- function(df, filename){write.table(df, file = filename, sep = "\t", na = NULL, row.names = F, col.names = F, quote = F)}
write.gmt3 <- function(gmtlist, filename){
  if (file.exists(filename)){
    file.remove(filename)
  }
  gsname = names(gmtlist)
  for (i in seq_along(gmtlist)) {
    cat(gsname[i], paste0(opt$species, "_", opt$database), gmtlist[[i]], file = filename,
        append = TRUE, sep = "\t")
    cat("\n", append = TRUE, file = filename)
  }
  print(paste0("saving as ",filename))
}
#####################
### Main ############
print(paste("Trying to get gmt file for",species,"from",db,"..."))
write.gmt3(get_gmt(species, db), out)
print("Done.")


# spezien <- c("hsa","mmu","rno","dre")
# debes <- c("kegg","wikipathways","reactome","panglaodb")
# for (spezie in spezien) {
#   for(debe in debes) {
#     if (!file.exists(paste0(spezie,debe,".gmt"))){
#       try(write.gmt3(get_gmt(spezie,debe),paste0(spezie,debe,".gmt")))
#     }
#
#   }
#
# }
#####################
# write.gmt3(pathways.mapped.aslist, "gmt3wikipathways.gmt")
# test1 <- get_gmt()
# test2 <- get_gmt("hsa","reactome")
# write.gmt(ct, "Documents/panglao_ct.gmt")
# write.gmt2(test1, "gmt2.gmt")
# panglao_ct <- gmtPathways("Documents/panglao_ct.gmt")
# panglao_ct2 <- gmtPathways("Documents/panglao_ct2.gmt")
# abac1 <- unique(unlist(panglao_ct))
