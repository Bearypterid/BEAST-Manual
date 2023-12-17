message("Create a directory in your home directory for this script to save files to, then change the working directory in the next line")

setwd("~/RProject")#instead of /RProject it should be whatever your folder of interest is called
packages=c('biomaRt','BiocManager','plyr','dplyr','tibble','tidyselect')
install.packages(setdiff(packages, rownames(installed.packages())))
install.packages("devtools")
# devtools::install_version("dbplyr", version = "2.3.4") #This line is needed if the wrong dplyr version is installed and is interferring with biomart
library("biomaRt")
library("BiocManager")
library("plyr")
library("dplyr")
library("tibble")
library("tidyselect")
message("If biomart did not instal than uncomment and run the next two lines of code")
# BiocManager::install("biomaRt",force = TRUE)
# library("biomaRt")
biomaRt::listMarts()# These 3 lines are available for you to use and browse the other more niche databases
listEnsembl()
listEnsemblGenomes()

ens=useMart("ensembl")
anim=useEnsemblGenomes("metazoa_mart")
singleeuk=useEnsemblGenomes("protists_mart")
fungi=useEnsemblGenomes("fungi_mart")
plants=useEnsemblGenomes("plants_mart")

datasetsens=listDatasets(ens)
datasetsanim=listDatasets(anim)
datasetsprotist=listDatasets(singleeuk)
datasetsfungi=listDatasets(fungi)
datasetsplants=listDatasets(plants)

datasetused=datasetsens
datsetused=rownames(datasetused)=datasetused[,2]
message("Browse through yourdataset used to find your organisms of interest. Once you have found your ",
        "organisms of interest, copy and past the text from the information cell into the list function below. ",
        "The list function below will be used to search for and collect data on your species of interest, so make sure ",
        "the names in the list are exactly as written in the database you browsed, otehrwise there will be an error. ",
        "The first item in the list needs to be the database used, otherwise the script will not run.")

specieslist=list(ens,"Opossum","Wallaby","Platypus","Hyrax","Elephant genes","Koala","Tasmanian","Common wombat")

numberspecieslist=length(specieslist)
for (i in 1:numberspecieslist) {#This for loop creates the marts needed for the get bm function
  maht=paste("mart",i,sep = "")
  # speciesgen=paste("speciesgenes",i,sep = "")
  if (i==1) {
    assign(maht,specieslist[[i]])
  } else{
    searchspecies=specieslist[[i]]
    speciesdata=datasetused[c(paste(searchspecies)),]
    speciesdata=speciesdata[1,1]
    speciesgenes=useDataset(paste(speciesdata,sep = ""),mart=mart1)
    assign(maht,speciesgenes)
  }
}
message("Use speciesgenes0 to browse for possible genes of interest")

speciesgenes0=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name"),mart=mart2)#Use to find Genes

message("This is the gene of interest the rest of the code will use.")
filtervalue="PAOX"#Gene of interest

for (i in 2:numberspecieslist) {
  speciesgen=paste("speciesgenes",i,sep = "")
  speciesgenepresent=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name"),mart=get(paste("mart",i,sep = "")),filters = "external_gene_name", values = filtervalue)
  assign(speciesgen,speciesgenepresent)
  if (is.na(speciesgenepresent[1,1])) {
    message("The gene you have selected is not present in all species of interest.",
            " Please select a different gene and try again")
    {break}
  }
}
message("If your gene of interest was found in all species of interest, proceed.")
for (i in 2:numberspecieslist) {
  searchspecies=specieslist[[i]]
  speciesdata=datasetused[c(paste(searchspecies)),]
  geneseq=getSequence(id = filtervalue,type ="external_gene_name",seqType = "cdna",mart = get(paste("mart",i,sep = "")))
  geneid=get(paste("speciesgenes",i,sep = ""))
  geneid=geneid[1,1]
  Organism=speciesdata[1,2]
  geneseq[geneseq == filtervalue]=paste(geneid,Organism,filtervalue,"cDNA",sep = ".")
  if (i==2) {
    Fastaseq=geneseq
  }else{
  Fastaseq=rbind.data.frame(Fastaseq,geneseq)
  }
}
message("Make sure to delete any duplicate fasta files before the next line is run.",
        " Otherwise it will simply add your new sequences to the old one including any duplicates")
exportFASTA(Fastaseq,file = paste(filtervalue,"fasta",sep = "."))
