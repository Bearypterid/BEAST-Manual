setwd("~/RProject")
packages=c('dplyr','biomaRt','BiocManager','plyr','ggplot2')
install.packages(setdiff(packages, rownames(installed.packages())))
library ("biomaRt")
library("BiocManager")
library("plyr")
library("dplyr")
library("ggplot2")
ens=useMart("ensembl") #Using code from Rbasics file
listDatasets(ens)
datasets=listDatasets(ens)
opossum=useDataset("mdomestica_gene_ensembl",mart=ens)
atr=listAttributes(opossum)
listFilters(opossum)
wallaby=useDataset("neugenii_gene_ensembl",mart=ens)
atr2=listAttributes(wallaby)
platypus=useDataset("oanatinus_gene_ensembl",mart=ens)
hyrax=useDataset("pcapensis_gene_ensembl",mart=ens)
koala=useDataset("pcinereus_gene_ensembl",mart=ens)
tasdevil=useDataset("sharrisii_gene_ensembl",mart=ens)
wombat=useDataset("vursinus_gene_ensembl",mart=ens)
elephant=useDataset("lafricana_gene_ensembl",mart=ens)
OGenes=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name"),mart=opossum)
filtervalue="PAOX"
OGenes=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name","entrezgene_id","entrezgene_accession","entrezgene_description"),mart=opossum,filters = "external_gene_name", values = filtervalue)
WGenes=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name"),mart=wallaby,filters = "external_gene_name", values = filtervalue)
PGenes=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name"),mart=platypus,filters = "external_gene_name", values = filtervalue)
HGenes=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name"),mart=hyrax,filters = "external_gene_name", values = filtervalue)
KGenes=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name"),mart=koala,filters = "external_gene_name", values = filtervalue)
TDGenes=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name"),mart=tasdevil,filters = "external_gene_name", values = filtervalue)
WomGenes=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name"),mart=wombat,filters = "external_gene_name", values = filtervalue)
EGenes=getBM(attributes=c("ensembl_gene_id","name_1006","external_gene_name"),mart=elephant,filters = "external_gene_name", values = filtervalue)
oseq=getSequence(id = filtervalue,type ="external_gene_name",seqType = "cdna",mart = opossum)
id=OGenes[1,1]
Organism="gray_short-tailed_opossum_(Monodelphis_domestica)"
oseq[oseq == filtervalue]=paste(id,Organism,filtervalue,"cDNA",sep = ".")
wseq=getSequence(id = filtervalue,type ="external_gene_name",seqType = "cdna",mart = wallaby)
id=WGenes[1,1]
Organism="tammar_wallaby_(Notamacropus_eugenii)"
wseq[wseq == filtervalue]=paste(id,Organism,filtervalue,"cDNA",sep = ".")
pseq=getSequence(id = filtervalue,type ="external_gene_name",seqType = "cdna",mart = platypus)
id=PGenes[1,1]
Organism="platypus_(Ornithorhynchus_anatinus)"
pseq[pseq == filtervalue]=paste(id,Organism,filtervalue,"cDNA",sep = ".")
hseq=getSequence(id = filtervalue,type ="external_gene_name",seqType = "cdna",mart = hyrax)
id=HGenes[1,1]
Organism="rock_hyrax_(P._capensis)"
hseq[hseq == filtervalue]=paste(id,Organism,filtervalue,"cDNA",sep = ".")
kseq=getSequence(id = filtervalue,type ="external_gene_name",seqType = "cdna",mart = koala)
id=KGenes[1,1]
Organism="koala_(Phascolarctos_cinereus)"
kseq[kseq == filtervalue]=paste(id,Organism,filtervalue,"cDNA",sep = ".")
tdseq=getSequence(id = filtervalue,type ="external_gene_name",seqType = "cdna",mart = tasdevil)
id=TDGenes[1,1]
Organism="Tasmanian_devil_(Sarcophilus_harrisii)"
tdseq[tdseq == filtervalue]=paste(id,Organism,filtervalue,"cDNA",sep = ".")
wmseq=getSequence(id = filtervalue,type ="external_gene_name",seqType = "cdna",mart = wombat)
id=WomGenes[1,1]
Organism="common_wombat_(Vombatus_ursinus)"
wmseq[wmseq == filtervalue]=paste(id,Organism,filtervalue,"cDNA",sep = ".")
eseq=getSequence(id = filtervalue,type ="external_gene_name",seqType = "cdna",mart = elephant)
id=EGenes[1,1]
Organism="African_bush_elephant_(Loxodonta_africana)"
eseq[eseq == filtervalue]=paste(id,Organism,filtervalue,"cDNA",sep = ".")
Fastaseq=rbind.data.frame(oseq,wseq,pseq,hseq,kseq,tdseq,wmseq,eseq)
exportFASTA(Fastaseq,file = paste(filtervalue,"fasta",sep = "."))
# exportFASTA(oseq,file = "test.fasta")
# exportFASTA(oseq,file = paste("Opossum",filtervalue,"fasta",sep = "."))
# exportFASTA(wseq,file = paste("Wallaby",filtervalue,"fasta",sep = "."))
# exportFASTA(pseq,file = paste("Platypus",filtervalue,"fasta",sep = "."))
# exportFASTA(hseq,file = paste("Hyrax",filtervalue,"fasta",sep = "."))
# exportFASTA(kseq,file = paste("Koala",filtervalue,"fasta",sep = "."))
# exportFASTA(tdseq,file = paste("Tasmanian_Devil",filtervalue,"fasta",sep = "."))
# exportFASTA(wmseq,file = paste("Common_Wombat",filtervalue,"fasta",sep = "."))
# exportFASTA(eseq,file = paste("African_Elephant",filtervalue,"fasta",sep = "."))
# show(oseq)
# show(wseq)
# show(pseq)
# show(hseq)
# show(kseq)
# show(tdseq)
# show(wmseq)
# show(eseq)
# Opossum genes (ASM229v1)mdomestica_gene_ensembl
# Wallaby genes (Meug_1.0)neugenii_gene_ensembl
# Platypus genes (mOrnAna1.p.v1)oanatinus_gene_ensembl
# Hyrax genes (proCap1)pcapensis_gene_ensembl
# Koala genes (phaCin_unsw_v4.1)pcinereus_gene_ensembl
# Tasmanian devil genes (mSarHar1.11)sharrisii_gene_ensembl
# Common wombat genes (bare-nosed_wombat_genome_assembly)vursinus_gene_ensembl
# Elephant genes (Lox8afr3.0)lafricana_gene_ensembl
# go_id
# name_1006
# external_gene_name
# ensembl_gene_id