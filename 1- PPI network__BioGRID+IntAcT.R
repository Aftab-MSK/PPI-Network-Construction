#  21-12-2022
######################################################
## Loading library
######################################################
library(plyr) 
library(splitstackshape)

library("data.table")
library(tidyverse)

######################################################
## Loading data
######################################################

## Loading BioGrid (H.Sapiens)
# (https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-3.4.126/BIOGRID-ORGANISM-3.4.126.tab2.zip

BioGrid.raw <- read.csv ("./input/1. BIOGRID-ORGANISM-Homo_sapiens-4.4.216.tab2.txt", colClasses=as.vector (sapply("character", function (x) rep(x,24))), header=TRUE, sep="\t", stringsAsFactors = F)

# Loading Entrez_Gene to UniProtAC
uniprot <-  read.csv ("./input/4. uniprot_Short_20-12-2022.txt", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)
mapping_uniprot2GeneID <- concat.split.multiple(data = uniprot, split.cols = "Cross.reference..GeneID.", seps=";", direction = "long") ## Ignore the warning
mapping_uniprot2GeneID <- subset (mapping_uniprot2GeneID, select=c("Entry","Cross.reference..GeneID."))

uniprot.current <- read.csv ("./input/3. uniprot-_Full_20-12-2022.txt", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)
mapping_uniprot2GeneID <- merge (mapping_uniprot2GeneID, uniprot.current[,c("Entry","Entry.Name")],all.x=TRUE,sort=F)   # Entry.name =Entry.Name (i chnaged as per file)
######################################################
## Mapping BioGrid (adding UniProt ID)
######################################################

BioGrid.UniProtID <- merge (BioGrid.raw, mapping_uniprot2GeneID, by.x="Entrez.Gene.Interactor.A", by.y="Cross.reference..GeneID.", all.x=TRUE,sort=F)
names(BioGrid.UniProtID)[names(BioGrid.UniProtID)=="Entry"] <- "PARTICIPANT_A_Entry"
names(BioGrid.UniProtID)[names(BioGrid.UniProtID)=="Entry.Name"] <- "PARTICIPANT_A_Entry.name"     # Entry.name =Entry.Name (i chnaged as per file)

BioGrid.UniProtID <- merge (BioGrid.UniProtID, mapping_uniprot2GeneID, by.x="Entrez.Gene.Interactor.B", by.y="Cross.reference..GeneID.", all.x=TRUE,sort=F)
names(BioGrid.UniProtID)[names(BioGrid.UniProtID)=="Entry"] <- "PARTICIPANT_B_Entry"
names(BioGrid.UniProtID)[names(BioGrid.UniProtID)=="Entry.Name"] <- "PARTICIPANT_B_Entry.name"      # Entry.name =Entry.Name (i chnaged as per file)

BioGrid.UniProtID.unmapped <- BioGrid.UniProtID[(is.na(BioGrid.UniProtID$PARTICIPANT_A_Entry) | is.na(BioGrid.UniProtID$PARTICIPANT_B_Entry)),] # Record un-UniProt Entries
BioGrid.UniProtID <- BioGrid.UniProtID[!((is.na(BioGrid.UniProtID$PARTICIPANT_A_Entry) | is.na(BioGrid.UniProtID$PARTICIPANT_B_Entry))),] # Remove un-UniProt Entries

# write.csv (BioGrid.UniProtID.unmapped,"BioGrid.UniProtID.unmapped.csv")
BioGrid.slim <- data.frame (PARTICIPANT_A_Entry=BioGrid.UniProtID$PARTICIPANT_A_Entry, PARTICIPANT_B_Entry=BioGrid.UniProtID$PARTICIPANT_B_Entry, PARTICIPANT_A_Entry.name=BioGrid.UniProtID$PARTICIPANT_A_Entry.name, PARTICIPANT_B_Entry.name=BioGrid.UniProtID$PARTICIPANT_B_Entry.name, Source="BioGrid", Experiment=BioGrid.UniProtID$Experimental.System, Literature=BioGrid.UniProtID$Pubmed.ID, Interaction=BioGrid.UniProtID$Experimental.System.Type, stringsAsFactors = F)
write.csv (BioGrid.slim,"BioGrid_30.csv")

# BioGrid.slim <-BioGrid.slim[, -3]    ## remove unnessary colum before merging
# write.csv(mapping_uniprot2GeneID, "Uniprot+Entrez+UniPKB.csv")  ##IDS

######################################################
# Loading IntAct
######################################################
## The raw data of IntAct of can be downloaded from  (http://ftp.ebi.ac.uk/pub/databases/intact/2022-07-11/ 
# Loading IntAct
######################################################

# system ("curl -O ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip")
# system ("curl -O ftp://ftp.ebi.ac.uk/pub/databases/intact/2019-05-02/psi25/species/human.zip")
#IntAct.raw <- fread("unzip -p intact.zip",sep="\t", header=T, stringsAsFactors = F, showProgress = T, select=c(1,2,7,9,10,11,12), nThread=8, quote="")

IntAct.raw <- fread("./intact___22-11-2022/intact.txt",sep="\t", header=T, stringsAsFactors = F, showProgress = T, select=c(1,2,7,9,10,11,12), nThread=8, quote="")

#IntAct.raw <- read.csv (unz("intact.zip", "intact.txt"), header=F, sep="\t",stringsAsFactors = F)

IntAct.raw.tmp <- IntAct.raw[grepl ("9606",IntAct.raw$`Taxid interactor A`) & grepl ("9606",IntAct.raw$`Taxid interactor B`),] # Retrieve  Human-Human PPIs
IntAct.raw.tmp <- IntAct.raw.tmp[!(grepl ("intact|ensembl",IntAct.raw.tmp$`#ID(s) interactor A`) | grepl ("intact|ensembl", IntAct.raw.tmp$`ID(s) interactor B`)),] # Remove un-UniProt Entries

ac1 <- gsub ("uniprotkb:","",IntAct.raw.tmp$`#ID(s) interactor A`)
ac1 <- gsub ("-[A-Za-z0-9_]*","",ac1)
ac2 <- gsub ("uniprotkb:","",IntAct.raw.tmp$`ID(s) interactor B`)
ac2 <- gsub ("-[A-Za-z0-9_]*","",ac2)

IntAct.raw.tmp <- data.frame (PARTICIPANT_A_Entry=ac1, PARTICIPANT_B_Entry=ac2, Source="IntAct",
                              Experiment=IntAct.raw.tmp$`Interaction detection method(s)`, 
                              Literature=IntAct.raw.tmp$`Publication Identifier(s)`,
                              Interaction=IntAct.raw.tmp$`Interaction type(s)`, stringsAsFactors = F)

IntAct.raw.tmp <- IntAct.raw.tmp [!IntAct.raw.tmp$PARTICIPANT_A_Entry=="",]
IntAct.raw.tmp <- IntAct.raw.tmp [!IntAct.raw.tmp$PARTICIPANT_B_Entry=="",]

#IntAct.raw3 <- merge (IntAct.raw3, mapping_uniprot2GeneID, by.x="PARTICIPANT_A_Entry", by.y="Entry", all.x=TRUE)
#names(IntAct.raw3)[names(IntAct.raw3)=="Entry.name"] <- "PARTICIPANT_A_Entry.name"
#IntAct.raw3 <- merge (IntAct.raw3, mapping_uniprot2GeneID, by.x="PARTICIPANT_B_Entry", by.y="Entry", all.x=TRUE)
#names(IntAct.raw3)[names(IntAct.raw3)=="Entry.name"] <- "PARTICIPANT_B_Entry.name"

IntAct.slim <- IntAct.raw.tmp ; remove (IntAct.raw.tmp, ac1, ac2)



save (IntAct.slim, file="IntAct.Rdata")
write.csv (IntAct.slim, "IntAct_30.csv", row.names = F)


######################################################
# Merging BioGrid and IntAct
######################################################
interactions <- rbind (BioGrid.slim, IntAct.slim)
unique (interactions$Interaction)

###### Filter #######
interactions <- interactions[!interactions$Interaction == "genetic",] ## Discard "genetic" interactions
interactions <- interactions[!interactions$Interaction == "psi-mi:MI:0208(genetic interaction)",] ## Discard "genetic" interactions
sel.delete.experiment <- c("Co-localization","genetic interference","Synthetic Rescue","Synthetic Growth Defect","Synthetic Lethality")
sel.delete.interaction <- c("colocalization","genetic interaction")
interactions <- interactions [! grepl (paste0(sel.delete.experiment,collapse = "|"),interactions$Experiment), ]
interactions <- interactions [! grepl (paste0 (sel.delete.interaction, collapse="|"), interactions$Interaction),]
ppi <- interactions
str (ppi)
save (ppi,file="ppi_2022_1227.Rdata")

write.csv(ppi, file = "PPI-2023.csv")
