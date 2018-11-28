# Load Packages, Import Data, Create & Clean Physeq Objects ----------------------------------------------

library("phyloseq")
packageVersion("phyloseq")
library("ggplot2")
packageVersion("ggplot2")
theme_set(theme_bw())
library("reshape2")
packageVersion("reshape2")
sessionInfo()
set.seed(1)


#CORN
#Import OTU table ('#' in cell A1 has been removed)
Corn_OTU_table <- read.table("./Corn Input Files/Clean_OTU_table.txt", sep="\t", row.names=1, header=TRUE)
#Convert to matrix, then to phyloseq format.
Corn_OTU_table <- data.matrix(Corn_OTU_table)
Corn_OTU_table = otu_table(Corn_OTU_table, taxa_are_rows = TRUE)
#Import taxonomy (has been formatted in Excel to be tab-delimited)
Corn_TAX <- read.table("./Corn Input Files/Clean_otus_tax_assignments.txt", sep="\t", row.names=1, stringsAsFactors = FALSE, header=FALSE, na.strings="")
colnames(Corn_TAX) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
#convert to matrix, then phyloseq format.
Corn_TAX <- as.matrix(Corn_TAX)
Corn_TAX = tax_table(Corn_TAX)
#Import sample data
Corn_METADATA <- read.csv("./Corn Input Files/Metadata.csv", sep=",", row.names=1, header=TRUE, na.strings="")
#convert to phyloseq format.
Corn_METADATA<- sample_data(Corn_METADATA)
#Create phyloseq object:
Corn_physeq = phyloseq(Corn_OTU_table, Corn_TAX, Corn_METADATA)
Corn_physeq
#Check for singletons and remove if any.
any(taxa_sums(Corn_physeq) <= 1) #FALSE
#Remove taxa with Domain designated as "Eukaryota" (Chloroplast and mitochondria have already been removed on HPCC).
Corn_physeq <- subset_taxa(Corn_physeq, Domain!= "Eukaryota") #12695 taxa remaining
Corn_physeq
sample_sums(Corn_physeq)
sorted <- sort(sample_sums(Corn_physeq))
sorted <- as.data.frame(sorted)
sorted
#Actually do the rarefying here.
set.seed(1)
Corn.rfy <- rarefy_even_depth(Corn_physeq, sample.size = 22556, replace = FALSE)
#1873 OTUs were removed: no longer present in any sample after subsampling.


#BEAN
#Import OTU table ('#' in cell A1 has been removed)
Bean_OTU_table <- read.table("./Bean Input Files/Clean_OTU_table.txt", sep="\t", row.names=1, header=TRUE)
#Convert to matrix, then to phyloseq format.
Bean_OTU_table <- data.matrix(Bean_OTU_table)
Bean_OTU_table = otu_table(Bean_OTU_table, taxa_are_rows = TRUE)
#Import taxonomy (has been formatted in Excel to be tab-delimited)
Bean_TAX <- read.table("./Bean Input Files/Clean_otus_tax_assignments.txt", sep="\t", row.names=1, stringsAsFactors = FALSE, header=FALSE, na.strings="")
colnames(Bean_TAX) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
#convert to matrix, then phyloseq format.
Bean_TAX <- as.matrix(Bean_TAX)
Bean_TAX = tax_table(Bean_TAX)
#Import sample data
Bean_METADATA <- read.csv("./Bean Input Files/bean_seed_rhizo_map.txt", sep="\t", row.names=1, header=TRUE, na.strings="")
#convert to phyloseq format.
Bean_METADATA<- sample_data(Bean_METADATA)
#Create phyloseq object:
Bean_physeq = phyloseq(Bean_OTU_table, Bean_TAX, Bean_METADATA)
Bean_physeq
#Check for singletons and remove if any.
any(taxa_sums(Bean_physeq) <= 1) #FALSE
#Remove taxa with Domain designated as "Eukaryota" (Chloroplast and mitochondria have already been removed on HPCC).
Bean_physeq <- subset_taxa(Bean_physeq, Domain!= "Eukaryota") #8148 taxa remaining
Bean_physeq
sample_sums(Bean_physeq)
sorted <- sort(sample_sums(Bean_physeq))
sorted <- as.data.frame(sorted)
sorted
#Actually do the rarefying here.
set.seed(1)
Bean.rfy <- rarefy_even_depth(Bean_physeq, sample.size = 37815, replace = FALSE)
#52 OTUs were removed: no longer present in any sample after subsampling.





# EXPLORE OCCURRENCE OF RNA1DNA0 AND INFLUENCE OF SUBSAMPLING  --------

#function
Corn_RNA1DNA0_calc<-function(rarefy_level){
  #First rarefy.
  set.seed(1)
  Corn_OTU.rfy <- rarefy_even_depth(Corn_physeq, sample.size = rarefy_level, replace = FALSE)
  Corn_OTU.rfy <- otu_table(Corn_OTU.rfy)
  Corn_OTU.rfy <- as.data.frame(Corn_OTU.rfy)
  #Make OTU table containing only 16sRNA
  Corn_OTUsRNA <- Corn_OTU.rfy[,grep('RNA',names(Corn_OTU.rfy))]
  #Make OTU table containing only 16sDNA
  Corn_OTUsDNA <- Corn_OTU.rfy[,grep('DNA',names(Corn_OTU.rfy))]
  #REMOVE SAMPLES THAT DO NOT HAVE A MATE IN THE OTHER TABLE.
  names(Corn_OTUsRNA) = gsub(pattern = "RNA", replacement = "", x = names(Corn_OTUsRNA))
  names(Corn_OTUsDNA) = gsub(pattern = "DNA", replacement = "", x = names(Corn_OTUsDNA))
  cols_to_keep <- intersect(colnames(Corn_OTUsRNA),colnames(Corn_OTUsDNA))
  Corn_OTUsRNA <- Corn_OTUsRNA[,cols_to_keep, drop=FALSE]
  Corn_OTUsDNA <- Corn_OTUsDNA[,cols_to_keep, drop=FALSE]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Corn_totalOTUs <- data.frame(Corn_OTUsRNA+Corn_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Corn_totalOTUs <- colSums(Corn_totalOTUs !=0)
    #in RNA table, replace everything that is not a 1 with 0.
  #This will retain RNA=1 as 1, and everything else is 0.
  Corn_OTUsRNA[Corn_OTUsRNA<1] = 0
  Corn_OTUsRNA[Corn_OTUsRNA>1] = 0
    #In DNA table, replace all 0's with 1, and everything else with 0.
  #First, make things that are currently 1 into a higher number.
  Corn_OTUsDNA[Corn_OTUsDNA>=1] = 2
  #then, make everything that's currently 0 into a 1.
  Corn_OTUsDNA[Corn_OTUsDNA<1] = 1
  #then, make everything greater than 1 into 0.
  Corn_OTUsDNA[Corn_OTUsDNA>1] = 0
  #sum the two tables
  Corn_RNA1DNA0 <- data.frame(Corn_OTUsRNA+Corn_OTUsDNA)
  #Count the number of 2's for each sample (RNA=1 coded as 1 and DNA=0 coded as 1)
  Corn_RNA1DNA0 <- colSums(Corn_RNA1DNA0 ==2)
  #calculate percent for each sample
  PercentRNA1DNA0 <- Corn_RNA1DNA0/Corn_totalOTUs*100
  PercentRNA1DNA0 <- as.data.frame(PercentRNA1DNA0)
  #make new vector with rarefy level
  subsamp_level<- rep(rarefy_level, nrow(PercentRNA1DNA0))
  #make new vector with sample names
  sample<-row.names(PercentRNA1DNA0)
  as.data.frame(cbind(PercentRNA1DNA0,subsamp_level,sample))
  }

Corn_RNA1DNA0_5000<-Corn_RNA1DNA0_calc(5000)   #5110 OTUs removed
Corn_RNA1DNA0_10000<-Corn_RNA1DNA0_calc(10000) #3573 OTUs removed
Corn_RNA1DNA0_15000<-Corn_RNA1DNA0_calc(15000) #2648 OTUs removed
Corn_RNA1DNA0_20000<-Corn_RNA1DNA0_calc(20000) #2070 OTUs removed
Corn_RNA1DNA0_25000<-Corn_RNA1DNA0_calc(25000) #1652 OTUs removed
Corn_RNA1DNA0_30000<-Corn_RNA1DNA0_calc(30000) #1404 OTUs removed
Corn_RNA1DNA0_35000<-Corn_RNA1DNA0_calc(35000) #1096 OTUs removed
Corn_RNA1DNA0_40000<-Corn_RNA1DNA0_calc(40000) #909 OTUs removed
Corn_RNA1DNA0_45000<-Corn_RNA1DNA0_calc(45000) #755 OTUs removed
Corn_RNA1DNA0_50000<-Corn_RNA1DNA0_calc(50000) #612 OTUs removed
Corn_RNA1DNA0_55000<-Corn_RNA1DNA0_calc(55000) #492 OTUs removed
Corn_RNA1DNA0_60000<-Corn_RNA1DNA0_calc(60000) #406 OTUs removed
Corn_RNA1DNA0_65000<-Corn_RNA1DNA0_calc(65000) #315 OTUs removed
Corn_RNA1DNA0_70000<-Corn_RNA1DNA0_calc(70000) #302 OTUs removed
Corn_RNA1DNA0_75000<-Corn_RNA1DNA0_calc(75000) #259 OTUs removed
Corn_RNA1DNA0_80000<-Corn_RNA1DNA0_calc(80000) #230 OTUs removed
Corn_RNA1DNA0_85000<-Corn_RNA1DNA0_calc(85000) #265 OTUs removed
Corn_RNA1DNA0_90000<-Corn_RNA1DNA0_calc(90000) #290 OTUs removed
Corn_RNA1DNA0_95000<-Corn_RNA1DNA0_calc(95000) #269 OTUs removed
Corn_RNA1DNA0_100000<-Corn_RNA1DNA0_calc(100000) #303 OTUs removed

Corn_RNA1DNA0 <- rbind(Corn_RNA1DNA0_5000,Corn_RNA1DNA0_10000,Corn_RNA1DNA0_15000,Corn_RNA1DNA0_20000,Corn_RNA1DNA0_25000,Corn_RNA1DNA0_30000,Corn_RNA1DNA0_35000,Corn_RNA1DNA0_40000,Corn_RNA1DNA0_45000,Corn_RNA1DNA0_50000,Corn_RNA1DNA0_55000,Corn_RNA1DNA0_60000,Corn_RNA1DNA0_65000,Corn_RNA1DNA0_70000,Corn_RNA1DNA0_75000,Corn_RNA1DNA0_80000,Corn_RNA1DNA0_85000,Corn_RNA1DNA0_90000,Corn_RNA1DNA0_95000,Corn_RNA1DNA0_100000)
treatmentnames=(as(sample_data(Corn_METADATA),"matrix"))
treatmentnames<-as.data.frame(treatmentnames)
treatmentnames <- subset(treatmentnames, Type=='RNA')
row.names(treatmentnames) = gsub(pattern = "RNA", replacement = "", x = row.names(treatmentnames))
treatmentnames<-subset(treatmentnames,select=-c(Type))
Corn_RNA1DNA0 <- merge(Corn_RNA1DNA0,treatmentnames,by.x='sample',by.y='row.names')

#function
Bean_RNA1DNA0_calc<-function(rarefy_level){
  #First rarefy.
  set.seed(1)
  Bean_OTU.rfy <- rarefy_even_depth(Bean_physeq, sample.size = rarefy_level, replace = FALSE)
  Bean_OTU.rfy <- otu_table(Bean_OTU.rfy)
  Bean_OTU.rfy <- as.data.frame(Bean_OTU.rfy)
  names(Bean_OTU.rfy) = gsub(pattern = "cDNA", replacement = "RNA", x = names(Bean_OTU.rfy))
  #Make OTU table containing only 16sRNA
  Bean_OTUsRNA <- Bean_OTU.rfy[,grep('RNA',names(Bean_OTU.rfy))]
  #Make OTU table containing only 16sDNA
  Bean_OTUsDNA <- Bean_OTU.rfy[,grep('DNA',names(Bean_OTU.rfy))]
  #REMOVE SAMPLES THAT DO NOT HAVE A MATE IN THE OTHER TABLE.
  names(Bean_OTUsRNA) = gsub(pattern = "RNA", replacement = "", x = names(Bean_OTUsRNA))
  names(Bean_OTUsDNA) = gsub(pattern = "DNA", replacement = "", x = names(Bean_OTUsDNA))
  cols_to_keep <- intersect(colnames(Bean_OTUsRNA),colnames(Bean_OTUsDNA))
  Bean_OTUsRNA <- Bean_OTUsRNA[,cols_to_keep, drop=FALSE]
  Bean_OTUsDNA <- Bean_OTUsDNA[,cols_to_keep, drop=FALSE]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Bean_totalOTUs <- data.frame(Bean_OTUsRNA+Bean_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Bean_totalOTUs <- colSums(Bean_totalOTUs !=0)
  #in RNA table, replace everything that is not a 1 with 0.
  #This will retain RNA=1 as 1, and everything else is 0.
  Bean_OTUsRNA[Bean_OTUsRNA<1] = 0
  Bean_OTUsRNA[Bean_OTUsRNA>1] = 0
  #In DNA table, replace all 0's with 1, and everything else with 0.
  #First, make things that are currently 1 into a higher number.
  Bean_OTUsDNA[Bean_OTUsDNA>=1] = 2
  #then, make everything that's currently 0 into a 1.
  Bean_OTUsDNA[Bean_OTUsDNA<1] = 1
  #then, make everything greater than 1 into 0.
  Bean_OTUsDNA[Bean_OTUsDNA>1] = 0
  #sum the two tables
  Bean_RNA1DNA0 <- data.frame(Bean_OTUsRNA+Bean_OTUsDNA)
  #Count the number of 2's for each sample (RNA=1 coded as 1 and DNA=0 coded as 1)
  Bean_RNA1DNA0 <- colSums(Bean_RNA1DNA0 ==2)
  #calculate percent for each sample
  PercentRNA1DNA0 <- Bean_RNA1DNA0/Bean_totalOTUs*100
  PercentRNA1DNA0 <- as.data.frame(PercentRNA1DNA0)
  #make new vector with rarefy level
  subsamp_level<- rep(rarefy_level, nrow(PercentRNA1DNA0))
  #make new vector with sample names
  sample<-row.names(PercentRNA1DNA0)
  as.data.frame(cbind(PercentRNA1DNA0,subsamp_level,sample))
}

Bean_RNA1DNA0_5000<-Bean_RNA1DNA0_calc(5000)   #2219 OTUs removed
Bean_RNA1DNA0_10000<-Bean_RNA1DNA0_calc(10000) #1264 OTUs removed
Bean_RNA1DNA0_15000<-Bean_RNA1DNA0_calc(15000) #770 OTUs removed
Bean_RNA1DNA0_20000<-Bean_RNA1DNA0_calc(20000) #451 OTUs removed
Bean_RNA1DNA0_25000<-Bean_RNA1DNA0_calc(25000) #248 OTUs removed
Bean_RNA1DNA0_30000<-Bean_RNA1DNA0_calc(30000) #154 OTUs removed
Bean_RNA1DNA0_35000<-Bean_RNA1DNA0_calc(35000) #68 OTUs removed
Bean_RNA1DNA0_40000<-Bean_RNA1DNA0_calc(40000) #69 OTUs removed. 1 sample removed.
Bean_RNA1DNA0_45000<-Bean_RNA1DNA0_calc(45000) #702 OTUs removed. 16 saamples removed.
Bean_RNA1DNA0_50000<-Bean_RNA1DNA0_calc(50000) #2218 OTUs removed. 37 samples removed.

Bean_RNA1DNA0 <- rbind(Bean_RNA1DNA0_5000,Bean_RNA1DNA0_10000,Bean_RNA1DNA0_15000,Bean_RNA1DNA0_20000,Bean_RNA1DNA0_25000,Bean_RNA1DNA0_30000,Bean_RNA1DNA0_35000,Bean_RNA1DNA0_40000,Bean_RNA1DNA0_45000,Bean_RNA1DNA0_50000)
treatmentnames=(as(sample_data(Bean_METADATA),"matrix"))
treatmentnames<-as.data.frame(treatmentnames)
treatmentnames <- subset(treatmentnames, Type=='RNA')
row.names(treatmentnames) = gsub(pattern = "cDNA", replacement = "", x = row.names(treatmentnames))
treatmentnames<-subset(treatmentnames,select=-c(Type))
Bean_RNA1DNA0 <- merge(Bean_RNA1DNA0,treatmentnames,by.x='sample',by.y='row.names')


# EXPLORE OCCURRENCE OF RNA no DNA AND INFLUENCE OF SUBSAMPLING  --------
#function
Corn_RNAnoDNA_calc<-function(rarefy_level){
  #First rarefy.
  set.seed(1)
  Corn_OTU.rfy <- rarefy_even_depth(Corn_physeq, sample.size = rarefy_level, replace = FALSE)
  Corn_OTU.rfy <- otu_table(Corn_OTU.rfy)
  Corn_OTU.rfy <- as.data.frame(Corn_OTU.rfy)
  #Make OTU table containing only 16sRNA
  Corn_OTUsRNA <- Corn_OTU.rfy[,grep('RNA',names(Corn_OTU.rfy))]
  #Make OTU table containing only 16sDNA
  Corn_OTUsDNA <- Corn_OTU.rfy[,grep('DNA',names(Corn_OTU.rfy))]
  #REMOVE SAMPLES THAT DO NOT HAVE A MATE IN THE OTHER TABLE.
  names(Corn_OTUsRNA) = gsub(pattern = "RNA", replacement = "", x = names(Corn_OTUsRNA))
  names(Corn_OTUsDNA) = gsub(pattern = "DNA", replacement = "", x = names(Corn_OTUsDNA))
  cols_to_keep <- intersect(colnames(Corn_OTUsRNA),colnames(Corn_OTUsDNA))
  Corn_OTUsRNA <- Corn_OTUsRNA[,cols_to_keep, drop=FALSE]
  Corn_OTUsDNA <- Corn_OTUsDNA[,cols_to_keep, drop=FALSE]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Corn_totalOTUs <- data.frame(Corn_OTUsRNA+Corn_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Corn_totalOTUs <- colSums(Corn_totalOTUs !=0)
  #in RNA table, replace RNA>=1 as 1, and keep all zeros as zero.
  Corn_OTUsRNA[Corn_OTUsRNA>=1] = 1
  Corn_OTUsRNA[Corn_OTUsRNA<1] = 0
    #In DNA table, replace all 0's with 1, and everything else with 0.
  #First, make things that are currently >=1 into a higher number.
  Corn_OTUsDNA[Corn_OTUsDNA>=1] = 2
  #then, make everything that's currently 0 into a 1.
  Corn_OTUsDNA[Corn_OTUsDNA<1] = 1
  #then, make everything greater than 1 into 0.
  Corn_OTUsDNA[Corn_OTUsDNA>1] = 0
  #sum the two tables
  Corn_RNAnoDNA <- data.frame(Corn_OTUsRNA+Corn_OTUsDNA)
  #Count the number of 2's for each sample (RNA=1 coded as 1 and DNA=0 coded as 1)
  Corn_RNAnoDNA <- colSums(Corn_RNAnoDNA ==2)
  #calculate percent for each sample
  PercentRNAnoDNA <- Corn_RNAnoDNA/Corn_totalOTUs*100
  PercentRNAnoDNA <- as.data.frame(PercentRNAnoDNA)
  #make new vector with rarefy level
  subsamp_level<- rep(rarefy_level, nrow(PercentRNAnoDNA))
  #make new vector with sample names
  sample<-row.names(PercentRNAnoDNA)
  as.data.frame(cbind(PercentRNAnoDNA,subsamp_level,sample))
}

Corn_RNAnoDNA_5000<-Corn_RNAnoDNA_calc(5000)   #5110 OTUs removed
Corn_RNAnoDNA_10000<-Corn_RNAnoDNA_calc(10000) #3573 OTUs removed
Corn_RNAnoDNA_15000<-Corn_RNAnoDNA_calc(15000) #2648 OTUs removed
Corn_RNAnoDNA_20000<-Corn_RNAnoDNA_calc(20000) #2070 OTUs removed
Corn_RNAnoDNA_25000<-Corn_RNAnoDNA_calc(25000) #1652 OTUs removed. 2 samples removed.
Corn_RNAnoDNA_30000<-Corn_RNAnoDNA_calc(30000) #1404 OTUs removed. 2 samples removed.
Corn_RNAnoDNA_35000<-Corn_RNAnoDNA_calc(35000) #1096 OTUs removed 4 samples removed.
Corn_RNAnoDNA_40000<-Corn_RNAnoDNA_calc(40000) #909 OTUs removed 6 samples removed.
Corn_RNAnoDNA_45000<-Corn_RNAnoDNA_calc(45000) #755 OTUs removed. 8 samples removed.
Corn_RNAnoDNA_50000<-Corn_RNAnoDNA_calc(50000) #612 OTUs removed. 10 samples removed.
Corn_RNAnoDNA_55000<-Corn_RNAnoDNA_calc(55000) #492 OTUs removed. 10 samples removed.
Corn_RNAnoDNA_60000<-Corn_RNAnoDNA_calc(60000) #406 OTUs removed. 11 samples removed.
Corn_RNAnoDNA_65000<-Corn_RNAnoDNA_calc(65000) #315 OTUs removed. 14 samples removed.
Corn_RNAnoDNA_70000<-Corn_RNAnoDNA_calc(70000) #302 OTUs removed. 17 samples removed.
Corn_RNAnoDNA_75000<-Corn_RNAnoDNA_calc(75000) #259 OTUs removed. 22 samples removed.
Corn_RNAnoDNA_80000<-Corn_RNAnoDNA_calc(80000) #230 OTUs removed. 22 samples removed.
Corn_RNAnoDNA_85000<-Corn_RNAnoDNA_calc(85000) #265 OTUs removed. 28 samples removed.
Corn_RNAnoDNA_90000<-Corn_RNAnoDNA_calc(90000) #290 OTUs removed. 29 samples removed.
Corn_RNAnoDNA_95000<-Corn_RNAnoDNA_calc(95000) #269 OTUs removed. 31 samples removed.
Corn_RNAnoDNA_100000<-Corn_RNAnoDNA_calc(100000) #303 OTUs removed. 33 samples removed.

Corn_RNAnoDNA <- rbind(Corn_RNAnoDNA_5000,Corn_RNAnoDNA_10000,Corn_RNAnoDNA_15000,Corn_RNAnoDNA_20000,Corn_RNAnoDNA_25000,Corn_RNAnoDNA_30000,Corn_RNAnoDNA_35000,Corn_RNAnoDNA_40000,Corn_RNAnoDNA_45000,Corn_RNAnoDNA_50000,Corn_RNAnoDNA_55000,Corn_RNAnoDNA_60000,Corn_RNAnoDNA_65000,Corn_RNAnoDNA_70000,Corn_RNAnoDNA_75000,Corn_RNAnoDNA_80000,Corn_RNAnoDNA_85000,Corn_RNAnoDNA_90000,Corn_RNAnoDNA_95000,Corn_RNAnoDNA_100000)
treatmentnames=(as(sample_data(Corn_METADATA),"matrix"))
treatmentnames<-as.data.frame(treatmentnames)
treatmentnames <- subset(treatmentnames, Type=='RNA')
row.names(treatmentnames) = gsub(pattern = "RNA", replacement = "", x = row.names(treatmentnames))
treatmentnames<-subset(treatmentnames,select=-c(Type))
Corn_RNAnoDNA <- merge(Corn_RNAnoDNA,treatmentnames,by.x='sample',by.y='row.names')



#function
Bean_RNAnoDNA_calc<-function(rarefy_level){
  #First rarefy.
  set.seed(1)
  Bean_OTU.rfy <- rarefy_even_depth(Bean_physeq, sample.size = rarefy_level, replace = FALSE)
  Bean_OTU.rfy <- otu_table(Bean_OTU.rfy)
  Bean_OTU.rfy <- as.data.frame(Bean_OTU.rfy)
  names(Bean_OTU.rfy) = gsub(pattern = "cDNA", replacement = "RNA", x = names(Bean_OTU.rfy))
  #Make OTU table containing only 16sRNA
  Bean_OTUsRNA <- Bean_OTU.rfy[,grep('RNA',names(Bean_OTU.rfy))]
  #Make OTU table containing only 16sDNA
  Bean_OTUsDNA <- Bean_OTU.rfy[,grep('DNA',names(Bean_OTU.rfy))]
  #REMOVE SAMPLES THAT DO NOT HAVE A MATE IN THE OTHER TABLE.
  names(Bean_OTUsRNA) = gsub(pattern = "RNA", replacement = "", x = names(Bean_OTUsRNA))
  names(Bean_OTUsDNA) = gsub(pattern = "DNA", replacement = "", x = names(Bean_OTUsDNA))
  cols_to_keep <- intersect(colnames(Bean_OTUsRNA),colnames(Bean_OTUsDNA))
  Bean_OTUsRNA <- Bean_OTUsRNA[,cols_to_keep, drop=FALSE]
  Bean_OTUsDNA <- Bean_OTUsDNA[,cols_to_keep, drop=FALSE]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Bean_totalOTUs <- data.frame(Bean_OTUsRNA+Bean_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Bean_totalOTUs <- colSums(Bean_totalOTUs !=0)
  #in RNA table, replace RNA>=1 as 1, and keep all zeros as zero.
  Bean_OTUsRNA[Bean_OTUsRNA>=1] = 1
  Bean_OTUsRNA[Bean_OTUsRNA<1] = 0
  #In DNA table, replace all 0's with 1, and everything else with 0.
  #First, make things that are currently >=1 into a higher number.
  Bean_OTUsDNA[Bean_OTUsDNA>=1] = 2
  #then, make everything that's currently 0 into a 1.
  Bean_OTUsDNA[Bean_OTUsDNA<1] = 1
  #then, make everything greater than 1 into 0.
  Bean_OTUsDNA[Bean_OTUsDNA>1] = 0
  #sum the two tables
  Bean_RNAnoDNA <- data.frame(Bean_OTUsRNA+Bean_OTUsDNA)
  #Count the number of 2's for each sample (RNA=1 coded as 1 and DNA=0 coded as 1)
  Bean_RNAnoDNA <- colSums(Bean_RNAnoDNA ==2)
  #calculate percent for each sample
  PercentRNAnoDNA <- Bean_RNAnoDNA/Bean_totalOTUs*100
  PercentRNAnoDNA <- as.data.frame(PercentRNAnoDNA)
  #make new vector with rarefy level
  subsamp_level<- rep(rarefy_level, nrow(PercentRNAnoDNA))
  #make new vector with sample names
  sample<-row.names(PercentRNAnoDNA)
  as.data.frame(cbind(PercentRNAnoDNA,subsamp_level,sample))
}

Bean_RNAnoDNA_5000<-Bean_RNAnoDNA_calc(5000)   #2219 OTUs removed
Bean_RNAnoDNA_10000<-Bean_RNAnoDNA_calc(10000) #1264 OTUs removed
Bean_RNAnoDNA_15000<-Bean_RNAnoDNA_calc(15000) #770 OTUs removed
Bean_RNAnoDNA_20000<-Bean_RNAnoDNA_calc(20000) #451 OTUs removed
Bean_RNAnoDNA_25000<-Bean_RNAnoDNA_calc(25000) #248 OTUs removed
Bean_RNAnoDNA_30000<-Bean_RNAnoDNA_calc(30000) #154 OTUs removed
Bean_RNAnoDNA_35000<-Bean_RNAnoDNA_calc(35000) #68 OTUs removed
Bean_RNAnoDNA_40000<-Bean_RNAnoDNA_calc(40000) #69 OTUs removed. 1 sample removed.
Bean_RNAnoDNA_45000<-Bean_RNAnoDNA_calc(45000) #702 OTUs removed. 16 saamples removed.
Bean_RNAnoDNA_50000<-Bean_RNAnoDNA_calc(50000) #2218 OTUs removed. 37 samples removed.

Bean_RNAnoDNA <- rbind(Bean_RNAnoDNA_5000,Bean_RNAnoDNA_10000,Bean_RNAnoDNA_15000,Bean_RNAnoDNA_20000,Bean_RNAnoDNA_25000,Bean_RNAnoDNA_30000,Bean_RNAnoDNA_35000,Bean_RNAnoDNA_40000,Bean_RNAnoDNA_45000,Bean_RNAnoDNA_50000)
treatmentnames=(as(sample_data(Bean_METADATA),"matrix"))
treatmentnames<-as.data.frame(treatmentnames)
treatmentnames <- subset(treatmentnames, Type=='RNA')
row.names(treatmentnames) = gsub(pattern = "cDNA", replacement = "", x = row.names(treatmentnames))
treatmentnames<-subset(treatmentnames,select=-c(Type))
Bean_RNAnoDNA <- merge(Bean_RNAnoDNA,treatmentnames,by.x='sample',by.y='row.names')




# PREPARE SUPPLEMENTARY FIGURES 1 and 2 --------

#Plot all 4 graphs together
Corn_RNAnoDNA_Fig<- ggplot(Corn_RNAnoDNA, aes(subsamp_level, PercentRNAnoDNA, colour=Treatment))+
  geom_smooth()+
  #Add "geom_point()" below for Supp Figure 2.
  #geom_point()+
  ylab("OTUs with RNA > 0 and DNA = 0 (%)")+
  coord_cartesian(ylim=c(0,65))+
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())+
  ggtitle("(A)")+
  xlab("Level of Subsampling")+
  scale_color_manual(limits=c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"),values=c("chocolate4","magenta4","hotpink","blue","red","yellow","orange","green4"))
Corn_RNAnoDNA_Fig

Corn_RNA1DNA0_Fig<- ggplot(Corn_RNA1DNA0, aes(subsamp_level, PercentRNA1DNA0, colour=Treatment))+
  geom_smooth()+
  #Add "geom_point()" below for Supp Figure 2.
  #geom_point()+
  ylab("OTUs with RNA = 1 and DNA = 0 (%)")+
  coord_cartesian(ylim=c(0,65))+
  ggtitle("(C)")+
  theme_bw()+
  theme(legend.position="bottom",legend.title=element_text(size=8))+
  theme(legend.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 2,byrow=TRUE))+
  #ggtitle("Corn: Percent of Total OTUs with RNA1, DNA0")+
  xlab("Level of Subsampling")+
  scale_color_manual(limits=c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"),values=c("chocolate4","magenta4","hotpink","blue","red","yellow","orange","green4"))
Corn_RNA1DNA0_Fig

Bean_RNAnoDNA_Fig<- ggplot(Bean_RNAnoDNA, aes(subsamp_level, PercentRNAnoDNA, colour=Treatment))+
  geom_smooth()+
  #Add "geom_point()" below for Supp Figure 2.
  #geom_point()+
  coord_cartesian(ylim=c(0,65))+
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+
  ggtitle("(B)")+
  #xlab("Level of Subsampling")
  scale_color_manual(limits=c("Control", "Drought", "Nutrients"),values=c("blue","red","green4"))
Bean_RNAnoDNA_Fig

Bean_RNA1DNA0_Fig<- ggplot(Bean_RNA1DNA0, aes(subsamp_level, PercentRNA1DNA0, colour=Treatment))+
  geom_smooth()+
  #Add "geom_point()" below for Supp Figure 2.
  #geom_point()+
  ylab("OTUs with RNA=1, DNA=0 (%)")+
  coord_cartesian(ylim=c(0,65))+
  ggtitle("(D)")+
  theme_bw()+
  theme(legend.position="bottom",axis.title.y=element_blank(),axis.text.y=element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  #ggtitle("Bean: Percent of Total OTUs with RNA1, DNA0")+
  xlab("Level of Subsampling")+
  scale_color_manual(limits=c("Control", "Drought", "Nutrients"),values=c("blue","red","green4"))
Bean_RNA1DNA0_Fig

library("gridExtra")
SuppFigure1 <- grid.arrange(Corn_RNAnoDNA_Fig,Bean_RNAnoDNA_Fig,Corn_RNA1DNA0_Fig,Bean_RNA1DNA0_Fig, ncol = 2,nrow=2, widths=c(2,1.8),heights = c(2, 3))
SuppFigure1
ggsave("SuppFig1.pdf",SuppFigure1,width=174,height=174, units="mm")

#Before attempting to make SuppFig2, be sure to add 'geom_point' in the ggplots above.
library("gridExtra")
SuppFigure2 <- grid.arrange(Corn_RNAnoDNA_Fig,Bean_RNAnoDNA_Fig,Corn_RNA1DNA0_Fig,Bean_RNA1DNA0_Fig, ncol = 2,nrow=2, widths=c(2,1.8),heights = c(2, 3))
SuppFigure2
ggsave("SuppFig2.pdf",SuppFigure2,width=174,height=174, units="mm")



# CORN:  COMPARING 4 METHODS FOR DEALING WITH PHANTOMS --------------------

#function
Corn_AB_method_calc<-function(threshold){
  #Make OTU table all the OTUs in the DNA&RNA dataset.
  Corn_OTUs <- otu_table(Corn.rfy)
  Corn_OTUs <- as.data.frame(Corn_OTUs)
  #Make OTU table containing only 16sRNA
  Corn_OTUsRNA <- Corn_OTUs[,grep('RNA',names(Corn_OTUs))]
  #Make OTU table containing only 16sDNA
  Corn_OTUsDNA <- Corn_OTUs[,grep('DNA',names(Corn_OTUs))]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Corn_totalOTUs <- data.frame(Corn_OTUsRNA+Corn_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Corn_totalOTUs <- colSums(Corn_totalOTUs !=0)
  #Make table of 16s ratios (RNA/DNA) and rename column headers.
  Corn_OTUs_16sRatio <- data.frame(Corn_OTUsRNA/Corn_OTUsDNA)
  #Correct the values in the table as follows:
  #1. When RNA=0 and DNA>0, the value is '0'. We want these cases to remain as '0'.
  #2. When RNA=0 and DNA=0, the value is 'NaN'. Must change these to 'NA'.
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  Corn_OTUs_16sRatio[is.nan(Corn_OTUs_16sRatio)] <- NA
  #3. When RNA>0 and DNA=0, the value is 'Inf'. Must change these to 100 (to indicate 'active', regardless of what threshold is chosen as 'active')
  is.infinite.data.frame <- function(x)
    do.call(cbind, lapply(x, is.infinite))
  Corn_OTUs_16sRatio[is.infinite(Corn_OTUs_16sRatio)] <- 100 
  #Change 16s ratios to 1 or 0 based on 16sRatio greater than or equal to a given threshold.
  Corn_OTUs_16sRatio[Corn_OTUs_16sRatio<=threshold] = 0
  Corn_OTUs_16sRatio[Corn_OTUs_16sRatio>threshold] = 1
  #Calculatethe percent of OTUs with 16sRatio greater than or equal to the given threshold.
  #This gives us the Percent of OTUs that are active.
  Corn_PercActive <- colSums(Corn_OTUs_16sRatio, na.rm=TRUE)/Corn_totalOTUs*100
  Corn_PercActive <- as.data.frame(Corn_PercActive)
  #make new vector with method
  method<- rep('Method1', nrow(Corn_PercActive))
  #make new vector with sample names
  sample<-row.names(Corn_PercActive)
  #make new vector with threshold
  threshold_level <- rep(threshold,nrow(Corn_PercActive))
    as.data.frame(cbind(Corn_PercActive,method,sample,threshold_level))
}

Corn_ABmethod_threshold_1<-Corn_AB_method_calc(1)
Corn_ABmethod_threshold_2<-Corn_AB_method_calc(2)
Corn_ABmethod_threshold_5<-Corn_AB_method_calc(5)



#function
Corn_PK_method_calc<-function(threshold){
  #Make OTU table all the OTUs in the DNA&RNA dataset.
  Corn_OTUs <- otu_table(Corn.rfy)
  Corn_OTUs <- as.data.frame(Corn_OTUs)
  #Make OTU table containing only 16sRNA
  Corn_OTUsRNA <- Corn_OTUs[,grep('RNA',names(Corn_OTUs))]
  #Make OTU table containing only 16sDNA
  Corn_OTUsDNA <- Corn_OTUs[,grep('DNA',names(Corn_OTUs))]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Corn_totalOTUs <- data.frame(Corn_OTUsRNA+Corn_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Corn_totalOTUs <- colSums(Corn_totalOTUs !=0)
  #add 1 to each DNA OTU that is currently a zero, so no 0's in denominator
  Corn_OTUsDNA<-replace(Corn_OTUsDNA, Corn_OTUsDNA == 0, 1)
  #Make table of 16s ratios (RNA/DNA).
  Corn_OTUs_16sRatio <- data.frame(Corn_OTUsRNA/Corn_OTUsDNA)
  #Change 16s ratios to 1 or 0 based on 16sRatio greater than or equal to a given threshold.
  Corn_OTUs_16sRatio[Corn_OTUs_16sRatio<=threshold] = 0
  Corn_OTUs_16sRatio[Corn_OTUs_16sRatio>threshold] = 1
  #Calculatethe percent of OTUs with 16sRatio greater than or equal to the given threshold.
  #This gives us the Percent of OTUs that are active.
  Corn_PercActive <- colSums(Corn_OTUs_16sRatio, na.rm=TRUE)/Corn_totalOTUs*100
  Corn_PercActive <- as.data.frame(Corn_PercActive)
  #make new vector with method
  method<- rep('Method2', nrow(Corn_PercActive))
  #make new vector with sample names
  sample<-row.names(Corn_PercActive)
  #make new vector with threshold
  threshold_level <- rep(threshold,nrow(Corn_PercActive))
  as.data.frame(cbind(Corn_PercActive,method,sample,threshold_level))
}

Corn_PKmethod_threshold_1<-Corn_PK_method_calc(1)
Corn_PKmethod_threshold_2<-Corn_PK_method_calc(2)
Corn_PKmethod_threshold_5<-Corn_PK_method_calc(5)


#function
Corn_PK2016_method_calc<-function(threshold){
  #Make OTU table all the OTUs in the DNA&RNA dataset.
  Corn_OTUs <- otu_table(Corn.rfy)
  Corn_OTUs <- as.data.frame(Corn_OTUs)
  #Make OTU table containing only 16sRNA
  Corn_OTUsRNA <- Corn_OTUs[,grep('RNA',names(Corn_OTUs))]
  #Make OTU table containing only 16sDNA
  Corn_OTUsDNA <- Corn_OTUs[,grep('DNA',names(Corn_OTUs))]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Corn_totalOTUs <- data.frame(Corn_OTUsRNA+Corn_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Corn_totalOTUs <- colSums(Corn_totalOTUs !=0)
  #add 1 to each DNA OTU
  Corn_OTUsDNA<-Corn_OTUsDNA+1
  #Make table of 16s ratios (RNA/DNA).
  Corn_OTUs_16sRatio <- data.frame(Corn_OTUsRNA/Corn_OTUsDNA)
  #Change 16s ratios to 1 or 0 based on 16sRatio greater than or equal to a given threshold.
  Corn_OTUs_16sRatio[Corn_OTUs_16sRatio<=threshold] = 0
  Corn_OTUs_16sRatio[Corn_OTUs_16sRatio>threshold] = 1
  #Calculatethe percent of OTUs with 16sRatio greater than or equal to the given threshold.
  #This gives us the Percent of OTUs that are active.
  Corn_PercActive <- colSums(Corn_OTUs_16sRatio, na.rm=TRUE)/Corn_totalOTUs*100
  Corn_PercActive <- as.data.frame(Corn_PercActive)
  #make new vector with method
  method<- rep('Method3', nrow(Corn_PercActive))
  #make new vector with sample names
  sample<-row.names(Corn_PercActive)
  #make new vector with threshold
  threshold_level <- rep(threshold,nrow(Corn_PercActive))
  as.data.frame(cbind(Corn_PercActive,method,sample,threshold_level))
}

Corn_PK2016method_threshold_1<-Corn_PK2016_method_calc(1)
Corn_PK2016method_threshold_2<-Corn_PK2016_method_calc(2)
Corn_PK2016method_threshold_5<-Corn_PK2016_method_calc(5)


#function
Corn_Denef_method_calc<-function(threshold){
  #Make OTU table all the OTUs in the DNA&RNA dataset.
  Corn_OTUs <- otu_table(Corn.rfy)
  Corn_OTUs <- as.data.frame(Corn_OTUs)
  #Make OTU table containing only 16sRNA
  Corn_OTUsRNA <- Corn_OTUs[,grep('RNA',names(Corn_OTUs))]
  #Make OTU table containing only 16sDNA
  Corn_OTUsDNA <- Corn_OTUs[,grep('DNA',names(Corn_OTUs))]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Corn_totalOTUs <- data.frame(Corn_OTUsRNA+Corn_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Corn_totalOTUs <- colSums(Corn_totalOTUs !=0)
  #add 1 to every cell, treating every cell the same.
  Corn_OTUsDNA<-Corn_OTUsDNA+1
  Corn_OTUsRNA<-Corn_OTUsRNA+1
  #Make table of 16s ratios (RNA/DNA) and rename column headers.
  Corn_OTUs_16sRatio <- data.frame(Corn_OTUsRNA/Corn_OTUsDNA)
  #Change 16s ratios to 1 or 0 based on 16sRatio greater than or equal to a given threshold.
  Corn_OTUs_16sRatio[Corn_OTUs_16sRatio<=threshold] = 0
  Corn_OTUs_16sRatio[Corn_OTUs_16sRatio>threshold] = 1
  #Calculatethe percent of OTUs with 16sRatio greater than or equal to the given threshold.
  #This gives us the Percent of OTUs that are active.
  Corn_PercActive <- colSums(Corn_OTUs_16sRatio, na.rm=TRUE)/Corn_totalOTUs*100
  Corn_PercActive <- as.data.frame(Corn_PercActive)
  #make new vector with method
  method<- rep('Method4', nrow(Corn_PercActive))
  #make new vector with sample names
  sample<-row.names(Corn_PercActive)
  #make new vector with threshold
  threshold_level <- rep(threshold,nrow(Corn_PercActive))
  as.data.frame(cbind(Corn_PercActive,method,sample,threshold_level))
}

Corn_Denefmethod_threshold_1<-Corn_Denef_method_calc(1)
Corn_Denefmethod_threshold_2<-Corn_Denef_method_calc(2)
Corn_Denefmethod_threshold_5<-Corn_Denef_method_calc(5)



#Comparing the 4 methods:

Corn_PercActive_compare <-rbind(Corn_ABmethod_threshold_1,Corn_ABmethod_threshold_2,Corn_ABmethod_threshold_5,Corn_PKmethod_threshold_1,Corn_PKmethod_threshold_2,Corn_PKmethod_threshold_5,Corn_PK2016method_threshold_1,Corn_PK2016method_threshold_2,Corn_PK2016method_threshold_5,Corn_Denefmethod_threshold_1,Corn_Denefmethod_threshold_2,Corn_Denefmethod_threshold_5)

treatmentnames=(as(sample_data(Corn_METADATA),"matrix"))
treatmentnames<-as.data.frame(treatmentnames)
treatmentnames <- subset(treatmentnames, Type=='RNA')
treatmentnames<-subset(treatmentnames,select=-c(Type))
Corn_PercActive_compare <- merge(Corn_PercActive_compare,treatmentnames,by.x='sample',by.y='row.names')

labels <- c('1'="Threshold=1", '2' = "Threshold=2",'5'="Threshold=5")

#Reorder treatment 'factors' for the figure (so they are no longer alphabetical, but are in the order I want)
Corn_PercActive_compare$Treatment = factor(Corn_PercActive_compare$Treatment,c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"))



# BEAN: COMPARE 4 Methods FOR DEALING WITH PHANTOMS -----------------------------------------------------------

#function
Bean_AB_method_calc<-function(threshold){
  #Make OTU table all the OTUs in the DNA&RNA dataset.
  Bean_OTUs <- otu_table(Bean.rfy)
  Bean_OTUs <- as.data.frame(Bean_OTUs)
  names(Bean_OTUs) = gsub(pattern = "cDNA", replacement = "RNA", x = names(Bean_OTUs))
  #Make OTU table containing only 16sRNA
  Bean_OTUsRNA <- Bean_OTUs[,grep('RNA',names(Bean_OTUs))]
  #Make OTU table containing only 16sDNA
  Bean_OTUsDNA <- Bean_OTUs[,grep('DNA',names(Bean_OTUs))]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Bean_totalOTUs <- data.frame(Bean_OTUsRNA+Bean_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Bean_totalOTUs <- colSums(Bean_totalOTUs !=0)
  #Make table of 16s ratios (RNA/DNA) and rename column headers.
  Bean_OTUs_16sRatio <- data.frame(Bean_OTUsRNA/Bean_OTUsDNA)
  #Correct the values in the table as follows:
  #1. When RNA=0 and DNA>0, the value is '0'. We want these cases to remain as '0'.
  #2. When RNA=0 and DNA=0, the value is 'NaN'. Must change these to 'NA'.
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  Bean_OTUs_16sRatio[is.nan(Bean_OTUs_16sRatio)] <- NA
  #3. When RNA>0 and DNA=0, the value is 'Inf'. Must change these to 100 (to indicate 'active', regardless of what threshold is chosen as 'active')
  is.infinite.data.frame <- function(x)
    do.call(cbind, lapply(x, is.infinite))
  Bean_OTUs_16sRatio[is.infinite(Bean_OTUs_16sRatio)] <- 100 
  #Change 16s ratios to 1 or 0 based on 16sRatio greater than or equal to a given threshold.
  Bean_OTUs_16sRatio[Bean_OTUs_16sRatio<=threshold] = 0
  Bean_OTUs_16sRatio[Bean_OTUs_16sRatio>threshold] = 1
  #Calculatethe percent of OTUs with 16sRatio greater than or equal to the given threshold.
  #This gives us the Percent of OTUs that are active.
  Bean_PercActive <- colSums(Bean_OTUs_16sRatio, na.rm=TRUE)/Bean_totalOTUs*100
  Bean_PercActive <- as.data.frame(Bean_PercActive)
  #make new vector with method
  method<- rep('Method1', nrow(Bean_PercActive))
  #make new vector with sample names
  sample<-row.names(Bean_PercActive)
  #make new vector with threshold
  threshold_level <- rep(threshold,nrow(Bean_PercActive))
  as.data.frame(cbind(Bean_PercActive,method,sample,threshold_level))
}

Bean_ABmethod_threshold_1<-Bean_AB_method_calc(1)
Bean_ABmethod_threshold_2<-Bean_AB_method_calc(2)
Bean_ABmethod_threshold_5<-Bean_AB_method_calc(5)


#function
Bean_PK_method_calc<-function(threshold){
  #Make OTU table all the OTUs in the DNA&RNA dataset.
  Bean_OTUs <- otu_table(Bean.rfy)
  Bean_OTUs <- as.data.frame(Bean_OTUs)
  names(Bean_OTUs) = gsub(pattern = "cDNA", replacement = "RNA", x = names(Bean_OTUs))
  #Make OTU table containing only 16sRNA
  Bean_OTUsRNA <- Bean_OTUs[,grep('RNA',names(Bean_OTUs))]
  #Make OTU table containing only 16sDNA
  Bean_OTUsDNA <- Bean_OTUs[,grep('DNA',names(Bean_OTUs))]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Bean_totalOTUs <- data.frame(Bean_OTUsRNA+Bean_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Bean_totalOTUs <- colSums(Bean_totalOTUs !=0)
  #add 1 to each DNA OTU that is currently a zero, so no 0's in denominator
  Bean_OTUsDNA<-replace(Bean_OTUsDNA, Bean_OTUsDNA == 0, 1)
  #Make table of 16s ratios (RNA/DNA) and rename column headers.
  Bean_OTUs_16sRatio <- data.frame(Bean_OTUsRNA/Bean_OTUsDNA)
  #Change 16s ratios to 1 or 0 based on 16sRatio greater than or equal to a given threshold.
  Bean_OTUs_16sRatio[Bean_OTUs_16sRatio<=threshold] = 0
  Bean_OTUs_16sRatio[Bean_OTUs_16sRatio>threshold] = 1
  #Calculatethe percent of OTUs with 16sRatio greater than or equal to the given threshold.
  #This gives us the Percent of OTUs that are active.
  Bean_PercActive <- colSums(Bean_OTUs_16sRatio, na.rm=TRUE)/Bean_totalOTUs*100
  Bean_PercActive <- as.data.frame(Bean_PercActive)
  #make new vector with method
  method<- rep('Method2', nrow(Bean_PercActive))
  #make new vector with sample names
  sample<-row.names(Bean_PercActive)
  #make new vector with threshold
  threshold_level <- rep(threshold,nrow(Bean_PercActive))
  as.data.frame(cbind(Bean_PercActive,method,sample,threshold_level))
}

Bean_PKmethod_threshold_1<-Bean_PK_method_calc(1)
Bean_PKmethod_threshold_2<-Bean_PK_method_calc(2)
Bean_PKmethod_threshold_5<-Bean_PK_method_calc(5)


#function
Bean_PK2016_method_calc<-function(threshold){
  #Make OTU table all the OTUs in the DNA&RNA dataset.
  Bean_OTUs <- otu_table(Bean.rfy)
  Bean_OTUs <- as.data.frame(Bean_OTUs)
  names(Bean_OTUs) = gsub(pattern = "cDNA", replacement = "RNA", x = names(Bean_OTUs))
  #Make OTU table containing only 16sRNA
  Bean_OTUsRNA <- Bean_OTUs[,grep('RNA',names(Bean_OTUs))]
  #Make OTU table containing only 16sDNA
  Bean_OTUsDNA <- Bean_OTUs[,grep('DNA',names(Bean_OTUs))]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Bean_totalOTUs <- data.frame(Bean_OTUsRNA+Bean_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Bean_totalOTUs <- colSums(Bean_totalOTUs !=0)
  #add 1 to each DNA OTU
  Bean_OTUsDNA<-Bean_OTUsDNA+1
  #Make table of 16s ratios (RNA/DNA) and rename column headers.
  Bean_OTUs_16sRatio <- data.frame(Bean_OTUsRNA/Bean_OTUsDNA)
  #Change 16s ratios to 1 or 0 based on 16sRatio greater than or equal to a given threshold.
  Bean_OTUs_16sRatio[Bean_OTUs_16sRatio<=threshold] = 0
  Bean_OTUs_16sRatio[Bean_OTUs_16sRatio>threshold] = 1
  #Calculatethe percent of OTUs with 16sRatio greater than or equal to the given threshold.
  #This gives us the Percent of OTUs that are active.
  Bean_PercActive <- colSums(Bean_OTUs_16sRatio, na.rm=TRUE)/Bean_totalOTUs*100
  Bean_PercActive <- as.data.frame(Bean_PercActive)
  #make new vector with method
  method<- rep('Method3', nrow(Bean_PercActive))
  #make new vector with sample names
  sample<-row.names(Bean_PercActive)
  #make new vector with threshold
  threshold_level <- rep(threshold,nrow(Bean_PercActive))
  as.data.frame(cbind(Bean_PercActive,method,sample,threshold_level))
}

Bean_PK2016method_threshold_1<-Bean_PK2016_method_calc(1)
Bean_PK2016method_threshold_2<-Bean_PK2016_method_calc(2)
Bean_PK2016method_threshold_5<-Bean_PK2016_method_calc(5)



#function
Bean_Denef_method_calc<-function(threshold){
  #Make OTU table all the OTUs in the DNA&RNA dataset.
  Bean_OTUs <- otu_table(Bean.rfy)
  Bean_OTUs <- as.data.frame(Bean_OTUs)
  names(Bean_OTUs) = gsub(pattern = "cDNA", replacement = "RNA", x = names(Bean_OTUs))
  #Make OTU table containing only 16sRNA
  Bean_OTUsRNA <- Bean_OTUs[,grep('RNA',names(Bean_OTUs))]
  #Make OTU table containing only 16sDNA
  Bean_OTUsDNA <- Bean_OTUs[,grep('DNA',names(Bean_OTUs))]
  #Count the totalnumber of OTUs present in each sample (i.e. the # of OTUs that have either DNA or RNA present, or both)
  #First, add the RNA and DNA datasets together for each sample.
  Bean_totalOTUs <- data.frame(Bean_OTUsRNA+Bean_OTUsDNA)
  #Lastly, count the number of cells with a non-zero value (in other words, cells that have either RNA reads, or DNA reads, or both)
  Bean_totalOTUs <- colSums(Bean_totalOTUs !=0)
  #add 1 to every cell, treating every cell the same.
  Bean_OTUsDNA<-Bean_OTUsDNA+1
  Bean_OTUsRNA<-Bean_OTUsRNA+1
  #Make table of 16s ratios (RNA/DNA) and rename column headers.
  Bean_OTUs_16sRatio <- data.frame(Bean_OTUsRNA/Bean_OTUsDNA)
  #Change 16s ratios to 1 or 0 based on 16sRatio greater than or equal to a given threshold.
  Bean_OTUs_16sRatio[Bean_OTUs_16sRatio<=threshold] = 0
  Bean_OTUs_16sRatio[Bean_OTUs_16sRatio>threshold] = 1
  #Calculatethe percent of OTUs with 16sRatio greater than or equal to the given threshold.
  #This gives us the Percent of OTUs that are active.
  Bean_PercActive <- colSums(Bean_OTUs_16sRatio, na.rm=TRUE)/Bean_totalOTUs*100
  Bean_PercActive <- as.data.frame(Bean_PercActive)
  #make new vector with method
  method<- rep('Method4', nrow(Bean_PercActive))
  #make new vector with sample names
  sample<-row.names(Bean_PercActive)
  #make new vector with threshold
  threshold_level <- rep(threshold,nrow(Bean_PercActive))
  as.data.frame(cbind(Bean_PercActive,method,sample,threshold_level))
}

Bean_Denefmethod_threshold_1<-Bean_Denef_method_calc(1)
Bean_Denefmethod_threshold_2<-Bean_Denef_method_calc(2)
Bean_Denefmethod_threshold_5<-Bean_Denef_method_calc(5)




#Comparing the 4 methods:
Bean_PercActive_compare <-rbind(Bean_ABmethod_threshold_1,Bean_ABmethod_threshold_2,Bean_ABmethod_threshold_5,Bean_PKmethod_threshold_1,Bean_PKmethod_threshold_2,Bean_PKmethod_threshold_5,Bean_PK2016method_threshold_1,Bean_PK2016method_threshold_2,Bean_PK2016method_threshold_5,Bean_Denefmethod_threshold_1,Bean_Denefmethod_threshold_2,Bean_Denefmethod_threshold_5)
treatmentnames=(as(sample_data(Bean_METADATA),"matrix"))
treatmentnames<-as.data.frame(treatmentnames)
treatmentnames <- subset(treatmentnames, Type=='RNA')
row.names(treatmentnames) = gsub(pattern = "cDNA", replacement = "RNA", x = row.names(treatmentnames))
treatmentnames<-subset(treatmentnames,select=-c(Type))
Bean_PercActive_compare <- merge(Bean_PercActive_compare,treatmentnames,by.x='sample',by.y='row.names')
labels <- c('1'="Threshold=1", '2' = "Threshold=2",'5'="Threshold=5")



# PREPARE FIGURE 2 ----------------------------------------------------------

Corn_fig2 <-ggplot(data=Corn_PercActive_compare, aes(x=method, y=Corn_PercActive, fill=Treatment))+
  #need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
  geom_boxplot(outlier.shape=NA) + 
  scale_size(guide=FALSE)+
  scale_x_discrete(name="Method for Calculating 16S Ratios")+
  scale_y_continuous(name="Percent Active OTUs",limits=c(0,70))+
  theme_bw(base_size=12)+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(size=7))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  #scale_x_discrete(limits=c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"))+
  scale_fill_manual(name="Treatment", breaks=c("Pre-dry", "Post-dry", "Post-water", "Control", "ABA", "IAA", "JA", "SA"), values=c("chocolate4", "magenta4","hotpink", "blue", "red","yellow", "orange", "green4"))+
  ggtitle("(A)")+
  facet_grid(threshold_level ~ ., labeller=labeller(threshold_level=labels))
Corn_fig2

Bean_fig2 <- ggplot(data=Bean_PercActive_compare, aes(x=method, y=Bean_PercActive, fill=Treatment))+
  #need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
  geom_boxplot(outlier.shape=NA) + 
  scale_size(guide=FALSE)+
  scale_x_discrete(name="Method for Calculating 16S Ratios")+
  scale_y_continuous(name="Percent Active",limits=c(0,70))+
  theme_bw(base_size=12)+
  scale_fill_manual(limits=c("Control","Drought","Nutrients"),values=c("blue","red","green4"))+
  theme(legend.position="bottom",axis.title.y=element_blank(),axis.text.y=element_blank())+
  theme(legend.text = element_text(size=7))+
  #scale_color_manual(limits=c("Control", "Drought", "Nutrients"),values=c("blue","red","green4"))+
  ggtitle("(B)")+
  facet_grid(threshold_level ~ ., labeller=labeller(threshold_level=labels))
Bean_fig2



library("grid")
pdf("Figure2.pdf",width=7, height=6,paper = "letter")
grid.newpage()
grid.draw(cbind(ggplotGrob(Corn_fig2), ggplotGrob(Bean_fig2), size="last"))
dev.off()



# PREPARE FIGURE 3 (rrnDB operons) --------------------------------------------------------------

#Need to get rrnDB data for mean 16s operon count for each phylum.
#My first step is to check how similar the NCBI and RDP taxonomy data are in rrnDB.

#First for Corn:

#get list of phylum names in rarefied OTU table
#Corn_taxalist <- tax_table(Corn.rfy)
#get only unique phylum names
#Corn_taxalist <- subset(Corn_taxalist, select=c("Phylum"))
#Corn_taxalist <- unique(Corn_taxalist)
#Corn_taxalist <- as.data.frame(Corn_taxalist)

#download from rrnDB (NCBI) and save as CSV file
#rrnDB_NCBI <- read.csv("./Corn Input Files/rrnDB-5.4_pantaxa_stats_NCBI.csv", sep=",", row.names=1, header=TRUE, na.strings="")
#remove all but rank, name,and mean data
#rrnDB_NCBI <- subset(rrnDB_NCBI, select=c("rank","name","mean"))
#subset to include only level of "phylum"
#rrnDB_NCBI <- rrnDB_NCBI[(rrnDB_NCBI$rank=="phylum"),]
#merge with list of phyla in the Corn dataset.
#newdataset1 <- merge(Corn_taxalist, rrnDB_NCBI, by.x="Phylum", by.y="name")
#rename column header
#colnames(newdataset1)[colnames(newdataset1) == 'mean'] <- 'NCBIMean'

#download from rrnDB (RDP) and save as CSV file
#rrnDB_RDP <- read.csv("./Corn Input Files/rrnDB-5.4_pantaxa_stats_RDP.csv", sep=",", row.names=1, header=TRUE, na.strings="")
#remove all but rank, name,and mean data
#rrnDB_RDP <- subset(rrnDB_RDP, select=c("rank","name","mean"))
#subset to include only level of "phylum"
#rrnDB_RDP <- rrnDB_RDP[(rrnDB_RDP$rank=="phylum"),]
#merge with list of phyla in the Corn dataset.
#newdataset2 <- merge(Corn_taxalist, rrnDB_RDP, by.x="Phylum", by.y="name")
#rename column header
#colnames(newdataset2)[colnames(newdataset2) == 'mean'] <- 'RDPMean'

#newdataset3 <- merge(newdataset1, newdataset2, by="Phylum")

#cor.test(newdataset3$NCBIMean, newdataset3$RDPMean)
#very high correlation. p-value of 2.2*10^-16.
#confident both taxonomies will give similar results.



#Now for Bean:

#get list of phylum names in rarefied OTU table
#Bean_taxalist <- tax_table(Bean.rfy)
#get only unique phylum names
#Bean_taxalist <- subset(Bean_taxalist, select=c("Phylum"))
#Bean_taxalist <- unique(Bean_taxalist)
#Bean_taxalist <- as.data.frame(Bean_taxalist)

#download from rrnDB (NCBI) and save as CSV file
#rrnDB_NCBI <- read.csv("./Corn Input Files/rrnDB-5.4_pantaxa_stats_NCBI.csv", sep=",", row.names=1, header=TRUE, na.strings="")
#remove all but rank, name,and mean data
#rrnDB_NCBI <- subset(rrnDB_NCBI, select=c("rank","name","mean"))
#subset to include only level of "phylum"
#rrnDB_NCBI <- rrnDB_NCBI[(rrnDB_NCBI$rank=="phylum"),]
#merge with list of phyla in the Corn dataset.
#newdataset1 <- merge(Bean_taxalist, rrnDB_NCBI, by.x="Phylum", by.y="name")
#rename column header
#colnames(newdataset1)[colnames(newdataset1) == 'mean'] <- 'NCBIMean'

#download from rrnDB (RDP) and save as CSV file
#rrnDB_RDP <- read.csv("./Corn Input Files/rrnDB-5.4_pantaxa_stats_RDP.csv", sep=",", row.names=1, header=TRUE, na.strings="")
#remove all but rank, name,and mean data
#rrnDB_RDP <- subset(rrnDB_RDP, select=c("rank","name","mean"))
#subset to include only level of "phylum"
#rrnDB_RDP <- rrnDB_RDP[(rrnDB_RDP$rank=="phylum"),]
#merge with list of phyla in the Corn dataset.
#newdataset2 <- merge(Bean_taxalist, rrnDB_RDP, by.x="Phylum", by.y="name")
#rename column header
#colnames(newdataset2)[colnames(newdataset2) == 'mean'] <- 'RDPMean'
#newdataset3 <- merge(newdataset1, newdataset2, by="Phylum")
#cor.test(newdataset3$NCBIMean, newdataset3$RDPMean)
#very high correlation. p-value of 2.2*10^-16.
#confident both taxonomies will give similar results.



###It is clear that the NCBI and RDP taxonomies give extremely highly correlated
###values for mean operon numbers for each phylum.
###I will stick with the RDP taxonomies for this study.






rrnDB_RDP <- read.csv("./Corn Input Files/rrnDB-5.4_pantaxa_stats_RDP.csv", sep=",", row.names=1, header=TRUE, na.strings="")
#remove all but rank, name,and mean data
rrnDB_RDP <- subset(rrnDB_RDP, select=c("rank","name","mean"))
#subset to include only level of "phylum"
rrnDB_RDP <- rrnDB_RDP[(rrnDB_RDP$rank=="phylum"),]

#rename Cyanobacteria/Chloroplast to read "Cyanobacteria"
rrnDB_RDP$name<-gsub('Cyanobacteria/Chloroplast', 'Cyanobacteria', rrnDB_RDP$name)



#merge with 16s ratio data.
#Make OTU table all the OTUs in the DNA&RNA dataset.
Corn_OTUs <- otu_table(Corn.rfy)
Corn_OTUs <- as.data.frame(Corn_OTUs)
#Make OTU table containing only 16sRNA
Corn_OTUsRNA <- Corn_OTUs[,grep('RNA',names(Corn_OTUs))]
#Make OTU table containing only 16sDNA
Corn_OTUsDNA <- Corn_OTUs[,grep('DNA',names(Corn_OTUs))]
#add 1 to each DNA OTU (Kearns et al. 2016 method)
Corn_OTUsDNA<-Corn_OTUsDNA+1
#Make table of 16s ratios (RNA/DNA).
Corn_OTUs_16sRatio <- data.frame(Corn_OTUsRNA/Corn_OTUsDNA)
#Get taxonomy.
Corn_taxonomy <- tax_table(Corn.rfy)
Corn_taxonomy <- subset(Corn_taxonomy, select = c("Phylum"))
#Now build the data.frame holding 16s ratio and taxonomy data.
Corn_OTUs_16sRatio <- merge(Corn_OTUs_16sRatio, Corn_taxonomy, by="row.names")
row.names(Corn_OTUs_16sRatio) <- Corn_OTUs_16sRatio$Row.names
colnames(Corn_OTUs_16sRatio)[colnames(Corn_OTUs_16sRatio) == 'Row.names'] <- 'OTU ID'
#melt data
Corn_OTUs_16sRatio <- melt(Corn_OTUs_16sRatio, id=c("OTU ID", "Phylum"))

Corn_OTUs_16sRatio <- merge(Corn_OTUs_16sRatio, rrnDB_RDP, by.x="Phylum", by.y="name")

#log transform the 16S ratio values for plotting
Corn_OTUs_16sRatio$value<-log(Corn_OTUs_16sRatio$value)

#Remove rows with -Inf
Corn_OTUs_16sRatio<-Corn_OTUs_16sRatio[!(Corn_OTUs_16sRatio$value=="-Inf"),]

pal26<-c("#89C5DA","#DA5724","#ffa500","#74D944","#DA5724","#508578","#CE50CA","#3F4921","#652926","#D7C1B1","#AD6F3B","#689030","#CD9BCD","#C0717C","#CBD588","#D14285","#5F7FC7","#673770","#8A7C64")

#plot it with ggplot
h1 <-ggplot(Corn_OTUs_16sRatio, aes(mean, value, colour=Phylum))+
  geom_point()+
  xlab("Average Copy Number at Phylum level")+
  ylab("Log10 16S ratio")+
  coord_cartesian(ylim=c(-10,10))+
  theme_bw()+
  theme(legend.position="bottom",legend.title=element_text(size=8))+
  theme(legend.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 6,byrow=TRUE))+
  scale_colour_manual(values=pal26)+
  ggtitle("(A)")+
  theme(legend.text = element_text(margin = margin(r = 24, unit = "pt")))
h1

#pearson correlation
cor.test(Corn_OTUs_16sRatio$mean, Corn_OTUs_16sRatio$value)
#r=-0.069, p=2.2*10^-16



#merge with 16s ratio data.
#Make OTU table all the OTUs in the DNA&RNA dataset.
Bean_OTUs <- otu_table(Bean.rfy)
Bean_OTUs <- as.data.frame(Bean_OTUs)
names(Bean_OTUs) = gsub(pattern = "cDNA", replacement = "RNA", x = names(Bean_OTUs))
#Make OTU table containing only 16sRNA
Bean_OTUsRNA <- Bean_OTUs[,grep('RNA',names(Bean_OTUs))]
#Make OTU table containing only 16sDNA
Bean_OTUsDNA <- Bean_OTUs[,grep('DNA',names(Bean_OTUs))]
#add 1 to each DNA OTU
Bean_OTUsDNA<-Bean_OTUsDNA+1
#Make table of 16s ratios (RNA/DNA).
Bean_OTUs_16sRatio <- data.frame(Bean_OTUsRNA/Bean_OTUsDNA)
#Get taxonomy.
Bean_taxonomy <- tax_table(Bean.rfy)
Bean_taxonomy <- subset(Bean_taxonomy, select = c("Phylum"))
#Now build the data.frame holding 16s ratio and taxonomy data.
Bean_OTUs_16sRatio <- merge(Bean_OTUs_16sRatio, Bean_taxonomy, by="row.names")
row.names(Bean_OTUs_16sRatio) <- Bean_OTUs_16sRatio$Row.names
colnames(Bean_OTUs_16sRatio)[colnames(Bean_OTUs_16sRatio) == 'Row.names'] <- 'OTU ID'
#melt data
Bean_OTUs_16sRatio <- melt(Bean_OTUs_16sRatio, id=c("OTU ID", "Phylum"))

Bean_OTUs_16sRatio <- merge(Bean_OTUs_16sRatio, rrnDB_RDP, by.x="Phylum",by.y="name")

#log transform the 16S ratio values for plotting
Bean_OTUs_16sRatio$value<-log(Bean_OTUs_16sRatio$value)

#Remove rows with -Inf
Bean_OTUs_16sRatio<-Bean_OTUs_16sRatio[!(Bean_OTUs_16sRatio$value=="-Inf"),]

pal26<-c("#89C5DA","#DA5724","#ffa500","#74D944","#DA5724","#508578","#CE50CA","#3F4921","#652926","#D7C1B1","#AD6F3B","#689030","#CD9BCD","#C0717C","#CBD588","#D14285","#673770","#8A7C64")


#plot it with ggplot
h2 <-ggplot(Bean_OTUs_16sRatio, aes(mean, value, colour=Phylum))+
  geom_point()+
  xlab("Average Copy Number at Phylum level")+
  coord_cartesian(ylim=c(-10,10))+
  theme_bw()+
  theme(legend.position="bottom",legend.title=element_text(size=8),axis.title.y=element_blank(),axis.text.y=element_blank())+
  theme(legend.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 5,byrow=TRUE))+
  scale_colour_manual(values=pal26)+
  ggtitle("(B)")
h2

#pearson correlation
cor.test(Bean_OTUs_16sRatio$mean, Bean_OTUs_16sRatio$value)
#r=-0.017, p=0.00012




#Plotting the two together
library("grid")
pdf("Figure3.pdf", height = 6, width = 14)
grid.newpage()
grid.draw(cbind(ggplotGrob(h1), ggplotGrob(h2), size="last"))
dev.off()




# PREPARE FIGURE 4 -----------------------------------------------------------


#CORN

Corn_CellStainData <- read.csv("./Corn Input Files/CellCountData.csv", sep=",",header=TRUE, na.strings="")

#Find mean of technical reps for each within each column by the Taxonomy of Interest.
library("dplyr")
Corn_CellStainData <-Corn_CellStainData %>% group_by(Treatment,Bio_Replicate) %>% summarise_all(funs(mean))
Corn_CellStainData <- as.data.frame(Corn_CellStainData)

#Prepare data for figure
treatmentnames=(as(sample_data(Corn_METADATA),"matrix"))
treatmentnames<-as.data.frame(treatmentnames)
treatmentnames <- subset(treatmentnames, Type=='RNA')
treatmentnames<-subset(treatmentnames,select=-c(Type))
Corn_RatioData <- merge(Corn_PK2016method_threshold_1,treatmentnames,by.x='sample',by.y='row.names')

Corn_RatioData<-subset(Corn_RatioData,select=-c(sample,threshold_level,method))

#Reorder treatment 'factors' for the figure (so they are no longer alphabetical, but are in the order I want)
Corn_RatioData$Treatment = factor(Corn_RatioData$Treatment,c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"))

j1<-ggplot(data=Corn_RatioData, aes(x=Treatment, y=Corn_PercActive, fill=Treatment))+
  #need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
  geom_boxplot(outlier.shape=NA) + 
  scale_size(guide=FALSE)+
  scale_x_discrete(name="Treatment")+
  scale_y_continuous(name="16s RNA:DNA Active OTUs (%)",limits=c(0,75))+
  theme_bw(base_size=12)+
  theme(legend.position="none")+
  scale_fill_manual(name="Treatment", breaks=c("Pre-dry", "Post-dry", "Post-water", "Control", "ABA", "IAA", "JA", "SA"), values=c("chocolate4", "magenta4","hotpink", "blue", "red","yellow", "orange", "green4"))+
  ggtitle("(C)")
j1

#Reorder treatment 'factors' for the figure (so they are no longer alphabetical, but are in the order I want)
Corn_CellStainData$Treatment = factor(Corn_CellStainData$Treatment,c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"))

j2<-ggplot(data=Corn_CellStainData, aes(x=Treatment, y=Perc_Active, fill=Treatment))+
  #need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
  geom_boxplot(outlier.shape=NA) + 
  scale_size(guide=FALSE)+
  scale_y_continuous(name="CTC Stain: Active Cells (%)",limits=c(0,75))+
  theme_bw(base_size=12)+
  scale_x_discrete(limits=c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"))+
  #scale_fill_manual(limits=c("CTC-Stain"),values=c("red"))+
  ggtitle("(A)")+
  theme(legend.position="none")+
  scale_fill_manual(name="Treatment", breaks=c("Pre-dry", "Post-dry", "Post-water", "Control", "ABA", "IAA", "JA", "SA"), values=c("chocolate4", "magenta4","hotpink", "blue", "red","yellow", "orange", "green4"))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
j2






#BEAN

Bean_CellStainData <- read.csv("./Bean Input Files/CellCountData.csv", sep=",",header=TRUE, na.strings="")
Bean_CellStainData <- subset(Bean_CellStainData, select = -c(Type, Bio_Replicate,Tech_Replicate))
#Find mean of bio/technical reps within each sample.
library("dplyr")
Bean_CellStainData <-Bean_CellStainData %>% group_by(Sample) %>% summarise_all(funs(mean))
Bean_CellStainData <- as.data.frame(Bean_CellStainData)


#Prepare data for figure
treatmentnames=(as(sample_data(Bean_METADATA),"matrix"))
treatmentnames<-as.data.frame(treatmentnames)
treatmentnames <- subset(treatmentnames, Type=='RNA')
treatmentnames<-subset(treatmentnames,select=-c(Type))
row.names(treatmentnames) = gsub(pattern = "cDNA", replacement = "RNA", x = row.names(treatmentnames))
Bean_RatioData <- merge(Bean_PK2016method_threshold_1,treatmentnames,by='row.names')

Bean_CellStainData$Treatment <- Bean_RatioData$Treatment

j3<-ggplot(data=Bean_RatioData, aes(x=Treatment, y=Bean_PercActive,fill=Treatment))+
  #need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
  geom_boxplot(outlier.shape=NA) + 
  scale_size(guide=FALSE)+
  scale_x_discrete(name="Treatment")+
  scale_y_continuous(name="16s RNA:DNA Active OTUs (%)",limits=c(0,75))+
  theme_bw(base_size=12)+
  ggtitle("(D)")+
  theme(legend.position="none")+
  scale_fill_manual(limits=c("Control","Drought","Nutrients"),values=c("blue","red","green4"))+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank())
j3

j4<-ggplot(data=Bean_CellStainData, aes(x=Treatment, y=Perc_Active, fill=Treatment))+
  #need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
  geom_boxplot(outlier.shape=NA) + 
  scale_size(guide=FALSE)+
  scale_y_continuous(name="CTC Stain: Active Cells (%)",limits=c(0,75))+
  scale_x_discrete(limits=c("Control","Drought","Nutrients"))+
  theme_bw(base_size=12)+
  ggtitle("(B)")+
  theme(legend.position="none")+
  scale_fill_manual(limits=c("Control","Drought","Nutrients"),values=c("blue","red","green4"))+
  
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank())
j4


#Making Figure4
library("gridExtra")
Figure4 <- grid.arrange(j2,j4,j1,j3, ncol = 2,nrow=2, widths=c(3,2.8),heights = c(2.5, 3))
Figure4
ggsave("Figure4.pdf",Figure4,width=300,height=174, units="mm")



#Do ANOVA to compare percent of OTUs that are active.
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Corn_16sRatio_fit <- aov(Corn_PercActive ~ Treatment, data=Corn_RatioData, na.action=na.exclude)
drop1(Corn_16sRatio_fit,~.,test="F") # type III SS and F Tests
kruskal.test(Corn_PercActive ~ Treatment, data=Corn_RatioData, na.action=na.exclude)
#TESTING ASSUMPTIONS
Corn_16sRatio__resids <- residuals(Corn_16sRatio_fit)
Corn_16sRatio__preds <- predict(Corn_16sRatio_fit)
plot(Corn_16sRatio__resids ~ Corn_16sRatio__preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(Corn_PercActive ~ Treatment, data=Corn_RatioData, na.action=na.exclude)
library("car")
leveneTest(Corn_PercActive ~ Treatment, data=Corn_RatioData, na.action=na.exclude)
#Test for normality
shapiro.test(Corn_16sRatio_resids)
qqnorm(Corn_16sRatio_resids)
qqline(Corn_16sRatio_resids)
TukeyHSD(Corn_16sRatio_fit)

#find mean percent activity within each treatment
#allows one to calculate average difference between 'post-water' and phytohormone treatments
library("plyr")
corn_table <- ddply(Corn_RatioData,~Treatment,summarise,mean=mean(Corn_PercActive))
corn_table

#Do ANOVA to compare percent of OTUs that are active.
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Bean_16sRatio_fit <- aov(Bean_PercActive ~ Treatment, data=Bean_RatioData, na.action=na.exclude)
drop1(Bean_16sRatio_fit,~.,test="F") # type III SS and F Tests
kruskal.test(Bean_PercActive ~ Treatment, data=Bean_RatioData, na.action=na.exclude)
#TESTING ASSUMPTIONS
Bean_16sRatio_resids <- residuals(Bean_16sRatio_fit)
Bean_16sRatio_preds <- predict(Bean_16sRatio_fit)
plot(Bean_16sRatio_resids ~ Bean_16sRatio_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(Bean_PercActive ~ Treatment, data=Bean_RatioData, na.action=na.exclude)
library("car")
leveneTest(Bean_PercActive ~ Treatment, data=Bean_RatioData, na.action=na.exclude)
#Test for normality
shapiro.test(Bean_16sRatio_resids)
qqnorm(Bean_16sRatio_resids)
qqline(Bean_16sRatio_resids)
TukeyHSD(Bean_16sRatio_fit)


#remove unnecessary columns
bean_table<-subset(Bean_RatioData,select=c(Bean_PercActive, Treatment))
#quick way to give us one row per treatment and subsampling level.
bean_table<-bean_table %>% group_by(Treatment) %>% summarise_all(funs(mean))
bean_table



#Do ANOVA to compare percent of cells that are active.
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Corn_CellStain_fit <- aov(Perc_Active ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
drop1(Corn_CellStain_fit,~.,test="F") # type III SS and F Tests
kruskal.test(Perc_Active ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
#TESTING ASSUMPTIONS
Corn_CellStain_resids <- residuals(Corn_CellStain_fit)
Corn_CellStain_preds <- predict(Corn_CellStain_fit)
plot(Corn_CellStain_resids ~ Corn_CellStain_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(Perc_Active ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
library("car")
leveneTest(Perc_Active ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
#Test for normality
shapiro.test(Corn_CellStain_resids)
qqnorm(Corn_CellStain_resids)
qqline(Corn_CellStain_resids)
TukeyHSD(Corn_CellStain_fit)

#find mean percent activity within each treatment
#allows one to calculate average difference between 'post-water' and phytohormone treatments
library("plyr")
corn_table <- ddply(Corn_CellStainData,~Treatment,summarise,mean=mean(Perc_Active))
corn_table

#Do ANOVA to compare percent of cells that are active.
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Bean_CellStain_fit <- aov(Perc_Active ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
drop1(Bean_CellStain_fit,~.,test="F") # type III SS and F Tests
kruskal.test(Perc_Active ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
#TESTING ASSUMPTIONS
Bean_CellStain_resids <- residuals(Bean_CellStain_fit)
Bean_CellStain_preds <- predict(Bean_CellStain_fit)
plot(Bean_CellStain_resids ~ Bean_CellStain_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(Perc_Active ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
library("car")
leveneTest(Perc_Active ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
#Test for normality
shapiro.test(Bean_CellStain_resids)
qqnorm(Bean_CellStain_resids)
qqline(Bean_CellStain_resids)
TukeyHSD(Bean_CellStain_fit)




#remove unnecessary columns
bean_table<-subset(Bean_CellStainData,select=c(Perc_Active, Treatment))
#quick way to give us one row per treatment and subsampling level.
bean_table<-bean_table %>% group_by(Treatment) %>% summarise_all(funs(mean))
bean_table











# PREPARE SUPPLEMENTARY FIGURE 3 -----------------------------------------


#Prepare data for figure
library("dplyr")

#First, find nmber of replicate samples per treatment and subsampling level.
Corn_RNAnoDNA <- add_count(Corn_RNAnoDNA, subsamp_level,Treatment)
#remove unnecessary columns
Corn_RNAnoDNA<-subset(Corn_RNAnoDNA,select=-c(sample,PercentRNAnoDNA))
#find mean of column "n"...
#quick way to give us one row per treatment and subsampling level.
Corn_RNAnoDNA<-Corn_RNAnoDNA %>% group_by(subsamp_level,Treatment) %>% summarise_all(funs(mean))

#First, find nmber of replicate samples per treatment and subsampling level.
Bean_RNAnoDNA <- add_count(Bean_RNAnoDNA, subsamp_level,Treatment)
#remove unnecessary columns
Bean_RNAnoDNA<-subset(Bean_RNAnoDNA,select=-c(sample,PercentRNAnoDNA,Rep))
#find mean of column "n"...
#quick way to give us one row per treatment and subsampling level.
Bean_RNAnoDNA<-Bean_RNAnoDNA %>% group_by(subsamp_level,Treatment) %>% summarise_all(funs(mean))

#Plot the two graphs together
Corn_RNAnoDNA_N_Fig<- ggplot(Corn_RNAnoDNA, aes(subsamp_level, n, colour=Treatment))+
  geom_point(position=position_jitterdodge(jitter.height = .1,dodge.width = 0))+
  geom_line()+
  ylab("Number of samples (n)")+
  coord_cartesian(ylim=c(1,8))+
  ggtitle("(A)")+
  theme_bw()+
  theme(legend.position="bottom",legend.title=element_text(size=8))+
  theme(legend.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 2,byrow=TRUE))+
  xlab("Level of Subsampling")+
  scale_color_manual(limits=c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"),values=c("chocolate4","magenta4","hotpink","blue","red","yellow","orange","green4"))
Corn_RNAnoDNA_N_Fig

Bean_RNAnoDNA_N_Fig<- ggplot(Bean_RNAnoDNA, aes(subsamp_level, n, colour=Treatment))+
  geom_point(position=position_jitterdodge(jitter.height = .1,dodge.width = 0))+
  geom_line()+
  coord_cartesian(ylim=c(1,8))+
  ggtitle("(B)")+
  theme_bw()+
  theme(legend.position="bottom",axis.title.y=element_blank(),axis.text.y=element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  #ggtitle("Bean: Percent of Total OTUs with RNA1, DNA0")+
  xlab("Level of Subsampling")+
  scale_color_manual(limits=c("Control", "Drought", "Nutrients"),values=c("blue","red","green4"))
Bean_RNAnoDNA_N_Fig

library("gridExtra")
SuppFig3 <- grid.arrange(Corn_RNAnoDNA_N_Fig,Bean_RNAnoDNA_N_Fig,ncol = 2, widths=c(2,1.8),heights = c(2))
SuppFig3
ggsave("SuppFig3.pdf",SuppFig3,width=174,height=100, units="mm")



# PREPARE SUPPLEMENTARY FIGURE 4 ------------------------------------------

#Reorder treatment 'factors' for the figure (so they are no longer alphabetical, but are in the order I want)
Corn_CellStainData$Treatment = factor(Corn_CellStainData$Treatment,c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"))
require(scales)
k1<-ggplot(data=Corn_CellStainData, aes(x=Treatment, y=Syto_Counts_per_gram, fill=Treatment))+
  #need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
  geom_boxplot(outlier.shape=NA) + 
  scale_size(guide=FALSE)+
  scale_x_discrete(name="Treatment")+
  scale_x_discrete(limits=c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"))+
  scale_y_continuous(name="Syto Counts Per Gram Soil", limits=c(0,400000),labels=comma)+
  theme_bw(base_size=12)+
  theme(legend.position="none")+
  scale_fill_manual(name="Treatment", breaks=c("Pre-dry", "Post-dry", "Post-water", "Control", "ABA", "IAA", "JA", "SA"), values=c("chocolate4", "magenta4","hotpink", "blue", "red","yellow", "orange", "green4"))+
  ggtitle("(C)")
k1

k2<-ggplot(data=Corn_CellStainData, aes(x=Treatment, y=CTC_Counts_per_gram, fill=Treatment))+
  #need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
  geom_boxplot(outlier.shape=NA) + 
  scale_size(guide=FALSE)+
  scale_y_continuous(name="CTC Counts Per Gram Soil", limits=c(0,400000),labels=comma)+
  theme_bw(base_size=12)+
  scale_x_discrete(limits=c("Pre-dry","Post-dry","Post-water","Control","ABA","IAA","JA","SA"))+
  #scale_fill_manual(limits=c("CTC-Stain"),values=c("red"))+
  ggtitle("(A)")+
  theme(legend.position="none")+
  scale_fill_manual(name="Treatment", breaks=c("Pre-dry", "Post-dry", "Post-water", "Control", "ABA", "IAA", "JA", "SA"), values=c("chocolate4", "magenta4","hotpink", "blue", "red","yellow", "orange", "green4"))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
k2



k3<-ggplot(data=Bean_CellStainData, aes(x=Treatment, y=Syto_Counts_per_gram,fill=Treatment))+
  #need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
  geom_boxplot(outlier.shape=NA) + 
  scale_size(guide=FALSE)+
  scale_x_discrete(name="Treatment")+
  scale_y_continuous(name="Syto Counts Per Gram Soil", limits=c(0,1000000),labels=comma)+
  theme_bw(base_size=12)+
  ggtitle("(D)")+
  theme(legend.position="none")+
  scale_fill_manual(limits=c("Control","Drought","Nutrients"),values=c("blue","red","green4"))+
  theme(axis.title.y=element_blank())
k3

k4<-ggplot(data=Bean_CellStainData, aes(x=Treatment, y=CTC_Counts_per_gram, fill=Treatment))+
  #need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
  geom_boxplot(outlier.shape=NA) + 
  scale_size(guide=FALSE)+
  scale_y_continuous(name="CTC Counts Per Gram Soil", limits=c(0,1000000),labels=comma)+
  scale_x_discrete(limits=c("Control","Drought","Nutrients"))+
  theme_bw(base_size=12)+
  ggtitle("(B)")+
  theme(legend.position="none")+
  scale_fill_manual(limits=c("Control","Drought","Nutrients"),values=c("blue","red","green4"))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(axis.title.y=element_blank())
k4


#Making Supp Fig4
library("gridExtra")
SuppFig4 <- grid.arrange(k2,k4,k1,k3, ncol = 2,nrow=2, widths=c(3,2.8),heights = c(2.5, 3))
SuppFig4
ggsave("SuppFig4.pdf",SuppFig4,width=300,height=174, units="mm")





#Do ANOVA to compare Syto Counts

options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Corn_CellStain_fit <- aov(Syto_Counts_per_gram ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
drop1(Corn_CellStain_fit,~.,test="F") # type III SS and F Tests
kruskal.test(Syto_Counts_per_gram ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
#TESTING ASSUMPTIONS
Corn_CellStain_resids <- residuals(Corn_CellStain_fit)
Corn_CellStain_preds <- predict(Corn_CellStain_fit)
plot(Corn_CellStain_resids ~ Corn_CellStain_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(Perc_Active ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
library("car")
leveneTest(Perc_Active ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
#Test for normality
shapiro.test(Corn_CellStain_resids)
qqnorm(Corn_CellStain_resids)
qqline(Corn_CellStain_resids)
TukeyHSD(Corn_CellStain_fit)

options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Bean_CellStain_fit <- aov(Syto_Counts_per_gram ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
drop1(Bean_CellStain_fit,~.,test="F") # type III SS and F Tests
kruskal.test(Syto_Counts_per_gram ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
#TESTING ASSUMPTIONS
Bean_CellStain_resids <- residuals(Bean_CellStain_fit)
Bean_CellStain_preds <- predict(Bean_CellStain_fit)
plot(Bean_CellStain_resids ~ Bean_CellStain_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(Perc_Active ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
library("car")
leveneTest(Perc_Active ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
#Test for normality
shapiro.test(Bean_CellStain_resids)
qqnorm(Bean_CellStain_resids)
qqline(Bean_CellStain_resids)
TukeyHSD(Bean_CellStain_fit)

#calculate average difference between 'post-water' and phytohormone treatments
#allows one to calculate average difference between 'post-water' and phytohormone treatments
library("plyr")
corn_table <- ddply(Corn_CellStainData,~Treatment,summarise,mean=mean(Syto_Counts_per_gram))
corn_table



#Do ANOVA to compare CTC Counts

options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Corn_CellStain_fit <- aov(CTC_Counts_per_gram ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
drop1(Corn_CellStain_fit,~.,test="F") # type III SS and F Tests
kruskal.test(CTC_Counts_per_gram ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
#TESTING ASSUMPTIONS
Corn_CellStain_resids <- residuals(Corn_CellStain_fit)
Corn_CellStain_preds <- predict(Corn_CellStain_fit)
plot(Corn_CellStain_resids ~ Corn_CellStain_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(CTC_Counts_per_gram ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
library("car")
leveneTest(CTC_Counts_per_gram ~ Treatment, data=Corn_CellStainData, na.action=na.exclude)
#Test for normality
shapiro.test(Corn_CellStain_resids)
qqnorm(Corn_CellStain_resids)
qqline(Corn_CellStain_resids)
TukeyHSD(Corn_CellStain_fit)

options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Bean_CellStain_fit <- aov(CTC_Counts_per_gram ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
drop1(Bean_CellStain_fit,~.,test="F") # type III SS and F Tests
kruskal.test(CTC_Counts_per_gram ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
#TESTING ASSUMPTIONS
Bean_CellStain_resids <- residuals(Bean_CellStain_fit)
Bean_CellStain_preds <- predict(Bean_CellStain_fit)
plot(Bean_CellStain_resids ~ Bean_CellStain_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
bartlett.test(CTC_Counts_per_gram ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
library("car")
leveneTest(CTC_Counts_per_gram ~ Treatment, data=Bean_CellStainData, na.action=na.exclude)
#Test for normality
shapiro.test(Bean_CellStain_resids)
qqnorm(Bean_CellStain_resids)
qqline(Bean_CellStain_resids)
TukeyHSD(Bean_CellStain_fit)


#find mean percent activity within each treatment
#allows one to calculate average difference between 'post-water' and phytohormone treatments
library("plyr")
corn_table <- ddply(Corn_CellStainData,~Treatment,summarise,mean=mean(CTC_Counts_per_gram))
corn_table

