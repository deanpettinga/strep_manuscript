##
#  title: "Uncovering the hidden diversity and functional roles of root endophytic Streptomyces in Sorghum bicolor under drought stress"
#author: "Daniel Caddell"
#Last edit on 02/25/2019
#
###################################################################
## Analysis of Primer bias between Illumina MiSeq and PacBio CCS ##
###################################################################

### Phylum analysis to describe primer bias using Silva TestPrime
library(dada2);packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(phyloseq) 
library("RColorBrewer")

####################
## Illumina MiSeq ##
####################
setwd("~working_directory/")

biom_file <- "otu.biom"
map_file <- "Metadata_Epicon_CSP2017.txt"

otutax <- import_biom(biom_file, parseFunction = parse_taxonomy_greengenes) #this is the biom file which is converted from itagger output http://biom-format.org/documentation/biom_conversion.html
sam <- import_qiime_sample_data(map_file) #this is your metadata for your project
all <- merge_phyloseq(otutax, sam) #also can add: ,tree)
all
colnames(tax_table(all)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
colnames(tax_table(all))
rank_names(all)

SampleOrderRZ <-  scale_x_discrete(limits = c(#"DC7_fastq","DC9_fastq","DC12_fastq", #TP8, control, soil
  #"DC25_fastq","DC27_fastq","DC36_fastq", #TP8_24h, control, soil                              
  #"DC10_fastq", "DC11_fastq", "DC8_fastq", #TP8, drought, soil
  #"DC26_fastq","DC28_fastq","DC35_fastq", #TP8_24h, drought, soil                             
  "DC2_fastq","DC3_fastq","DC6_fastq",#TP8, control, rhizosphere
  "DC19_fastq","DC21_fastq","DC24_fastq", #TP8_24h, control, rhizosphere                              
  "DC1_fastq","DC4_fastq","DC5_fastq",#TP8, drought, rhizosphere
  "DC20_fastq","DC22_fastq","DC23_fastq", #TP8_24h, drought, rhizosphere                              
  "DC13_fastq","DC15_fastq","DC18_fastq",#TP8, control, root
  "DC29_fastq","DC31_fastq","DC34_fastq",#TP8_24h, control, root                              
  "DC14_fastq","DC16_fastq","DC17_fastq", #TP8, drought, root
  "DC30_fastq","DC32_fastq","DC33_fastq"), #TP8_24h, drought, root
  labels = c(#"TP8_Control_Soil","TP8_Control_Soil","TP8_Control_Soil",
    #"24h_Control_Soil","24h_Control_Soil","24h_Control_Soil",                               
    #"TP8_Drought_Soil","TP8_Drought_Soil","TP8_Drought_Soil",
    #"24h_Watered_Soil","24h_Watered_Soil","24h_Watered_Soil",
    "TP8_Control_Rhizosphere","TP8_Control_Rhizosphere","TP8_Control_Rhizosphere",
    "24h_Control_Rhizosphere","24h_Control_Rhizosphere","24h_Control_Rhizosphere",                              
    "TP8_Drought_Rhizosphere","TP8_Drought_Rhizosphere","TP8_Drought_Rhizosphere",
    "24h_Watered_Rhizosphere","24h_Watered_Rhizosphere","24h_Watered_Rhizosphere",
    "TP8_Control_Root","TP8_Control_Root","TP8_Control_Root",
    "24h_Control_Root","24h_Control_Root","24h_Control_Root",                              
    "TP8_Drought_Root","TP8_Drought_Root","TP8_Drought_Root",
    "24h_Watered_Root","24h_Watered_Root","24h_Watered_Root"))

all.tf <- transform_sample_counts(all, function(OTU) 100*OTU/sum(OTU))
all_tf.df<-psmelt(all.tf)

phylum.list <- c('Acidobacteria', 'Actinobacteria', 'Bacteroidetes','Firmicutes','Proteobacteria')
'%ni%' <- Negate('%in%')
all_tf.df$Phylum <- as.character(all_tf.df$Phylum)
all_tf.df[all_tf.df$Phylum %ni% phylum.list,]$Phylum <-'Other'
all_tf.df$Phylum <- as.factor(all_tf.df$Phylum)

tempV3V4 <- c(labels=c('Acidobacteria', 'Actinobacteria', 'Bacteroidetes', 'Firmicutes', 'Proteobacteria', 'Other'))
all_tf.df$Phylum <- factor(all_tf.df$Phylum, levels = c(tempV3V4))

ggplot(all_tf.df[order(all_tf.df$Phylum),], aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(stat = "identity", position ="stack", width = 0.95)+
  SampleOrderRZ +
  scale_fill_manual(labels=c('Acidobacteria', 'Actinobacteria', 'Bacteroidetes', 'Firmicutes', 'Proteobacteria', 'Other'),
                    values=c('Acidobacteria'='dark red','Actinobacteria'='#c03728','Bacteroidetes'='#fd8f24', 
                             'Firmicutes'='#f5c04a', 'Proteobacteria'='#919c4c', 'Other'='#4f5157'))+
  ggtitle("Illumina MiSeq OTUs") +     
  theme(axis.text.x=element_text(angle = 270, hjust=0, vjust=0.5),
        plot.title = element_text(lineheight=.8, face="bold"),
        legend.position="none") #Export size 5''x 4.5''

#################
## PacBio CCS ##
################
setwd("~working_directory/")

###First need to replace Devin's R data objects into new global environment
st2 <- readRDS ("all_st2.rds")
tax2 <- readRDS ("all_tax2_Silva128.RDS")



ft2 <- sweep(st2, 1, rowSums(st2), "/")
#read in metadata file
df <- read.table("metadata copy.txt", header=TRUE) #metadata copy.txt has long names for treatment and sampletype, metadata shorthand.

#make the rownames <- SampleID column
rownames(df) <- df$SampleID 
df <- df[sample.names2,-1]
#reorder the metadata rows
df<-df[c("RD1", "RD2", "RD3", "RW1", "RW2", "RW3", "ZD1", "ZD2", "ZD3", "ZW1", "ZW2", "ZW3", "RD1-2", "RD2-2", "RD3-2", "RW1-2", "RW2-2", "RW3-2", "ZD1-2", "ZD2-2", "ZD3-2", "ZW1-2", "ZW2-2", "ZW3-2"),]
df$Timepoint <- as.factor(df$Timepoint)
df$Rep <- as.factor(df$Rep)
df$SampleType <- as.factor(df$SampleType)
df$Treatment <- as.factor(df$Treatment)
df$SampleID <- as.factor(rownames(df))
df$SampleID <- factor(df$SampleID,levels=c("RD1", "RD2", "RD3", "RW1", "RW2", "RW3", "ZD1", "ZD2", "ZD3", "ZW1", "ZW2", "ZW3", "RD1-2", "RD2-2", "RD3-2", "RW1-2", "RW2-2", "RW3-2", "ZD1-2", "ZD2-2", "ZD3-2", "ZW1-2", "ZW2-2", "ZW3-2"))

SampleOrder <-   scale_x_discrete(limits = c(#"ZW1","ZW2","ZW3",#TP8, control, rhizosphere
                                             #"ZW1-2","ZW2-2","ZW3-2", #TP8_24h, control, rhizosphere                              
                                             #"ZD1","ZD2","ZD3",#TP8, drought, rhizosphere
                                             #"ZD1-2","ZD2-2","ZD3-2", #TP8_24h, drought, rhizosphere                              
                                             "RW1","RW2","RW3",#TP8, control, root
                                             #"RW1-2","RW2-2","RW3-2",#TP8_24h, control, root                              
                                             "RD1","RD2","RD3"), #TP8, drought, root
                                             #"RD1-2","RD2-2","RD3-2"), #TP8_24h, drought, root
                                  labels = c(#"TP8_Control_Rhizosphere","TP8_Control_Rhizosphere","TP8_Control_Rhizosphere",
                                             #"24h_Control_Rhizosphere","24h_Control_Rhizosphere","24h_Control_Rhizosphere",                              
                                             #"TP8_Drought_Rhizosphere","TP8_Drought_Rhizosphere","TP8_Drought_Rhizosphere",
                                             #"24h_Watered_Rhizosphere","24h_Watered_Rhizosphere","24h_Watered_Rhizosphere",
                                             "TP8_Control_Root","TP8_Control_Root","TP8_Control_Root",
                                             #"24h_Control_Root","24h_Control_Root","24h_Control_Root",                              
                                             "TP8_Drought_Root","TP8_Drought_Root","TP8_Drought_Root"))
                                             #"24h_Watered_Root","24h_Watered_Root","24h_Watered_Root"))
ReplicateMerge <- scale_x_discrete(limits = c( "TP8_Control_Root",
                                               "TP8_Drought_Root" ),                        
                                               
                                   labels = c( "TP8_Control_Root",
                                               "TP8_Drought_Root"
                                   ))

ps2 <- phyloseq(otu_table(ft2, taxa_are_rows=FALSE), sample_data(df))
tax_table(ps2)<-tax2

ps2.t <- transform_sample_counts(ps2, function(OTU) 100*OTU/sum(OTU))
ps2.t.df<-psmelt(ps2.t)

write.csv(ps2.t@tax_table,"ps2.t@tax_table.csv")
write.csv(ps2.t@otu_table,"ps2.t@otu_table.csv")

phylum.list <- c('Acidobacteria', 'Actinobacteria', 'Bacteroidetes','Firmicutes','Proteobacteria')
'%ni%' <- Negate('%in%')
ps2.t.df$Phylum <- as.character(ps2.t.df$Phylum)
ps2.t.df[ps2.t.df$Phylum %ni% phylum.list,]$Phylum <-'Other'
ps2.t.df$Phylum <- as.factor(ps2.t.df$Phylum)


tempPB <- c(labels=c('Acidobacteria', 'Actinobacteria', 'Bacteroidetes', 'Firmicutes', 'Proteobacteria', 'Other'))
ps2.t.df$Phylum <- factor(ps2.t.df$Phylum, levels = c(tempPB))

ggplot(ps2.t.df[order(ps2.t.df$Phylum),], aes(x=Sample, y=Abundance, fill=Phylum))+
  geom_bar(stat = "identity", position ="stack", width = 0.95)+
  ReplicateMerge +
  scale_fill_manual(labels=c('Acidobacteria', 'Actinobacteria', 'Bacteroidetes', 'Firmicutes', 'Proteobacteria', 'Other'),
                    values=c('Acidobacteria'='dark red','Actinobacteria'='#c03728','Bacteroidetes'='#fd8f24', 
                             'Firmicutes'='#f5c04a', 'Proteobacteria'='#919c4c', 'Other'='#4f5157'))+
  ggtitle("PacBio CCS ASVs") +     
  theme(axis.text.x=element_text(angle = 270, hjust=0, vjust=0.5),
        plot.title = element_text(lineheight=.8, face="bold"),
        legend.position="none") #export size 5''x 5''

##############################
## Shared legend for phylum ##
##############################
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)+
  legend(0.2,1, legend =c('Acidobacteria', 'Actinobacteria', 'Bacteroidetes', 'Firmicutes', 'Proteobacteria', 'Other'), 
         pch=15, pt.cex=3, cex=1, bty='n',
         col = c('Acidobacteria'='dark red','Actinobacteria'='#c03728','Bacteroidetes'='#fd8f24', 
                 'Firmicutes'='#f5c04a', 'Proteobacteria'='#919c4c', 'Other'='#4f5157'))
mtext("Phylum", line=0, adj=0.27, padj=0.75, cex=1.5) #export size 6''x 4'' or 6"x 4.5"     




##########
#Shared tables with ASV from Daniel
#########
otu.table <- read.csv ("C:/Users/fonse/Documents/Mis documentos/ESTANCIA POSDOCTORAL/PosdocDevin/Strep project/Figure 1A_files/ps2.t@otu_table.csv")
tax.table <- read.csv ("C:/Users/fonse/Documents/Mis documentos/ESTANCIA POSDOCTORAL/PosdocDevin/Strep project/Figure 1A_files/ps2.t@tax_table.csv")
DC.table <- read.csv ("C:/Users/fonse/Documents/Mis documentos/ESTANCIA POSDOCTORAL/PosdocDevin/Strep project/Figure 1A_files/ASV_Ac_DC.csv")

#merge otu table with taxa and Ac groups
otu.tax.table <- merge(otu.table, tax.table, sort = FALSE)
otu.DC.table <- merge(otu.tax.table, DC.table, sort = FALSE)
otu.DC.table <- as.matrix(otu.DC.table)
write.table (otu.DC.table, "ASV.DC.table.csv")

###exporting read counts from Daniel's data
st2.t <- t(st2)
write.csv(st2.t,"otu.readcounts_table_full.csv")

#merge previous tables with read counts
otu.DC.table.ed <- read.csv ("C:/Users/fonse/Documents/Mis documentos/ESTANCIA POSDOCTORAL/PosdocDevin/Strep project/Figure 1A_files/ASV.DC.table_full_edited.csv")
readcounts.table <- read.csv ("C:/Users/fonse/Documents/Mis documentos/ESTANCIA POSDOCTORAL/PosdocDevin/Strep project/Figure 1A_files/otu.readcounts_table_full.csv")

readcounts.ASV.table <- merge(otu.DC.table.ed, readcounts.table, sort = FALSE)
write.table (otu.DC.table, "ASV.DC.table.csv")


