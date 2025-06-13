#  title: "Uncovering the hidden diversity and functional roles of root endophytic Streptomyces in Sorghum bicolor under drought stress"
#author: "Citlali Fonseca-Garcia"
#Last edit on 06/11/2025

library(ggplot2); packageVersion("ggplot2")
library(patchwork);packageVersion("patchwork")
library(agricolae);packageVersion("agricolae")
library(vegan); packageVersion("vegan")
library(mctoolsr); packageVersion("mctoolsr")

setwd("~working_directory/")


####Non-Polar data E
Strep12_neg <- read.csv("Strep12_merge-root-swap_NonPolar_E.csv")
Strep12_neg.test <- Strep12_neg[1:48,] #modify the last column if you want to exclude the control samples
#Make a PCA plot
pc <- Strep12_neg
com <- pc[,9:ncol(pc)]
m_com <- as.matrix(com)
# log transform pc
log.com <- log(com+1)
log.m_com <- as.matrix(log.com)

# PCA plot
Log_Strep12_neg.pca <- prcomp(log.com[1:48,]) #modify the last column if you want to exclude the control samples
Log_Strep12_neg.pca
summary(Log_Strep12_neg.pca)


pca.plot <- ggplot2::autoplot(Log_Strep12_neg.pca, x=1, y=2, data = Strep12_neg.test, colour = "Ac_group")+
  geom_point(size = 14, aes(colour = Ac_group))+ 
  #geom_vline(xintercept = 0,linetype="longdash",color="gray")+
  #geom_hline(yintercept = 0,linetype="longdash",color="gray")+
  geom_text(size = 4, color = "black", aes(label=Isolate_num), show.legend = F, hjust=0.5, vjust=0.5)+
  theme_set(theme_bw(base_size = 20)) + theme( plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values=c("#9933FF","#6699FF", "#CCCC99","#669933","#CC9966"))


pca.plot #export size 10x6
ggsave("TWYE_NonPolar_noControlPC12.png", plot = pca.plot, width=12, height=6.5, dpi=600)

####Polar data E
Strep12_neg <- read.csv("Strep12_merge-root-swap_Polar_E_metadata.csv")
Strep12_neg.test <- Strep12_neg[1:48,]#modify the last column if you want to exclude the control samples
#Make a PCA plot
pc <- Strep12_neg
com <- pc[,9:ncol(pc)]
m_com <- as.matrix(com)
# log transform pc
log.com <- log(com+1)
log.m_com <- as.matrix(log.com)

# PCA plot
Log_Strep12_neg.pca <- prcomp(log.com[1:48,])#modify the last column if you want to exclude the control samples
Log_Strep12_neg.pca
summary(Log_Strep12_neg.pca)


pca.plot <- ggplot2::autoplot(Log_Strep12_neg.pca, x=1, y=2, data = Strep12_neg.test, colour = "Ac_group")+
  geom_point(size = 14, aes(colour = Ac_group))+ 
  #geom_vline(xintercept = 0,linetype="longdash",color="gray")+
  #geom_hline(yintercept = 0,linetype="longdash",color="gray")+
  geom_text(size = 4, color = "black", aes(label=Isolate_num), show.legend = F, hjust=0.5, vjust=0.5)+
  theme_set(theme_bw(base_size = 20)) + theme( plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values=c("#9933FF","#6699FF", "#CCCC99","#669933","#CC9966"))


pca.plot #export size 10x6

ggsave("TWYE_Polar_noControl_PC12.png", plot = pca.plot, width=12, height=6.5, dpi=600)

####Non polar data Root-D and Root-W 
Strep12_neg <- read.csv("Strep12_merge-root-swap_NonPolar_metadata.csv")
Strep12_neg.test <- Strep12_neg[1:104,]#modify the last column if you want to exclude the control samples
#Make a PCA plot
pc <- Strep12_neg
com <- pc[,9:ncol(pc)]
m_com <- as.matrix(com)
# log transform pc
log.com <- log(com+1)
log.m_com <- as.matrix(log.com)

# PCA plot
Log_Strep12_neg.pca <- prcomp(log.com[1:104,])#modify the last column if you want to exclude the control samples
Log_Strep12_neg.pca
summary(Log_Strep12_neg.pca)


pca.plot <- ggplot2::autoplot(Log_Strep12_neg.pca, x=1, y=2, data = Strep12_neg.test, colour = "Ac_group")+
  geom_point(size = 12, aes(colour = Ac_group, shape = Treatment))+ 
  #geom_vline(xintercept = 0,linetype="longdash",color="gray")+
  #geom_hline(yintercept = 0,linetype="longdash",color="gray")+
  geom_text(size = 3, color = "black", aes(label=Isolate_num), show.legend = F, hjust=0.5, vjust=0.5)+
  theme_set(theme_bw(base_size = 18)) + theme( plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values=c("#9933FF","#6699FF", "#CCCC99","#669933","#CC9966"))


pca.plot #export size 10x6

ggsave("Root_NonPolar.png", plot = pca.plot, width=12, height=6.5, dpi=600)


###PERMANOVA from osmotic stress and siderophores phenotyping data 
DataBase <- read.table("DC-Strains_phenotypes_O_S.txt", header=TRUE) #update path info
Metadata <- read.table("DC-Strains_phenotypes_metadata_O_S.txt", header=TRUE) #update path info

rownames(DataBase) <- DataBase[,1]
DataBase <- DataBase[,-1 ] 
DataBase <- as.data.frame(DataBase)

rownames(Metadata) <- Metadata[,1]
Metadata <- Metadata[,-1]
Metadata <- as.data.frame(Metadata)

typeof(DataBase)
#convert the df in numeric
i <- c(1:48)    # Specify columns you want to change
DataBase[ , i] <- apply(DataBase[ , i], 2,            # Specify own function within apply
                        function(x) as.numeric(as.character(x)))

matrix <- vegdist (DataBase, method  = "euclidean", na.rm = TRUE)
matrix <- as.matrix (matrix) 

test <- adonis(matrix ~ Isolate, data = Metadata)
test[["aov.tab"]]
write.table(test[["aov.tab"]], file = "aov.tab_GlobalIsolate_O_S.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = TRUE)

test <- adonis(matrix ~ Ac_group, data = Metadata)
test[["aov.tab"]]
write.table(test[["aov.tab"]], file = "aov.tab_GlobalAc_group_O_S.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = TRUE)

#adonis posthoc test: Pairwise permanova

Metadata$Ac_group <- as.factor(Metadata$Ac_group)
Metadata$Isolate <- as.factor(Metadata$Isolate)

testT <- calc_pairwise_permanovas(matrix, Metadata, "Isolate")
testT
write.table(testT, file = "aov.tab_pairwiseIsolate_O_S.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = TRUE)

testG <- calc_pairwise_permanovas(matrix, Metadata, "Ac_group")
testG
write.table(testG, file = "aov.tab_pairwiseAc_group_O_S.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = TRUE)
