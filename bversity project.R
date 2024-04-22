library(Biobase)
library(limma)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(geneplotter)
library(pheatmap)
library(enrichplot)
library(tidyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(GEOquery)
library(umap)
library(org.Hs.eg.db)
library(AnnotationDbi)

#BiocManager::install("org.Hs.eg.db",force=TRUE)
 setwd("C:/Users/rajas/OneDrive/Documents/capstoneproject")
#phenodata

gse <- getGEO(filename = "C:/Users/rajas/OneDrive/Documents/capstoneproject/GSE200306_series_matrix.txt",GSEMatrix = TRUE)
phenoData <- pData(gse)

#modify rownames as geo_accession
rownames(phenoData) <- phenoData$geo_accession
# to add column with sample_id unique to each sample mentioned in expressiondata columns
phenoData_title <- dplyr::select(phenoData, title)

# Extract the desired information using pattern matching and substitution:
phenoData_extracted <- sub(".*\\[(.*)\\].*", "\\1", phenoData_title$title)
phenoData_extracted<- data.frame(title = phenoData_extracted, row.names = rownames(phenoData_title))
phenoData_extracted<-phenoData_extracted %>% dplyr::rename(sample_id=title)


## Extract patient information using regular expressions:
patient_info <- sub("^(.*?),.*", "\\1",phenoData_title$title)
##Print the extracted patient information:
patient_info<- data.frame(title = patient_info, row.names = rownames(phenoData_title))
patient_info<-patient_info %>% dplyr::rename(patient_info=title)


## Check if any values are duplicated:

#no_duplicates <- !any(duplicated(phenoData_extracted$sample_id))

##Print the result (TRUE for no duplicates, FALSE otherwise):
print(no_duplicates)

##to see both tissue extractd from same patient
#no_duplicates <- !any(duplicated(patient_info$patient_info))
## Print the result (TRUE for no duplicates, FALSE otherwise):
#print(no_duplicates)

#library(dplyr)
#add sample_id from phenoData_extracted to phenoData datframe 

phenoData$sample_id <- phenoData_extracted$sample_id
phenoData$patient_info <- patient_info$patient_info
phenoData <-relocate(phenoData,geo_accession, .before = 1)
phenoData <- relocate(phenoData, sample_id, .before = 2)
phenoData <- relocate(phenoData, characteristics_ch1, .before = 3)


library(S4Vectors)
phenoData <- head(phenoData, 79)
phenoData_mixed<-phenoData %>% filter(grepl("ln class: Mixed \\(III/IV\\+V\\)", characteristics_ch1))

#rownames(phenoData_mixed)

phenoData<- phenoData %>% filter(!grepl("ln class: Mixed \\(III/IV\\+V\\)", characteristics_ch1))

phenoData<-phenoData[ ,c(1:4,42)]

phenoData <- mutate(phenoData,
                    stages = case_when(
                      characteristics_ch1 == "ln class: Class IV" ~ "stage4",
                      characteristics_ch1 == "ln class: Class III" ~ "stage3",
                      TRUE ~ "control"  # Default label for other values
                    )
)

#normalised.expressDatagl
#has 158 columns of sample and 3 other columnswith gene-name in one such column
# modify to make gene_name as rows and samples in column

normalised.expressDatagl <- read.delim("C:/Users/rajas/OneDrive/Documents/capstoneproject/GSE200306_normalized_Glom_expression_data.txt")

first_column <- normalised.expressDatagl$Name

rownames(normalised.expressDatagl) <- first_column

normalised.expressDatagl<- normalised.expressDatagl %>%
  dplyr::select(-Name,-Accession, -Code.Class)


exprsData<- normalised.expressDatagl

#remove column of mixed and other normalisation colums from normalised.expressionDatagl 
exprsData<- select(exprsData,
                   -X7a_7 , -X12A_8, -X15A_5 ,-X23A_9,-X32A,-X10A,-X7b_7,-X15B_5,-X10B ,-X24B , -X4B_7 , -X16b_9,  
                    -CR1_mean, -CR2_mean,-NORM_mean,-NR1_mean,-NR2_mean,-PR1_mean,-PR2_mean,-CR1_vs_Normal_est,
                    -NR1_vs_Normal_est,  - PR1_vs_Normal_est,-CR2_vs_Normal_est,-NR2_vs_Normal_est,-PR2_vs_Normal_est,      
                    -CR2_vs_CR1_est,-NR2_vs_NR1_est,-PR2_vs_PR1_est,-CR1_vs_NR1_est,         
                    -CR2_vs_NR2_est,-CR1_vs_PR1_est,-CR2_vs_PR2_est,-PR1_vs_NR1_est,         
                    -PR2_vs_NR2_est,-CR2_1_vs_NR2_1_est,-CR2_1_vs_PR2_1_est,-NR2_1_vs_PR2_1_est,    
                    -fold.CR1_vs_Normal_est,-p.CR1_vs_Normal_est,-fold.NR1_vs_Normal_est,-p.NR1_vs_Normal_est,    
                    -fold.PR1_vs_Normal_est,-p.PR1_vs_Normal_est,-fold.CR2_vs_Normal_est,-p.CR2_vs_Normal_est,    
                    -fold.NR2_vs_Normal_est,-p.NR2_vs_Normal_est,-fold.PR2_vs_Normal_est,-p.PR2_vs_Normal_est,    
                    -fold.CR2_vs_CR1_est,-p.CR2_vs_CR1_est,-fold.NR2_vs_NR1_est,-p.NR2_vs_NR1_est,       
                    -fold.PR2_vs_PR1_est,-p.PR2_vs_PR1_est,-fold.CR1_vs_NR1_est,-p.CR1_vs_NR1_est,       
                    -fold.CR2_vs_NR2_est,-p.CR2_vs_NR2_est,-fold.CR1_vs_PR1_est,-p.CR1_vs_PR1_est,       
                    -fold.CR2_vs_PR2_est,-p.CR2_vs_PR2_est,-fold.PR1_vs_NR1_est,-p.PR1_vs_NR1_est,      
                    -fold.PR2_vs_NR2_est,-p.PR2_vs_NR2_est,-fold.CR2_1_vs_NR2_1_est,-p.CR2_1_vs_NR2_1_est,
                    -fold.CR2_1_vs_PR2_1_est,-p.CR2_1_vs_PR2_1_est,-fold.NR2_1_vs_PR2_1_est,-p.NR2_1_vs_PR2_1_est,
                    -flag.CR1_vs_Normal,-flag.NR1_vs_Normal,-flag.PR1_vs_Normal,-flag.CR2_vs_Normal,     
                    -flag.NR2_vs_Normal,-flag.PR2_vs_Normal,-flag.CR2_vs_CR1,-flag.NR2_vs_NR1,        
                    -flag.PR2_vs_PR1,-flag.CR1_vs_NR1,-flag.CR2_vs_NR2,-flag.CR1_vs_PR1,        
                    -flag.CR2_vs_PR2,-flag.PR1_vs_NR1,-flag.PR2_vs_NR2,-flag.CR2_1_vs_NR2_1, 
                    -flag.CR2_1_vs_PR2_1,-flag.NR2_1_vs_PR2_1)

#tochange column names in esprsData from sample_id to respective geo_accession number

old_names <- c("X1A.G_3" ,"X02A_5" ,"X06A_1" ,"X11A_5" ,"X20a_8","X26A_2","X30A","X37A_5","X39A_1","X48a_9",
               "X49A_10","X50A_10","X51A_7","X55A_9","X58a_8","X5a_7","X13A","X17a_7","X24A","X27A_9",
               "X38A.G_3","X14a_8","X18A_2","X21A_1","X22A_10","X28A_1","X35A_2","X56A_10","X1B.G_3","X02B_5",
               "X11B_5","X19B_2","X20b_8","X23B_9","X26B_2","X30B","X32B","X37B_5","X39B_1","X48b_9",
                "X49B_10","X50B_10","X51B_7","X53B_1","X55B_9","X58b_8","X05B","X13B","X17B","X27B",
                "X31B_8","X38B.G_3","X14b_8","X18B_2","X21B_1","X22B_10","X28B_1","X56B_10","X1XN_10","X2XN_10",
                 "X3XNG_1","X4XNG_1","X5XNG1_2","X7XNG_2","X10Y","X15Y","X16Y")

new_names<- c("GSM6030648","GSM6030649","GSM6030650","GSM6030652","GSM6030655","GSM6030657",
              "GSM6030658","GSM6030660","GSM6030661","GSM6030662","GSM6030663","GSM6030664",
              "GSM6030665","GSM6030666","GSM6030667","GSM6030668","GSM6030670","GSM6030671",
              "GSM6030672","GSM6030673","GSM6030674","GSM6030675","GSM6030676","GSM6030677",
              "GSM6030678","GSM6030679","GSM6030680","GSM6030681","GSM6030682","GSM6030683",
              "GSM6030685","GSM6030687","GSM6030688","GSM6030689","GSM6030690","GSM6030691",
              "GSM6030692","GSM6030693","GSM6030694","GSM6030695","GSM6030696","GSM6030697",
              "GSM6030698","GSM6030699","GSM6030700","GSM6030701","GSM6030702","GSM6030704",
              "GSM6030705","GSM6030707","GSM6030708","GSM6030709","GSM6030711","GSM6030713",
              "GSM6030714","GSM6030715","GSM6030716","GSM6030717","GSM6030718","GSM6030719",
              "GSM6030720","GSM6030721","GSM6030722","GSM6030723","GSM6030724","GSM6030725","GSM6030726")


# Check if entries are present in column names
#present_cols <- old_names %in% names(exprsData)

# Print results
#if (any(present_cols)) {
 # cat("Entries present as column names:", old_names[present_cols], sep = ", ")
#} else {
# print("No entries from the list are present as column names.")
#}

# TO GET matching AND not matchimg ones

#matching_cols <- old_names[old_names %in% names(exprsData)]

#matching_cols <- old_names[!old_names %in% names(exprsData)]

#Arrange exprsData columns in order like in old names
#head(exprsData)

exprsData <- exprsData[, old_names]


### WORKED CODE TO RENAME ULTIPLE COLUMN NAMES

colnames(exprsData)<-  c("GSM6030648","GSM6030649","GSM6030650","GSM6030652","GSM6030655","GSM6030657",
                         "GSM6030658","GSM6030660","GSM6030661","GSM6030662","GSM6030663","GSM6030664",
                         "GSM6030665","GSM6030666","GSM6030667","GSM6030668","GSM6030670","GSM6030671",
                         "GSM6030672","GSM6030673","GSM6030674","GSM6030675","GSM6030676","GSM6030677",
                         "GSM6030678","GSM6030679","GSM6030680","GSM6030681","GSM6030682","GSM6030683",
                         "GSM6030685","GSM6030687","GSM6030688","GSM6030689","GSM6030690","GSM6030691",
                         "GSM6030692","GSM6030693","GSM6030694","GSM6030695","GSM6030696","GSM6030697",
                         "GSM6030698","GSM6030699","GSM6030700","GSM6030701","GSM6030702","GSM6030704",
                         "GSM6030705","GSM6030707","GSM6030708","GSM6030709","GSM6030711","GSM6030713",
                         "GSM6030714","GSM6030715","GSM6030716","GSM6030717","GSM6030718","GSM6030719",
                         "GSM6030720","GSM6030721","GSM6030722","GSM6030723","GSM6030724","GSM6030725","GSM6030726")

# Assuming exprsData is your data frame

# Sort columns alphabetically (ascending order)
exprsData <- exprsData[, sort(colnames(exprsData))]

#featuredata
gse <- getGEO( "GSE200306", GSEMatrix = TRUE)

#Get feature data

featureData <- gse$GSE200306_series_matrix.txt.gz@featureData@data

genes_to_remove <- c("BLNK", "C3", "C8B", "DPP4", "KIT")

featureData <- featureData[!rownames(featureData) %in% genes_to_remove, ]

#View(featureData)
#View(exprsData)## counts_data
#View(phenoData)##sample_info



# Create the ExpressionSet object

gset<-ExpressionSet(as.matrix(exprsData))
pData(gset)<-phenoData
featureData(gset) <- as(featureData,"AnnotatedDataFrame")

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }


# assign samples to groups and set up design matrix
gs <- factor(phenoData$stages)
groups <- make.names(c("stage3","stage4","control"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] 

fit <- lmFit(gset, design)  # fit linear model


# Contrast 1: Genes differentially expressed between Class III and Class IV
# Assuming "stagefactor" column names match model fit coefficients
# Assuming stage3 level is 2 and stage4 level is 3 (verify your actual levels)
stage3_vs_stage4 <- makeContrasts(contrasts=c(paste(groups[2],"-",groups[3],sep="")), levels=design)
fit_2 <- contrasts.fit(fit, stage3_vs_stage4)

fit2 <- contrasts.fit(fit, stage3_vs_stage4)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","logFC","GB_ACC","GO.Annotation","Official.Full.Name"))

sg34 <- tT[tT$adj.P.Val < 0.05, ]
sg34UP<- sg34[sg34$logFC > 1, ]
sg34DOWN <- sg34[sg34$logFC < -1, ]
# Contrast 2: Genes differentially expressed between Class IV and Healthy Control
stage4_vs_control <- makeContrasts(contrasts=c(paste(groups[1],"-",groups[3],sep="")), levels=design)

fit4 <- contrasts.fit(fit,stage4_vs_control)
fit4 <- eBayes(fit4, 0.01)
tT4 <- topTable(fit4, adjust="fdr", sort.by="B", number=250)
tT4<- subset(tT4, select=c("ID","adj.P.Val","P.Value","logFC","GB_ACC","GO.Annotation","Official.Full.Name"))

sg4c<- tT4[tT4$adj.P.Val < 0.05, ]
sg4cUP<- sg4c[sg4c$logFC > 1, ]
 sg4cDOWN <- sg4c[sg4c$logFC < -1, ]
 
# Contrast 3: Genes differentially expressed between Class III and Healthy Control
stage3_vs_control <- makeContrasts(contrasts=c(paste(groups[1],"-",groups[2],sep="")), levels=design)

fit3<- contrasts.fit(fit,stage3_vs_control)
fit3<- eBayes(fit3, 0.01)
tT3 <- topTable(fit3, adjust="fdr", sort.by="B", number=250)
tT3<- subset(tT3, select=c("ID","adj.P.Val","P.Value","logFC","GB_ACC","GO.Annotation","Official.Full.Name"))

sg3c<- tT3[tT3$adj.P.Val < 0.05, ]
sg3cUP<- sg3c[sg3c$logFC > 1, ]
sg3cDOWN <- sg3c[sg3c$logFC < -1, ]

siggenes <- rbind(sg34,sg3c,sg4c)
#to print in console
sapply(rownames(sg34), function(x) cat(x, "\n"))

#GO ONTOLOGY

gene_upregulated<-rownames(siggenes[siggenes$logFC>1,])

gene_downregulated<-rownames(siggenes[siggenes$logFC<-1])

#BP=BIOLOGICAL PROCESS,MP=MOLECULAR FUNCTION,C=CELLULAR COMPONENT,

GO_UP<- enrichGO(gene = gene_upregulated,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont="BP")

GO_DOWN<- enrichGO(gene=gene_downregulated,OrgDb = "org.Hs.eg.db",keyType="SYMBOL",ont="BP")

as.data.frame(GO_UP)
GOfit_U<-plot(barplot(GO_UP, showCategory = 20))

as.data.frame(GO_DOWN)
GOfit_D<-plot(barplot(GO_DOWN, showCategory = 20))


###visualisation


#Exploratory Graph

#Before applying mulitiple testing on the data we should examine exploratory graphs like PCA and heatmaps to assesss our data.
fig <- function(width, heigth){
  options(repr.plot.width = width, repr.plot.height = heigth)
}

#Create a PCA plot examine the variation in the data by the phenotype variable of interest.
# Creates PCA Plots
fig(12,8)
GSE200306_exprs <- Biobase::exprs(gset)
PCA <- prcomp(t(GSE200306_exprs), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     
                     Phenotype =
                       Biobase::pData(gset)$stages)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Phenotype)) +
  ggtitle("PCA plot of the GSE200306") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5,size=25, face='bold'),
        axis.text.x = element_text(size=12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=18, face='bold'),
        axis.title.y = element_text(size=18, face='bold'),
        legend.title = element_text(size=18, face='bold'),
        legend.text = element_text(size=18) ) +
  scale_color_manual(values = c("hotpink", "deepskyblue","red"))



# The following will produce basic volcano plot using limma function:

#volcanoplot(fit2, coef=1L, main=colnames(fit2)[1L], pch=20,
#  highlight=length(which(dT[,1L]!=0)), names=rep('+', nrow(fit2)))


#volcanoplot(fit3, coef=1L, main=colnames(fit3)[1L], pch=20,
  #          highlight=length(which(dT[,1L]!=0)), names=rep('+', nrow(fit3)))

#volcanoplot(fit4, coef=1L, main=colnames(fit4)[1L], pch=20,
 #           highlight=length(which(dT[,1L]!=0)), names=rep('+', nrow(fit4)))

#Volcano Plots

volcano_names <- ifelse(abs(fit2$coefficients)>=1,
                        as.character(fit2$ID), NA)

volcanoplot(fit2, coef = 1L, style = "p-value", highlight = 100,
            names = volcano_names,
            xlab = "logFC", ylab = NULL, pch=16, cex=0.35)


EnhancedVolcano(tT,
                lab = as.character(tT$ID),
                x = 'logFC',
                title=" stage3 vs stage4",
                y = 'adj.P.Val')


#stage3 vs control

volcano_names3 <- ifelse(abs(fit3$coefficients)>=1,
                        as.character(fit3$ID), NA)

volcanoplot(fit3, coef = 1L, style = "p-value", highlight = 100,
            names = volcano_names3,
            xlab = "logFC", ylab = NULL, pch=16, cex=0.35)


EnhancedVolcano(tT3,
                lab = as.character(tT3$ID),
                x = 'logFC',
                title=" stage3 vs control",
                y = 'adj.P.Val')



#stage4 vs control

volcano_names4 <- ifelse(abs(fit4$coefficients)>=1,
                         as.character(fit4$ID), NA)

volcanoplot(fit4, coef = 1L, style = "p-value", highlight = 100,
            names = volcano_names4,
            xlab = "logFC", ylab = NULL, pch=16, cex=0.35)


EnhancedVolcano(tT4,
                lab = as.character(tT4$ID),
                x = 'logFC',
                title=" stage4 vs control",
                y = 'adj.P.Val')

#Plotting a heatmap to examine the sample to sample distances and to see how well the samples cluster to stages
annotation_for_heatmap <- data.frame(Phenotype = Biobase::pData(gset)$stages)
row.names(annotation_for_heatmap) <- row.names(pData(gset))

dists <- as.matrix(dist(t(GSE200306_exprs), method = "manhattan"))
rownames(dists) <- row.names(pData(gset))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA
ann_colors <- list(
  Phenotype = c(stage3 = "hotpink", stage4 = "deepskyblue",control="red")
)
pheatmap(dists, col = (hmcol),
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE),
                           max(dists, na.rm = TRUE)),
         
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the GSE200306 samples")

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE200306", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()



# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE200306", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)


# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, stage3 vs stage4")

plotSA(fit3, main="Mean Variance trend,stage3 vs control")

plotSA(fit4,main="Mean Variance trend,stag4 vs control")
