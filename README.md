Hi 
The repository contain files related to capstone project on Identification of potential drug target for Lupus Nephritis such as my workplan documents , code source , results and conclusions obtained from analysis


1. Disease & Population
Lupus Nephritis is a systemic autoimmune disease. Systemic means that it can potentially affect 
all body organs and tissues. Autoimmune means that the immune system instead of fighting 
infections and tumours, it starts to attack body tissues. Genetic problems, sunlight exposure, 
hormonal reasons, and infections may be the reason for the immune system dysfunction.

Autoimmune B-lymphocytes produce autoantibodies such as antinuclear antibodies (ANA), 
anti-dsDNA, anti-Smith and antiphospholipid antibodies. These are often used by doctors to 
help make the diagnosis of lupus. 

Autoantibodies bind self-particles like DNA, and cell proteins, The complexes they form 
(immune complexes) travel in the blood circulation and go to different organs. There they cause 
inflammation by different mechanisms including the activation of "complement factors" C4 
and C3. 

Systematic lupus erythematosus, or SLE, is an autoimmune disease that can affect nearly every 
organ system and has a wide range of clinical presentations. One of the most dangerous side 
effects of SLE is kidney involvement. There is still much to learn about the mechanisms of 
action underlying SLE and its complex etiopathogenesis. The risk of SLE is an excessive 
generation of radioactive debris as a result of aberrant and massive apoptotic events that are 
mistaken for foreign bodies. This abnormal antigen presentation eventually results in a loss of 
tolerance to B and T cells. Due to this loss of tolerance, T cells become hyperactive, which in 
turn causes the production of inflammatory cytokines. B cells also become hyperactive, which 
results in the massive production of autoantibodies and the development of immune response


Classification of Lupus Nephritis :

• Class I. Minimal mesangial lupus nephritis (LN) 

• Class II. Mesangial proliferative LN 

• Class III. Focal LN (less than 50% of glomeruli are involved) 

• Class IV. Diffuse LN (more than 50% of glomeruli are involved)

• Class V. Membranous LN 

• Class VI. Advanced sclerosing LN

A patient might have either: 

• Class III/IV LN, or Class V LN or 

• both Class III/IV+Class V 

• Class III or IV nephritis (with or without Class V) is worse than Class V LN because 
it often affects kidney function. It requires more aggressive treatment 

• Class V nephritis (membranous) without Class III/IV LN causes proteinuria without 
affecting kidney function in the beginning. 

• However, Class V can cause severe (nephrotic) proteinuria, with often very bad leg 
edema, and high cholesterol, and may cause a blood clot. This needs to be treated 
more aggressively. 


DATA ACQUISTION

"Data for this project was obtained from the NCBI Gene Expression Omnibus (GEO) public 
repository. We searched GEO using keywords ‘lupus nephritis’, ‘chronic kidney disease’ , 
‘Homo Sapiens’and certain criteria such as a sample size of above 50 and gene expression 
profile. The selected dataset (GEO Series accession number: GSE200306) profiled gene 
expression in glomerular tissues and controls. The series matrix file was downloaded for 
analysis and normalised glomeruli expression data was also obtained from database.

The overall design used in the dataset was 58 paired kidney biopsies from patients with 
proliferative lupus nephritis (class III, class IV or class III/IV+V) that were initially 
evaluated. This is a well-defined, well-curated dataset with clinical and histologic data 
available for all patients. Glomeruli and tubulointerstitium (TI) were isolated separately and 
samples with adequate RNA were submitted for nanostring analysis. Overall, 70 glomerular 
samples and 92 TI samples were submitted for nanostring analysis (35 glomerular and 46 TI 
pairs). Pre-implantation donor kidney biopsies (n=10) served as healthy controls.

Data processing:For nanostring data, raw counts were normalized to the positive spike-in 
controls and then log2 transformed. To reduce technical bias, genes with an expression level 
below the mean plus two standard deviations of the negative controls were filtered out. Quantile 
normalization was applied to the remaining transcripts (522 for glomeruli and 502 for TI). 
Linear mixed-effect models were used to identify differentially expressed genes by taking into 
account the correlation between repeated measures before and after treatment. To improve the 
stability of variance estimation, moderated t-tests were employed for differential expression 
detection.

PROJECT MANAGEMENT:

The data chosen has 181 samples out of which 79 glomeruli tissue samples were taken for this 
analysis 28 class 3 Lupus Nephritis and 31 class4 Lupus Nephritis and 9 control samples of 
glomeruli tissues were further considered for Data preprocessing and exploration to find 
potential targets for Lupus Nephritis by finding genes up or downregulated in the pathway 
leading to calcineurin activation. A phosphatase called calcineurin is involved in the synthesis 
of inflammatory mediators by T cells as well as the upkeep and appropriate operation of the 
glomerular filtration membrane, both of which are critical for the development of lupus 
nephritis. As a result, calcineurin inhibitors (CNIs) are a promising tactic for blocking T cells, 
and in recent years, CNI-based research and clinical trials have grown dramatically. The use 
and mechanisms of action of CNIs as a therapeutic approach for LN are the main topics of this 
review, along with the role of calcineurin in the pathogenesis of SLE. 

1. Data preprocessing
   
▪ Data Manipulation:The following R packages were used:

• Biobase 

• dplyr 

• tidyr 

• GEOquery.

▪ Data Normalization:

• limma 

• Biobase 

▪ Data Visualization:

• ggplot2 

• RColorBrewer 

• geneplotter 

• pheatmap.

• enrichplot 

• EnhancedVolcano 

• umap 

Others used:

clusterProfiler (analysis): Designed for functional enrichment analysis of gene 
sets.

org.Hs.eg.db (annotation): Provides gene annotations for the human genome.
AnnotationDbi (database access): Facilitates access to annotation databases.

 Exploratory data analysis:
   
• Volcano Plots

• Pheatmaps

• Functional enrichment analysis

Differential gene expression analysis:
   
• linear model statistical testing

• multiple testing through FDR control
