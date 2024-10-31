###linear mixed effects models are appropriate for the study of non-independent data points
library(lme4)
library(lmerTest) # gets p-values in a straightforward way
library(ggplot2)
library(dplyr) 
library(corrplot)
library(gridExtra)
library(ggcorrplot)
library(ggfortify)
library(factoextra)
library(MuMIn) # gets Pseudo-R-squared values from GLMM 
library(broom)





####################################### 1: Load Amino acid properties, conduct analyses #######################################


aa <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')

### Hydropathy index
KyteDoolittle <- c(1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3)
### a more modern hydrophobicity index
GenAlgorithm3 <- c(7.2,-3.8,-55.1,-35.2,1.9,-1.3,-23.0,9.9,-41.2,17.3,4.6,-29.1,-47.0,-30.2,-74.4,-3.5,-2.5,8.7,3.6,-18.4)

### A standard. Taken from AAindex database
MolWeights <- c(89.1,121.2,133.1,147.1,165.2,75.1,155.2,131.2,146.2,131.2,149.2,132.1,115.1,146.2,174.2,105.1,119.1,117.1,204.2,181.2)

### Volume https://www.genome.jp/entry/aaindex:KRIW790103
### Volume and MolWeights are very strongly correlated
Volume <- c(27.5,44.6,40.0,62.0,115.5,0.0,79.0,93.5,100.0,93.5,94.1,58.7,41.9,80.7,105.0,29.3,51.3,71.5,145.5,117.3)

### Tien, M. Z.; Meyer, A. G.; Sydykova, D. K.; Spielman, S. J.; Wilke, C. O. (2013). "Maximum allowed solvent accessibilites of residues in proteins
### Theoretical values, of MaxASA
MaxASA <- c(129,167,193,225,240,104,224,197,236,201,224,195,159,223,274,155,172,174,285,263)

### Composition of amino acids in nuclear proteins (percent) (Cedano et al., 1997)
Comp <- c(8.3,1.6,4.7,6.5,2.7,6.3,2.1,3.7,7.9,7.4,2.3,3.7,6.9,4.7,8.7,8.8,5.1,5.3,0.7,2.4)

### Metabolic costs of amino acid biosynthesis in E. coli (Akashi et al. 2002 https://www.pnas.org/doi/10.1073/pnas.062526999)
Cost <- c(11.7,24.7,12.7,15.3,52,11.7,38.3,32.3,30.3,27.3,34.3,14.7,20.3,16.3,27.3,11.7,18.7,23.3,74.3,50)

### Normalized flexibility parameters (B-values), average (Vihinen et al., 1994) VINM940101
Flexibility <- c(0.984,0.906,1.068,1.094,0.915,1.031,0.95,0.927,1.102,0.935,0.952,1.048,1.049,1.037,1.008,1.046,0.997,0.931,0.904,0.929)

### Optimized relative partition energies - method D (Miyazawa-Jernigan, 1999) MIYS990105
ContactEnergy <- c(-0.02,-0.32,0.19,0.21,-0.33,-0.02,-0.02,-0.28,0.3,-0.32,-0.25,0.1,0.11,0.15,0.08,0.11,0.05,-0.23,-0.27,-0.23)

### Stats taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC11001767/#bib17, SERT-StructNet- a protein structure prediction tool 
### acid-base properties of molecules with acidic protons (pKa1), acid-base properties of molecules with basic protons (pKb2), 
### isoelectric point PH (pl4), and hydrophobicity (H). Taken from AAindex database by the authors.
pKa1 <- c(0.62,0.29,-0.9,-0.74,1.19,0.48,-0.4,1.38,-1.5,1.06,0.64,-0.78,0.12,-0.85,-2.53,-0.18,-0.05,1.08,0.81,0.26)
pKb2 <- c(2.34,1.96,1.88,2.19,1.83,2.34,1.82,2.36,2.18,2.36,2.28,2.02,1.99,2.17,2.17,2.21,2.09,2.32,2.83,2.2)
pl4 <- c(9.69,10.28,9.6,9.67,9.13,9.6,9.17,9.6,8.95,9.6,9.21,8.8,10.6,9.13,9.04,9.15,9.1,9.62,9.39,9.11)
H <- c(6,5.07,3.65,4.25,5.48,5.79,7.59,5.97,9.74,5.98,5.74,5.41,6.3,5.65,10.76,5.68,5.6,5.96,5.89,5.66)


Property_rawvalues <- as.data.frame(cbind(GenAlgorithm3, KyteDoolittle, H, Cost, MolWeights, MaxASA, Volume, Flexibility, ContactEnergy, pKa1, pKb2, pl4))
PropertyList <- list(GenAlgorithm3, KyteDoolittle, H, Cost, MolWeights, MaxASA, Volume, Flexibility, ContactEnergy, pKa1, pKb2, pl4)
names(PropertyList) <- c('GenAlgorithm3', 'KyteDoolittle', 'H', 'Cost', 'MolWeights', 'MaxASA', 'Volume', 'Flexibility', 'ContactEnergy', 'pKa1', 'pKb2', 'pl4')
colnames(Property_rawvalues) <- c('Hydrophobicity GA3', 'Hydrophobicity KD', 'Hydrophobicity H', 'Metabolic Cost', 'Molecular Weights', 'Max ASA', 'Volume', 'Flexibility', 'Contact Energy', 'pKa1 acidic', 'pKb2 basic', 'pl4 isoelectic point')



####################################### 1.1 Generate matrices of differences in aa properties, and create overall dataframe #######################################

PropertyDF <- data.frame()

for (i in 1:length(PropertyList)){
  property_name <- names(PropertyList)[i]
  property <- PropertyList[[i]]
  
  PropertyDif <- outer(property, property, '-')
  colnames(PropertyDif) <- aa
  rownames(PropertyDif) <- aa
  PropertyDif <- data.frame(PropertyDif)
  
  PropertyDifLabels <- list()
  PropertyDifList <- list()
  
  for (x in rownames(PropertyDif)){
    
    for (y in colnames(PropertyDif)){
      print(paste(x,y, sep = ''))
      PropertyDifLabels[[length(PropertyDifLabels)+1]] = paste(x,y, sep = '')
    }
    
    PropertyDif[x]
    PropertyDifList <- append(PropertyDifList, PropertyDif[x])
  }
  
  PropertyDifList <- unlist(PropertyDifList)
  MoWeightDifLabels <- unlist(PropertyDifLabels)
  
  summary(PropertyDifList)
  DF <- do.call(rbind.data.frame, Map('c', PropertyDifLabels, PropertyDifList))
  print(property_name)
  colnames(DF) <- c('PropertyDifLabels', paste(property_name, 'List', sep = ''))
  
  PropertyDF <- append(PropertyDF, DF)
}  

PropertyDF <- data.frame(PropertyDF) 
PropertyDF <- PropertyDF  %>% dplyr::select(-contains("."))
PropertyDF <- PropertyDF %>% mutate(across(-c(PropertyDifLabels), as.numeric))
### removing same to same amino acid contrasts:
PropertyDF <- PropertyDF[which(PropertyDF$GenAlgorithm3List != 0),]



####################################### 1.2 Generate correlation matrix results on raw values #######################################

### data = Property_rawvalues only

### calculate correlation p-values
p.mat <- model.matrix(~0+., data= Property_rawvalues) %>% 
  cor_pmat(use="pairwise.complete.obs")

model.matrix(~0+., data= Property_rawvalues) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag=FALSE, type="upper",
             hc.order = TRUE,
             lab=TRUE, lab_size=4, 
             p.mat = p.mat,
             insig = "pch",
             colors = c("#6D9EC1", "white", "#E46726"),
             ggtheme = ggplot2::theme_bw())



####################################### 1.3 PCA analysis on amino acid differences in properties, extract PC1 and 2 for analyses  #######################################


summary(PropertyDF[2:length(PropertyDF)])
### Pretty names for plotting
colnames(PropertyDF) <- c('PropertyDifLabels', 'Hydrophobicity GA3', 'Hydrophobicity KD', 'Hydrophobicity H', 'Metabolic Cost', 'Molecular Weights', 'Max ASA', 'Volume', 'Flexibility', 'Contact Energy', 'pKa1 acidic', 'pKb2 basic', 'pl4 isoelectic point')
                          
Property_pca <- prcomp(PropertyDF[2:length(PropertyDF)], scale = TRUE)
print(Property_pca)
summary(Property_pca)
get_eigenvalue(Property_pca)

fviz_pca_var(Property_pca, 
             col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#autoplot(Property_pca, loadings.label=TRUE)

var <- get_pca_var(Property_pca)
var$cos2
corrplot(var$cos2, is.corr = FALSE)

fviz_eig(Property_pca, col.var="blue")
fviz_cos2(Property_pca, choice = "var", axes = 1)
fviz_cos2(Property_pca, choice = "var", axes = 2)

PCA_DF <- as.data.frame(cbind(PropertyDF$PropertyDifLabels, Property_pca$x[,1:2]))

colnames(PCA_DF) <- c('PropertyDifLabels', 'PC1, Charge', 'PC2, Size')

### Finally, we replace our full Property dataframe with the Charge and Size principle components
PropertyDF <- PCA_DF





####################################### 2: Load Mutational parameters #######################################


### the standard 20 amino acids
### this is in a different order than in the rest of R code, be careful this list is replaced
aa <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')

### calculate mean GC content of all amino acid codons
AAs <- 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
Base1 <- strsplit('TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG', split = '')
Base2 <- strsplit('TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG', split = "")
Base3 <- strsplit('TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG', split = '')

GC_ratio <- list()

for (a in aa){
  positions <- unlist(gregexpr(a, AAs))
  total = 0
  GC = 0
  for (pos in positions){
    total = total + 3
    print(pos)
    if (Base1[[1]][pos] == 'C' || Base1[[1]][pos] == 'G'){
      GC = GC + 1 
    }
    if (Base2[[1]][pos] == 'C' || Base2[[1]][pos] == 'G'){
      GC = GC + 1 
    }
    if (Base3[[1]][pos] == 'C' || Base3[[1]][pos] == 'G'){
      GC = GC + 1 
    }
  }
  GC_ratio[a] <- GC/total
}
GC_ratio_df <- t(as.data.frame(GC_ratio))
colnames(GC_ratio_df) <- c('GC_ratio')

GCDif <- outer(GC_ratio_df, GC_ratio_df, '-')
GCDif <- data.frame(GCDif)
colnames(GCDif) <- aa

GCDifLabels <- list()
GCDifList <- list()

for (x in rownames(GCDif)){
  
  for (y in colnames(GCDif)){
    print(paste(x,y, sep = ''))
    GCDifLabels[[length(GCDifLabels)+1]] = paste(x,y, sep = '')
  }
  
  GCDif[x]
  GCDifList <- append(GCDifList, GCDif[x])
}
GCDifList <- unlist(GCDifList)
MoWeightDifLabels <- unlist(GCDifLabels)

summary(GCDifList)
GC_DF <- do.call(rbind.data.frame, Map('c', GCDifList, GCDifLabels))
colnames(GC_DF) <- c('GCDifList', 'GCDifLabels')
GC_DF$GCDifList <- as.numeric(GC_DF$GCDifList)



### calculate transition:transversions between amino acids
AAs <- 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
PurPyr1 <- strsplit('YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU', split = '')
PurPyr2 <- strsplit('YYYYYYYYUUUUUUUUYYYYYYYYUUUUUUUUYYYYYYYYUUUUUUUUYYYYYYYYUUUUUUUU', split = "")
PurPyr3 <- strsplit('YYUUYYUUYYUUYYUUYYUUYYUUYYUUYYUUYYUUYYUUYYUUYYUUYYUUYYUUYYUUYYUU', split = '')

pyr_to_total <- list()

for (a in aa){
  positions <- unlist(gregexpr(a, AAs))
  total = 0
  pyr = 0
  for (pos in positions){
    print(pos)
    total = total + 3
    if (PurPyr1[[1]][pos] == 'Y'){
      pyr = pyr + 1 }
    if (PurPyr2[[1]][pos] == 'Y'){
      pyr = pyr + 1}   
    if (PurPyr3[[1]][pos] == 'Y'){
      pyr = pyr + 1   
    }}
  print(a)
  print(pyr/total)
  pyr_to_total[a] <- pyr/total
}

pyr_to_total_df <- t(as.data.frame(pyr_to_total))
colnames(pyr_to_total_df) <- c('pyr_to_total_ratio')

transition_to_transversion <- outer(pyr_to_total_df, pyr_to_total_df, '-')
transition_to_transversion <- data.frame(transition_to_transversion)
colnames(transition_to_transversion) <- aa


TransitionDifLabels <- list()
TransitionDifList <- list()

for (x in rownames(transition_to_transversion)){
  
  for (y in colnames(transition_to_transversion)){
    print(paste(x,y, sep = ''))
    TransitionDifLabels[[length(TransitionDifLabels)+1]] = paste(x,y, sep = '')
  }
  
  transition_to_transversion[x]
  TransitionDifList <- append(TransitionDifList, transition_to_transversion[x])
}
TransitionDifList <- unlist(TransitionDifList)
MolWeightDifLabels <- unlist(TransitionDifLabels)

summary(TransitionDifList)
Trans_DF <- do.call(rbind.data.frame, Map('c', TransitionDifList, TransitionDifLabels))
colnames(Trans_DF) <- c('TransitionDifList', 'TransitionDifLabels')
Trans_DF$TransitionDifList <- as.numeric(Trans_DF$TransitionDifList)



### the G matrix is the minimum number of mutational steps between pairs of amino acids. 
G_matrix <- read.delim("~/Downloads/Matrices-normalized/universal_code_reformat.dat")
rownames(G_matrix) <- G_matrix[,1]
G_matrix<- G_matrix[,-1]
G_matrix<- G_matrix[,-20]
G_matrix<- G_matrix[-1,]


G_matrixLabels <- list()
G_matrixList <- list()

for (x in colnames(G_matrix)){
  
  for (y in rownames(G_matrix)){
    print(paste(x,y, sep = ''))
    G_matrixLabels <- append(G_matrixLabels, paste(x,y, sep = ''))
  }
  
  print(G_matrix[x,])
  G_matrixList <- append(G_matrixList, G_matrix[x])
  
}

for (x in colnames(G_matrix)){
  
  for (y in rownames(G_matrix)){
    print(paste(x,y, sep = ''))
    G_matrixLabels <- append(G_matrixLabels, paste(y,x, sep = ''))
  }
  
  print(G_matrix[x,])
  G_matrixList <- append(G_matrixList, G_matrix[x])
  
}

G_matrixList <- unlist(G_matrixList)
G_matrixLabels <- unlist(G_matrixLabels)

G_matrix_DF <- do.call(rbind.data.frame, Map('c', G_matrixList, G_matrixLabels))
colnames(G_matrix_DF) <- c('G_matrixList', 'G_matrixLabels')

G_matrix_DF <- G_matrix_DF[which(G_matrix_DF$G_matrixList != 'NA'),]
G_matrix_DF$G_matrixList <- as.numeric(G_matrix_DF$G_matrixList)





####################################### 3: Load substitution matrix #######################################
### Which we calculate by multiplying the entries of the exchangeability matrix by the frequency of amino acids

### We reset our amino acids here- to the standard order, rather than alphabetical by single letter codes
aa <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

matrices <- c('bird', 'insect', 'mammal', 'pfam', 'plant', 'yeast')

Q.pfam <- read.table("~/Downloads/Matrices-normalized/Q.pfam", fill = TRUE, col.names = aa )

if (is.na(Q.pfam[1,20])){
  Q.pfam <- rbind(NA, Q.pfam)
}

AA_freqs <- Q.pfam[21,]
Q.pfam.rates <- Q.pfam[1:20,]

Q.pfam.subs <- Q.pfam.rates*AA_freqs[col(Q.pfam.rates)]
rownames(Q.pfam.subs) <- aa

#### This is the part that fills in the reverse side of the reversible Q matrices
if (is.na(Q.pfam.subs[1,20])){
  Q.pfam.subs_reverse <- as.data.frame(t(Q.pfam.rates*AA_freqs[row(Q.pfam.rates)]))
  colnames(Q.pfam.subs_reverse) <- aa
  
  Q.pfam.subs[is.na(Q.pfam.subs)] <- 0
  Q.pfam.subs_reverse[is.na(Q.pfam.subs_reverse)] <- 0
  
  Q.pfam.subs.df <- cbind(Q.pfam.subs + Q.pfam.subs_reverse)
} else {
  Q.pfam.subs.df <- Q.pfam.subs
}


rownames(Q.pfam.subs.df)

Q.pfam.subs.dfLabels <- list()
Q.pfam.subs.dfList <- list()

for (x in rownames(Q.pfam.subs.df)){
  
  for (y in colnames(Q.pfam.subs.df)){
    print(paste(x,y, sep = ''))
    Q.pfam.subs.dfLabels[[length(Q.pfam.subs.dfLabels)+1]] = paste(x,y, sep = '')
  }
  
  print(Q.pfam.subs.df[x,])
  Q.pfam.subs.dfList <- append(Q.pfam.subs.dfList, Q.pfam.subs.df[x,])
  
}

Q.pfam.subs.dfList <- unlist(Q.pfam.subs.dfList)
Q.pfam.subs.dfLabels <- unlist(Q.pfam.subs.dfLabels)

summary(Q.pfam.subs.dfList)
Q_subs_DF <- do.call(rbind.data.frame, Map('c', Q.pfam.subs.dfList, Q.pfam.subs.dfLabels))
colnames(Q_subs_DF) <- c('Q.pfam.subs.dfList', 'Q.pfam.subs.dfLabels')
Q_subs_DF$Q.pfam.subs.dfList <- as.numeric(Q_subs_DF$Q.pfam.subs.dfList)





####################################### 4: Putting everything togather, running models #######################################


DF <- merge(PropertyDF, Q_subs_DF, by.x = 'PropertyDifLabels', by.y = 'Q.pfam.subs.dfLabels')
DF <- merge(DF, Trans_DF, by.x = 'PropertyDifLabels', by.y = 'TransitionDifLabels')
DF <- merge(DF, G_matrix_DF, by.x = 'PropertyDifLabels', by.y = 'G_matrixLabels')
DF <- merge(DF, GC_DF, by.x = 'PropertyDifLabels', by.y = 'GCDifLabels')
DF$AminoAcidGroup1 <- substring(DF$PropertyDifLabels, 1, 1)
DF$AminoAcidGroup2 <- substring(DF$PropertyDifLabels, 2, 2)

DF$'PC1, Charge' <- as.numeric(DF$'PC1, Charge')
DF$'PC2, Size' <- as.numeric(DF$'PC2, Size')





####################################### 4.1 Suppl. Fig. 2- check distribution of substitution values  #######################################


plot_1 <- ggplot(DF, aes(x=Q.pfam.subs.dfList))+geom_histogram(color="black", fill="white") + 
  xlab("Amino acid exchangeabilities, Q Pfam") +
  theme_bw()

plot_2 <- ggplot(DF, aes(x=log10(Q.pfam.subs.dfList)))+geom_histogram(color="black", fill="white") + 
  xlab("Amino acid exchangeabilities, Q Pfam, log transformed") +
  theme_bw()

### Distribution of exchangeability values
grid.arrange(plot_1, plot_2, ncol = 2)





####################################### 4.2 Main modelling #######################################


model_all <- lmer(log10(DF$Q.pfam.subs.dfList) ~ abs(DF$'PC1, Charge') + abs(DF$'PC2, Size') + abs(DF$GCDifList) + abs(DF$G_matrixList) + abs(DF$TransitionDifList) + (1|DF$AminoAcidGroup1) + (1|DF$AminoAcidGroup2))
summary(model_all)

model_all_interact <- lmer(log10(DF$Q.pfam.subs.dfList) ~ abs(DF$'PC1, Charge') * abs(DF$'PC2, Size') + abs(DF$GCDifList) + abs(DF$G_matrixList) + abs(DF$TransitionDifList) + (1|DF$AminoAcidGroup1) + (1|DF$AminoAcidGroup2))
summary(model_all_interact)

### model is improved with interaction term.
anova(model_all, model_all_interact)

### GC not significant in model, take it out and assess
model_noGC <- lmer(log10(DF$Q.pfam.subs.dfList) ~ abs(DF$'PC1, Charge') * abs(DF$'PC2, Size') + abs(DF$G_matrixList) + abs(DF$TransitionDifList) + (1|DF$AminoAcidGroup1) + (1|DF$AminoAcidGroup2))

###comparing the model to a model without GC, GC doesn't significantly improve model fit
anova(model_all_interact, model_noGC)

###Best Pfam Q matrix model
summary(model_noGC)
model_noGC_ranmut <- lmer(log10(DF$Q.pfam.subs.dfList) ~ abs(DF$'PC1, Charge') * abs(DF$'PC2, Size') + (1|abs(DF$G_matrixList)) + abs(DF$TransitionDifList) + (1|DF$AminoAcidGroup1) + (1|DF$AminoAcidGroup2))

r.squaredGLMM(model_noGC)
r.squaredGLMM(model_noGC_ranmut)





####################################### 4.3 Mutation steps modelling #######################################


### Try analyses on only amino acids that require a single mutational step
### Keeping GC out, as it continues to have no effect in models 

DF_1step <- DF[which(DF$G_matrixList == 1),]
model_1step <- lmer(log10(DF_1step$Q.pfam.subs.dfList) ~ abs(DF_1step$'PC1, Charge') * abs(DF_1step$'PC2, Size') + abs(DF_1step$TransitionDifList) + abs(DF_1step$GCDifList)+ (1|DF_1step$AminoAcidGroup1) + (1|DF_1step$AminoAcidGroup2))
model_1step_random <- lmer(log10(DF$Q.pfam.subs.dfList) ~ abs(DF$'PC1, Charge') * abs(DF$'PC2, Size') + abs(DF$TransitionDifList) + (1|DF$G_matrixList) + (1|DF$AminoAcidGroup1) + (1|DF$AminoAcidGroup2))

###better to include mutational steps as a fixed effect
anova(model_1step_random, model_all_interact)

summary(model_1step)


DF_2step <- DF[which(DF$G_matrixList == 2),]
model_2step <- lmer(log10(DF_2step$Q.pfam.subs.dfList) ~ abs(DF_2step$'PC1, Charge') * abs(DF_2step$'PC2, Size') + abs(DF_2step$TransitionDifList) + abs(DF_2step$GCDifList) + (1|DF_2step$AminoAcidGroup1) + (1|DF_2step$AminoAcidGroup2))
summary(model_2step)

###too small for complex models
DF_3step <- DF[which(DF$G_matrixList == 3),]
model_3step <- lmer(log10(DF_3step$Q.pfam.subs.dfList) ~ abs(DF_3step$'PC1, Charge') * abs(DF_3step$'PC2, Size') + abs(DF_3step$TransitionDifList) + abs(DF_3step$GCDifList) + (1|DF_3step$AminoAcidGroup1) + (1|DF_3step$AminoAcidGroup2))
summary(model_3step)




### removal of sites with more than 1 mutational step removes the influence of transition:transversion differences
summary(model_1step)
model_1step_nomut <- lmer(log10(DF_1step$Q.pfam.subs.dfList) ~ abs(DF_1step$'PC1, Charge') * abs(DF_1step$'PC2, Size') + (1|DF_1step$AminoAcidGroup1) + (1|DF_1step$AminoAcidGroup2))
### including GC and transition:transversion does not result in a significantly better model fit
summary(model_1step_nomut)
anova(model_1step, model_1step_nomut)

### removal of sites with more than 1 mutational step removes the influence of GC and transition:transversion differences
model_2step_nomut <- lmer(log10(DF_2step$Q.pfam.subs.dfList) ~ abs(DF_2step$'PC1, Charge') * abs(DF_2step$'PC2, Size') + (1|DF_2step$AminoAcidGroup1) + (1|DF_2step$AminoAcidGroup2))
### including GC and transition:transversion does not result in a significantly better model fit
summary(model_2step_nomut)
anova(model_2step, model_2step_nomut)

### removal of sites with more than 1 mutational step removes the influence of GC and transition:transversion differences
model_3step_nomut <- lmer(log10(DF_3step$Q.pfam.subs.dfList) ~ abs(DF_3step$'PC1, Charge') * abs(DF_3step$'PC2, Size') + (1|DF_3step$AminoAcidGroup1) + (1|DF_3step$AminoAcidGroup2))
### including GC and transition:transversion does not result in a significantly better model fit
summary(model_3step_nomut)
anova(model_3step, model_3step_nomut)

model_3step_simple <- lmer(log10(DF_3step$Q.pfam.subs.dfList) ~ abs(DF_3step$'PC1, Charge') + (1|DF_3step$AminoAcidGroup1) + (1|DF_3step$AminoAcidGroup2))


### Facet wrap plots for different mutational steps
facet_1 <-
  ggplot(DF,aes(x=abs(DF$'PC1, Charge'),y=log10(Q.pfam.subs.dfList))) + #can set color aes#,col=AminoAcidGroup1
  geom_jitter() + 
  geom_smooth(method = "lm", col = "black") +
  geom_point(col = '#73648A') + 
  facet_wrap(~G_matrixList) +
  scale_y_continuous(breaks = log(pretty(exp(log10(DF$Q.pfam.subs.dfList)))), labels = pretty(exp(log10(DF$Q.pfam.subs.dfList)))) +
  theme_bw() + xlab("PC1, Difference in charge-associated metrics") + ylab("Amino acid exchangeabilities, Q Pfam")


facet_2 <- 
  ggplot(DF,aes(x=abs(DF$'PC2, Size'),y=log10(Q.pfam.subs.dfList))) + #can set color aes#,col=AminoAcidGroup1
  geom_jitter() + 
  geom_smooth(method = "lm", col = "black") +
  geom_point(col = '#E96C1E') + 
  facet_wrap(~G_matrixList) +
  scale_y_continuous(breaks = log(pretty(exp(log10(DF$Q.pfam.subs.dfList)))), labels = pretty(exp(log10(DF$Q.pfam.subs.dfList)))) +
  theme_bw() + xlab("PC2, Difference in size-associated metrics") + ylab("Amino acid exchangeabilities, Q Pfam")

##### Supplementary figure 2
grid.arrange(facet_1, facet_2, nrow = 2)


### Neat plots
plot_2 <- ggplot(DF) +
  aes(abs(DF$'PC2, Size'), log10(Q.pfam.subs.dfList)) +
  geom_smooth(method = "lm", col = "black") +
  geom_point(aes(col = as.factor(abs(G_matrixList)))) +
  scale_color_manual(values= c("#EAC072", "#E96C1E", "#A10702")) +  
  scale_y_continuous(breaks = log(pretty(exp(log10(DF$Q.pfam.subs.dfList)))), labels = pretty(exp(log10(DF$Q.pfam.subs.dfList)))) +
  theme_bw() + 
  ylab("") + xlab("PC2, Difference in size-associated metrics") +  labs(col =  "Mutational\nSteps") 


plot_1 <- ggplot(DF) +
  aes(abs(DF$'PC1, Charge'), log10(Q.pfam.subs.dfList)) +
  geom_smooth(method = "lm", col = "black") +
  geom_point(aes(col = as.factor(abs(G_matrixList)))) +
  scale_color_manual(values= c("#D8C9F3", "#73648A", "#130E19")) +
  scale_y_continuous(breaks = log(pretty(exp(log10(DF$Q.pfam.subs.dfList)))), labels = pretty(exp(log10(DF$Q.pfam.subs.dfList)))) +
  theme_bw() + 
  ylab("Amino acid exchangeabilities, Q Pfam") + xlab("PC1, Difference in charge-associated metrics") +  labs(col =  "Mutational\nSteps") 

grid.arrange(plot_1, plot_2, ncol = 2)

### just for visualising- removing amino acid pairs separated by more than 1 mutational step.
ggplot(DF[which(DF$G_matrixList == 1),]) +
  aes('PC1, Charge', Q.PfamList) +
  geom_point(aes(col = abs(CompList))) +
  geom_smooth(method = "loess", col = "black", span = 0.8) +
  #  geom_smooth(method = "lm", formula = y ~ poly(x, 5, raw = TRUE), col = "black") +
  theme_bw()





####################################### 5. Across taxa comparisons #######################################


####################################### 5.1 Load all non-reversible matrices #######################################

aa <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

PfamNQ <- read.table("~/Downloads/Matrices-normalized/non_reversible_matrices/NQ.pfam", fill = TRUE, col.names = aa )
BirdNQ <- read.table("~/Downloads/Matrices-normalized/non_reversible_matrices/NQ.bird", fill = TRUE, col.names = aa )
InsectNQ <- read.table("~/Downloads/Matrices-normalized/non_reversible_matrices/NQ.insect", fill = TRUE, col.names = aa )
MammalNQ <- read.table("~/Downloads/Matrices-normalized/non_reversible_matrices/NQ.mammal", fill = TRUE, col.names = aa )
PlantNQ <- read.table("~/Downloads/Matrices-normalized/non_reversible_matrices/NQ.plant", fill = TRUE, col.names = aa )
YeastNQ <- read.table("~/Downloads/Matrices-normalized/non_reversible_matrices/NQ.yeast", fill = TRUE, col.names = aa )


matrices <- c('PfamNQ', 'BirdNQ', 'InsectNQ', 'MammalNQ', 'PlantNQ', 'YeastNQ')


for (matrixgroup in matrices){
  Q.pfam <- get(matrixgroup)
  outfile <- paste(matrixgroup, 'subs.df', sep = '')
  print(outfile)
  
  if (is.na(Q.pfam[1,20])){
    Q.pfam <- rbind(NA, Q.pfam)
  }
  
  AA_freqs <- Q.pfam[21,]
  Q.pfam.rates <- Q.pfam[1:20,]
  
  Q.pfam.subs <- Q.pfam.rates*AA_freqs[col(Q.pfam.rates)]
  rownames(Q.pfam.subs) <- aa
  
  #### This is the part that fills in the reverse side of the reversible Q matrices
  if (is.na(Q.pfam.subs[1,20])){
    Q.pfam.subs_reverse <- as.data.frame(t(Q.pfam.rates*AA_freqs[row(Q.pfam.rates)]))
    colnames(Q.pfam.subs_reverse) <- aa
    
    Q.pfam.subs[is.na(Q.pfam.subs)] <- 0
    Q.pfam.subs_reverse[is.na(Q.pfam.subs_reverse)] <- 0
    
    Q.pfam.subs.df <- cbind(Q.pfam.subs + Q.pfam.subs_reverse)
  } else {
    Q.pfam.subs.df <- Q.pfam.subs
  }
  assign(outfile, Q.pfam.subs.df)
}


for (matrixgroup in matrices){
  Q.pfam <- get(matrixgroup)
  outfile <- paste(matrixgroup, 'exchange.df', sep = '')
  print(outfile)
  
  if (is.na(Q.pfam[1,20])){
    Q.pfam <- rbind(NA, Q.pfam)
  }
  
  AA_freqs <- Q.pfam[21,]
  Q.pfam.rates <- as.data.frame(Q.pfam[1:20,])
  assign(outfile, Q.pfam.rates)
}


rownames(PfamNQexchange.df) <- aa
rownames(BirdNQexchange.df) <- aa
rownames(InsectNQexchange.df) <- aa
rownames(MammalNQexchange.df) <- aa
rownames(PlantNQexchange.df) <- aa
rownames(YeastNQexchange.df) <- aa


####################################### 5.2: Combine substitution matrices across groups #######################################

matrices <- c('PfamNQsubs.df', 'YeastNQsubs.df', 'InsectNQsubs.df', 'PlantNQsubs.df', 'BirdNQsubs.df', 'MammalNQsubs.df')

for (matrixgroup in matrices){
  matrix_compare <- get(matrixgroup)
  
  outfile <- paste(matrixgroup, '.vector', sep = '')
  
  rownames(matrix_compare)
  
  matrix_compareLabels <- list()
  matrix_compareList <- list()
  
  for (x in rownames(matrix_compare)){
    
    for (y in colnames(matrix_compare)){
      print(paste(x,y, sep = ''))
      matrix_compareLabels[[length(matrix_compareLabels)+1]] = paste(x,y, sep = '')
    }
    
    print(matrix_compare[x,])
    matrix_compareList <- append(matrix_compareList, matrix_compare[x,])
    
  }
  matrix_compareList <- unlist(matrix_compareList)
  matrix_compareLabels <- unlist(matrix_compareLabels)
  
  summary(matrix_compareList)
  Vector_DF <- do.call(rbind.data.frame, Map('c', matrix_compareList, matrix_compareLabels))
  colnames(Vector_DF) <- c('matrix_compareList', 'matrix_compareLabels')
  Vector_DF$matrix_compareList <- as.numeric(Vector_DF$matrix_compareList)
  
  assign(outfile, Vector_DF)
}

matrices <- c('PfamNQsubs.df', 'YeastNQsubs.df', 'InsectNQsubs.df', 'PlantNQsubs.df', 'BirdNQsubs.df', 'MammalNQsubs.df')

X = cbind(PfamNQsubs.df.vector$matrix_compareList, 
          YeastNQsubs.df.vector$matrix_compareList,
          InsectNQsubs.df.vector$matrix_compareList,
          PlantNQsubs.df.vector$matrix_compareList,
          BirdNQsubs.df.vector$matrix_compareList, 
          MammalNQsubs.df.vector$matrix_compareList)

NQsubs.df.vector_all <- cbind(YeastNQsubs.df.vector$matrix_compareLabels, X)
colnames(NQsubs.df.vector_all) <- c('AA', 'Yeast', 'Pfam', 'Insect', 'Plant', 'Bird', 'Mammal')
NQsubs.df.vector_all <- data.frame(NQsubs.df.vector_all)
NQsubs.df.vector <- as.data.frame(lapply(NQsubs.df.vector_all[2:7],as.numeric))
NQsubs.df.vector <- cbind(NQsubs.df.vector_all$AA, NQsubs.df.vector)
NQsubs.df.vector <- NQsubs.df.vector[which(NQsubs.df.vector$Pfam >= 0),]
colnames(NQsubs.df.vector) <- c('AA', 'Yeast', 'Pfam', 'Insect', 'Plant', 'Bird', 'Mammal')

#########  comparative plots- pairwise with whichever taxonomic groups are of interest
plot(NQsubs.df.vector$Mammal, NQsubs.df.vector$Yeast, pch = 20, abline(a=0, b=1))
text(NQsubs.df.vector$Mammal+0.001, NQsubs.df.vector$Yeast, labels =  NQsubs.df.vector$AA, cex = 0.7)


#######################################  5.3: Combine substitution matrices with property and mutation statistics #######################################

TaxaDF <- merge(PropertyDF, NQsubs.df.vector, by.x = 'PropertyDifLabels', by.y = 'AA')
TaxaDF <- merge(TaxaDF, Trans_DF, by.x = 'PropertyDifLabels', by.y = 'TransitionDifLabels')
TaxaDF <- merge(TaxaDF, G_matrix_DF, by.x = 'PropertyDifLabels', by.y = 'G_matrixLabels')
TaxaDF <- merge(TaxaDF, GC_DF, by.x = 'PropertyDifLabels', by.y = 'GCDifLabels')
TaxaDF$AminoAcidGroup1 <- substring(TaxaDF$PropertyDifLabels, 1, 1)
TaxaDF$AminoAcidGroup2 <- substring(TaxaDF$PropertyDifLabels, 2, 2)

TaxaDF$'PC1, Charge' <- as.numeric(DF$'PC1, Charge')
TaxaDF$'PC2, Size' <- as.numeric(DF$'PC2, Size')


#######################################  5.4: Calculate statistics over all matrices  #######################################

### Run models on every matrix
for (i in colnames(TaxaDF)){
  if (i %in% c('Pfam', 'Yeast', 'Insect', 'Plant',  'Bird', 'Mammal') ) {
    
    Taxamodel_all <- lmer(log10(TaxaDF[[i]]) ~ abs(TaxaDF$'PC1, Charge') * abs(TaxaDF$'PC2, Size') + abs(TaxaDF$GCDifList) + abs(TaxaDF$G_matrixList) + abs(TaxaDF$TransitionDifList)
                      + (1|TaxaDF$AminoAcidGroup1) + (1|TaxaDF$AminoAcidGroup2))
    
    Taxamodel_all_noGC <- lmer(log10(TaxaDF[[i]]) ~ abs(TaxaDF$'PC1, Charge') * abs(TaxaDF$'PC2, Size') + abs(TaxaDF$G_matrixList) + abs(TaxaDF$TransitionDifList)
                           + (1|TaxaDF$AminoAcidGroup1) + (1|TaxaDF$AminoAcidGroup2))
    
    print(anova(Taxamodel_all, Taxamodel_all_noGC)) ### Including GC does not improve model fit, for any dataset
    
    Taxamodel_ranmut <- lmer(log10(TaxaDF[[i]]) ~ abs(TaxaDF$'PC1, Charge') * abs(TaxaDF$'PC2, Size') + (1|abs(TaxaDF$G_matrixList)) + abs(TaxaDF$TransitionDifList) 
                         + (1|TaxaDF$AminoAcidGroup1) + (1|TaxaDF$AminoAcidGroup2))
    
    print(summary(Taxamodel_all_noGC))
    print(r.squaredGLMM(Taxamodel_all_noGC))
    print(r.squaredGLMM(Taxamodel_ranmut))
  }
}

