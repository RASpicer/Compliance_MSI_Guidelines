# MSI Compliance Analysis
Analysis by Rachel Spicer, github:RASpicer  
2017/08/05  



# Data Processing
This RMarkdown contains the code used for analysis for the 2017 paper <b>Compliance with Minimum Information Guidelines in Public Metabolomics Repositories</b>.

## Functions

This section includes code for extracting only the minimal or best practice/ optional reporting standards.


```r
# Calculate the percentage compliance with the minimal reporting standards for clinical data
MinPerClin <- function(x) {
  # Calculate totals
  Totals = as.data.frame(colSums(x, na.rm = TRUE))
  # Extract only minimal reporting standards
  MinRep = c(Totals[1:7,], Totals[14:22,], Totals[42:44,], Totals[46,], Totals[48,], Totals[55,])
  # Calculate percentage
  Per = MinRep/(nrow(x)) * 100
  Per = round(Per, digits = 2)
  # Get column names
  Names = c(colnames(x[,1:7]), colnames(x[,14:22]), colnames(x[,42:44]), colnames(x[46]), colnames(x[48]), colnames(x[55]))
  MinMeta = as.data.frame(cbind(Names, Per))
  colnames(MinMeta) = c("Metadata", "Percentage")
  MinMeta$Percentage = as.numeric(as.character(MinMeta$Percentage))
  # Return percentage
  return(MinMeta)
}

# Calculate the percentage compliance with the minimal reporting standards for preclinical data
MinPerPreClin <- function(x,y) {
  # Calculate totals
  Totals = as.data.frame(colSums(x, na.rm = TRUE))
  # Extract only minimal reporting standards
  MinRep = c(Totals[1,], Totals[3:7,], Totals[14:17,], Totals[22:23,], Totals[25:30,], Totals[32:34,], Totals[36:37,], Totals[39:45,])
  # Calculate percentage
  Per = MinRep/(nrow(x)) * 100
  Per = round(Per, digits = 2)
  # Get column names
  Names = c(colnames(x[1]), colnames(x[,3:7]), colnames(x[,14:17]), colnames(x[,22:23]), colnames(x[,25:30]), colnames(x[,32:34]), colnames(x[,36:37]), colnames(x[,39:45]))
  MinMeta = as.data.frame(cbind(Names, Per))
  colnames(MinMeta) = c("Metadata", "Percentage")
  MinMeta$Percentage = as.numeric(as.character(MinMeta$Percentage))
  # Return percentage
  return(MinMeta)
}

# Calculate the percentage compliance with the minimal reporting standards for microbial/cell line data
MinPerCellLine <- function(x) {
  # Calculate totals
  Totals = as.data.frame(colSums(x, na.rm = TRUE))
  # Extract only minimal reporting standards
  MinRepSt = Totals[1:15,]
  # Calculate percentage
  Per = MinRepSt/(nrow(x)) * 100
  Per = round(Per, digits = 2)
  # Get column names
  Names = c(colnames(x[,1:15]))
  MinMeta = as.data.frame(cbind(Names, Per))
  colnames(MinMeta) = c("Metadata", "Percentage")
  MinMeta$Percentage = as.numeric(as.character(MinMeta$Percentage))
  return(MinMeta)
}

# Calculate the percentage compliance with the minimal reporting standards for plant studies
MinPerPlant <- function(x) {
  # Calculate totals
  Totals = as.data.frame(colSums(x, na.rm = TRUE))
  # Calculate percentage
  Per = Totals/(nrow(x)) * 100
  Per = round(Per, digits = 2)
  # Get column names
  Names = colnames(x)
  MinMeta = as.data.frame(cbind(Names, Per))
  colnames(MinMeta) = c("Metadata", "Percentage")
  MinMeta$Percentage = as.numeric(as.character(MinMeta$Percentage))
  # Return percentage
  return(MinMeta)
}

# Calculate the percentage compliance with the optional reporting standards for microbial/cell line data
OpPerCellLine <- function(x) {
  # Calculate totals
  Totals = as.data.frame(colSums(x, na.rm = TRUE))
  # Extract only optional reporting standards
  OpRep = Totals[16:54,]
  # Calculate percentage
  Per = OpRep/(nrow(x)) * 100
  Per = round(Per, digits = 2)
  # Get column names
  Names = c(colnames(x[,16:54]))
  OpMeta = as.data.frame(cbind(Names, Per))
  colnames(OpMeta) = c("Metadata", "Percentage")
  OpMeta$Percentage = as.numeric(as.character(OpMeta$Percentage))
  # Return percentage
  return(OpMeta)
}

# Calculate the percentage compliance with the optional reporting standards for clinical data
OpPerClin <- function(x) {
  # Calculate totals
  Totals = as.data.frame(colSums(x, na.rm = TRUE))
  OpRep = c(Totals[8:13,], Totals[23:41,], Totals[45,], Totals[47,], Totals[49:54,])
  # Calculate percentage
  Per = OpRep/(nrow(x)) * 100
  Per = round(Per, digits = 2)
  # Get column names
  Names = c(colnames(x[,8:13]), colnames(x[,23:41]), colnames(x[45]), colnames(x[47]), colnames(x[,49:54]))
  OpMeta = as.data.frame(cbind(Names, Per))
  colnames(OpMeta) = c("Metadata", "Percentage")
  OpMeta$Percentage = as.numeric(as.character(OpMeta$Percentage))
  # Return percentage
  return(OpMeta)
}

# Calculate the percentage compliance with the optional reporting standards for preclinical data
OpPerPreClin <- function(x) {
  # Calculate totals
  Totals = as.data.frame(colSums(x, na.rm = TRUE))
  # Extract only minimal reporting standards
  OpRep = c(Totals[2,], Totals[8:13,], Totals[18:21,], Totals[24,], Totals[31,], Totals[35,], Totals[38,], Totals[46,])
  # Calculate percentage
  Per = OpRep/(nrow(x)) * 100
    Per = round(Per, digits = 2)
  # Get column names
  Names = c(colnames(x[2]), colnames(x[,8:13]), colnames(x[,18:21]),colnames(x[24]), colnames(x[31]), colnames(x[35]), colnames(x[38]), colnames(x[46]))
  OpMeta = as.data.frame(cbind(Names, Per))
  colnames(OpMeta) = c("Metadata", "Percentage")
  OpMeta$Percentage = as.numeric(as.character(OpMeta$Percentage))
  # Return percentage
  return(OpMeta)
}

# ggplot theme
theme_adjbw <- function(base_size = 12, base_family = ""){
  theme_bw() +
    theme(
      axis.text = element_text(colour = "black"),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5),
      # Remove gridlines and borders
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
}
```

## Data Processing

Data for the compliance to one of the four sets of reporting standards is first split by repository. The percentage of studies that comply with each reporting standard is then calculated. Reporting standards are then divided into minimal and best practice/optional.


```r
# read compliance data for each study type
Mouse <- read.xlsx("../MousePreclinicalRawData4Analysis.xlsx", sheetIndex = 1, 
    check.names = FALSE)
Plant <- read.xlsx("../ArabidopsisPlantRawData4Analysis.xlsx", sheetIndex = 1, 
    check.names = FALSE)
HumanClin <- read.xlsx("../HumanClinicalRawData4Analysis.xlsx", sheetIndex = 1, 
    check.names = FALSE)
HumanCell <- read.xlsx("../HumanCellRawData4Analysis.xlsx", sheetIndex = 1, 
    check.names = FALSE)

# Extract only MetaboLights studies
MLMouse <- Mouse[grep("MetaboLights", Mouse$`Study Repository `), ]
MLPlant <- Plant[grep("MetaboLights", Plant$`Study Repository `), ]
MLHumanClin <- HumanClin[grep("MetaboLights", HumanClin$`Study Repository`), 
    ]
MLHumanCell <- HumanCell[grep("MetaboLights", HumanCell$`Study Repository `), 
    ]

# Extract only Metabolomics Workbench studies
MWMouse <- Mouse[grep("Metabolomics Workbench", Mouse$`Study Repository `), 
    ]
MWPlant <- Plant[grep("Metabolomics Workbench", Plant$`Study Repository `), 
    ]
MWHumanClin <- HumanClin[grep("Metabolomics Workbench", HumanClin$`Study Repository`), 
    ]
MWHumanCell <- HumanCell[grep("Metabolomics Workbench", HumanCell$`Study Repository `), 
    ]

# Extract only MetaPhen studies
MPPlant <- Plant[grep("MetaPhen", Plant$`Study Repository `), ]

# Extract only MeRy-B studies
MeRPlant <- Plant[grep("MeRy-B", Plant$`Study Repository `), ]

# Minimal - MetaboLights
MLHumanCellP <- MinPerCellLine(MLHumanCell[, -(1:2)])
MLHumanCellP$Order <- rep(1)
MLHumanClinP <- MinPerClin(MLHumanClin[, -(1:2)])
MLHumanClinP$Order <- rep(2)
MLMouseP <- MinPerPreClin(MLMouse[, -(1:2)])
MLMouseP$Order <- rep(3)
MLPlantP <- MinPerPlant(MLPlant[, -(1:2)])
MLPlantP$Order <- rep(4)

# Minimal - Metabolomics Workbench
MWHumanCellP <- MinPerCellLine(MWHumanCell[, -(1:2)])
MWHumanCellP$Order <- rep(1)
MWHumanClinP <- MinPerClin(MWHumanClin[, -(1:2)])
MWHumanClinP$Order <- rep(2)
MWMouseP <- MinPerPreClin(MWMouse[, -(1:2)])
MWMouseP$Order <- rep(3)
MWPlantP <- MinPerPlant(MWPlant[, -(1:2)])
MWPlantP$Order <- rep(4)

# Minimal - MetaPhen
MPPlantP <- MinPerPlant(MPPlant[, -(1:2)])
MPPlantP$Order <- rep(1)

# Minimal - MeRy-B
MeRPlantP <- MinPerPlant(MeRPlant[, -(1:2)])
MeRPlantP$Order <- rep(4)

# Best Practice - MetaboLights
MLHumanCellOP <- OpPerCellLine(MLHumanCell[, -(1:2)])
MLHumanCellOP$Order <- rep(1)
MLHumanClinOP <- OpPerClin(MLHumanClin[, -(1:2)])
MLHumanClinOP$Order <- rep(2)
MLMouseOP <- OpPerPreClin(MLMouse[, -(1:2)])
MLMouseOP$Order <- rep(3)

# Best Practice - Metabolomics Workbench
MWHumanCellOP <- OpPerCellLine(MWHumanCell[, -(1:2)])
MWHumanCellOP$Order <- rep(1)
MWHumanClinOP <- OpPerClin(MWHumanClin[, -(1:2)])
MWHumanClinOP$Order <- rep(2)
MWMouseOP <- OpPerPreClin(MWMouse[, -(1:2)])
MWMouseOP$Order <- rep(3)
```

# Figures
Code that was used to generate the raw figures from the paper. The raw figures were then edited in Adobe Illustrator to add letters to indicate statistical significance and to change colours to patterns in figure 2.

## Figure 1a. Compliance to Minimal MSI Reporting Standards in MetaboLights
Combined box-and-whisker and dot plot showing the percentage compliance with the MSI minimum reporting standards within the MetaboLights repository. and b) Metabolomics Workbench repositories. Each “dot” represents compliance to a single reporting standard.


```r
ML <- rbind(MLHumanCellP, MLHumanClinP, MLMouseP, MLPlantP)
ML$Order <- as.factor(ML$Order)
ggplot(ML, aes(x = Order, y = Percentage)) + geom_boxplot(fill = "#7fc97f", 
    outlier.shape = 3, outlier.color = "red") + geom_jitter() + scale_x_discrete(labels = c(expression(italic("In vitro")), 
    "Clinical", "Pre-clinical", "Plant")) + scale_y_continuous(expand = c(0, 
    0), limits = c(-1, 101)) + theme_adjbw() + xlab("MSI Guideline") + ylab("Compliance (%)")
```

<img src="figs/MLMin-1.png" style="display: block; margin: auto;" />

## Figure 1b. Compliance to Minimal MSI Reporting Standards in Metabolomics Workbench
Combined box-and-whisker and dot plot showing the percentage compliance with the MSI minimum reporting standards within the Metabolomics Workbench repository. Red ‘+’ indicate outliers, each “dot” represents compliance to a single reporting standard.


```r
MW <- rbind(MWHumanCellP, MWHumanClinP, MWMouseP, MWPlantP)
MW$Order <- as.factor(MW$Order)
ggplot(MW, aes(x = Order, y = Percentage)) + geom_boxplot(fill = "#beaed4", 
    outlier.shape = 3, outlier.color = "red") + geom_jitter() + scale_x_discrete(labels = c(expression(italic("In vitro")), 
    "Clinical", "Pre-clinical", "Plant")) + theme_adjbw() + scale_y_continuous(expand = c(0, 
    0), limits = c(-1, 101)) + xlab("MSI Guideline") + ylab("Compliance (%)")
```

<img src="figs/MWMin-1.png" style="display: block; margin: auto;" />

## Figure 2a. Comparison of Compliance to Minimal and Best Practice Reporting Standards in MetaboLights
Box-and-whisker plot showing the percentage compliance to the minimal and optional/best practice reporting standards in MetaboLights. Minimal reporting standards are indicated in green, optional/best practice reporting standards are blue and red ‘+’ indicate outliers.


```r
# Compare minimal and optional reporting standards within databases
# MetaboLights
MLMin <- rbind(MLHumanCellP, MLHumanClinP, MLMouseP)
MLO <- rbind(MLHumanCellOP, MLHumanClinOP, MLMouseOP)
MLMin$Req <- rep("Minimal")
MLO$Req <- rep("Optional")
MLComp <- rbind(MLMin, MLO)
MLComp$Order <- as.factor(MLComp$Order)
# Plot box and Whisker diagram
ggplot(MLComp, aes(x=Order, y=Percentage, fill = Req))  + 
  geom_boxplot(outlier.shape = 3, outlier.color = "red")  +
  scale_x_discrete(labels = c(expression(italic("In vitro")),
                              "Clinical",
                              "Pre-clinical",
                              "Plant")) +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    # Remove gridlines and borders
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.position="top") +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 100)) +
  xlab("MSI Guideline") +
  ylab("Compliance (%)") +
  scale_fill_manual(values = c("#7fc97f", "#a6cee3"))
```

<img src="figs/MLCom-1.png" style="display: block; margin: auto;" />

## Figure 2b.Comparison of Compliance to Minimal and Best Practice Reporting Standards in Metabolomics Workbench
Box-and-whisker plot showing the percentage compliance to the minimal and optional/best practice reporting standards in Metabolomics Workbench. Minimal reporting standards are indicated in purple, optional/best practice reporting standards are pink and red ‘+’ indicate outliers.


```r
# Compare minimal and optional reporting standards within databases
# MetaboLights
MWMin <- rbind(MWHumanCellP, MWHumanClinP, MWMouseP)
MWO <- rbind(MWHumanCellOP, MWHumanClinOP, MWMouseOP)
MWMin$Req <- rep("Minimal")
MWO$Req <- rep("Optional")
MWComp <- rbind(MWMin, MWO)
MWComp$Order <- as.factor(MWComp$Order)
# Plot box and Whisker diagram
ggplot(MWComp, aes(x=Order, y=Percentage, fill = Req))  + 
  geom_boxplot(outlier.shape = 3, outlier.color = "red")  +
  scale_x_discrete(labels = c(expression(italic("In vitro")),
                              "Clinical",
                              "Pre-clinical",
                              "Plant")) +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    # Remove gridlines and borders
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.position="top") +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 100)) +
  xlab("MSI Guideline") +
  ylab("Compliance (%)") +
  scale_fill_manual(values = c("#beaed4", "#fb9a99"))
```

<img src="figs/MWCom-1.png" style="display: block; margin: auto;" />

## Supplementary Figure 1a. Compliance to Best Practice MSI Reporting Standards in MetaboLights
Combined box-and-whisker and dot plot showing the percentage compliance with the MSI best practice/optional reporting standards within the MetaboLights repository. Red `+' indicate outliers. 


```r
MLO <- rbind(MLHumanCellOP, MLHumanClinOP, MLMouseOP)
MLO$Order <- as.factor(MLO$Order)
ggplot(MLO, aes(x = Order, y = Percentage)) + geom_boxplot(fill = "#7fc97f", 
    outlier.shape = 3, outlier.color = "red") + geom_jitter() + scale_x_discrete(labels = c(expression(italic("In vitro")), 
    "Clinical", "Pre-clinical", "Plant")) + scale_y_continuous(expand = c(0, 
    0), limits = c(-1, 101)) + theme_adjbw() + xlab("MSI Guideline") + ylab("Compliance (%)")
```

<img src="figs/MLOp-1.png" style="display: block; margin: auto;" />

## Supplementary Figure 1b. Compliance to Best Practice MSI Reporting Standards in Metabolomics Workbench
Combined box-and-whisker and dot plot showing the percentage compliance with the MSI best practice/optional reporting standards within the Metabolomics Workbench Repository. Red `+' indicate outliers.


```r
MWO <- rbind(MWHumanCellOP, MWHumanClinOP, MWMouseOP)
MWO$Order <- as.factor(MWO$Order)
ggplot(MWO, aes(x = Order, y = Percentage)) + geom_boxplot(fill = "#beaed4", 
    outlier.shape = 3, outlier.color = "red") + geom_jitter() + scale_x_discrete(labels = c(expression(italic("In vitro")), 
    "Clinical", "Pre-clinical", "Plant")) + theme_adjbw() + scale_y_continuous(expand = c(0, 
    0), limits = c(-1, 101)) + xlab("MSI Guideline") + ylab("Compliance (%)")
```

<img src="figs/MWOp-1.png" style="display: block; margin: auto;" />

## Supplementary Figure 2. Species per Repository

The frequency of studies including different species in the metabolomics data repositories: MeRy-B, MetaboLights, Metabolomics Workbench and MetaPhen. For a species to be plotted as an individual band it must have been found in a minimum of nine studies across the repositories. Species found in less than nine studies across the four repositories are reported as Other Species.

```r
# Load csv
SpeciesbyRepository <- read.csv("../SpeciesperRepository.csv")
# Replace all spaces with line graps for a cleaner plot
levels(SpeciesbyRepository$Repository) <- gsub(" ", "\n", levels(SpeciesbyRepository$Repository))
ggplot(SpeciesbyRepository, aes(x=Repository, y=Frequency, fill=Species))  +
  geom_bar(colour="black", stat="identity", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 400)) +
  scale_fill_manual(values = c('#80b1d3','#bebada','#b3de69',"#8dd3c7",
                                          '#fb8072',           
                                          '#ffffb3',
                                          '#fdb462',
                                          '#fccde5',
                                          '#d9d9d9'),
                               labels = c(expression(italic("Arabidopsis thaliana")), 
                                          expression(italic("Escherichia coli")), 
                                          expression(italic("Homo sapiens")),
                                          expression(italic("Lycopersicon esculentum")),
                                          expression(italic("Mus musculus")),
                                          expression(italic("Oryza sativas")),
                                          "Other Species",
                                          expression(italic("Rattus norvegicus")),
                                          expression(italic("Solanum lycopersicum"))))  +
  ylab("Number of Studies") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"),
        legend.text.align = 0,
        #legend.text = element_text(face = "italic"),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
```

<img src="figs/SpeciesperRep-1.png" style="display: block; margin: auto;" />

# Statistical Analysis
For analysing differences within repositories (ML and MW) Kruskal Wallis tests, followed by Dunn post-hoc tests with Benjamini-Hochberg correction were used. A Kruskal Wallis test was also used for comparing differences between repositories, both between ML and MW and between all four repositories for the plant guidelines. Mann-Whitney U tests were used to analyse difference between minimal and best practice reporting standards. 


```r
# Between MetaboLights minimal
KruskalML <- kruskal.test(Percentage ~ Order, data = ML)
DunnML <- dunnTest(Percentage ~ Order, data = ML, method = "bh")

# Between Metabolomics Workbench minimal
KruskalMW <- kruskal.test(Percentage ~ Order, data = MW)
DunnMW <- dunnTest(Percentage ~ Order, data = MW, method = "bh")

# Between MetaboLights best practice/optional
KruskalMLO <- kruskal.test(Percentage ~ Order, data = MLO)
DunnMLO <- dunnTest(Percentage ~ Order, data = MLO, method = "bh")

# Between Metabolomics Workbench best practice/optional
KruskalMWO <- kruskal.test(Percentage ~ Order, data = MWO)
DunnMWO <- dunnTest(Percentage ~ Order, data = MWO, method = "bh")

# Between repositories for plant guidelines only Add repostitory to plant
# data
MLPlantP$Repository <- rep("MetaboLights")
MWPlantP$Repository <- rep("Metabolomics Workbench")
MPPlantP$Repository <- rep("MetaPhen")
MeRPlantP$Repository <- rep("MeRy-B")
Plant <- rbind(MLPlantP, MWPlantP, MPPlantP, MeRPlantP)
Plant$Repository <- as.factor(Plant$Repository)
KruskalPlant <- kruskal.test(Percentage ~ Repository, data = Plant)
DunnPlant <- dunnTest(Percentage ~ Repository, data = Plant, method = "bh")

# Compare Between repositories - MetaboLights and Metabolomics Workbench
# Kruskal Wallis - Minimal
ML$Repository <- rep("MetaboLights")
MW$Repository <- rep("Metablomics Workbench")
MLMW <- rbind(ML, MW)
MLMW$Repository <- as.factor(MLMW$Repository)
KruskalMLMW <- kruskal.test(Percentage ~ Repository, data = MLMW)

# Kruskal Wallis - Optional
MLO$Repository <- rep("MetaboLights")
MWO$Repository <- rep("Metablomics Workbench")
MLMWO <- rbind(MLO, MWO)
MLMWO$Repository <- as.factor(MLMWO$Repository)
KruskalMLMWO <- kruskal.test(Percentage ~ Repository, data = MLMWO)

# Between Minimum and Optional MetaboLights Human - Cell Line
MLHumanCellP$MinBP <- rep("Min")
MLHumanCellOP$MinBP <- rep("BP")
MLHumanCellCom <- rbind(MLHumanCellP, MLHumanCellOP)
MLHumanCellCom$MinBP <- as.factor(MLHumanCellCom$MinBP)
MannWhitMLHumanCell <- wilcox.test(Percentage ~ MinBP, data = MLHumanCellCom)
# Human - Clinical
MLHumanClinP$MinBP <- rep("Min")
MLHumanClinOP$MinBP <- rep("BP")
MLHumanClinCom <- rbind(MLHumanClinP, MLHumanClinOP)
MLHumanClinCom$MinBP <- as.factor(MLHumanClinCom$MinBP)
MannWhitMLHumanClin <- wilcox.test(Percentage ~ MinBP, data = MLHumanClinCom)
# Mouse - Pre-clinical
MLMouseP$MinBP <- rep("Min")
MLMouseOP$MinBP <- rep("BP")
MLMouseCom <- rbind(MLMouseP, MLMouseOP)
MLMouseCom$MinBP <- as.factor(MLMouseCom$MinBP)
MannWhitMLMouse <- wilcox.test(Percentage ~ MinBP, data = MLMouseCom)

# Metabolomics Workbench Human - Cell Line
MWHumanCellP$MinBP <- rep("Min")
MWHumanCellOP$MinBP <- rep("BP")
MWHumanCellCom <- rbind(MWHumanCellP, MWHumanCellOP)
MWHumanCellCom$MinBP <- as.factor(MWHumanCellCom$MinBP)
MannWhitMWHumanCell <- wilcox.test(Percentage ~ MinBP, data = MWHumanCellCom)
# Human - Clinical
MWHumanClinP$MinBP <- rep("Min")
MWHumanClinOP$MinBP <- rep("BP")
MWHumanClinCom <- rbind(MWHumanClinP, MWHumanClinOP)
MWHumanClinCom$MinBP <- as.factor(MWHumanClinCom$MinBP)
MannWhitMWHumanClin <- wilcox.test(Percentage ~ MinBP, data = MWHumanClinCom)
# Mouse - Pre-clinical
MWMouseP$MinBP <- rep("Min")
MWMouseOP$MinBP <- rep("BP")
MWMouseCom <- rbind(MWMouseP, MWMouseOP)
MWMouseCom$MinBP <- as.factor(MWMouseCom$MinBP)
MannWhitMWMouse <- wilcox.test(Percentage ~ MinBP, data = MWMouseCom)
```

# Supplementary Tables

## Supplementary Table 1. Compliance of Human Clinical Studies to Minimal MSI Reporting Standards 


```r
HumanClinP <- cbind(MLHumanClinP[, 1:2], MWHumanClinP[, 2])
rownames(HumanClinP) <- HumanClinP$Metadata
colnames(HumanClinP) <- c("Metadata", "MetaboLights", "Metabolomics Workbench")
HumanClinP$total <- rowSums(HumanClinP[2:3])
HumanClinP <- HumanClinP[order(-HumanClinP$total), ]
kable(HumanClinP[, 2:3], caption = "The percentage of *Homo sapiens* studies in each repository that comply with each mammalian clinical trials and human studies minimal reporting standard. ")
```



Table: The percentage of *Homo sapiens* studies in each repository that comply with each mammalian clinical trials and human studies minimal reporting standard. 

                                      MetaboLights   Metabolomics Workbench
-----------------------------------  -------------  -----------------------
Biofluid or Tissue                           98.28                    90.91
Number of Groups                             86.21                    93.94
Disease Status                               60.34                    61.62
Sample Storage Temperature                   58.62                    53.54
Ethical Approval                             87.93                    15.15
Age Range                                    79.31                    19.19
Weight range and Height and/or BMI           46.55                    32.32
Fasting Status                               31.03                    21.21
Gender                                       46.55                     3.03
Inclusion Criteria                           41.38                     5.05
Exclusion Criteria                           37.93                     7.07
Volume or Quantity of Collection             27.59                     5.05
Anticoagulant                                10.34                    10.10
Treatment                                     3.45                    11.11
Ethnicity                                     6.90                     6.06
Trial Type                                   10.34                     1.01
Treatment Dose                                3.45                     6.06
Location of Collection                        6.90                     2.02
Treatment Duration                            1.72                     5.05
Treatment Route                               1.72                     2.02
Treatment Vehicle                             0.00                     0.00
Bacteriostatic Agent                          0.00                     0.00

## Supplementary Table 2. Compliance of Human Cell Line Studies to Minimal MSI Reporting Standards


```r
HumanCellP <- cbind(MLHumanCellP[, 1:2], MWHumanCellP[, 2])
rownames(HumanCellP) <- HumanCellP$Metadata
colnames(HumanCellP) <- c("Metadata", "MetaboLights", "Metabolomics Workbench")
HumanCellP$total <- rowSums(HumanCellP[2:3])
HumanCellP <- HumanCellP[order(-HumanCellP$total), ]
kable(HumanCellP[, 2:3], caption = "The percentage of *Homo sapiens* studies in each repository that comply with each microbial and *in vitro* minimal reporting standard.")
```



Table: The percentage of *Homo sapiens* studies in each repository that comply with each microbial and *in vitro* minimal reporting standard.

                                           MetaboLights   Metabolomics Workbench
----------------------------------------  -------------  -----------------------
Metabolism Quenching Method                       83.33                    62.22
Harvesting Method                                 83.33                    48.89
Metabolite Extraction                             77.78                    51.11
Sample Storage                                    72.22                    37.78
Normalisation by Cell Number                       5.56                     2.22
Cell Integrity                                     5.56                     0.00
Stability                                          5.56                     0.00
Temperature from Sampling to Quenching             0.00                     4.44
Time Until Quenching                               0.00                     2.22
Extracellular Metabolites Discriminated            0.00                     2.22
Recovery from Extraction                           0.00                     2.22
Sample Clean-up                                    0.00                     0.00
Sample Storage Duration                            0.00                     0.00
Quality Control                                    0.00                     0.00
Detection Limit                                    0.00                     0.00

## Supplementary Table 3. Compliance of Mouse Pre-Clinical Studies to Minimal MSI Reporting Standards


```r
MouseP <- cbind(MLMouseP[, 1:2], MWMouseP[, 2])
rownames(MouseP) <- MouseP$Metadata
colnames(MouseP) <- c("Metadata", "MetaboLights", "Metabolomics Workbench")
MouseP$total <- rowSums(MouseP[2:3])
MouseP <- MouseP[order(-MouseP$total), ]
kable(MouseP[, 2:3], caption = "The percentage of *Mus musculus* studies in each repository that comply with each pre-clinical minimal reporting standard.")
```



Table: The percentage of *Mus musculus* studies in each repository that comply with each pre-clinical minimal reporting standard.

                                           MetaboLights   Metabolomics Workbench
----------------------------------------  -------------  -----------------------
Biofluid or Tissue                               100.00                    89.01
Number of Groups                                  75.86                    91.21
Strain                                           100.00                    61.54
Treatment                                         75.86                    78.02
Sex                                               93.10                    46.15
Collection Time                                   89.66                    40.66
Age at Collection or Euthanization                68.97                    46.15
Sample Storage Temperature                        65.52                    35.16
Age at Study Start                                75.86                    24.18
Treatment Duration                                72.41                    27.47
Animal Supplier                                   82.76                    10.99
Treatment Dose                                    72.41                    19.78
Diet                                              44.83                    40.66
Treatment Route                                   68.97                    14.29
ab lib or Restricted Diet                         68.97                     6.59
Light Cycle                                       55.17                     9.89
Treatment Vehicle                                 55.17                     6.59
Euthanasia Method                                 41.38                    18.68
Tissue Processing                                 27.59                    25.27
Group or Individual Housing                       31.03                     5.49
Fasting Status                                    27.59                     8.79
Weight Range                                      34.48                     0.00
Volume or Quantity of Sample Collection            6.90                    12.09
Tap or Purified Water                              6.90                     0.00
Anticoagulant                                      6.90                     0.00
Collection Method                                  3.45                     1.10
Collection Frequency                               3.45                     1.10
Location of Sample Collection                      0.00                     2.20
Collection Duration                                0.00                     0.00
Bacteriostatic Agent                               0.00                     0.00

## Supplementary Table 4. Compliance of Arabidopis Studies to Minimal MSI Reporting Standards


```r
PlantP <- cbind(MeRPlantP[, 1:2], MLPlantP[, 2], MWPlantP[, 2], MPPlantP[, 2])
colnames(PlantP) <- c("Metadata", "MeRy-B", "MetaboLights", "Metabolomics Workbench", 
    "MetaPhen")
PlantP$total <- rowSums(PlantP[2:5])
PlantP <- PlantP[order(-PlantP$total), ]
kable(PlantP[, 2:5], caption = "The percentage of *Arabidopsis thaliana* studies in each repository that comply with each plant minimal reporting standard.")
```



Table: The percentage of *Arabidopsis thaliana* studies in each repository that comply with each plant minimal reporting standard.

                                  MeRy-B   MetaboLights   Metabolomics Workbench   MetaPhen
-------------------------------  -------  -------------  -----------------------  ---------
Organ or Cell Type                   100            100                   100.00      94.59
Plant Growth Stage                   100             65                   100.00      78.38
Genotype                             100             70                    66.67     100.00
Growth Support                        50             95                    66.67      97.30
Light                                 50             90                    66.67      97.30
Metabolism Quenching Method          100             30                   100.00      64.86
Date(s) of Plant Establishment        50             80                    66.67      89.19
Temperature                           50             70                    66.67      78.38
Harvest Time, Date                     0             65                   100.00      78.38
Nutrients Regime                      50             40                    66.67      64.86
Harvest Method                        50             30                   100.00      35.14
Biosource Amount                      50             65                    66.67      13.51
Sample Storage                        50             40                    66.67      32.43
Humidity                               0             20                    66.67      32.43
Treatment                              0             25                    33.33      59.46
Treatment Time                         0             25                    33.33      56.76
Treatment Dose                         0             20                    33.33      59.46
Plot Design                            0              5                    66.67       5.41
Growth Location                        0             15                    33.33       0.00
Watering Regime                        0              0                    33.33       2.70

## Supplementary Table 5. Compliance of Human Clinical Studies to Best Practice MSI Reporting Standards   


```r
HumanClinOP <- cbind(MLHumanClinOP[, 1:2], MWHumanClinOP[, 2])
rownames(HumanClinOP) <- HumanClinOP$Metadata
colnames(HumanClinOP) <- c("Metadata", "MetaboLights", "Metabolomics Workbench")
HumanClinOP$total <- rowSums(HumanClinOP[2:3])
HumanClinOP <- HumanClinOP[order(-HumanClinOP$total), ]
kable(HumanClinOP[, 2:3], caption = "The percentage of *Homo sapiens* studies in each repository that comply with each mammalian clinical trials and human studies optional reporting standard.")
```



Table: The percentage of *Homo sapiens* studies in each repository that comply with each mammalian clinical trials and human studies optional reporting standard.

                                    MetaboLights   Metabolomics Workbench
---------------------------------  -------------  -----------------------
Arterial or Venous Blood                   31.03                     4.04
Speed of Centrifugation                    25.86                     6.06
Temperature of Centrifugation              22.41                     6.06
Time of Centrifugation                     22.41                     4.04
Smoking Status                              6.90                    12.12
Time from Collection to Freezing            6.90                     2.02
Drug Consumption                            3.45                     1.01
Alcohol Consumption                         3.45                     0.00
Hemoglobin                                  3.45                     0.00
Platelets                                   3.45                     0.00
White Blood Count                           3.45                     0.00
Diet                                        1.72                     0.00
Malnutrition                                1.72                     0.00
Creatinine                                  1.72                     0.00
Hemocrit                                    1.72                     0.00
Sodium                                      1.72                     0.00
Potassium                                   1.72                     0.00
Sample Storage Duration                     1.72                     0.00
Mid Flow or Total Urine                     1.72                     0.00
Metal Exposure                              0.00                     0.00
Urea                                        0.00                     0.00
Glucose                                     0.00                     0.00
Total Cholesterol                           0.00                     0.00
HDL Cholesterol                             0.00                     0.00
LDL Cholesterol                             0.00                     0.00
Triglycerides                               0.00                     0.00
Total Protein                               0.00                     0.00
Albumin                                     0.00                     0.00
Bilirubin                                   0.00                     0.00
ALT                                         0.00                     0.00
ALP                                         0.00                     0.00
GT                                          0.00                     0.00
Hemolysis                                   0.00                     0.00

## Supplementary Table 6. Compliance of Human Cell Line Studies to Best Practice MSI Reporting Standards  


```r
HumanCellOP <- cbind(MLHumanCellOP[, 1:2], MWHumanCellOP[, 2])
rownames(HumanCellOP) <- HumanCellOP$Metadata
colnames(HumanCellOP) <- c("Metadata", "MetaboLights", "Metabolomics Workbench")
HumanCellOP$total <- rowSums(HumanCellOP[2:3])
HumanCellOP <- HumanCellOP[order(-HumanCellOP$total), ]
kable(HumanCellOP[, 2:3], caption = "The percentage of *Homo sapiens* studies in each repository that comply with each microbial and *in vitro* best practice reporting standard")
```



Table: The percentage of *Homo sapiens* studies in each repository that comply with each microbial and *in vitro* best practice reporting standard

                                          MetaboLights   Metabolomics Workbench
---------------------------------------  -------------  -----------------------
Cell Type                                       100.00                    77.78
Treatment                                        88.89                    60.00
Treatment Dose                                   88.89                    44.44
Medium or Substrate                              94.44                    37.78
Medium or Substrate Concentration                88.89                    31.11
Treatment Time                                   77.78                    33.33
Treatment Vehicle                                83.33                    26.67
Growth Container                                 61.11                    24.44
Cell Supplier                                    66.67                    17.78
Medium or Substrate Supplier                     66.67                    17.78
CO2                                              50.00                    31.11
Temperature                                      44.44                    22.22
Harvesting Time                                  22.22                    31.11
Isotopic Labelling                               11.11                    31.11
Replicates                                       33.33                     6.67
Subculturing and Splitting Protocols             33.33                     0.00
Inoculation Size                                 11.11                     2.22
Immortalized or Transformed                      11.11                     0.00
Growth Support                                   11.11                     0.00
pO2                                               5.56                     4.44
pH                                                5.56                     2.22
Additional -omics Datasets                        5.56                     0.00
Growth Container Supplier                         5.56                     0.00
Growth Support Supplier                           5.56                     0.00
Humidity                                          0.00                     2.22
Growth Configuration                              0.00                     0.00
Gas Composition                                   0.00                     0.00
Stirrer Speed                                     0.00                     0.00
Evaporation                                       0.00                     0.00
Growth Rate                                       0.00                     0.00
Pretreatment                                      0.00                     0.00
Pretreatment Time                                 0.00                     0.00
Harvesting Growth Phase                           0.00                     0.00
Number of Generations Until Harvesting            0.00                     0.00
Stabilization Time                                0.00                     0.00
Number of Culture Passages                        0.00                     0.00
Marker of Differentiated Stage                    0.00                     0.00
Harvesting Cell Density                           0.00                     0.00
Harvesting Depletion of Nutrients                 0.00                     0.00

## Supplementary Table 7. Compliance of Mouse Pre-Clinical Studies to Best Practice MSI Reporting Standards 


```r
MouseOP <- cbind(MLMouseOP[, 1:2], MWMouseOP[, 2])
rownames(MouseOP) <- MouseOP$Metadata
colnames(MouseOP) <- c("Metadata", "MetaboLights", "Metabolomics Workbench")
MouseOP$total <- rowSums(MouseOP[2:3])
MouseOP <- MouseOP[order(-MouseOP$total), ]
kable(MouseOP[, 2:3], caption = "The percentage of *Mus musculus* studies in each repository that comply with each pre-clinical optional reporting standard.")
```



Table: The percentage of *Mus musculus* studies in each repository that comply with each pre-clinical optional reporting standard.

                                                 MetaboLights   Metabolomics Workbench
----------------------------------------------  -------------  -----------------------
Use of Anesthesia                                       41.38                     7.69
Environmental Control: Temperature                      31.03                     2.20
Acclimation Duration to Experimental Facility           24.14                     3.30
Germ-free or Conventional Housing                       24.14                     2.20
Fasting Duration                                        17.24                     5.49
Bedding Type                                             6.90                    15.38
Environmental Control: Humidity                         13.79                     2.20
Anesthesia Time                                         13.79                     0.00
Anesthesia Dose                                          6.90                     2.20
Cage Cleaning Frequency                                  3.45                     0.00
Cage Type                                                3.45                     0.00
Inclusion Criteria                                       3.45                     0.00
Additional Phenotypic Model                              0.00                     3.30
Sample Storage Duration                                  0.00                     2.20
Temperature of Collection Tube                           0.00                     1.10
Body Weights or Food Consumption                         0.00                     0.00

