##### Single nucleus RNA-seq Data Processing #####

## This script is stuctured such that raw data and reference files should be in the following directory formats:
## 1 Working Directory
## 1 ./sampleInfo.txt
## 1 ./AGI_Synonyms_AtPlastidGenes.csv
## 1 ./Romanowski.cyclinggenes.csv
## 1 ./features.tsv
## 1 ./JTK_PseudoBulk _NO_SOUP_0.05adjp.xlsx
## -- 2 ./{ID}cr15k/
## -- -- 3 ./outs/
## -- -- -- 4 ./filtered_feature_bc_matrix/
## -- -- -- 4 ./raw_feature_bc_matrix/

##### Read in data #####
sampleInfo <- read_csv("./sampleInfo.txt", col_names = F)
AtPlastidGenes <- read.csv("./AGI_Synonyms_AtPlastidGenes.csv", col.names = c("gid","synonym"))
AtMitoGenes <- AtPlastidGenes[88:207,] #%>% dplyr::select('synonym')
AtChloroplastGenes <- AtPlastidGenes[1:87,] #%>% dplyr::select('synonym')
AtChloroplastGenes[88,1:2] <- c("ATCG00020","psbA")
AtPlastidGenes <- AtPlastidGenes #%>% dplyr::select('synonym')
AtPlastidGenes[208,1:2] <- c("ATCG00020","psbA")

AtCycling <- readxl::read_excel("./Romanowski.cyclinggenes.csv")
AtCycling <- AtCycling %>% dplyr::select(Gene_ID)
colnames(AtCycling) <- c("gid")
AtAllSynonyms <- read_tsv("./features.tsv", col_names = c("gid","synonym","X3"))
AtCycling <- inner_join(AtCycling, AtAllSynonyms) #%>% dplyr::select(synonym)
AlveenaCyclingGenes <- readxl::read_excel("/home/greenham/myers998/snTC/JTK_PseudoBulk _NO_SOUP_0.05adjp.xlsx") %>% dplyr::select(CycID)
colnames(AlveenaCyclingGenes) <- c("gid")

AtCycling <- full_join(AtCycling, AlveenaCyclingGenes)

##### NO SOUP CORRECTIONS #####
#### ZT00C ####
i <- sampleInfo[3,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT00C", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "Two"

#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT00C.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)

#### ZT00D ####
i <- sampleInfo[4,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT00d", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts)) 
soupOutMAD[["day"]] <- "Two"

#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT00D.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)


#### ZT04A ####
i <- sampleInfo[5,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT04A", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT04A.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)


#### ZT04B ####
i <- sampleInfo[6,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT04B", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT04B.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)


#### ZT08A ####
i <- sampleInfo[7,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT08A", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT08A.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)

#### ZT08B ####
i <- sampleInfo[8,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT08B", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT08B.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)


#### ZT12A ####
i <- sampleInfo[9,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT12A", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT12A.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)

#### ZT12B ####
i <- sampleInfo[10,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT12B", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT12B.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)

#### ZT16A ####
i <- sampleInfo[11,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT16A", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT16A.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)

#### ZT16B ####
i <- sampleInfo[12,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT16B", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT16B.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)

#### ZT20A ####
i <- sampleInfo[13,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT20A", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT20A.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)

#### ZT20B ####
i <- sampleInfo[14,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT20B", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT20B.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)

#### ZT24A ####
i <- sampleInfo[15,1]

## Read in data
filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)

# scDblFinder
scDbls <- scDblFinder(filt.matrix)

# Create Seurat object and add contamination metadata
soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT24A", min.cells = 1, min.features = 100)
table(soupOut$orig.ident)

tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
colnames(tmpAll) <- c("synonym")
tmpMito <- inner_join(tmpAll, AtMitoGenes)
tmpMito <- as.list(tmpMito)
tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
tmpChloro <- as.list(tmpChloro)

tmpCycling <- inner_join(tmpAll, AtCycling)
tmpCycling <- as.list(tmpCycling)
soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)

soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)


soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]

VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
soupOut <- subset(soupOut, multiplet == "singlet")

soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)

counts <- GetAssayData(soupOutMAD, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
soupOutMAD[["day"]] <- "One"
#write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# scTransform
ZT24A.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))

rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)

# #### ZT24B ####
# i <- sampleInfo[16,1]
# 
# ## Read in data
# filt.matrix <- Read10X(paste("./", i, "cr15k/outs/filtered_feature_bc_matrix/", sep = ""), gene.column = 1)
# raw.matrix <- Read10X(paste("./", i, "cr15k/outs/raw_feature_bc_matrix/", sep = ""), gene.column = 1)
# 
# # scDblFinder
# scDbls <- scDblFinder(filt.matrix)
# 
# # Create Seurat object and add contamination metadata
# soupOut <- CreateSeuratObject(counts = filt.matrix, project = "ZT24B", min.cells = 1, min.features = 100)
# table(soupOut$orig.ident)
# 
# tmpAll <- as.data.frame(rownames(soupOut@assays$RNA@features@.Data))
# colnames(tmpAll) <- c("synonym")
# tmpMito <- inner_join(tmpAll, AtMitoGenes)
# tmpMito <- as.list(tmpMito)
# tmpChloro <- inner_join(tmpAll, AtChloroplastGenes)
# tmpChloro <- as.list(tmpChloro)
# 
# tmpCycling <- inner_join(tmpAll, AtCycling)
# tmpCycling <- as.list(tmpCycling)
# soupOut[["percent.cyc"]] <- PercentageFeatureSet(soupOut, features = tmpCycling$synonym)
# 
# soupOut[["percent.mt"]] <- PercentageFeatureSet(soupOut, features = tmpMito$synonym)
# soupOut[["percent.cp"]] <- PercentageFeatureSet(soupOut, features = tmpChloro$synonym)
# 
# 
# soupOut[["multiplet"]] <- scDbls$scDblFinder.class[match(rownames(soupOut@meta.data), rownames(scDbls@colData))]
# 
# VlnPlot(soupOut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.cp"), ncol = 4, pt.size = 0)
# soupOut <- subset(soupOut, multiplet == "singlet")
# 
# soupOutMAD <- filterbyMAD(soupOut, tmpMito$gid, tmpChloro$gid)
# 
# counts <- GetAssayData(soupOutMAD, assay = "RNA")
# counts <- counts[-(which(rownames(counts) %in% AtPlastidGenes$gid)),]
# soupOutMAD <- subset(soupOutMAD, features = rownames(counts))
# soupOutMAD[["day"]] <- "One"
# #write.csv(counts, paste("20240506.RawCounts.Sample",i,".NoSoup.csv", sep = ""))
# # scTransform
# ZT24B.sct <- soupOutMAD #%>% SCTransform(vars.to.regress = c("percent.cp","percent.mt"))
# 
# rm(i, scDbls, soupOut, soupOutMAD, tmpAll, tmpChloro, tmpMito, tmpCycling, srat, raw.matrix, filt.matrix, filt.matrix, counts, meta, soup.channel, umap)

##### Merging samples #####
snTClistQCdNoSoup <- list(ZT00C.sct, ZT00D.sct, ZT04A.sct, ZT04B.sct, ZT08A.sct, ZT08B.sct, ZT12A.sct, ZT12B.sct, ZT16A.sct, ZT16B.sct, ZT20A.sct, ZT20B.sct, ZT24A.sct, ZT24B.sct)
saveRDS(snTClistQCdNoSoup, "20240507.snTClist.SoupComp.WithoutSoup.BulkCycUpdated.RDS")

merged_seurat <- merge(x = snTClistQCdNoSoup[[1]],
                       y = snTClistQCdNoSoup[2:13])

merged_seurat15k <- merged_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 15000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30)

tmp <- JoinLayers(merged_seurat15k)
tmp2 <- rownames(tmp@assays$RNA$data)
tmp3 <- as.data.frame(tmp2)
colnames(tmp3) <- "gid"
tmp4 <- anti_join(tmp3, AlveenaCyclingGenes)
integrated_seurat <- IntegrateLayers(
  object = merged_seurat15k,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = F,
  features = tmp4$gid
)
integrated_seurat <- JoinLayers(integrated_seurat)
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:30, reduction = "integrated.rpca")
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.3)
integrated_seurat15k <- RunUMAP(integrated_seurat, reduction = "integrated.rpca", dims = 1:30)
p <- DimPlot_scCustom(integrated_seurat15k)

##### Subclustering #####
tmp <- FindSubCluster(integrated_seurat15k, c(8), graph.name = "RNA_snn", resolution = 0.1)
Idents(tmp) <- "sub.cluster"
tmp <- FindSubCluster(tmp, c(12), graph.name = "RNA_snn", resolution = 0.1)
Idents(tmp) <- "sub.cluster"
DefaultAssay(tmp) <- "RNA"

saveRDS(object = tmp, file = "20241215.snTC.DROP24B.Subclustered.RDS")

tmpMarkers <- FindAllMarkers(tmp, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "LR")
write_tsv(tmpMarkers, "20241215.snTC.DROP24B.Subclustered.Markers.tsv")


