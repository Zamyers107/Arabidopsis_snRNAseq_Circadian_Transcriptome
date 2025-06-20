##### MAD Filtering of a Seurat Object #####
##### ZM 1/29/2024 #####

## Inputs
# 
#
## Outputs
#
#

filterbyMAD <- function(soupSeuratObject,mGIDs,cpGIDs){
  
  ##Cell quality control##
  #Set ggplot theme as classic
  theme_set(theme_classic())
  
  ##Load the dataset and calculate QC metrics
  # Load the raw filtered_gene_bc_matrix outputed by Cell Ranger v2.1.1
  #Countdata <- Read10X(data.dir = "/Users/alveenazulfiqar/Desktop/snTC_Data_Files/snTC/soupX_ZT00Ccr15k_filt/")
  # SO = Seurat Object
  SO <- soupSeuratObject
  SO@meta.data$Barcodes <- rownames(SO@meta.data)
  
  AtMitoGenes <- mGIDs
  AtChloroplastGenes <- cpGIDs

  ##Low quality cell filtering
  Cell.QC.Stat <- SO@meta.data
  
  #Filtering cells based on percentage of mitochondrial transcripts
  #max.mito.thr <- median(Cell.QC.Stat$percent.mt) + 3*mad(Cell.QC.Stat$percent.mt)
  #min.mito.thr <- median(Cell.QC.Stat$percent.mt) - 3*mad(Cell.QC.Stat$percent.mt)
  #p4 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mt)) +
  #  geom_point() +
  #  geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
  #  geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
  #  annotate(geom = "text", label = paste0(as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr | Cell.QC.Stat$percent.mt < min.mito.thr)[2])," cells removed\n",
  #                                         as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr | Cell.QC.Stat$percent.mt < min.mito.thr)[1])," cells remain"), x = 6000, y = 0.1)
  #ggMarginal(p4, type = "histogram", fill="lightgrey", bins=100)
  # Filter cells based on these thresholds
  #Cell.QC.Stat <- Cell.QC.Stat %>% filter(percent.mt < max.mito.thr) %>% filter(percent.mt > min.mito.thr)
  
  # Calculate max.cp.thr using median and adjusted MAD
  max.cp.thr <- median(Cell.QC.Stat$percent.cp) + 3*mad(Cell.QC.Stat$percent.cp)
  min.cp.thr <- median(Cell.QC.Stat$percent.cp) - 3*mad(Cell.QC.Stat$percent.cp)
  p5 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.cp)) +
    geom_point() +
    geom_hline(aes(yintercept = max.cp.thr), colour = "red", linetype = 2) +
    geom_hline(aes(yintercept = min.cp.thr), colour = "red", linetype = 2) +
    annotate(geom = "text", label = paste0(as.numeric(table(Cell.QC.Stat$percent.cp > max.cp.thr | Cell.QC.Stat$percent.cp < min.cp.thr)[2])," cells removed\n",
                                           as.numeric(table(Cell.QC.Stat$percent.cp > max.cp.thr | Cell.QC.Stat$percent.cp < min.cp.thr)[1])," cells remain"), x = 6000, y = 0.1)
  ggMarginal(p5, type = "histogram", fill="lightgrey", bins=100)
  # Filter cells based on these thresholds
  
  Cell.QC.Stat <- Cell.QC.Stat %>% filter(percent.cp < max.cp.thr) %>% filter(percent.cp > min.cp.thr)
  #Filtering cells based on number of genes and transcripts detected:Remove cells with to few gene detected or with to many UMI counts
  # Set low and hight thresholds on the number of detected genes
  min.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) - 3*mad(log10(Cell.QC.Stat$nFeature_RNA))
  max.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) + 3*mad(log10(Cell.QC.Stat$nFeature_RNA))
  # Set hight threshold on the number of transcripts
  max.nUMI.thr <- median(log10(Cell.QC.Stat$nCount_RNA)) + 3*mad(log10(Cell.QC.Stat$nCount_RNA))
  # Gene/UMI scatter plot before filtering
  p1 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.Genes.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.Genes.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = max.nUMI.thr), colour = "red", linetype = 2)
  ggMarginal(p1, type = "histogram", fill="lightgrey")
  # Filter cells base on both metrics
  Cell.QC.Stat <- Cell.QC.Stat %>% filter(log10(nFeature_RNA) > min.Genes.thr) %>% filter(log10(nCount_RNA) < max.nUMI.thr)
  lm.model <- lm(data = Cell.QC.Stat, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
  p2 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.Genes.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.Genes.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = max.nUMI.thr), colour = "red", linetype = 2) +
    geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
    annotate(geom = "text", label = paste0(dim(Cell.QC.Stat)[1], " QC passed cells"), x = 4, y = 3.8)
  ggMarginal(p2, type = "histogram", fill="lightgrey")
  ###Filter cells below the main population nUMI/nGene relationship
  Cell.QC.Stat$valideCells <- log10(Cell.QC.Stat$nFeature_RNA) > (log10(Cell.QC.Stat$nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - 0.09))
  p3 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point(aes(colour = valideCells)) +
    geom_smooth(method="lm") +
    geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
    theme(legend.position="none") +
    annotate(geom = "text", label = paste0(as.numeric(table(Cell.QC.Stat$valideCells)[2]), " QC passed cells\n",
                                           as.numeric(table(Cell.QC.Stat$valideCells)[1]), " QC filtered"), x = 4, y = 3.8)
  ggMarginal(p3, type = "histogram", fill="lightgrey")
  # Remove unvalid cells
  Cell.QC.Stat <- Cell.QC.Stat %>% filter(valideCells)
  SO_Filtered <- subset(SO, cells = rownames(Cell.QC.Stat))
  
  
  SO_Filtered
  
}
