libraryList <- c("flowCore", "limma", "edgeR", "FlowSOM", "flowType", "ncdfFlow", "openCyto", "flowStats", "tidyverse", "gplots", "pvclust")

lapply(libraryList, require, quietly = T, character.only = TRUE)

load("../Cox-et-al-meat-allergy-data/flowTypeResultsAllEventsAllFilesMeatAllergyProjectUpdate.Rdata") # loads up ResList from flowType run

countTable <- ResList %>%
  lapply(function(x){return(x@CellFreqs)}) %>%
  bind_cols() %>%
  setNames(names(ResList))

cell_count_total <- lapply(countTable, max) %>%
  data.frame(FileNames = names(.),
             `n()` = unlist(.)) %>%
  select(FileNames:`n..`)

load("../Cox-et-al-meat-allergy-data/FlowTypeProportionsAndStatsTable.Rdata") # loads up statsTable which contains decimal values for each immunophenotype
row.names(countTable) <- row.names(statsTable)

propData <- DGEList(countTable,
                    lib.size = cell_count_total$n..[cell_count_total$FileNames %in% colnames(countTable)])

new_targets <- read_csv("../Cox-et-al-meat-allergy-data/targets_data_frame.csv")

exprDesign <- new_targets$Condition
exprDesign <- factor(exprDesign, levels =c("HD", "AD"))
gender <- factor(new_targets$gender, levels = c("male", "female"))
age <- new_targets$age

design <- model.matrix(~0 + exprDesign + age)

# estimate glm robust dispersion, quasi-likelihood framework F test
fit <- estimateGLMRobustDisp(propData[,!is.na(age)], design)

fit_glm_robust <- fit

fit_fit <- glmQLFit(fit_glm_robust, design, robust=TRUE)

res <- glmQLFTest(fit_fit, contrast = c(-1,1,0))

topTable <- topTags(res, n = Inf)$table %>%
  rownames_to_column("Immunophenotype") %>%
  mutate(Observation = rep("Count", nrow(.)),
         Case = rep("Allergy", nrow(.)),
         Baseline = rep("Healthy", nrow(.)))

heat_map_data_frame <- topTable %>%
  filter(FDR < 0.05,
         logFC > 0) %>%
  left_join(cpm(propData, log = T) %>%
              as.data.frame() %>%
              rownames_to_column("Immunophenotype"))

heat_map_matrix <- as.matrix(heat_map_data_frame[,c(10:48)])
row.names(heat_map_matrix) <- heat_map_data_frame$Immunophenotype

column_labels <- data.frame(SampleID = colnames(heat_map_matrix)) %>%
  left_join(new_targets)

dendrogram <- pvclust(t(heat_map_matrix),
                      method.dist = "euclidean",
                      nboot = 1)

tiff(filename = "../Cox-et-al-meat-allergy-data/log-cpm-heat-map.tiff",
     width = 8,
     height = 8,
     units = "in",
     res = 600)
heatmap.2(heat_map_matrix,
          Colv = F,
          Rowv = as.dendrogram(dendrogram$hclust),
          dendrogram = "row",
          col = bluered(25),
          trace = "none",
          scale = "none",
          labCol = column_labels$SampleID,
          density.info = "none",
          labRow = F)
dev.off()

marker_list <- c("CD196","HLA\\-DR", "IgD", "CD25", "CD43", "CD45", "CD279\\-FITC", "CD38", "CD24", "CXCR5", "IgM", "CXCR4")

immunophenotype_matrix <- matrix(nrow = length(heat_map_data_frame$Immunophenotype),
                                 ncol = length(marker_list))
for( immunophenotype in heat_map_data_frame$Immunophenotype){
  for( marker in marker_list){
    if(!str_detect(immunophenotype, marker)){
      marker_val <- 0
    }else if(str_detect(immunophenotype, paste0(marker, "\\+"))){
      marker_val <- 1
    }else if(str_detect(immunophenotype, paste0(marker, "\\-"))){
      marker_val <- -1
    }else{
      marker_val <- NA
    }
    immunophenotype_matrix[which(heat_map_data_frame$Immunophenotype %in% immunophenotype), which(marker_list %in% marker)] <- marker_val
  }
}

colnames(immunophenotype_matrix) <- gsub("\\-FITC|\\\\", "", marker_list)

tiff(filename = "../Cox-et-al-meat-allergy-data/labels-heat-map.tiff",
     width = 4,
     height = 8,
     units = "in",
     res = 600)
heatmap.2(immunophenotype_matrix[,-6],
          Colv = T,
          Rowv = as.dendrogram(dendrogram$hclust),
          dendrogram = "none",
          col = bluered(3),
          trace = "none",
          scale = "none",
          key = F,
          labRow = F)
dev.off()