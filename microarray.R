# GSE31684: Male vs. Female
#---- Load libraries------------------------------------------------------------
library(GEOquery)
library(limma)
library(ggplot2)
library(dplyr)
#---- Parameters ---------------------------------------------------------------
geo_id <- "GSE31684"
#---- Input and output directories ---------------------------------------------
# Define directory paths
dir_data <- './data/'
dir_output <- './output/'
# Check if the input directory exists, if not create it
if (!dir.exists(dir_data)) {
  dir.create(dir_data)
}
# Check if the output directory exists, if not create it
if (!dir.exists(dir_output)) {
  dir.create(dir_output)
}
#---- Raw data -----------------------------------------------------------------
# Download the GEO dataset
gds <- getGEO(geo_id)
eset <- gds[[1]]
# Save expression matrix
exp_matrix <- exprs(eset)
write.csv(exp_matrix, file.path(dir_data, 'exp.csv'))
# Save metadata
metadata <- pData(phenoData(eset))
write.csv(metadata, file.path(dir_data, 'metadata.csv'), row.names = FALSE)
# Save probe information
probe <- pData(featureData(eset))
write.csv(probe, file.path(dir_data, 'probe.csv'), row.names = FALSE)
#---- PCA ----------------------------------------------------------------------
pca_result <- prcomp(t(exp_matrix), center = TRUE, scale. = FALSE)
pca_data <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], 
                       Label = metadata$`gender:ch1`)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Label)) +
  geom_point() +
  theme_bw() +
  labs(x = "Principal Component 1", 
       y = "Principal Component 2", 
       color = "Gender") +
  theme(axis.title = element_text(size = 14), # Axis title
        axis.text = element_text(size = 14), # Axis text
        legend.title = element_text(size = 14), # Legend title
        legend.text = element_text(size = 14) # Legend text
        )
ggsave(file.path(dir_output, 'PCA.png'), width = 8, height = 5)
#---- Differential analysis ----------------------------------------------------
# Design matrix
metadata$`gender:ch1` <- factor(metadata$`gender:ch1`, 
                                levels = c("male", "female"))
design <- model.matrix(~0+metadata$`gender:ch1`)
colnames(design) <- c('male', 'female')
# Fit the lm model for differential gene expression analysis
fit <- lmFit(eset, design)
cont.matrix <- makeContrasts(malevsfemale=male-female, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
# Save the result
res <- topTable(fit2, adjust="BH", n=Inf) %>%
       select(ID, Gene.Symbol, ENTREZ_GENE_ID, logFC, AveExpr, t, P.Value, adj.P.Val)
write.csv(res, file.path(dir_output, 'differential_analysis.csv'), row.names = FALSE)
#---- GSEA ---------------------------------------------------------------------
rnk <- res[, c("Gene.Symbol", "t")]
idx <- grep("///", rnk$Gene.Symbol)
# Entries with single genes
rnk.sing <- rnk[-idx, ]
# Entries with multiple genes
rnk.mult <- rnk[idx, ]
ss <- strsplit(rnk.mult$Gene.Symbol, " /// ")
rnk.mult.expanded <- do.call(rbind, lapply(1:length(ss), function(i) {
  data.frame(
    Gene.Symbol = ss[[i]],
    t = rep(rnk.mult$t[i], length(ss[[i]]))
  )
}))
# Generate ranked genes
rnk.c <- rbind(rnk.sing, rnk.mult.expanded)
# Remove rows with empty gene symbol or test stats
rnk.c <- rnk.c[rnk.c$Gene.Symbol != "", ]
rnk.c <- rnk.c[!is.na(rnk.c$t), ]
# Take an average for genes with multiple probes
rnk.c <- rnk.c %>%
  group_by(Gene.Symbol) %>%
  summarise(t = mean(t, na.rm = TRUE)) %>%
  ungroup()
rnk.c <- rnk.c[order(rnk.c$t, decreasing=TRUE), ]
write.table(rnk.c, file.path(dir_output, 'gsea.rnk'), row.names = FALSE, 
            col.names = FALSE, sep = '\t')