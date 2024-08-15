library(edgeR)
library(dplyr)
library(here)
# open counts

# crear nombre de muestras
samples <- c("AB33_1", "AB33_2", "AB33_3",
             "AB30_1", "AB30_2", "AB30_3",
             "PAU33_1", "PAU33_2", "PAU33_3",
             "PAU30_1", "PAU30_2", "PAU30_3",
             "CTL_1", "CTL_2", "CTL_3")

abundancia <- read.table("test_data/RSEM.isoform.counts.matrix",
                         sep = "\t", header = TRUE, row.names = 1)
abundancia <- as.matrix(abundancia)


head(abundancia)
colnames(abundancia) <- samples


condicion <- gsub("*_.", "", samples)
replica <- gsub(".*_", "", samples)

colData <- data.frame(samples = samples,
                      condicion = condicion,
                      replica = replica)

rownames(colData) <- samples


# Crear el objeto DGEList
dge <- DGEList(counts=abundancia, group=colData$condicion)

# Filtrado de genes con baja expresión
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]


#guardar nueva table

write.table(dge$counts, file = here("test_data/RSEM_filter_counts.txt"),
            row.names = TRUE, 
            sep = "\t",
            quote = FALSE)



# Normalización de los datos
dge <- calcNormFactors(dge)



# Diseño experimental
design <- model.matrix(~0 + as.factor(colData$condicion))
colnames(design) <- levels(as.factor(condicion))

# Estimar la dispersión
dge <- estimateDisp(dge, design)

# Ajuste del modelo lineal generalizado
fit <- glmQLFit(dge, design)

# Comparaciones específicas
comparisons <- list(
  AB33_vs_CTL = makeContrasts(AB33 - CTL, levels=design),
  AB30_vs_CTL = makeContrasts(AB30 - CTL, levels=design),
  PAU33_vs_CTL = makeContrasts(PAU33 - CTL, levels=design),
  PAU30_vs_CTL = makeContrasts(PAU30 - CTL, levels=design)
)




# Initialize an empty list to store results
all_results <- list()

# Loop through each contrast, perform LRT, and store results
for (contrast_name in names(comparisons)) {
  glrt <- glmLRT(fit, contrast = comparisons[[contrast_name]])
  top_genes <- topTags(glrt, n=Inf)$table
  top_genes <- top_genes |> 
    tibble::rownames_to_column("transcript")
  
  # Add a column with the name of the contrast
  top_genes <- top_genes |> 
    mutate(contrast = contrast_name)
  
  # add up or downreg colum
  top_genes <- top_genes |> 
    mutate(expression = ifelse(logFC < 1, "down", "up"))
  
  # Store the result in the list
  all_results[[contrast_name]] <- top_genes
}

# Combine all results into a single data frame
final_results <- bind_rows(all_results)

# View the combined table
head(final_results)


write.table(final_results, file = here("test_data/DET_test.txt"), sep = "\t",
            row.names = FALSE,quote = FALSE )





