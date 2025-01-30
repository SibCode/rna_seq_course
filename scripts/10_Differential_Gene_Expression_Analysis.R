# Install and load the required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("cachem")
BiocManager::install("SummarizedExperiment", force = TRUE)
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("org.Mm.eg.db")
install.packages("ggrepel")
install.packages("dplyr")

library(SummarizedExperiment)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Set working directory
setwd("~/OneDrive - Universitaet Bern/MScBInf/1. Semester/RNA Seq/rna_seq_course")

# Load the summary data
summary_data <- read.table("outputs/09_run_feature_count/feature_counts.txt.summary",
                           header=TRUE, row.names=1)

# Trim the sample names to remove the path and .bam ending
colnames(summary_data) <- gsub(".*/|\\.bam", "", colnames(summary_data))

# Trim periods instead of slashes (different OS, just to make sure)
colnames(summary_data) <- gsub(".*\\.|\\.bam", "", colnames(summary_data))

# Calculate the proportion of reads overlapping with annotated genes
proportion_overlapping <- summary_data["Assigned", ] / colSums(summary_data)

# Calculate the number of unassigned reads due to ambiguity convert to numeric
unassigned_ambiguity <- summary_data["Unassigned_Ambiguity", ]
unassigned_ambiguity <- as.numeric(as.character(unassigned_ambiguity))

# Calculate the average number of unassigned reads due to ambiguity
average_unassigned_ambiguity <- mean(unassigned_ambiguity)

cat("Proportion of reads overlapping with annotated genes in each sample:\n")
print(proportion_overlapping)
cat("\n")
cat("Average number of unassigned reads due to ambiguity:\n")
print(average_unassigned_ambiguity)

# Plot the proportion overlapping genes in a bar plot
# Convert proportion_overlapping to a data frame for plotting
proportion_overlapping_df <- data.frame(
  Sample = names(proportion_overlapping),
  Proportion =  as.numeric(proportion_overlapping * 100) # Convert to percentage
)
colnames(proportion_overlapping_df) <- c("Sample", "Proportion")

# Create the bar plot
ggplot(proportion_overlapping_df, aes(x = Sample, y = Proportion)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  xlab("Sample") +
  ylab("Percentage of Overlapping Reads") +
  ggtitle("Percentage of Reads Overlapping with Annotated Genes in Each Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 100) 


# Read the featureCounts file
countData <- read.table("outputs/09_run_feature_count/feature_counts.txt",
                        header=TRUE, row.names=1, skip=1)

# Remove columns containing Chr, Start, End, Strand
countData <- countData[ , !(colnames(countData) %in% c("Chr",
                                                       "Start",
                                                       "End",
                                                       "Strand"))]

# Trim the sample names to remove the path and .bam ending
colnames(countData) <- gsub(".*/|\\.bam", "", colnames(countData))

# Trim periods instead of slashes (different OS, just to make sure)
colnames(countData) <- gsub(".*\\.|\\.bam", "", colnames(countData))

# Reorder countData columns to match samples in order:
sampleOrder <- c("SRR7821921", "SRR7821922", "SRR7821918", "SRR7821919",
                 "SRR7821920", "SRR7821937", "SRR7821938", "SRR7821939",
                 "SRR7821949", "SRR7821950", "SRR7821951", "SRR7821952",
                 "SRR7821953", "SRR7821968", "SRR7821969", "SRR7821970")
countData <- countData[, sampleOrder]
                       
# Create a sample data frame with conditions according to sample / condition
sampleData <- data.frame(row.names = sampleOrder,
                         condition = factor(c(rep("Lung_WT_Case", 5),
                                              rep("Lung_WT_Control", 3),
                                              rep("Blood_WT_Case", 5),
                                              rep("Blood_WT_Control", 3))))
                          
# Create DESeqDataSet 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleData,
                              design = ~ condition)

# Normalize the counts
dds <- DESeq(dds)

# Remove the dependence of the variance on the mean
# by variance stabilizing transformation (vst)
vsd <- vst(dds, blind = TRUE)

### Exploratory data analysis

# Plot PCA
plotPCA(vsd, intgroup = "condition")

# Extract PCA data for ggplot2 (more beautiful)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Determine the range of PC1 and PC2 for expanding the scales
range_PC1 <- range(pcaData$PC1)
range_PC2 <- range(pcaData$PC2)

# Create the PCA plot with customized margins
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  theme(plot.margin = margin(25, 25, 25, 25)) +
  theme(legend.position = c(0.925, 0.85)) +  
  guides(color = guide_legend(override.aes = list(size = 4.0))) +
  scale_x_continuous(limits = c(range_PC1[1] - 15, range_PC1[2] + 15)) +
  scale_y_continuous(limits = c(range_PC2[1] - 15, range_PC2[2] + 15))


# Extract results for lung cases (Lung_WT_Case vs Lung_WT_Control)
lung_res <- results(dds, contrast = c("condition", "Lung_WT_Case", "Lung_WT_Control"))

# Extract results for blood cases (Blood_WT_Case vs Blood_WT_Control)
blood_res <- results(dds, contrast = c("condition", "Blood_WT_Case", "Blood_WT_Control"))

# Interesting Genes:
# Irf7 is NOT expressed in blood, but higher in lung tissue
# Furthermore Irf7 is supressed in the diseased samples in blood, but not lung
# Isg15, Rsad2, Apol9a is supressed in blood, and expressed higher in lung
# Define interesting genes and their labels
interesting_genes <- c("ENSMUSG00000057346", "ENSMUSG00000020641", "ENSMUSG00000035692", "ENSMUSG00000025498")
gene_labels <- c("Apol9a", "Rsad2", "Isg15", "Irf7")

# Create a combined data frame for the volcano plot
lung_volcano_data <- data.frame(
  log2FoldChange = lung_res$log2FoldChange,
  negLog10padj = -log10(lung_res$padj),
  Gene = rownames(lung_res),
  Tissue = "Lung"
)

blood_volcano_data <- data.frame(
  log2FoldChange = blood_res$log2FoldChange,
  negLog10padj = -log10(blood_res$padj),
  Gene = rownames(blood_res),
  Tissue = "Blood"
)

# Combine the data
combined_volcano_data <- rbind(lung_volcano_data, blood_volcano_data)

# Initialize Label column with empty strings
combined_volcano_data$Label <- ""

# Assign labels to interesting genes
for (i in 1:length(interesting_genes)) {
  combined_volcano_data$Label[combined_volcano_data$Gene == interesting_genes[i]] <- gene_labels[i]
}

# Assign labels to interesting genes with the respective tissue type
for (i in 1:length(interesting_genes)) {
  lung_mask <- combined_volcano_data$Gene == interesting_genes[i] & combined_volcano_data$Tissue == "Lung"
  blood_mask <- combined_volcano_data$Gene == interesting_genes[i] & combined_volcano_data$Tissue == "Blood"
  
  combined_volcano_data$Label[lung_mask] <- paste(gene_labels[i], "(L)")
  combined_volcano_data$Label[blood_mask] <- paste(gene_labels[i], "(B)")
}

# Volcano plot with highlights and labels, closing in on the interesting genes
ggplot(combined_volcano_data, aes(x = log2FoldChange, y = negLog10padj, color = Tissue)) +
  geom_point(alpha = 0.4, size = 1.75) +
  geom_point(data = subset(combined_volcano_data, Label != ""), aes(x = log2FoldChange, y = negLog10padj), size = 1.5, fill = NA, color = "black", stroke = 1) +
  geom_text_repel(data = subset(combined_volcano_data, Label != ""), aes(label = Label), vjust = 1.5, color = "black", nudge_y = 0.5) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-Value") +
  ggtitle("Volcano Plot of Differentially Expressed Genes in Lung and Blood Tissues") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "blue", linetype = "dashed") +
  xlim(-10, 10) +
  ylim(0, 100) +
  theme(legend.position = c(0.925, 0.85)) +  
  guides(color = guide_legend(override.aes = list(size = 4.0))) +
  annotate("text", x = -9, y = 95, label = "(B): Blood", hjust = 0, vjust = 1) +
  annotate("text", x = -9, y = 90, label = "(L): Lung", hjust = 0, vjust = 1)


### Differential expression analysis

# Extract the results for lung cases (Lung_WT_Case vs Lung_WT_Control)
lung_res <- results(dds, contrast = c("condition",
                                      "Lung_WT_Case",
                                      "Lung_WT_Control"))

# Extract the results for blood cases (Blood_WT_Case vs Blood_WT_Control)
blood_res <- results(dds, contrast = c("condition",
                                       "Blood_WT_Case",
                                       "Blood_WT_Control"))

# Filter significant genes (padj < 0.05) 
lung_sig <- lung_res[!is.na(lung_res$padj) & lung_res$padj < 0.05, ]
blood_sig <- blood_res[!is.na(blood_res$padj) & blood_res$padj < 0.05, ]

# Count number of DE genes for lung and blood cases
num_lung_DE_genes <- nrow(lung_sig)
num_blood_DE_genes <- nrow(blood_sig)

# Count up-regulated and down-regulated genes for lung and blood cases 
num_lung_upregulated <- sum(lung_sig$log2FoldChange > 0)
num_lung_downregulated <- sum(lung_sig$log2FoldChange < 0)
num_blood_upregulated <- sum(blood_sig$log2FoldChange > 0)
num_blood_downregulated <- sum(blood_sig$log2FoldChange < 0)

# Output the results
cat("Lung cases:\n")
cat("Number of DE genes: ", num_lung_DE_genes, "\n")
cat("Up-regulated genes: ", num_lung_upregulated, "\n")
cat("Down-regulated genes: ", num_lung_downregulated, "\n")
cat("\n")
cat("Blood cases:\n")
cat("Number of DE genes: ", num_blood_DE_genes, "\n")
cat("Up-regulated genes: ", num_blood_upregulated, "\n")
cat("Down-regulated genes: ", num_blood_downregulated, "\n")
cat("\n")


### Overrepresentation analysis

# Extract Ensembl IDs of DE genes for lung cases (padj < 0.05)
lung_DE_genes <- rownames(lung_sig)
# Extract Ensembl IDs of DE genes for blood cases (padj < 0.05)
blood_DE_genes <- rownames(blood_sig)

# Extract Ensembl IDs of all measured genes
all_genes <- rownames(countData)

# Perform overrepresentation analysis for lung cases 
lung_enrich <- enrichGO(gene = lung_DE_genes,
                        universe = all_genes,
                        OrgDb = org.Mm.eg.db,
                        ont = "ALL",
                        keyType = "ENSEMBL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

# Perform overrepresentation analysis for blood cases
blood_enrich <- enrichGO(gene = blood_DE_genes,
                         universe = all_genes,
                         OrgDb = org.Mm.eg.db,
                         ont = "ALL", 
                         keyType = "ENSEMBL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)

# Display top GO terms for lung cases 
print(head(lung_enrich))

# Display top GO terms for blood cases
print(head(blood_enrich))

# Convert enrichGO results to data frame and add Tissue column
lung_enrich_df <- as.data.frame(lung_enrich)
lung_enrich_df$Tissue <- "Lung"
blood_enrich_df <- as.data.frame(blood_enrich)
blood_enrich_df$Tissue <- "Blood"

# Select top 10 GO terms for each tissue
lung_top10 <- lung_enrich_df %>% top_n(10, wt = -p.adjust)
blood_top10 <- blood_enrich_df %>% top_n(10, wt = -p.adjust)

# Combine the top GO terms for both tissues
combined_enrich_df <- rbind(lung_top10, blood_top10)

# Create the combined bar plot using ggplot2
ggplot(combined_enrich_df, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = Tissue)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top GO Terms for Lung and Blood Cases",
       x = "GO Term",
       y = "-Log10 Adjusted P-Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

