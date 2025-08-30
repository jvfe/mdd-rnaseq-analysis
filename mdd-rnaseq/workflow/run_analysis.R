library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

metadata_path <- "data/metadata.csv"
if (!file.exists(metadata_path)) {
    stop("Error: Metadata file not found at '", metadata_path, "'.")
}
full_meta <- read.csv(metadata_path)

meta <- full_meta[full_meta$condition %in% c("CTRL_male", "MDD_male"), ]
meta$condition <- factor(meta$condition) # Drop unused factor levels

print("Filtered metadata for analysis:")
print(meta)

gene_counts_path <- "results/tximport/gene_counts.tsv"
if (!file.exists(gene_counts_path)) {
    stop("Error: Gene counts file not found at '", gene_counts_path, "'.")
}
gene_counts <- read.table(gene_counts_path, header = TRUE, row.names = 1)
colnames(gene_counts) <- gsub("_T1", "", colnames(gene_counts))
gene_counts <- gene_counts[, meta$sample]

print("Gene counts matrix dimensions after filtering:")
print(dim(gene_counts))

dge_gene <- DGEList(counts = gene_counts, group = meta$condition)

min_samples <- min(table(meta$condition))
keep_gene <- rowSums(cpm(dge_gene) > 1) >= min_samples
dge_gene <- dge_gene[keep_gene, , keep.lib.sizes = FALSE]

dge_gene <- calcNormFactors(dge_gene)

design <- model.matrix(~0 + group, data = dge_gene$samples)
colnames(design) <- levels(dge_gene$samples$group)
dge_gene <- estimateDisp(dge_gene, design)

contrast <- makeContrasts(MDD_male - CTRL_male, levels = design)
fit_gene <- glmQLFit(dge_gene, design)
qlf_gene <- glmQLFTest(fit_gene, contrast = contrast)

top_degs <- topTags(qlf_gene, n = Inf)$table
top_degs$significant <- top_degs$FDR < 0.05
results_file_gene <- "results/edgeR_MDDmale_vs_CTRLmale_GENE_results.csv"
write.csv(top_degs, file = results_file_gene)
message("Gene-level DE results saved to '", results_file_gene, "'")


tx_counts_path <- "results/tximport/transcript_counts.tsv"
if (!file.exists(tx_counts_path)) {
    stop("Error: Transcript counts file not found at '", tx_counts_path, "'.")
}
tx_counts <- read.table(tx_counts_path, header = TRUE, row.names = 1)
colnames(tx_counts) <- gsub("_T1", "", colnames(tx_counts))
tx_counts <- tx_counts[, meta$sample]

print("Transcript counts matrix dimensions after filtering:")
print(dim(tx_counts))

dge_tx <- DGEList(counts = tx_counts, group = meta$condition)

keep_tx <- rowSums(cpm(dge_tx) > 1) >= min_samples
dge_tx <- dge_tx[keep_tx, , keep.lib.sizes = FALSE]

dge_tx <- calcNormFactors(dge_tx)

dge_tx <- estimateDisp(dge_tx, design)

fit_tx <- glmQLFit(dge_tx, design)
qlf_tx <- glmQLFTest(fit_tx, contrast = contrast)

top_dets <- topTags(qlf_tx, n = Inf)$table
top_dets$significant <- top_dets$FDR < 0.05
results_file_tx <- "results/edgeR_MDDmale_vs_CTRLmale_TRANSCRIPT_results.csv"
write.csv(top_dets, file = results_file_tx)
message("Transcript-level DE results saved to '", results_file_tx, "'")



sig_genes <- rownames(top_degs[top_degs$significant, ])

length(sig_genes) # q3

if (length(sig_genes) > 0) {
    message(paste("Found", length(sig_genes), "significant genes for enrichment analysis."))
    entrez_ids_gene <- bitr(sig_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    go_enrich_gene <- enrichGO(gene = entrez_ids_gene$ENTREZID, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    
    if (!is.null(go_enrich_gene) && nrow(go_enrich_gene@result) > 0) {
        go_res_gene <- go_enrich_gene@result
        go_res_gene$Description[which.min(go_res_gene$p.adjust)] #q5
        write.csv(go_res_gene, "results/GO_enrichment_GENE_results.csv")
        go_plot_gene <- dotplot(go_enrich_gene, showCategory = 20) + ggtitle("GO Enrichment (Gene Level)")
        ggsave("results/GO_enrichment_GENE_dotplot.png", plot = go_plot_gene, width = 10, height = 8)
    } else {
        message("No significant GO terms found for gene-level results.")
    }
} else {
    message("No significant genes found. Skipping gene-level enrichment.")
}


sig_transcripts <- rownames(top_dets[top_dets$significant, ])
sig_transcripts_no_version <- gsub("\\..*$", "", sig_transcripts)

length(sig_transcripts_no_version) #q4

if (length(sig_transcripts) > 0) {
    message(paste("Found", length(sig_transcripts), "significant transcripts for enrichment analysis."))
    entrez_ids_tx <- bitr(sig_transcripts_no_version, fromType = "ENSEMBLTRANS", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    go_enrich_tx <- enrichGO(gene = entrez_ids_tx$ENTREZID, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    
    if (!is.null(go_enrich_tx) && nrow(go_enrich_tx@result) > 0) {
        write.csv(go_enrich_tx@result, "results/GO_enrichment_TRANSCRIPT_results.csv")
    } else {
        message("No significant GO terms found for transcript-level results.")
    }
} else {
    message("No significant transcripts found. Skipping transcript-level enrichment.")
}

message("Script finished successfully!")
