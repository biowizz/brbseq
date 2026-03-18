#!/usr/bin/env Rscript
# ===============================
# Step 7: Differential Expression
# DESeq2
# Comparisons:
#   1. RPRD1B: KD4h vs ctrl
#   2. RPRD2:  KD4h vs ctrl
# ===============================
# Run interactively or:
#   Rscript /scratch/sandeepm/pawel/7_deseq2.R
# ===============================

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(tibble)
    library(pheatmap)
    library(RColorBrewer)
})

set.seed(42)

# ---- Paths ----
counts_file <- "/scratch/sandeepm/pawel/results/counts/counts_raw.tsv"
out_dir     <- "/scratch/sandeepm/pawel/results/deseq2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Sample metadata ----
meta <- data.frame(
    sample = c(
        "1B_ctrl_1_F03","1B_ctrl_2_G03","1B_ctrl_3_H03",
        "1B_KD4h_1_A04","1B_KD4h_2_B04","1B_KD4h_3_C04",
        "D2_ctrl_1_G04", "D2_ctrl_2_H04","D2_ctrl_3_A05",
        "D2_KD4h_1_B05", "D2_KD4h_2_C05","D2_KD4h_3_D05"
    ),
    gene_target = c(rep("RPRD1B",6), rep("RPRD2",6)),
    condition   = c(rep("ctrl",3), rep("KD4h",3), rep("ctrl",3), rep("KD4h",3)),
    stringsAsFactors = FALSE
)

# ---- Load count matrix ----
counts <- read.table(counts_file, header=TRUE, sep="\t",
                     row.names=1, check.names=FALSE)
counts <- counts[, meta$sample]  # ensure column order matches metadata
cat("Count matrix loaded:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ==============================================================
# Helper: DESeq2 for one gene target
# ==============================================================
run_deseq2 <- function(target, counts_mat, metadata) {

    cat("\n--------------------------------------------\n")
    cat(" DESeq2:", target, "— KD4h vs ctrl\n")
    cat("--------------------------------------------\n")

    idx       <- metadata$gene_target == target
    meta_sub  <- metadata[idx, ]
    cnt_sub   <- counts_mat[, idx]

    meta_sub$condition <- factor(meta_sub$condition, levels=c("ctrl","KD4h"))

    dds <- DESeqDataSetFromMatrix(
        countData = round(cnt_sub),
        colData   = meta_sub,
        design    = ~ condition
    )

    # Pre-filter: at least 10 counts total
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    cat("Genes after pre-filter (>=10 counts):", nrow(dds), "\n")

    dds <- DESeq(dds)

    # LFC shrinkage with apeglm
    res <- lfcShrink(dds, coef="condition_KD4h_vs_ctrl",
                     type="apeglm", quiet=TRUE)

    res_df <- as.data.frame(res) %>%
        rownames_to_column("gene_id") %>%
        arrange(padj, desc(abs(log2FoldChange))) %>%
        mutate(
            sig       = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
            direction = case_when(
                sig & log2FoldChange >  0 ~ "UP",
                sig & log2FoldChange <  0 ~ "DOWN",
                TRUE                      ~ "NS"
            )
        )

    cat("DEGs (padj<0.05, |LFC|>1):",
        sum(res_df$direction=="UP",   na.rm=TRUE), "UP |",
        sum(res_df$direction=="DOWN", na.rm=TRUE), "DOWN\n")

    return(list(dds=dds, results=res_df))
}

# ==============================================================
# Run both comparisons
# ==============================================================
res_1B <- run_deseq2("RPRD1B", counts, meta)
res_D2 <- run_deseq2("RPRD2",  counts, meta)

# ---- Save DEG tables ----
write.csv(res_1B$results,
          file.path(out_dir, "DEGs_RPRD1B_KD4h_vs_ctrl.csv"),
          row.names=FALSE, quote=FALSE)

write.csv(res_D2$results,
          file.path(out_dir, "DEGs_RPRD2_KD4h_vs_ctrl.csv"),
          row.names=FALSE, quote=FALSE)

cat("\nDEG tables saved to:", out_dir, "\n")

# ==============================================================
# PCA — all 12 samples
# ==============================================================
cat("\nGenerating PCA plot...\n")

dds_all <- DESeqDataSetFromMatrix(
    countData = round(counts[rowSums(counts) >= 10, ]),
    colData   = meta %>%
                    mutate(condition   = factor(condition,   levels=c("ctrl","KD4h")),
                           gene_target = factor(gene_target)),
    design    = ~ gene_target + condition
)
dds_all <- estimateSizeFactors(dds_all)
vsd     <- vst(dds_all, blind=TRUE)

pca_data <- plotPCA(vsd, intgroup=c("gene_target","condition"), returnData=TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data,
                aes(x=PC1, y=PC2, color=gene_target, shape=condition, label=name)) +
    geom_point(size=4) +
    geom_text_repel(size=3, show.legend=FALSE) +
    scale_color_manual(values=c("RPRD1B"="#185FA5","RPRD2"="#D85A30")) +
    labs(title  = "PCA — all 12 samples (VST)",
         x      = sprintf("PC1: %d%% variance", pct_var[1]),
         y      = sprintf("PC2: %d%% variance", pct_var[2]),
         color  = "Gene target", shape = "Condition") +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold"))

ggsave(file.path(out_dir, "PCA_all_samples.pdf"), p_pca, width=7, height=5, dpi=150)

# ==============================================================
# Volcano plots
# ==============================================================
plot_volcano <- function(res_df, title, out_path) {
    df <- res_df %>%
        filter(!is.na(padj)) %>%
        mutate(neg_log10_padj = pmin(-log10(padj), 50))

    top <- df %>% filter(sig==TRUE) %>% arrange(padj) %>% head(15)

    p <- ggplot(df, aes(x=log2FoldChange, y=neg_log10_padj, color=direction)) +
        geom_point(size=0.8, alpha=0.7) +
        geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
        geom_vline(xintercept=c(-1,1),      linetype="dashed", color="grey40") +
        scale_color_manual(values=c("UP"="#D85A30","DOWN"="#185FA5","NS"="grey70")) +
        geom_text_repel(data=top, aes(label=gene_id),
                        size=2.5, max.overlaps=20, box.padding=0.3) +
        labs(title    = title,
             x        = "log2 Fold Change (KD4h / ctrl)",
             y        = "-log10(adjusted p-value)",
             color    = NULL,
             subtitle = sprintf("%d UP  |  %d DOWN",
                                sum(df$direction=="UP"),
                                sum(df$direction=="DOWN"))) +
        theme_bw(base_size=12) +
        theme(plot.title=element_text(face="bold"), legend.position="top")

    ggsave(out_path, p, width=7, height=6, dpi=150)
    cat("Saved:", out_path, "\n")
}

plot_volcano(res_1B$results,
             "RPRD1B — 4h KD vs ctrl",
             file.path(out_dir, "volcano_RPRD1B.pdf"))

plot_volcano(res_D2$results,
             "RPRD2 — 4h KD vs ctrl",
             file.path(out_dir, "volcano_RPRD2.pdf"))

# ==============================================================
# Shared DEGs — genes responding to BOTH knockdowns
# ==============================================================
cat("\n=== Shared DEGs ===\n")

deg_1B <- res_1B$results %>% filter(sig==TRUE) %>% pull(gene_id)
deg_D2 <- res_D2$results %>% filter(sig==TRUE) %>% pull(gene_id)
shared <- intersect(deg_1B, deg_D2)

cat("RPRD1B DEGs:", length(deg_1B), "\n")
cat("RPRD2  DEGs:", length(deg_D2), "\n")
cat("Shared     :", length(shared), "\n")

if (length(shared) > 0) {
    shared_df <- res_1B$results %>%
        filter(gene_id %in% shared) %>%
        select(gene_id, log2FoldChange, padj) %>%
        rename(LFC_RPRD1B=log2FoldChange, padj_RPRD1B=padj) %>%
        left_join(
            res_D2$results %>%
                filter(gene_id %in% shared) %>%
                select(gene_id, log2FoldChange, padj) %>%
                rename(LFC_RPRD2=log2FoldChange, padj_RPRD2=padj),
            by="gene_id"
        ) %>%
        arrange(padj_RPRD1B)

    write.csv(shared_df,
              file.path(out_dir, "DEGs_shared_both_KDs.csv"),
              row.names=FALSE, quote=FALSE)
    cat("Saved:", file.path(out_dir, "DEGs_shared_both_KDs.csv"), "\n")
}

cat("\nStep 7 complete. Outputs in:", out_dir, "\n")
