#!/usr/bin/env Rscript
# ===============================
# Step 8: Functional Enrichment
# GO Biological Process + KEGG
# clusterProfiler
# Runs on:
#   - RPRD1B UP / DOWN / ALL DEGs
#   - RPRD2  UP / DOWN / ALL DEGs
#   - Shared DEGs (both KDs)
# ===============================
# Run interactively or:
#   Rscript /scratch/sandeepm/pawel/8_enrichment.R
# ===============================

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(enrichplot)
    library(dplyr)
    library(readr)
})

# ---- Paths ----
deseq2_dir <- "/scratch/sandeepm/pawel/results/deseq2"
out_dir    <- "/scratch/sandeepm/pawel/results/enrichment"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load DEG results ----
res_1B <- read.csv(file.path(deseq2_dir, "DEGs_RPRD1B_KD4h_vs_ctrl.csv"))
res_D2 <- read.csv(file.path(deseq2_dir, "DEGs_RPRD2_KD4h_vs_ctrl.csv"))

# ==============================================================
# Helper: ENSEMBL -> ENTREZ conversion
# GENCODE IDs have version suffix (ENSG000001.12) — strip it
# ==============================================================
to_entrez <- function(ensembl_ids) {
    clean <- sub("\\..*", "", ensembl_ids)
    bitr(clean, fromType="ENSEMBL", toType="ENTREZID",
         OrgDb=org.Hs.eg.db, drop=TRUE)
}

# ==============================================================
# Helper: run GO BP + KEGG, save results and plots
# ==============================================================
run_enrichment <- function(deg_df, label, direction="ALL") {

    cat("\n--------------------------------------------\n")
    cat(" Enrichment:", label, "| direction:", direction, "\n")
    cat("--------------------------------------------\n")

    # Select genes
    if (direction == "ALL") {
        genes_use <- deg_df %>% filter(sig==TRUE) %>% pull(gene_id)
    } else {
        genes_use <- deg_df %>% filter(sig==TRUE, direction==direction) %>% pull(gene_id)
    }

    bg_genes <- deg_df %>% pull(gene_id)

    if (length(genes_use) < 5) {
        cat("Too few DEGs for enrichment (<5). Skipping.\n")
        return(invisible(NULL))
    }
    cat("Input DEGs:", length(genes_use), "\n")

    # ID conversion
    gene_map <- to_entrez(genes_use)
    bg_map   <- to_entrez(bg_genes)
    cat("Mapped to ENTREZ:", nrow(gene_map), "/", length(genes_use), "\n")

    # Output subdirectory
    subdir <- file.path(out_dir, paste0(label, "_", direction))
    dir.create(subdir, showWarnings=FALSE)

    # ---- GO Biological Process ----
    ego <- enrichGO(
        gene          = gene_map$ENTREZID,
        universe      = bg_map$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        write.csv(as.data.frame(ego),
                  file.path(subdir, "GO_BP.csv"), row.names=FALSE)

        p1 <- dotplot(ego, showCategory=20, font.size=9) +
            ggtitle(paste("GO: Biological Process |", label, direction)) +
            theme(plot.title=element_text(size=10, face="bold"))
        ggsave(file.path(subdir, "GO_BP_dotplot.pdf"), p1, width=9, height=8, dpi=150)

        ego2 <- pairwise_termsim(ego)
        p2   <- emapplot(ego2, showCategory=30) +
            ggtitle(paste("GO BP enrichment map |", label, direction))
        ggsave(file.path(subdir, "GO_BP_emapplot.pdf"), p2, width=10, height=9, dpi=150)

        cat("GO BP significant terms:", nrow(as.data.frame(ego)), "\n")
    } else {
        cat("No significant GO BP terms.\n")
    }

    # ---- KEGG ----
    ekegg <- enrichKEGG(
        gene          = gene_map$ENTREZID,
        universe      = bg_map$ENTREZID,
        organism      = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2
    )

    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        write.csv(as.data.frame(ekegg),
                  file.path(subdir, "KEGG.csv"), row.names=FALSE)

        p3 <- dotplot(ekegg, showCategory=20, font.size=9) +
            ggtitle(paste("KEGG pathways |", label, direction)) +
            theme(plot.title=element_text(size=10, face="bold"))
        ggsave(file.path(subdir, "KEGG_dotplot.pdf"), p3, width=9, height=7, dpi=150)

        cat("KEGG significant pathways:", nrow(as.data.frame(ekegg)), "\n")
    } else {
        cat("No significant KEGG pathways.\n")
    }

    cat("Saved to:", subdir, "\n")
}

# ==============================================================
# Run all enrichments
# ==============================================================

# RPRD1B
run_enrichment(res_1B, label="RPRD1B", direction="UP")
run_enrichment(res_1B, label="RPRD1B", direction="DOWN")
run_enrichment(res_1B, label="RPRD1B", direction="ALL")

# RPRD2
run_enrichment(res_D2, label="RPRD2",  direction="UP")
run_enrichment(res_D2, label="RPRD2",  direction="DOWN")
run_enrichment(res_D2, label="RPRD2",  direction="ALL")

# Shared DEGs (both knockdowns)
shared_file <- file.path(deseq2_dir, "DEGs_shared_both_KDs.csv")
if (file.exists(shared_file)) {
    shared <- read.csv(shared_file) %>%
        mutate(
            sig       = TRUE,
            direction = ifelse(LFC_RPRD1B > 0, "UP", "DOWN"),
            gene_id   = gene_id
        )
    run_enrichment(shared, label="SHARED", direction="ALL")
    run_enrichment(shared, label="SHARED", direction="UP")
    run_enrichment(shared, label="SHARED", direction="DOWN")
}

cat("\nStep 8 complete. All enrichment results in:", out_dir, "\n")
