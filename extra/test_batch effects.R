# ==============================================================================
# Test different batch scenarios
# ==============================================================================

### Prepare example data
suppressPackageStartupMessages({
    library(dplyr)
    library(edgeR)
    library(magrittr)
    library(purrr)
    library(SingleCellExperiment)
    library(forcats)
    library(data.table)
    library(CellMixS)
    library(ggplot2)
})

#Add dummy batch to example data
data(sce)
x <- sce
x$batch_id <- fct_collapse(x$sample_id,
                             batch1 = c("ctrl101", "stim101"),
                             batch2 = "ctrl107",
                             batch3 = "stim107")
colData(x) <- colData(x)[,!names(colData(x)) %in% "group_id"]
x <- prepSim(x)

table(x$sample_id, x$batch_id)

# ==============================================================================
# Test different batch scenarios
# 1. Estimated batch from reference data with more cluster than in reference:
# As the reference data have a random batch effect, it needs to be increased by 
# rel_be to get visible.
# ------------------------------------------------------------------------------
sim <- simData(x, nc = 1000, nb = 2, nk = 3,
               p_dd = c(0.9, 0, 0.1, 0, 0, 0),
               ng = 1000, force = TRUE,
               rel_be = c(4, 2), lfc_be = 0,
               probs = list(NULL, NULL, c(1, 0)))

p <- visGroup(sim, "batch_id")
p[["data"]]$cluster_id <- colData(sim)[, "cluster_id"]
p + aes(shape = cluster_id) +
    scale_shape_manual(values=1:nlevels(sim$cluster_id))


# 2. Simulated batch by setting lfc_be != 0:
# lfc_be is representing sd of the normal distribution to sample lfcs from
# rel_be_c is used to change the batch effect relative between cluster.
# Cluster1 has the strongest batch effect, cluster2 should not have a batch effect
# ------------------------------------------------------------------------------
sim <- simData(x, nc = 1000, nb = 2, nk = 3, ns = 3,
               p_dd = c(0.9, 0, 0.1, 0, 0, 0),
               ng = 1000, force = TRUE, rel_be_c = c(2, 0, 1), 
               lfc_be = 1.5, probs = list(NULL, NULL, c(1, 0)))

p <- visGroup(sim, "batch_id")
p[["data"]]$cluster_id <- colData(sim)[, "cluster_id"]
p + aes(shape = cluster_id) +
    scale_shape_manual(values=1:nlevels(sim$cluster_id))


# 3. Estimated with to few batches in the reference and lfc_be == 0 (default):
# This should result in no batch effect and give a warning.
# If rel_be or rel_be_c are specified this will break 
# (maybe another stop message is needed here).
# ------------------------------------------------------------------------------
sim <- simData(x, nc = 500, nb = 4, nk = 3,
               p_dd = c(0.9, 0, 0.1, 0, 0, 0),
               ng = 1000, force = TRUE,
               probs = list(NULL, NULL, c(1, 0)))

p <- visGroup(sim, "batch_id")
p[["data"]]$cluster_id <- colData(sim)[, "cluster_id"]
p + aes(shape = cluster_id) +
    scale_shape_manual(values=1:nlevels(sim$cluster_id))

# 4. Estimated with no batch_id column in reference and lfc_be == 0 (default):
# This should result in no batch effect and give a warning.
# If rel_be or rel_be_c are specified this will break 
# (maybe another stop message is needed here).
# ------------------------------------------------------------------------------
sce <- prepSim(sce)
table(sce$sample_id, sce$cluster_id)
sim <- simData(sce, nc = 500, nb = 4, nk = 3,
               p_dd = c(0.9, 0, 0.1, 0, 0, 0),
               ng = 1000, force = TRUE,
               probs = list(NULL, NULL, c(1, 0)))

p <- visGroup(sim, "batch_id")
p[["data"]]$cluster_id <- colData(sim)[, "cluster_id"]
p + aes(shape = cluster_id) +
    scale_shape_manual(values=1:nlevels(sim$cluster_id))


# 5. Negative control:
# There are different options to simulate a dataset without any batch effect:
# Most reasonable is to specify only one batch, otherwise lfc_be = 0 and no 
# reference batch 
# ------------------------------------------------------------------------------
sim <- simData(x, nc = 500, nb = 1, nk = 3,
               p_dd = c(0.9, 0, 0.1, 0, 0, 0),
               ng = 1000, force = TRUE,
               probs = list(NULL, NULL, c(1, 0)))

p <- visGroup(sim, "batch_id")
p[["data"]]$cluster_id <- colData(sim)[, "cluster_id"]
p + aes(shape = cluster_id) +
    scale_shape_manual(values=1:nlevels(sim$cluster_id))




