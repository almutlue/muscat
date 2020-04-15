###----- Helper functions for batch effect 

# ==============================================================================
# Estimate genewise and clusterwise logFC of different batches. Runs edgeR's 
# glmQLFTest on each cluster comparing all batches against the same reference 
# batch. 
# ------------------------------------------------------------------------------
.get_batch_lfc <- function(sce){
    kids <- purrr::set_names(levels(sce$cluster_id))
    bids <- purrr::set_names(levels(sce$batch_id))
    ref <- bids[1]    # reference batch
    com <- bids[-1]   # comparison batches
    expr <- counts(sce)
    #Estimate clusterwise logFC between batches
    lfc_k <- lapply(kids, function (k) {
        cat(k, "..", sep = "")
        n <- sce$cluster_id == k
        es_tmp <- expr[, n]
        grp <- sce$batch_id[n]
        skip <- bids[table(grp) < 10]
        lfc_kb <- list()
        com <- com[!com %in% skip]
        #make sure enough cells from com and ref are in cluster  
        if (length(com) > 0 & !ref %in% skip) {
            dge <- DGEList(es_tmp, group = grp)
            dge <- calcNormFactors(dge)
            design <- model.matrix(~ 0 + grp)
            colnames(design) <- bids
            dge <- estimateDisp(dge, design = design)
            f <- glmQLFit(dge, design = design, robust = TRUE)
            lfc_kb <- lapply(com, function(cc){
                # defining contrast
                cont <- rep(0, length(bids))
                cont[bids %in% ref] <- 1
                cont[bids %in% cc] <- -1
                f <- glmQLFTest(f, contrast = cont)
                f <- f[rownames(sce),]
                res <- data.frame(logFC = f[["table"]]$logFC) %>% 
                    set_colnames(paste0(ref, "-", cc, "_logFC_", k))
                res
            }) %>% bind_cols()
        #Add mean lfc for cluster that were skipped due to cell number
        res_skip <- lapply(skip, function(cc){
            data.frame(logFC = rowMeans(lfc_kb)) %>% 
            set_colnames(paste0(ref, "-", cc, "_logFC_", k))
        }) %>% bind_cols()
        if (ncol(res_skip) > 0) {
            lfc_kb <- cbind(lfc_kb, res_skip)   
        }
        }else{
            # no batch effect in case not enough cells from the reference batch 
            # are in cluster (maybe add warning here)
            nam_exclude <- paste0(ref, "-", skip, "_logFC_", k)
            lfc_kb <- skip %>% map(., rep, x = 0, times = nrow(sce)) %>% 
                bind_cols() %>% 
                set_names(nam_exclude)
        }
        return(lfc_kb)
    }) %>% bind_cols()
    
    rd <- cbind(rowData(sce), lfc_k)
    rowData(sce)[,colnames(rd)] <- rd
    sce
}



# ==============================================================================
# sample lfc from a normal distribution
#  
#  
# ------------------------------------------------------------------------------
.sim_lfc_be <- function(lfc_be, gs){
    n <- length(gs)
    signs <- sample(c(-1, 1), n, TRUE)
    lfcs <- rnorm(n,
                  mean = 0,
                  sd = lfc_be) * signs
    names(lfcs) <- gs
    return(lfcs)
}


# ==============================================================================
# Multiply lfcs from estimated ones for more cluster.
# generate new estimates from the medians and standraddeviations of the 
# estimated lfcs to get a different estimate for all cluster.
# ------------------------------------------------------------------------------
.est_lfc_batch <- function(lfcb, kids){
    nk <- length(kids)
    lfc_b <- as.matrix(rowData(x)[, lfcb])
    if (length(lfcb) < nk) {
        med <- rowMedians(lfc_b)
        sds <- rowSds(lfc_b)
        lfcb_est <- lapply(seq_len(nk), function(k){
            signs <- sample(c(-1, 1), length(med), TRUE)
            lfc <- med + k * signs * sds
        }) %>% bind_cols()
        lfc_b <- cbind(lfc_b, lfcb_est)
    }
    lfc_b <- lfc_b[, seq_len(nk)] %>% set_names(paste0("logFC_", kids))
}

