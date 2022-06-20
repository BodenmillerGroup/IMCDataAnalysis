library(SingleCellExperiment)
library(stringr)
label_files <- list.files("data/gated_cells", 
                          full.names = TRUE, pattern = ".rds$")

# Read in SPE objects
spes <- lapply(label_files, readRDS)
names(spes) <- list.files("data/gated_cells", pattern = ".rds$")

# Read in current SPE object
spe <- readRDS("data/spe.rds")

new_spes <- lapply(spes, function(x){
        cur_spe <- spe[,spe$sample_id == unique(x$sample_id)]  
        cur_gates <- metadata(x)[grepl("cytomapper_gate", names(metadata(x)))]
        cur_gates <- cur_gates[order(as.numeric(str_split(names(cur_gates), "_", simplify = TRUE)[,3]), decreasing = FALSE)]
        
        cur_meta <- metadata(cur_spe)
        metadata(cur_spe) <- list()
        metadata(cur_spe)$metadata <- cur_meta
        
        for (i in 1:length(cur_gates)) {
            gate <- cur_gates[[i]]
            for (j in 1:nrow(gate$gate)){
                cur_val <- assay(cur_spe, gate$exprs_values)[rownames(gate$gate)[j],]
                cur_spe <- cur_spe[,cur_val > gate$gate[j,1] & cur_val < gate$gate[j,2]]
            }
            metadata(cur_spe)[[names(cur_gates)[i]]] <- gate
        }
        
        cur_spe$cytomapper_CellLabel <- unique(x$cytomapper_CellLabel)
        
        metadata(cur_spe)$cytomapper_SessionInfo <- metadata(x)$cytomapper_SessionInfo
        metadata(cur_spe)$cytomapper_GatingDate <- metadata(x)$cytomapper_GatingDate
        
        return(cur_spe)
})

lapply(1:length(new_spes), function(x){
    saveRDS(new_spes[[x]], file = paste0("data/gated_cells/", names(spes)[x]))
})
