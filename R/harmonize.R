#' @title harmonize
#' 
#' @description Perform batch correction on a Seurat object using Harmony
#'
#' @param seuratObj Object to perform batch correction on
#' @param batch.var Variable to use when grouping samples.  Typically a column in seuratObj@meta.data
#' @param reduction.use Dimensional reduction to use when performing batch correction. If 'pca' is selected 
#' but not yet been performed, RunPCA will be ran before Harmony; otherwise, an error will be thrown. 
#' Default: 'pca'
#' @param reduction.save Name to use when saving the Harmony embeddings
#' @param ... Additional parameters to pass to HarmonyMatrix
#'
#' @importFrom glue glue
#' @importFrom Seurat RunPCA GetDimReduction SetDimReduction FetchData
#'
#' @return A Seurat object with the results of HarmonyMatrix saved in seuratObj@dr$reduction.save
#' @export
#'
#' @examples
harmonize <- function(seuratObj, 
                      batch.var, 
                      reduction.use = "pca", 
                      reduction.save = "harmony", 
                      ...){
  if (reduction.use == "pca"){
    if (!hasName(seuratObj, "pca")){
      seuratObj <- RunPCA(object = seuratObj, pcs.print = FALSE, )
    }
  } else if (!hasName(seuratObj, reduction.use)){
    stop(glue("{reduction.use} has not been performed."))
  }
  
  cell_embeddings <- GetDimReduction(seuratObj, 
                                     reduction.type = reduction.use, 
                                     slot = "cell.embeddings")
  batch_vector <- FetchData(seuratObj, vars.all = batch.var)[,1]
  he <- HarmonyMatrix(cell_embeddings, batch_vector, ...)
  rownames(he) <- rownames(cell_embeddings)
  colnames(he) <- glue("harmony_{1:ncol(he)}")
  seuratObj <- SetDimReduction(seuratObj, reduction.type = reduction.save, slot = "cell.embeddings", new.data = he)
  seuratObj <- SetDimReduction(seuratObj, reduction.type = reduction.save, slot = "key", new.data = reduction.save)
  return(seuratObj)
}