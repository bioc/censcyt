#' Run \code{censcyt} pipeline
#' 
#' Wrapper function to run complete \code{censcyt} pipeline
#' 
#' This wrapper function runs the complete  \code{\link[diffcyt]{diffcyt}} analysis 
#' pipeline where the only difference is the analysis step which uses the functions
#' from \code{censcyt}
#' 
#' For more details about the functions for the individual steps, see 
#' \code{\link[diffcyt]{diffcyt}}, the \code{\link{diffcyt}} vignette, 
#' the \code{censcyt} package vignette and the function help pages. The following
#' is a slightly adapted summary from \code{\link[diffcyt]{diffcyt}}:
#' 
#' Running the individual functions may provide additional
#' flexibility, especially for complex analyses.
#' 
#' The input data can be provided as a \code{\link[flowCore]{flowSet}} or a list of
#' \code{\link[flowCore]{flowFrame}s}, \code{\link{DataFrame}s}, \code{data.frames}, or matrices
#' (one \code{flowFrame} or list item per sample). Alternatively, it is also possible to
#' provide the input as a \code{daFrame} object from the \code{CATALYST} Bioconductor
#' package (Chevrier, Crowell, Zanotelli et al., 2018). This can be useful when initial
#' exploratory analyses and clustering have been performed using \code{CATALYST}; the
#' \code{daFrame} object from \code{CATALYST} (containing cluster labels in the
#' \code{rowData}) can then be provided directly to the \code{censcyt} functions for
#' differential testing.
#' 
#' Minimum required arguments when not providing a \code{\link[flowCore]{flowSet}} or list of
#' \code{\link[flowCore]{flowFrame}s}, \code{\link{DataFrame}s}, \code{data.frames}, or matrices:
#' 
#' \itemize{
#' \item \code{d_input}
#' \item \code{experiment_info}
#' \item \code{marker_info}
#' \item either \code{design} or \code{formula} (depending on the differential testing
#' method used)
#' \item \code{contrast}
#' \item \code{analysis_type}
#' }
#' 
#' Minimum required arguments when providing a \code{CATALYST} \code{daFrame} object:
#' 
#' \itemize{
#' \item \code{d_input}
#' \item either \code{design} or \code{formula} (depending on the differential testing
#' method used)
#' \item \code{contrast}
#' \item \code{analysis_type}
#' }
#' 
#' 
#' @param d_input Input data. Must be either: (i) a \code{\link[flowCore]{flowSet}} or list of
#'   \code{\link[flowCore]{flowFrame}s}, \code{\link{DataFrame}s}, \code{data.frames}, or matrices
#'   as input (one \code{flowFrame} or list item per sample) (see
#'   \code{\link{prepareData}}); or (ii) a \code{CATALYST} \code{daFrame} (containing
#'   cluster labels in \code{rowData}; see vignette for an example).
#' 
#' @param experiment_info \code{data.frame}, \code{\link{DataFrame}}, or
#'   \code{\link{tbl_df}} of experiment information, for example sample IDs and group IDs.
#'   Must contain a column named \code{sample_id}. See \code{\link[diffcyt]{prepareData}}. (Not
#'   required when providing a \code{CATALYST} \code{daFrame} for \code{d_input}.)
#' 
#' @param marker_info \code{data.frame}, \code{\link{DataFrame}}, or \code{\link{tbl_df}}
#'   of marker information for each column of data. This should contain columns named
#'   \code{marker_name} and \code{marker_class}. The columns contain: (i) marker names
#'   (and any other column names); and (ii) a factor indicating the marker class for each
#'   column (with entries \code{"type"}, \code{"state"}, or \code{"none"}). See
#'   \code{\link[diffcyt]{prepareData}}. (Not required when providing a \code{CATALYST}
#'   \code{daFrame} for \code{d_input}.)
#' 
#' @param design Design matrix, created with \code{\link{createDesignMatrix}}. See
#'   \code{\link[diffcyt]{createDesignMatrix}}.
#' 
#' @param formula Model formula object, created with \code{\link{createFormula}}. See
#'   \code{\link{createFormula}}.
#' 
#' @param contrast Contrast matrix, created with \code{\link{createContrast}}. See
#'   \code{\link[diffcyt]{createContrast}}.
#' 
#' @param analysis_type Type of differential analysis to perform: differential abundance
#'   (DA) of cell populations. The only option at the moment is \code{"DA"}. 
#'   See \code{\link{testDA_censoredGLMM}}.
#' 
#' @param method_DA Method to use for calculating differential abundance (DA) tests.
#'   Currently the only option is \code{\link{testDA_censoredGLMM}}. 
#'   Default = \code{\link{testDA_censoredGLMM}}.
#' 
#' @param markers_to_test (Optional) Logical vector specifying which markers to test for
#'   differential expression (from the set of markers stored in the \code{assays} of
#'   \code{d_medians}; for method \code{\link{testDS_limma}} or \code{\link{testDS_LMM}}).
#'   Default = all 'cell state' markers, which are identified by the logical vector
#'   \code{id_state_markers} stored in the meta-data of \code{d_medians}. See
#'   \code{\link{testDS_limma}} or \code{\link{testDS_LMM}}.
#' 
#' @param clustering_to_use (Optional) Column name indicating which set of cluster labels
#'   to use for differential testing, when input data are provided as a \code{CATALYST}
#'   \code{daFrame} object containing multiple sets of cluster labels. (In this case, the
#'   \code{metadata} of the \code{daFrame} object is assumed to contain a data frame named
#'   \code{cluster_codes}; \code{clustering_to_use} is the column name of the selected
#'   column in \code{cluster_codes}. If \code{clustering_to_use} is provided, an
#'   identifier \code{clustering_name} to identify this column will also be saved in the
#'   \code{metadata} of the output object.) Default = NULL, in which case cluster labels
#'   stored in column named \code{cluster_id} in the \code{rowData} of the \code{daFrame}
#'   object are used.
#' 
#' @param cols_to_include Logical vector indicating which columns to include from the
#'   input data. Default = all columns. See \code{\link[diffcyt]{prepareData}}.
#' 
#' @param subsampling Whether to use random subsampling to select an equal number of cells
#'   from each sample. Default = FALSE. See \code{\link[diffcyt]{prepareData}}.
#' 
#' @param n_sub Number of cells to select from each sample by random subsampling, if
#'   \code{subsampling = TRUE}. Default = number of cells in smallest sample. See
#'   \code{\link[diffcyt]{prepareData}}.
#' 
#' @param seed_sub Random seed for subsampling. Set to an integer value to generate
#'   reproducible results. Default = \code{NULL}. See \code{\link[diffcyt]{prepareData}}.
#' 
#' @param transform Whether to apply 'arcsinh' transform. This may be set to FALSE if the
#'   input data has already been transformed. Default = TRUE. See
#'   \code{\link[diffcyt]{transformData}}.
#' 
#' @param cofactor Cofactor parameter for 'arcsinh' transform. Default = 5, which is
#'   appropriate for mass cytometry (CyTOF) data. For fluorescence flow cytometry, 
#'   cofactor = 150 is recommended instead. See \code{\link[diffcyt]{transformData}}.
#' 
#' @param cols_clustering Columns to use for clustering. Default = \code{NULL}, in which
#'   case markers identified as 'cell type' markers (with entries \code{"type"}) in the
#'   vector \code{marker_class} in the column meta-data of \code{d_se} will be used. See
#'   \code{\link[diffcyt]{generateClusters}}.
#' 
#' @param xdim Horizontal length of grid for self-organizing map for FlowSOM clustering
#'   (number of clusters = \code{xdim} * \code{ydim}). Default = 10 (i.e. 100 clusters).
#'   See \code{\link[diffcyt]{generateClusters}}.
#' 
#' @param ydim Vertical length of grid for self-organizing map for FlowSOM clustering
#'   (number of clusters = \code{xdim} * \code{ydim}). Default = 10 (i.e. 100 clusters).
#'   See \code{\link[diffcyt]{generateClusters}}.
#' 
#' @param meta_clustering Whether to include FlowSOM 'meta-clustering' step. Default =
#'   \code{FALSE}. See \code{\link[diffcyt]{generateClusters}}.
#' 
#' @param meta_k Number of meta-clusters for FlowSOM, if \code{meta-clustering = TRUE}.
#'   Default = 40. See \code{\link[diffcyt]{generateClusters}}.
#' 
#' @param seed_clustering Random seed for clustering. Set to an integer value to generate
#'   reproducible results. Default = \code{NULL}. See \code{\link[diffcyt]{generateClusters}}.
#' 
#' @param min_cells Filtering parameter. Default = 3. Clusters are kept for differential
#'   testing if they have at least \code{min_cells} cells in at least \code{min_samples}
#'   samples. See \code{\link{testDA_censoredGLMM}}.
#' 
#' @param min_samples Filtering parameter. Default = \code{number of samples / 2}, which
#'   is appropriate for two-group comparisons (of equal size). Clusters are kept for
#'   differential testing if they have at least \code{min_cells} cells in at least
#'   \code{min_samples} samples. See \code{\link{testDA_censoredGLMM}}.
#' 
#' @param normalize Whether to include optional normalization factors to adjust for
#'   composition effects. Default = FALSE. See \code{\link{testDA_censoredGLMM}}.
#' 
#' @param norm_factors Normalization factors to use, if \code{normalize = TRUE}. Default =
#'   \code{"TMM"}, in which case normalization factors are calculated automatically using
#'   the 'trimmed mean of M-values' (TMM) method from the \code{edgeR} package.
#'   Alternatively, a vector of values can be provided (the values should multiply to 1).
#'   See \code{\link{testDA_censoredGLMM}}.
#' 
#' @param verbose Whether to print status messages during each step of the pipeline.
#'   Default = TRUE.
#'   
#' @param mi_reps Number of imputations in multiple imputation. 
#' Default = 10. See \code{\link{testDA_censoredGLMM}}.
#'
#' @param imputation_method Method to be used in the imputation step. 
#' One of \code{km},\code{km_exp},\code{km_wei},\code{km_os}, \code{rs}, \code{mrl},
#' \code{cc}, \code{pmm}. See \code{\link{testDA_censoredGLMM}}.
#'  
#' @param BPPARAM Specification of parallelization option as one of
#'  \code{\link[BiocParallel]{BiocParallelParam}} if \code{BiocParallel} is available 
#'  otherwise no parallelization is used. e.g. \code{\link[BiocParallel]{MulticoreParam}}(workers=2) 
#'  for parallelization  with two cores. Default is \code{\link[BiocParallel]{SerialParam}}()
#'  (no parallelization). 
#' 
#' 
#' @return Returns a list containing the results object \code{res}, as well as the data
#'   objects \code{d_se}, \code{d_counts}, \code{d_medians},
#'   \code{d_medians_by_cluster_marker}, and \code{d_medians_by_sample_marker}. (If a
#'   \code{CATALYST} \code{daFrame} object was used as input, the output list contains
#'   objects \code{res}, \code{d_counts}, and \code{d_medians}.) 
#' 
#' 
#' @aliases censcyt-package
#' 
#' @importFrom SummarizedExperiment assays rowData colData 'colData<-' SummarizedExperiment
#' @importFrom S4Vectors metadata 'metadata<-'
#' @import diffcyt
#' @importFrom tidyr pivot_longer
#' 
#' @export
#' 
#' @examples
#' # Function to create random data (one sample)
#' fcs_sim <- function(n = 2000, mean = 0, sd = 1, ncol = 10, cofactor = 5) {
#'   d <- matrix(sinh(rnorm(n*ncol, mean, sd)) * cofactor,ncol=ncol)
#'   for(i in seq_len(ncol)){
#'     d[seq(n/ncol*(i-1)+1,n/ncol*(i)),i] <- sinh(rnorm(n/ncol, mean+5, sd)) * cofactor
#'   }
#'   colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
#'   d
#' }
#' 
#' # Create random data (without differential signal)
#' set.seed(123)
#' d_input <- lapply(1:50, function(i) fcs_sim())
#' 
#' # simulate survival time
#' d_surv <- simulate_singlecluster(50, formula(Y~Surv(X,I)))[c("X","I","TrVal")]
#' 
#' # Add differential abundance (DA) signal
#' for(i in 1:50){
#'   # number of cells in cluster 1 
#'   n_da <- round(sqrt(2000*d_surv$TrVal[i]))*9
#'   # set to no expression
#'   tmpd <- matrix(sinh(rnorm(n_da*10, 0, 1)) * 5, ncol=10)
#'   # increase expresion for cluster 1
#'   tmpd[ ,1] <- sinh(rnorm(n_da, 5, 1)) * 5
#'   d_input[[i]][seq_len(n_da), ] <- tmpd
#' }
#' 
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:50)),
#'   survival_time = d_surv$X,
#'   event_indicator= d_surv$I,
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   channel_name = paste0("channel", sprintf("%03d", 1:10)),
#'   marker_name = paste0("marker", sprintf("%02d", 1:10)),
#'   marker_class = factor(c(rep("type", 10)),
#'                         levels = c("type", "state", "none")),
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Create formula
#' da_formula <- createFormula(experiment_info, cols_fixed="survival_time",
#'                             cols_random = "sample_id",event_indicator = "event_indicator")
#' # Create contrast matrix
#' contrast <- diffcyt::createContrast(c(0, 1))
#' 
#' # Test for differential abundance (DA) of clusters
#' out_DA <- censcyt(d_input, experiment_info, marker_info,
#'                   formula = da_formula, contrast = contrast,
#'                   analysis_type = "DA", method_DA = "censcyt-DA-censored-GLMM",
#'                   seed_clustering = 123, verbose = FALSE, mi_reps = 3,
#'                   BPPARAM=BiocParallel::MulticoreParam(workers = 1),
#'                   imputation_method = "mrl",meta_clustering = TRUE, meta_k = 10)
#' 
#' # Display results for top DA clusters
#' diffcyt::topTable(out_DA, format_vals = TRUE)
#' 
#' 
#' # Plot heatmap for DA tests
#' diffcyt::plotHeatmap(out_DA, analysis_type = "DA")
#' 
censcyt <- function(d_input, experiment_info = NULL, marker_info = NULL, 
                    design = NULL, formula = NULL, contrast, 
                    analysis_type = c("DA"), 
                    method_DA = c("censcyt-DA-censored-GLMM"), 
                    markers_to_test = NULL, 
                    clustering_to_use = NULL, 
                    cols_to_include = NULL, 
                    subsampling = FALSE, n_sub = NULL, seed_sub = NULL, 
                    transform = TRUE, cofactor = 5, 
                    cols_clustering = NULL, xdim = 10, ydim = 10, 
                    meta_clustering = FALSE, meta_k = 40, seed_clustering = NULL, 
                    min_cells = 3, min_samples = NULL, 
                    normalize = FALSE, norm_factors = "TMM", 
                    verbose = TRUE, mi_reps = 10, 
                    imputation_method = c("km","km_exp","km_wei","km_os","rs","mrl","cc","pmm"),
                    BPPARAM=BiocParallel::SerialParam()) {
  # most of the following code is taken from the diffcyt wrapper function
  # get arguments
  analysis_type <- match.arg(analysis_type)
  method_DA <- match.arg(method_DA)
  imputation_method <- match.arg(imputation_method)
  
  # preliminary steps (if input is not a SingleCellExperiment object from CATALYST)
  if (!is(d_input, "SingleCellExperiment")) {
    if (is.null(experiment_info) | is.null(marker_info)) {
      stop("'experiment_info' and 'marker_info' must be provided (unless using a SingleCellExperiment ", 
           "object from CATALYST as input)")
    }
    # prepare data
    if (verbose) message("preparing data...")
    d_se <- prepareData(d_input, experiment_info, marker_info, cols_to_include, subsampling, n_sub, seed_sub)
    # transformation
    if (transform) {
      if (verbose) message("transforming data...")
      d_se <- transformData(d_se, cofactor)
    }
    # clustering
    if (verbose) message("generating clusters...")
    d_se <- generateClusters(d_se, cols_clustering, xdim, ydim, meta_clustering, meta_k, seed_clustering)
  }
  
  # alternatively, use SingleCellExperiment object from CATALYST (which already contains cluster labels) as input
  else if (is(d_input, "SingleCellExperiment")) {
    if (verbose) message("using SingleCellExperiment object from CATALYST as input")
    
    # select clustering to use
    
    if (is.null(clustering_to_use)) {
      stopifnot("cluster_id" %in% colnames(colData(d_input)))
      if (verbose) message("using cluster IDs stored in column named 'cluster_id' in 'colData' of ", 
                           "SingleCellExperiment object from CATALYST")
      # clustering identifier to store in metadata
      clustering_name <- colnames(metadata(d_input)$cluster_codes)[1]
      
    } else if (!is.null(clustering_to_use)) {
      stopifnot(as.character(clustering_to_use) %in% colnames(metadata(d_input)$cluster_codes))
      stopifnot("cluster_id" %in% colnames(colData(d_input)))
      if (verbose) message("using cluster IDs from clustering stored in column '", clustering_to_use, 
                           "' of 'cluster_codes' data frame in 'metadata' of SingleCellExperiment object from CATALYST")
      code_id <- colData(d_input)$cluster_id
      cluster_id <- metadata(d_input)$cluster_codes[, clustering_to_use][code_id]
      # store cluster labels in column 'cluster_id' in 'colData'; and add column 'code_id'
      # for original FlowSOM clustering codes
      stopifnot(length(cluster_id) == nrow(colData(d_input)), 
                length(code_id) == nrow(colData(d_input)))
      colData(d_input)$code_id <- code_id
      colData(d_input)$cluster_id <- cluster_id
      # clustering identifier to store in metadata
      clustering_name <- clustering_to_use
    }
    
    # unpack SingleCellExperiment (proteins x cells format) and create SummarizedExperiment (cells x proteins format)
    
    stopifnot("sample_id" %in% colnames(colData(d_input)))
    stopifnot("cluster_id" %in% colnames(colData(d_input)))
    stopifnot("cluster_codes" %in% names(metadata(d_input)))
    
    if (!("experiment_info" %in% names(metadata(d_input)))) {
      if (verbose) message("generating 'experiment_info' from input object")
      m <- match(levels(droplevels(factor(d_input$sample_id))), d_input$sample_id)
      experiment_info <- data.frame(colData(d_input)[m, ], check.names = FALSE, row.names = NULL)
      metadata <- as.list(c(metadata(d_input), experiment_info))
    } else {
      experiment_info <- metadata(d_input)$experiment_info
      metadata <- metadata(d_input)
    }
    
    # split cells by sample
    cs_by_s <- split(seq_len(ncol(d_input)), colData(d_input)$sample_id)
    # re-order according to experiment_info
    cs <- unlist(cs_by_s[as.character(experiment_info$sample_id)])
    es <- t(assays(d_input)[["exprs"]])[cs, , drop = FALSE]
    # create SummarizedExperiment (in transposed format compared to SingleCellExperiment)
    d_se <- SummarizedExperiment(
      assays = list(exprs = es), 
      rowData = colData(d_input)[cs, ], 
      colData = rowData(d_input), 
      metadata = metadata
    )
  }
  
  # calculate features
  if (verbose) message("calculating features...")
  d_counts <- calcCounts(d_se)
  
  # DA tests
  if (analysis_type == "DA" && method_DA == "censcyt-DA-censored-GLMM") {
    if (verbose) message("calculating DA tests using method 'censcyt-DA-censored-GLMM'...")
    res <- testDA_censoredGLMM(d_counts, formula, contrast, mi_reps, imputation_method, min_cells, min_samples, normalize, norm_factors, BPPARAM)
  }
  # return results and data objects
  if (!is(d_input, "SingleCellExperiment")) {
    return(list(
      res = res, 
      d_se = d_se, 
      d_counts = d_counts
    ))
  } else if (is(d_input, "SingleCellExperiment")) {
    # store clustering identifier in metadata
    metadata(res) <- as.list(c(metadata(res), clustering_name = clustering_name))
    # not returning input object, since it has been modified
    return(list(
      res = res, 
      d_counts = d_counts
    ))
  } 
}


