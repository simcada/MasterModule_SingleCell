#------------------------------------------------------------------
# subset of scCustomize utility functions 
#------------------------------------------------------------------

import::here("dplyr",  "mutate", .character_only=TRUE)
import::here("magrittr", "%>%", .character_only=TRUE)
import::here("glue", "glue_collapse", .character_only=TRUE)
import::here("purrr", "reduce", .character_only=TRUE)
import::here(
  "stringr",
  c("str_to_sentence", "str_to_upper"),
  .character_only=TRUE
)
import::here(
  "Seurat",
  c("Assays", "DefaultAssay", "GetAssayData"),
  .character_only=TRUE
)


#' Set a default value if an object is NOT null
#'
#' @param lhs An object to set if it's NOT null
#' @param rhs The value to provide if x is NOT null
#'
#' @return lhs if lhs is null, else rhs
#'
#' @author Hadley Wickham
#' @references https://adv-r.hadley.nz/functions.html#missing-arguments
#'
#' @noRd
#'
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}


#' Set a default value if an object is null
#'
#' @param lhs An object to set if it's null
#' @param rhs The value to provide if x is null
#'
#' @return rhs if lhs is null, else lhs
#'
#' @author Hadley Wickham
#' @references https://adv-r.hadley.nz/functions.html#missing-arguments
#'
#' @noRd
#'
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}


#' Check Seurat Object
#'
#' Checks if object is of class: Seurat and returns error message if not.
#'
#' @param seurat_object Seurat object name.
#'
#' @return stops function without error message
#'
#' @noRd
#'
Is_Seurat <- function(
    seurat_object
) {
  if (class(x = seurat_object)[[1]] != "Seurat") {
    stop("'seurat_object' provided is not an object of class: Seurat.")
  }
}


#' Automatically calculate a point size for ggplot2-based scatter plots
#
#' It happens to look good
#'
#' @param data a single value length vector corresponding to the number of cells.
#' @param raster If TRUE, point size is set to 1
#'
#' @return The "optimal" point size for visualizing these data
#'
#' @noRd
#'
#' @references This function and documentation text are modified versions of the `AutoPointSize` function
#' and documentation from Seurat (https://github.com/satijalab/seurat/blob/master/R/visualization.R) (Licence: GPL-3).
#' This version has been modified to take single value length input instead of data.frame input.
#'
AutoPointSize_scCustom <- function(data, raster = NULL) {
  # for single value
  if (is.null(x = nrow(x = data)) && length(x = data) == 1 && is.numeric(x = data)) {
    return(ifelse(
      test = isTRUE(x = raster),
      yes = 1,
      no = min(1583 / data, 1)
    ))
  } else {
    # for data frame/object based values (from Seurat, see documentation)
    return(ifelse(
      test = isTRUE(x = raster),
      yes = 1,
      no = min(1583 / nrow(x = data), 1)
    ))
  }
}


#' Check for assays present in object
#
#' Checks Seurat object for the presence of assays against list of specified assays
#'
#' @param seurat_object Seurat object name.
#' @param assay_list vector of genes to check.
#' @param print_msg logical. Whether message should be printed if all features are found.  Default is TRUE.
#' @param omit_warn logical. Whether to print message about features that are not found in current object.
#'
#' @return List of found vs not found assays.
#'
#' @noRd
#'
Assay_Present <- function(
    seurat_object,
    assay_list,
    print_msg = TRUE,
    omit_warn = TRUE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)
  
  # get all features
  possible_assays <- Assays(object = seurat_object)
  
  # If any features not found
  if (any(!assay_list %in% possible_assays)) {
    bad_assays <- assay_list[!assay_list %in% possible_assays]
    found_assays <- assay_list[assay_list %in% possible_assays]
    if (length(x = found_assays) == 0) {
      stop("No requested assays found.")
    }
    
    # Return message of assays not found
    if (length(x = bad_assays) > 0 && omit_warn) {
      warning("The following assays were omitted as they were not found",
              ": ", glue_collapse_scCustom(input_string = bad_assays, and = TRUE))
    }
    
    # Combine into list and return
    assay_list <- list(
      found_assays = found_assays,
      bad_assays = bad_assays
    )
    return(assay_list)
  }
  
  # Print all found message if TRUE
  if (print_msg) {
    message("All assays present.")
  }
  
  # Return full input gene list.
  # Combine into list and return
  assay_list <- list(
    found_assays = assay_list,
    bad_assays = NULL
  )
  return(assay_list)
}



#' Add percentage difference to DE results
#'
#' Adds new column labeled "pct_diff" to the data.frame output of \code{\link[Seurat]{FindMarkers}},  \code{\link[Seurat]{FindAllMarkers}}, or other DE test data.frames.
#'
#' @param marker_df data.frame containing the results of \code{\link[Seurat]{FindMarkers}},  \code{\link[Seurat]{FindAllMarkers}}, or other DE test data.frame.
#' @param pct.1_name the name of data.frame column corresponding to percent expressed in group 1.
#' Default is Seurat default "pct.1".
#' @param pct.2_name the name of data.frame column corresponding to percent expressed in group 2.
#' Default is Seurat default "pct.2".
#' @param overwrite logical.  If the `marker_df` already contains column named "pct_diff" whether to
#'  overwrite or return error message.  Default is FALSE.
#'
#' @importFrom dplyr mutate
#' @importFrom magrittr "%>%"
#'
#' @return Returns input `marker_df` with additional "pct_diff" column.
#'
#'
#' @examples
#' \dontrun{
#' marker_df <- FindAllMarkers(object = obj_name)
#' marker_df <- Add_Pct_Diff(marker_dataframe = marker_df)
#' # or piped with function
#' marker_df <- FindAllMarkers(object = obj_name) %>%
#'   Add_Pct_Diff()
#' }
#'
Add_Pct_Diff <- function(
    marker_df,
    pct.1_name = "pct.1",
    pct.2_name = "pct.2",
    overwrite = FALSE
) {
  # Check if percent difference exists already
  if ("pct_diff" %in% colnames(marker_dataframe)) {
    if (!overwrite) {
      stop("'pct_diff' column already present in `marker_dataframe: ", '"', deparse(expr = substitute(expr = marker_dataframe)), '"', ". To overwrite previous results set `overwrite = TRUE`.")
    } else {
      message("'pct_diff' column already present in `marker_dataframe: ", '"', deparse(expr = substitute(expr = marker_dataframe)), '"',
              "  Overwriting column as overwrite = TRUE.")
    }
  }
  # Add percentage difference
  pct_diff_df <- marker_dataframe %>%
    mutate(pct_diff = .data[[pct.1_name]] - .data[[pct.2_name]])
  return(pct_diff_df)
}


#' Check if genes/features are present
#'
#' Check if genes are present in object and return vector of found genes.  Return warning messages for
#' genes not found.
#'
#' @param data Name of input data.  Currently only data of classes: Seurat, liger, data.frame,
#' dgCMatrix, dgTMatrix, tibble are accepted.  Gene_IDs must be present in rownames of the data.
#' @param gene_list vector of genes to check.
#' @param case_check logical. Whether or not to check if features are found if the case is changed from the
#' input list (Sentence case to Upper and vice versa).  Default is TRUE.
#' @param case_check_msg logical. Whether to print message to console if alternate case features are found
#' in addition to inclusion in returned list.  Default is TRUE.
#' @param print_msg logical. Whether message should be printed if all features are found.  Default is TRUE.
#' @param omit_warn logical. Whether to print message about features that are not found in current object.
#'  Default is TRUE.
#' @param return_none logical. Whether list of found vs. bad features should still be returned if no
#' features are found.  Default is FALSE.
#' @param seurat_assay Name of assay to pull feature names from if `data` is Seurat Object.
#' Defaults to  `DefaultAssay(OBJ)` if NULL.
#'
#' @importFrom purrr reduce
#' @importFrom stringr str_to_upper, str_to_sentence
#'
#' @return A list of length 3 containing 1) found features, 2) not found features, 3) features found if
#' case was modified.
#'
#' @export
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' features <- Gene_Present(data = obj_name, gene_list = DEG_list, print_msg = TRUE, case_check = TRUE)
#' found_features <- features[[1]]
#' }
#'
Gene_Present <- function(
    data,
    gene_list,
    case_check = TRUE,
    case_check_msg = TRUE,
    print_msg = TRUE,
    omit_warn = TRUE,
    return_none = FALSE,
    seurat_assay = NULL
) {
  # Check object type
  # Seurat
  accepted_types <- c("data.frame", "dgCMatrix", "dgTMatrix", "tibble")
  if ((class(x = data)[[1]] == "Seurat")) {
    # set assay (if null set to active assay)
    assay <- seurat_assay %||% DefaultAssay(object = data)
    
    possible_features <- rownames(x = GetAssayData(object = data, assay = assay))
  } else if ((class(x = data)[[1]] == "liger")) {
    # get complete gene list
    length_liger <- length(data@raw.data)
    
    list_genes <- lapply(1:length_liger, function(x){
      rownames(x = data@raw.data[[x]])
    })
    
    possible_features <- reduce(list_genes, function(x, y) {
      union(x = x, y = y)})
  } else if ((class(x = data) %in% accepted_types)) {
    possible_features <- rownames(x = data)
  } else {
    all_accepted <- c(accepted_types, "Seurat", "liger")
    stop("Input data is currently accepted only in the following formats: \n",
         glue_collapse_scCustom(input_string = all_accepted, and = FALSE))
  }
  
  # If any features not found
  if (any(!gene_list %in% possible_features)) {
    bad_features <- gene_list[!gene_list %in% possible_features]
    found_features <- gene_list[gene_list %in% possible_features]
    if (length(x = found_features) == 0) {
      if (return_none) {
        # Combine into list and return
        feature_list <- list(
          found_features = NULL,
          bad_features = bad_features,
          wrong_case_found_features = NULL
        )
        return(feature_list)
      } else {
        stop("No requested features found.")
      }
    }
    
    # Return message of features not found
    if (length(x = bad_features) > 0 && omit_warn) {
      warning("The following features were omitted as they were not found",
              ": ", glue_collapse_scCustom(input_string = bad_features, and = TRUE))
    }
    
    # Check if features found if case is changed.
    if (case_check) {
      upper_bad_features <- str_to_upper(string = bad_features)
      upper_found_features <- upper_bad_features[upper_bad_features %in% possible_features]
      
      sentence_bad_features <- str_to_sentence(string = bad_features)
      sentence_found_features <- sentence_bad_features[sentence_bad_features %in% possible_features]
      
      # Combine case check
      wrong_case_found_features <- c(upper_found_features, sentence_found_features)
      
      # Additional messages if found.
      if (length(x = wrong_case_found_features) > 0) {
        if (case_check_msg) {
          warning("NOTE: However, the following features were found: ",
                  glue_collapse_scCustom(input_string = wrong_case_found_features, and = TRUE), ".\n",
                  "        Please check intended case of features provided.")
        }
        # Combine into list and return
        feature_list <- list(
          found_features = found_features,
          bad_features = bad_features,
          wrong_case_found_features = wrong_case_found_features
        )
        return(feature_list)
      }
    }
    # Combine into list and return
    feature_list <- list(
      found_features = found_features,
      bad_features = bad_features,
      wrong_case_found_features = "NA (check not performed.  Set 'case_check = TRUE' to perform check."
    )
    return(feature_list)
  }
  
  # Print all found message if TRUE
  if (print_msg) {
    message("All features present.")
  }
  
  # Return full input gene list.
  # Combine into list and return
  feature_list <- list(
    found_features = gene_list,
    bad_features = NULL,
    wrong_case_found_features = NULL
  )
  return(feature_list)
}


#' Custom glue collapse
#
#' Customized glue_collapse that is based on the number of items in input list
#'
#' @param input_string input string to be collapsed.
#' @param and logical.  Whether to use "and" or "or" in the collapsed string
#'
#' @return collapsed string
#'
#' @importFrom glue glue_collapse
#'
#' @noRd
#'
glue_collapse_scCustom <- function(
    input_string,
    and = TRUE
) {
  # Check length of input string
  input_length <- length(x = input_string)
  
  # set last seperator
  if (and) {
    last_sep <- " and "
  } else {
    last_sep <- " or "
  }
  
  if (input_length < 3) {
    glue_collapse(x = input_string, sep = ", ", last = last_sep)
  } else {
    glue_collapse(x = input_string, sep = ", ", last = paste0(",", last_sep))
  }
}


#' Calculate the percentage of a vector above some threshold
#'
#' @param x Vector of values
#' @param threshold Threshold to use when calculating percentage
#'
#' @return Returns the percentage of `x` values above the given threshold
#'
#' @author Satija Lab & all co-Authors of Seurat Package
#' @references See Utilities.R in source code of Seurat https://github.com/satijalab/seurat/blob/master/R/utilities.R  (Licence: GPL-3).
#'
#' @note modified from Seurat version to return a percentage instead of proportion/decimal as part of `Percent_Expressing` function.  To be replaced following next Seurat version update.
#'
#' @keywords internal
#'
#' @noRd
#'
PercentAbove_Seurat <- function(x, threshold) {
  return((length(x = x[x > threshold]) / length(x = x))*100)
}
