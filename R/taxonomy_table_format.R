#' taxonomy_table_format
#'
#' @description
#' strollur accepts two formats for classification data: `full` or `tidy`.
#'
#' @details
#' The `full` taxonomy table is a `data.frame` with two columns: `sequence_name`
#' or `bin_name`, and `taxonomy`. The `taxonomy` column contains the full
#' classification string.
#'
#' **The `full` taxonomy table format:**
#' ```
#' sequence_name
#' M00967_43_000000000-A3JHG_1_1102_26308_14496
#' M00967_43_000000000-A3JHG_1_1103_22654_26891
#'
#' taxonomy
#' Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);
#' Bacteria(100);Firmicutes(99);Clostridia(99);Clostridiales(99);
#' ```
#'
#' The `tidy` taxonomy table is a `data.frame` with 4 columns: `sequence_name`
#' or `bin_name`, `taxonomy`, `level`, and `confidence`. For example:
#'
#' **The `tidy` taxonomy table format:**
#' ```
#' sequence_name                                    taxonomy
#' M00967_43_000000000-A3JHG_1_1102_26308_14496     Bacteria
#' M00967_43_000000000-A3JHG_1_1102_26308_14496     "Bacteroidetes"
#' M00967_43_000000000-A3JHG_1_1102_26308_14496     "Bacteroidia"
#' M00967_43_000000000-A3JHG_1_1102_26308_14496     "Bacteroidales"
#' M00967_43_000000000-A3JHG_1_1103_22654_26891     Bacteria
#' M00967_43_000000000-A3JHG_1_1103_22654_26891     Firmicutes
#' M00967_43_000000000-A3JHG_1_1103_22654_26891     Clostridia
#' M00967_43_000000000-A3JHG_1_1103_22654_26891     Clostridiales
#'
#' level    confidence
#' 1        100
#' 2        100
#' 3        100
#' 4        100
#' 1        100
#' 2        99
#' 3        99
#' 4        99
#' ```
#'
#' @name taxonomy_table_format
#' @keywords internal
NULL
