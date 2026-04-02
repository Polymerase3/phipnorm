#' Flag peptides with anomalous stickiness
#'
#' After running \code{\link{brew}} with an input library (so that stickiness
#' coefficients \eqn{\delta_j} have been estimated), \code{flagSticky} identifies
#' peptides whose posterior \eqn{\delta} is in the extreme tails of the
#' log-scale distribution. These peptides are structurally unusual in the IP
#' process — over-represented or under-represented in the mockIP relative to
#' their input abundance — and enrichment calls for them should be interpreted
#' with caution.
#'
#' Stickiness is summarised across samples by taking the mean of the per-sample
#' posterior means on the log scale. This is equivalent to the geometric mean
#' of the per-sample posterior mean \eqn{\delta} values.
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object with a stickiness
#'   assay, as produced by \code{\link{brew}} when an input library is present.
#' @param assay.name character; name of the assay containing posterior mean
#'   \eqn{\delta} values. Default \code{"stickiness"}.
#' @param lower.quantile,upper.quantile numeric in (0, 1); quantiles of the
#'   log-\eqn{\delta} distribution used as flagging thresholds. Peptides outside
#'   \code{[lower.quantile, upper.quantile]} are flagged. Defaults are 0.01 and
#'   0.99, respectively.
#'
#' @return a logical vector of length \code{nrow(object)}, named by peptide.
#'   \code{TRUE} indicates a peptide flagged as anomalously sticky or slippery.
#'
#' @seealso \code{\link{brew}} for the fitting step that produces \eqn{\delta}
#'   estimates.
#'
#' @import PhIPData SummarizedExperiment
#' @importFrom stats quantile
#' @export
flagSticky <- function(object,
    assay.name = "stickiness",
    lower.quantile = 0.01,
    upper.quantile = 0.99) {

    if (!assay.name %in% assayNames(object)) {
        stop(
            "Assay '", assay.name, "' not found. ",
            "Run brew() with an input library group to estimate stickiness."
        )
    }
    if (lower.quantile <= 0 || lower.quantile >= upper.quantile ||
            upper.quantile >= 1) {
        stop("lower.quantile and upper.quantile must satisfy 0 < lower < upper < 1.")
    }

    delta_mat <- assay(object, assay.name)

    ## Summarise across samples on log scale (geometric mean of per-sample estimates)
    log_delta <- rowMeans(log(delta_mat), na.rm = TRUE)

    bounds <- quantile(log_delta, c(lower.quantile, upper.quantile), na.rm = TRUE)

    flag <- log_delta < bounds[[1]] | log_delta > bounds[[2]]
    flag[is.na(flag)] <- FALSE
    names(flag) <- rownames(object)

    flag
}
