#' Derive point estimates for c, pi, phi, and Z for a particular sample
#'
#' Posterior means are used as point estimates for \eqn{c}, \eqn{\pi},
#' \eqn{\phi}, and \eqn{Z}. As super-enriched peptides are tossed out before
#' MCMC sampling, super-enriched peptides return \code{NA} for the \eqn{\phi}
#' and \eqn{Z} point estimates. Indices corresponding to a particular peptide in
#' the MCMC sampler are mapped back to the original peptide names.
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object
#' @param file path to rds file
#' @param se.matrix logical matrix indicating which peptides were identified as
#' super-enriched peptides
#' @param burn.in number of iterations to be burned
#' @param post.thin thinning parameter
#'
#' @return list of point estimates for c, pi, phi and Z
summarizeRunOne <- function(object, file, se.matrix,
    burn.in = 0, post.thin = 1) {
    sample <- regmatches(file, regexec("/([^/]*)\\.rds", file))[[1]][2]

    rds <- readRDS(file)
    ## Support both the legacy format (plain mcmc object) and the wrapped format
    ## produced by brewOne (list with $mcmc and $sample_col).
    if (is.list(rds) && !is.null(rds$sample_col)) {
        mcmc_obj   <- rds$mcmc
        sample_col <- rds$sample_col
    } else {
        mcmc_obj   <- rds
        sample_col <- 1L
    }
    mcmc_matrix <- as.matrix(mcmc_obj)
    iter_ind <- seq(burn.in + 1, nrow(mcmc_matrix), by = post.thin)
    mcmc_matrix <- mcmc_matrix[iter_ind, , drop = FALSE]

    ## translate peptide indices
    pep_ind <- rep(NA, nrow(object))
    pep_ind[which(!se.matrix[, sample])] <- seq(sum(!se.matrix[, sample]))
    names(pep_ind) <- rownames(object)

    ## Extract parameter-specific columns using sample_col to handle the
    ## multi-column (input mode) case where c[j] and pi[j] are indexed by j.
    c_col  <- paste0("c[", sample_col, "]")
    pi_col <- paste0("pi[", sample_col, "]")
    samples_c  <- mcmc_matrix[, c_col,  drop = TRUE]
    samples_pi <- mcmc_matrix[, pi_col, drop = TRUE]

    z_pat   <- paste0("Z\\[.*,", sample_col, "\\]")
    phi_pat <- paste0("phi\\[.*,", sample_col, "\\]")
    samples_phi <- mcmc_matrix[, grepl(phi_pat, colnames(mcmc_matrix)), drop = FALSE]
    samples_Z   <- mcmc_matrix[, grepl(z_pat,   colnames(mcmc_matrix)), drop = FALSE]

    # summarize info
    point_c <- data.frame(
        parameter = "c",
        sample = sample,
        est_value = mean(samples_c)
    )
    point_pi <- data.frame(
        parameter = "pi",
        sample = sample,
        est_value = mean(samples_pi)
    )
    point_phi <- data.frame(
        parameter = "phi",
        sample = sample,
        peptide = rownames(object),
        est_value = unname(colMeans(samples_phi)[pep_ind]),
        est_enriched = unname((colSums(samples_phi * samples_Z) /
            colSums(samples_Z))
        [pep_ind])
    )
    point_Z <- data.frame(
        parameter = "Z",
        sample = sample,
        peptide = rownames(object),
        est_value = unname(colMeans(samples_Z)[pep_ind])
    )

    ## delta and tau_delta (input mode only)
    delta_cols <- grepl("^delta\\[", colnames(mcmc_matrix))
    point_delta <- if (any(delta_cols)) {
        samples_delta <- mcmc_matrix[, delta_cols, drop = FALSE]
        data.frame(
            parameter = "delta",
            sample = sample,
            peptide = rownames(object),
            est_value = unname(colMeans(samples_delta)[pep_ind])
        )
    } else {
        NULL
    }

    point_tau_delta <- if ("tau_delta" %in% colnames(mcmc_matrix)) {
        data.frame(
            parameter = "tau_delta",
            sample = sample,
            est_value = mean(mcmc_matrix[, "tau_delta"])
        )
    } else {
        NULL
    }

    list(
        point_c = point_c,
        point_pi = point_pi,
        point_phi = point_phi,
        point_Z = point_Z,
        point_delta = point_delta,
        point_tau_delta = point_tau_delta
    )
}

#' Summarize MCMC chain and return point estimates for BEER parameters
#'
#' Posterior means are used as point estimates for \eqn{c}, \eqn{\pi},
#' \eqn{\phi}, and \eqn{Z}. As super-enriched peptides are tossed out before
#' MCMC sampling, super-enriched peptides return \code{NA} for the \eqn{\phi}
#' and \eqn{Z} point estimates. Indices corresponding to a particular peptide in
#' the MCMC sampler are mapped back to the original peptide names.
#'
#' @param object a \code{\link[PhIPData]{PhIPData}} object
#' @param jags.files list of files containing MCMC sampling results
#' @param se.matrix logical matrix indicating which peptides were identified as
#' super-enriched peptides
#' @param burn.in number of iterations to be burned
#' @param post.thin thinning parameter
#' @param assay.names named vector of specifying where to store point estimates
#' @param BPPARAM \code{[BiocParallel::BiocParallelParam]} passed to
#' BiocParallel functions.
#'
#' @return PhIPData object with point estimates stored in the assays specified
#' by `assay.names`.
#'
#' @importFrom progressr handlers progressor
#' @importFrom BiocParallel bplapply
#' @import PhIPData SummarizedExperiment
summarizeRun <- function(object, jags.files, se.matrix,
    burn.in = 0, post.thin = 1,
    assay.names = c(
        phi = NULL, phi_Z = "logfc", Z = "prob",
        c = "sampleInfo", pi = "sampleInfo",
        delta = NA, tau_delta = "sampleInfo"
    ),
    BPPARAM = BiocParallel::bpparam()) {

    ## Check that all files are present
    if (!all(file.exists(jags.files))) {
        stop(paste0(
            "Cannot find the following files: ",
            paste0(jags.files[!file.exists(jags.files)],
                collapse = ", "
            )
        ))
    }

    samples <- vapply(
        regmatches(jags.files, regexec("/([^/]*)\\.rds", jags.files)),
        function(x) x[[2]],
        character(1)
    )
    names(jags.files) <- samples

    ## Pre-allocate containers
    point_c <- if (assay.names["c"] %in% colnames(sampleInfo(object))) {
        sampleInfo(object)[[assay.names["c"]]]
    } else {
        rep(NA, ncol(object))
    }
    point_pi <- if (assay.names["pi"] %in% colnames(sampleInfo(object))) {
        sampleInfo(object)[[assay.names["pi"]]]
    } else {
        rep(NA, ncol(object))
    }
    point_phi <- if (!is.null(assay.names["phi"]) &
        assay.names["phi"] %in% assayNames(object)) {
        assay(object, assay.names["phi"])
    } else {
        matrix(NA, nrow = nrow(object), ncol = ncol(object))
    }
    point_phi <- if (!is.null(assay.names["phi"]) &
        assay.names["phi"] %in% assayNames(object)) {
        assay(object, assay.names["phi"])
    } else {
        matrix(NA, nrow = nrow(object), ncol = ncol(object))
    }
    point_phi_Z <- if (!is.null(assay.names["phi_Z"]) &
        assay.names["phi_Z"] %in% assayNames(object)) {
        assay(object, assay.names["phi_Z"])
    } else {
        matrix(NA, nrow = nrow(object), ncol = ncol(object))
    }
    point_Z <- if (!is.null(assay.names["Z"]) &
        assay.names["Z"] %in% assayNames(object)) {
        assay(object, assay.names["Z"])
    } else {
        matrix(NA, nrow = nrow(object), ncol = ncol(object))
    }

    point_delta <- if (!is.na(assay.names["delta"]) &
        assay.names["delta"] %in% assayNames(object)) {
        assay(object, assay.names["delta"])
    } else {
        matrix(NA, nrow = nrow(object), ncol = ncol(object))
    }

    point_tau_delta <- if (!is.na(assay.names["tau_delta"]) &
        assay.names["tau_delta"] %in% colnames(sampleInfo(object))) {
        sampleInfo(object)[[assay.names["tau_delta"]]]
    } else {
        rep(NA, ncol(object))
    }

    names(point_c) <- names(point_pi) <- names(point_tau_delta) <-
        colnames(point_phi) <- colnames(point_phi_Z) <-
        colnames(point_Z) <- colnames(point_delta) <-
        colnames(object)
    rownames(point_phi) <- rownames(point_phi_Z) <- rownames(point_Z) <-
        rownames(point_delta) <- rownames(object)

    ## Summarize each file, use bplapply here to enable reading/import
    ## to be parallelized
    progressr::handlers("txtprogressbar")
    p <- progressr::progressor(along = jags.files)
    files_out <- bplapply(jags.files, function(file) {
        file_counter <- paste0(
            which(file == jags.files), " of ",
            length(jags.files)
        )
        p(file_counter, class = "sticky", amount = 1)
        summarizeRunOne(object, file, se.matrix, burn.in, post.thin)
    }, BPPARAM = BPPARAM)

    for (out in files_out) {
        sample <- out$point_c$sample

        point_c[sample] <- out$point_c$est_value
        point_pi[sample] <- out$point_pi$est_value
        point_phi[, sample] <- out$point_phi$est_value
        point_phi_Z[, sample] <- out$point_phi$est_enriched
        point_Z[, sample] <- out$point_Z$est_value
        if (!is.null(out$point_delta)) {
            point_delta[, sample] <- out$point_delta$est_value
        }
        if (!is.null(out$point_tau_delta)) {
            point_tau_delta[sample] <- out$point_tau_delta$est_value
        }
    }

    ## Assign sample-level scalars to sampleInfo
    if (!is.na(assay.names["c"])) object$c <- point_c
    if (!is.na(assay.names["pi"])) object$pi <- point_pi
    if (!is.na(assay.names["tau_delta"])) object$tau_delta <- point_tau_delta

    ## Assign phi, phi_Z, Z, and (when present) delta to assays
    assay <- c("phi", "phi_Z", "Z")[!is.na(assay.names[c("phi", "phi_Z", "Z")])]
    assays(object)[assay.names[assay]] <- list(
        phi = point_phi, phi_Z = point_phi_Z,
        Z = point_Z
    )[assay]
    if (!is.na(assay.names["delta"])) {
        assays(object)[[assay.names[["delta"]]]] <- point_delta
    }
    object
}
