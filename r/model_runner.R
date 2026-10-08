library(lme4)
library(lmerTest)
library(emmeans)

# Match uploaded sample tables; tables without IDs retain the positional API.
align_sample_tables <- function(metadata, data_matrix) {
    if (!all(vapply(list(metadata, data_matrix), function(table) {
        "SampleID" %in% names(table)
    }, logical(1)))) {
        if (nrow(metadata) != nrow(data_matrix)) {
            stop("Tables without SampleID must have the same number of rows.")
        }
        return(list(metadata = metadata, data_matrix = data_matrix))
    }
    for (table in list(metadata, data_matrix)) {
        ids <- as.character(table$SampleID)
        if (anyNA(ids) || any(!nzchar(ids))) stop("SampleID values must not be missing or blank.")
        if (anyDuplicated(ids)) stop("SampleID values must be unique in each table.")
    }
    metadata_ids <- as.character(metadata$SampleID)
    feature_ids <- as.character(data_matrix$SampleID)
    shared <- metadata_ids %in% feature_ids
    if (!any(shared)) stop("Metadata and feature matrix have no shared SampleID values.")
    if (!all(shared) || any(!feature_ids %in% metadata_ids)) {
        warning("Unmatched samples were removed before modeling.")
    }
    metadata <- metadata[shared, , drop = FALSE]
    data_matrix <- data_matrix[match(as.character(metadata$SampleID), feature_ids), , drop = FALSE]
    list(metadata = metadata, data_matrix = data_matrix)
}

converged_singular_fraction <- function(results) {
    valid <- results$converged %in% TRUE & !is.na(results$singular)
    if (any(valid)) mean(results$singular[valid]) else NA_real_
}

bind_feature_results <- function(results) {
    columns <- unique(unlist(lapply(results, names), use.names = FALSE))
    results <- lapply(results, function(result) {
        for (column in setdiff(columns, names(result))) result[[column]] <- NA
        result[, columns, drop = FALSE]
    })
    do.call(rbind, results)
}

validate_fixed_effects_rank <- function(formula, data) {
    design <- model.matrix(formula, data = data)
    if (qr(design)$rank < ncol(design)) {
        stop("Fixed-effects model is rank deficient: effects are not identifiable after covariate adjustment on the usable observations.")
    }
    invisible(TRUE)
}

# Build formulas with internal names, and use one complete-case population
# for all full/reduced fits belonging to a feature.
prepare_feature_model <- function(data, outcome, predictor, predictor_type,
                                  covariates, time = NULL, subject_id = NULL) {
    required <- unique(c(outcome, predictor, covariates, time, subject_id))
    if (!all(required %in% names(data))) stop("Model variables are missing from metadata.")
    categorical <- predictor_type %in% c("categorical", "ordered_categorical")
    expected_categories <- if (categorical) unique(as.character(na.omit(data[[predictor]]))) else NULL
    for (variable in unique(c(outcome, time,
                              if (predictor_type == "numeric") predictor))) {
        values <- data[[variable]]
        if (is.factor(values)) values <- as.character(values)
        converted <- suppressWarnings(as.numeric(values))
        if (any(!is.na(values) & is.na(converted))) {
            stop(paste0("Numeric variable '", variable, "' contains nonnumeric non-missing values."))
        }
        if (any(is.infinite(converted))) {
            stop(paste0("Numeric variable '", variable, "' contains infinite values."))
        }
        data[[variable]] <- converted
    }
    for (variable in covariates) {
        if (is.numeric(data[[variable]]) && any(is.infinite(data[[variable]]))) {
            stop(paste0("Numeric covariate '", variable, "' contains infinite values."))
        }
    }
    data <- data[complete.cases(data[, required, drop = FALSE]), , drop = FALSE]
    if (categorical) {
        missing_categories <- setdiff(expected_categories, as.character(data[[predictor]]))
        if (length(missing_categories)) {
            stop(paste("No usable observations for predictor categories:",
                       paste(missing_categories, collapse = ", ")))
        }
    }
    original <- names(data)
    # Reserve a prefix absent from original names and category labels, so
    # coefficient suffixes cannot be mistaken for internal variable names.
    labels <- c(original, unlist(lapply(data, function(values) {
        if (is.factor(values)) return(c(levels(values), as.character(values)))
        if (is.character(values)) return(values)
        character(0)
    }), use.names = FALSE))
    labels <- labels[!is.na(labels)]
    prefix <- "modelvar"
    while (any(grepl(prefix, labels, fixed = TRUE))) prefix <- paste0(prefix, "X")
    internal <- sprintf("%s%06dX", prefix, seq_along(original))
    names(data) <- internal
    mapping <- setNames(internal, original)
    list(data = data, mapping = mapping)
}

restore_model_labels <- function(result, mapping) {
    internal <- unname(mapping)
    originals <- names(mapping)
    # Matches are located in the untouched term, then replaced together.
    # Restored names are never scanned again, even if they contain tokens.
    restore_term <- function(value) {
        if (is.na(value)) return(NA_character_)
        matches <- gregexpr(paste(internal, collapse = "|"), value, perl = TRUE)
        tokens <- regmatches(value, matches)[[1]]
        if (length(tokens)) {
            regmatches(value, matches) <- list(originals[match(tokens, internal)])
        }
        value
    }
    for (table in c("overall", "pairwise")) {
        if (is.null(result[[table]])) next
        for (column in intersect(c("feature", "predictor", "term"),
                                 names(result[[table]]))) {
            values <- as.character(result[[table]][[column]])
            if (column == "term") {
                values <- vapply(values, restore_term, character(1), USE.NAMES = FALSE)
            } else {
                matched <- match(values, internal)
                found <- !is.na(matched)
                values[found] <- originals[matched[found]]
            }
            result[[table]][[column]] <- values
        }
    }
    result
}

run_cross_sectional_feature <- function(data, outcome, predictor, predictor_type,
                                       covariates = character(0), ordered_levels = NULL) {
    covariates <- as.character(unlist(covariates, use.names = FALSE))
    prepared <- prepare_feature_model(data, outcome, predictor, predictor_type, covariates)
    mapping <- prepared$mapping
    result <- run_cross_sectional_feature_impl(
        data = prepared$data, outcome = unname(mapping[outcome]),
        predictor = unname(mapping[predictor]), predictor_type = predictor_type,
        covariates = unname(mapping[covariates]), ordered_levels = ordered_levels)
    if (identical(result$overall$status, "estimate_only")) {
        result$overall$std_error <- NA_real_
        result$overall$p_value <- NA_real_
    }
    restore_model_labels(result, mapping)
}

run_longitudinal_feature <- function(data, outcome, predictor, predictor_type,
                                    time, subject_id,
                                    analysis_goal = c("overall_predictor_association",
                                                      "time_effect", "trajectory_difference"),
                                    covariates = character(0), ordered_levels = NULL,
                                    random_slope = TRUE) {
    analysis_goal <- match.arg(analysis_goal)
    covariates <- as.character(unlist(covariates, use.names = FALSE))
    prepared <- prepare_feature_model(data, outcome, predictor, predictor_type,
                                      covariates, time, subject_id)
    mapping <- prepared$mapping
    n_subjects <- length(unique(prepared$data[[mapping[[subject_id]]]]))
    if (random_slope && nrow(prepared$data) <= 2L * n_subjects) {
        stop(paste0("Random intercept + time slope is unsupported: ",
                    nrow(prepared$data), " usable observations for ", n_subjects,
                    " subjects (", 2L * n_subjects, " random effects). ",
                    "Select random intercept only and inspect again; the selected model was not changed."))
    }
    result <- run_longitudinal_feature_impl(
        data = prepared$data, outcome = unname(mapping[outcome]),
        predictor = unname(mapping[predictor]), predictor_type = predictor_type,
        time = unname(mapping[time]), subject_id = unname(mapping[subject_id]),
        analysis_goal = analysis_goal, covariates = unname(mapping[covariates]),
        ordered_levels = ordered_levels, random_slope = random_slope)
    if (!isTRUE(result$overall$converged)) {
        result$overall$error <- failure_message(result$diagnostics)
    }
    restore_model_labels(result, mapping)
}

# ============================================================
# Helper: prepare primary predictor
# ============================================================

prepare_primary_predictor <- function(
    data,
    predictor,
    predictor_type,
    ordered_levels = NULL
) {

    # --------------------------------------------------------
    # Numeric predictor
    # --------------------------------------------------------

    if (predictor_type == "numeric") {

        data[[predictor]] <- as.numeric(
            data[[predictor]]
        )

        return(
            list(
                data = data,
                model_predictor = predictor,
                strategy = "coefficient"
            )
        )
    }


    # --------------------------------------------------------
    # Ordered categorical predictor
    # --------------------------------------------------------

    if (predictor_type == "ordered_categorical") {

        if (
            is.null(ordered_levels) ||
            length(ordered_levels) < 2
        ) {

            stop(
                paste(
                    "Ordered levels must be supplied for",
                    predictor
                )
            )
        }


        # --------------------------------------------
        # Check observed levels
        # --------------------------------------------

        observed <- unique(
            as.character(
                data[[predictor]][
                    !is.na(data[[predictor]])
                ]
            )
        )

        missing_from_order <- setdiff(
            observed,
            ordered_levels
        )

        if (length(missing_from_order) > 0) {

            stop(
                paste(
                    "Values found in",
                    predictor,
                    "but missing from ordered_levels:",
                    paste(
                        missing_from_order,
                        collapse = ", "
                    )
                )
            )
        }


        # --------------------------------------------
        # Preserve specified biological order
        # --------------------------------------------

        data[[predictor]] <- factor(
            data[[predictor]],
            levels = ordered_levels,
            ordered = TRUE
        )


        # --------------------------------------------
        # Convert order to numeric trend score
        #
        # weak        -> 1
        # strong      -> 2
        # very strong -> 3
        # --------------------------------------------

        score_name <- paste0(
            predictor,
            "_trend"
        )

        data[[score_name]] <- as.numeric(
            data[[predictor]]
        )


        return(
            list(
                data = data,
                model_predictor = score_name,
                strategy = "ordered_trend"
            )
        )
    }


     # --------------------------------------------------------
    # Unordered categorical predictor
    # --------------------------------------------------------

    if (predictor_type == "categorical") {

        data[[predictor]] <- factor(
            data[[predictor]]
        )

        n_levels <- nlevels(
            data[[predictor]]
        )

        if (n_levels < 2) {
            stop(
                paste(
                    predictor,
                    "must contain at least 2 levels."
                )
            )
        }

        # Binary categorical
        if (n_levels == 2) {

            return(
                list(
                    data = data,
                    model_predictor = predictor,
                    strategy = "coefficient"
                )
            )
        }

        # Multi-level unordered categorical
        return(
            list(
                data = data,
                model_predictor = predictor,
                strategy = "overall_and_pairwise"
            )
        )
    }


    stop(
        paste(
            "Unknown predictor type:",
            predictor_type
        )
    )
}


# ============================================================
# Helper: build RHS
# ============================================================

build_rhs <- function(
    terms,
    covariates = character(0)
) {

    all_terms <- c(
        terms,
        covariates
    )

    all_terms <- all_terms[
        !is.na(all_terms)
    ]

    all_terms <- all_terms[
        all_terms != ""
    ]

    paste(
        all_terms,
        collapse = " + "
    )
}


# ============================================================
# Helper: extract coefficient from lm / lmer
# ============================================================

# ============================================================
# Extract primary coefficient
#
# strategy:
#   "coefficient"
#   "ordered_trend"
#   "overall_and_pairwise"
# ============================================================

# Select columns by their formula term, independently of factor-level text.
coefficient_for_variables <- function(fit, variables) {
    design <- model.matrix(fit)
    factors <- attr(terms(fit), "factors")
    term_ids <- which(vapply(seq_len(ncol(factors)), function(i) {
        setequal(rownames(factors)[factors[, i] != 0], variables)
    }, logical(1)))
    matched <- colnames(design)[attr(design, "assign") %in% term_ids]
    matched <- intersect(matched, rownames(summary(fit)$coefficients))
    if (length(matched) != 1L) {
        stop(paste("Expected one coefficient for", paste(variables, collapse = ":"),
                   "but found", length(matched)))
    }
    matched[[1]]
}

extract_primary_coefficient <- function(
    fit,
    predictor,
    strategy
) {

    # --------------------------------------------------------
    # Multi-category categorical predictor
    # does NOT have one primary coefficient
    # --------------------------------------------------------

    if (strategy == "overall_and_pairwise") {

        return(NULL)
    }


    # --------------------------------------------------------
    # Get coefficient table
    # --------------------------------------------------------

    sm <- summary(fit)$coefficients

    if (!strategy %in% c("coefficient", "ordered_trend")) {
        stop(paste("Unknown extraction strategy:", strategy))
    }
    term <- coefficient_for_variables(fit, predictor)

    # --------------------------------------------------------
    # Extract result
    # --------------------------------------------------------

    return(
        data.frame(

            term = term,

            estimate =
                sm[
                    term,
                    "Estimate"
                ],

            std_error =
                sm[
                    term,
                    "Std. Error"
                ],

            p_value =
                sm[
                    term,
                    "Pr(>|t|)"
                ],

            stringsAsFactors = FALSE
        )
    )
}

# ============================================================
# Cross-sectional analysis
# ============================================================

estimate_only_warning <- "Estimate only: insufficient residual degrees of freedom to calculate uncertainty or significance."

run_cross_sectional_feature_impl <- function(
    data,
    outcome,
    predictor,
    predictor_type,
    covariates = character(0),
    ordered_levels = NULL
) {

    # --------------------------------------------------------
    # Outcome must be numeric
    # --------------------------------------------------------

    data[[outcome]] <- as.numeric(
        data[[outcome]]
    )


    # --------------------------------------------------------
    # Prepare predictor and determine strategy
    # --------------------------------------------------------

    prepared <- prepare_primary_predictor(
        data = data,
        predictor = predictor,
        predictor_type = predictor_type,
        ordered_levels = ordered_levels
    )


    dat <- prepared$data

    model_predictor <- prepared$model_predictor

    strategy <- prepared$strategy


    # ========================================================
    # CASE 1
    #
    # Numeric predictor
    # Binary categorical predictor
    # Ordered categorical trend
    # ========================================================

    if (
        strategy %in%
        c(
            "coefficient",
            "ordered_trend"
        )
    ) {

        # ----------------------------------------------------
        # Build RHS
        # ----------------------------------------------------

        rhs <- build_rhs(
            terms = model_predictor,
            covariates = covariates
        )


        # ----------------------------------------------------
        # Build formula
        # ----------------------------------------------------

        formula <- as.formula(
            paste(
                outcome,
                "~",
                rhs
            )
        )


        # ----------------------------------------------------
        # Fit model
        # ----------------------------------------------------

        validate_fixed_effects_rank(formula, dat)
        fit <- lm(
            formula,
            data = dat
        )


        # ----------------------------------------------------
        # Extract primary predictor / trend coefficient
        # ----------------------------------------------------

        coef_result <- extract_primary_coefficient(
            fit = fit,
            predictor = model_predictor,
            strategy = strategy
        )


        # ----------------------------------------------------
        # Return main result
        # ----------------------------------------------------

        overall_result <- data.frame(

            feature = outcome,

            predictor = predictor,

            predictor_type = predictor_type,

            test_type = strategy,

            term = coef_result$term,

            estimate = coef_result$estimate,

            std_error = coef_result$std_error,

            p_value = coef_result$p_value,

            n_obs = nobs(fit),
            status = if (df.residual(fit) <= 0) "estimate_only" else "success",
            warning = if (df.residual(fit) <= 0) estimate_only_warning else NA_character_,
            stringsAsFactors = FALSE
        )


        return(
            list(
                overall = overall_result,
                pairwise = NULL
            )
        )
    }


    # ========================================================
    # CASE 2
    #
    # Unordered categorical predictor with >2 levels
    #
    # Return:
    #   overall association
    #   pairwise comparisons
    # ========================================================

    if (
        strategy ==
        "overall_and_pairwise"
    ) {

        # ----------------------------------------------------
        # Full model RHS
        #
        # predictor + covariates
        # ----------------------------------------------------

        full_rhs <- build_rhs(
            terms = model_predictor,
            covariates = covariates
        )


        full_formula <- as.formula(
            paste(
                outcome,
                "~",
                full_rhs
            )
        )


        # ----------------------------------------------------
        # Reduced model RHS
        #
        # covariates only
        # ----------------------------------------------------

        if (
            length(covariates) > 0
        ) {

            reduced_rhs <- build_rhs(
                terms = character(0),
                covariates = covariates
            )


            reduced_formula <- as.formula(
                paste(
                    outcome,
                    "~",
                    reduced_rhs
                )
            )

        } else {

            reduced_formula <- as.formula(
                paste(
                    outcome,
                    "~ 1"
                )
            )
        }


        # ----------------------------------------------------
        # Fit reduced and full models
        # ----------------------------------------------------

        validate_fixed_effects_rank(full_formula, dat)
        validate_fixed_effects_rank(reduced_formula, dat)
        fit0 <- lm(
            reduced_formula,
            data = dat
        )


        fit1 <- lm(
            full_formula,
            data = dat
        )


        # ----------------------------------------------------
        # Overall association
        # ----------------------------------------------------

        test <- anova(
            fit0,
            fit1
        )


        overall_p <- test$`Pr(>F)`[2]
        estimate_only <- df.residual(fit1) <= 0


        overall_result <- data.frame(

            feature = outcome,

            predictor = predictor,

            predictor_type = predictor_type,

            test_type = "overall_association",

            term = predictor,

            estimate = NA_real_,

            std_error = NA_real_,

            p_value = overall_p,

            n_obs = nobs(fit1),
            status = if (estimate_only) "estimate_only" else "success",
            warning = if (estimate_only) estimate_only_warning else NA_character_,
            stringsAsFactors = FALSE
        )


        # ----------------------------------------------------
        # Pairwise comparisons
        #
        # Always run
        # ----------------------------------------------------

        emm <- emmeans(
            fit1,
            specs = predictor
        )


        pw <- pairs(
            emm,
            adjust = "none"
        )


        pairwise_result <- as.data.frame(
            pw
        )


        pairwise_result$feature <- outcome

        pairwise_result$predictor <- predictor

        pairwise_result$predictor_type <- predictor_type

        pairwise_result$test_type <- "pairwise"


        pairwise_result <- pairwise_result[
            ,
            c(
                "feature",
                "predictor",
                "predictor_type",
                "test_type",
                "contrast",
                "estimate",
                "SE",
                "df",
                "t.ratio",
                "p.value"
            )
        ]


        if (estimate_only) {
            overall_result$p_value <- NA_real_
            pairwise_result$SE <- NA_real_
            pairwise_result$t.ratio <- NA_real_
            pairwise_result$p.value <- NA_real_
            pairwise_result$status <- "estimate_only"
            pairwise_result$warning <- estimate_only_warning
        }

        return(
            list(
                overall = overall_result,
                pairwise = pairwise_result
            )
        )
    }


    # ========================================================
    # Unknown strategy
    # ========================================================

    stop(
        paste(
            "Unknown predictor strategy:",
            strategy
        )
    )
}


run_cross_sectional_omics <- function(
    metadata,
    data_matrix,
    predictor,
    predictor_type,
    covariates = character(0),
    ordered_levels = NULL
) {

    aligned <- align_sample_tables(metadata, data_matrix)
    metadata <- aligned$metadata
    data_matrix <- aligned$data_matrix

    overall_results <- list()
    pairwise_results <- list()


    # ========================================================
    # Each column in data_matrix is one numeric outcome
    # ========================================================

    #features <- colnames(
    #    data_matrix
    #)
    features <- setdiff(
    colnames(data_matrix),
    "SampleID")
    
    # ========================================================
    # Loop through all features
    # ========================================================

    for (feature in features) {

        message(
            "Running: ",
            feature
        )


        # ----------------------------------------------------
        # Add current outcome to metadata
        #
        # Sample tables were aligned before entering the feature loop.
        # ----------------------------------------------------

        dat <- metadata

        outcome_column <- make.unique(c(names(dat), "analysis_outcome"))[ncol(dat) + 1L]

        dat[[outcome_column]] <- (
            data_matrix[[feature]]
        )


        # ----------------------------------------------------
        # Run one cross-sectional model
        #
        # If one feature fails, continue with the rest.
        # ----------------------------------------------------

        result <- tryCatch(

            {

                run_cross_sectional_feature(
                    data = dat,
                    outcome = outcome_column,
                    predictor = predictor,
                    predictor_type = predictor_type,
                    covariates = covariates,
                    ordered_levels = ordered_levels
                )
            },


            error = function(e) {

                list(

                    overall = data.frame(

                        feature =
                            feature,

                        predictor =
                            predictor,

                        predictor_type =
                            predictor_type,

                        test_type =
                            NA_character_,

                        term =
                            NA_character_,

                        estimate =
                            NA_real_,

                        std_error =
                            NA_real_,

                        p_value =
                            NA_real_,

                        n_obs =
                            NA_real_,

                        status = "failed",
                        error =
                            conditionMessage(e),

                        stringsAsFactors =
                            FALSE
                    ),

                    pairwise =
                        NULL
                )
            }
        )


        # ----------------------------------------------------
        # Store main result
        # ----------------------------------------------------

        result$overall$feature <- feature
        if (!is.null(result$pairwise)) result$pairwise$feature <- feature

        overall_results[[feature]] <- (
            result$overall
        )


        # ----------------------------------------------------
        # Store pairwise results if they exist
        #
        # Only multi-category unordered predictors generate
        # pairwise comparisons.
        # ----------------------------------------------------

        if (
            !is.null(
                result$pairwise
            )
        ) {

            pairwise_results[[feature]] <- (
                result$pairwise
            )
        }
    }


    # ========================================================
    # Combine main / overall results
    # ========================================================

    overall_results <- bind_feature_results(overall_results)

    rownames(
        overall_results
    ) <- NULL


    # ========================================================
    # BH FDR across outcomes
    #
    # Depending on predictor type, p_value represents:
    #
    # numeric:
    #   predictor coefficient
    #
    # binary categorical:
    #   group coefficient
    #
    # ordered categorical:
    #   ordered trend
    #
    # multi-category:
    #   overall F test
    # ========================================================

    overall_results$FDR <- p.adjust(
        overall_results$p_value,
        method = "BH"
    )


    # ========================================================
    # Combine pairwise results
    # ========================================================

    if (
        length(
            pairwise_results
        ) > 0
    ) {

        pairwise_results <- do.call(
            rbind,
            pairwise_results
        )

        rownames(
            pairwise_results
        ) <- NULL


        # ====================================================
        # Pairwise FDR
        #
        # At the feature level, use:
        #
        # pairs(
        #     emm,
        #     adjust = "none"
        # )
        #
        # so p.value contains raw pairwise p-values.
        #
        # BH correction is then applied across omics features
        # separately for each contrast.
        #
        # Example:
        #
        # A - B:
        # Gene1, Gene2, ..., Gene5000
        # -> BH correction
        #
        # A - C:
        # Gene1, Gene2, ..., Gene5000
        # -> separate BH correction
        # ====================================================

        pairwise_results$FDR <- ave(

            pairwise_results$p.value,

            pairwise_results$contrast,

            FUN = function(p) {

                p.adjust(
                    p,
                    method = "BH"
                )
            }
        )


    } else {

        pairwise_results <- NULL
    }


    # ========================================================
    # Return
    # ========================================================

    return(
        list(

            overall =
                overall_results,

            pairwise =
                pairwise_results
        )
    )
}



######longitudinal functions###########
# ============================================================
# Helper: identify predictor × time coefficient

# Find a single predictor × time interaction coefficient
#
# Works when predictor produces one model coefficient:
# - numeric predictor
# - binary categorical predictor
# - ordered trend predictor
#
# time is numeric.
# ============================================================

extract_interaction_coefficient <- function(
    fit,
    model_predictor,
    time
) {

    sm <- coef(
        summary(fit)
    )

    term <- coefficient_for_variables(fit, c(model_predictor, time))

    data.frame(
        term = term,

        estimate =
            sm[
                term,
                "Estimate"
            ],

        std_error =
            sm[
                term,
                "Std. Error"
            ],

        p_value =
            sm[
                term,
                "Pr(>|t|)"
            ],

        stringsAsFactors = FALSE
    )
}


# ============================================================
# Helper: extract numeric time coefficient
# ============================================================

extract_time_coefficient <- function(
    fit,
    time
) {

    sm <- summary(fit)$coefficients

    if (!time %in% rownames(sm)) {

        stop(
            paste(
                "Time coefficient",
                time,
                "was not found."
            )
        )
    }

    data.frame(
        term = time,

        estimate =
            sm[
                time,
                "Estimate"
            ],

        std_error =
            sm[
                time,
                "Std. Error"
            ],

        p_value =
            sm[
                time,
                "Pr(>|t|)"
            ],

        stringsAsFactors = FALSE
    )
}


# ============================================================
# 6. Build random-effect term
# ============================================================

build_random_term <- function(
    subject_id,
    time,
    random_slope
) {

    if (random_slope) {

        return(
            paste0(
                "(1 + ",
                time,
                " | ",
                subject_id,
                ")"
            )
        )
    }


    paste0(
        "(1 | ",
        subject_id,
        ")"
    )
}

# ============================================================
# 4. Check convergence
# ============================================================

check_lmer_convergence <- function(
    fit
) {

    messages <- (
        fit@optinfo$conv$lme4$messages
    )
    # A boundary/singular fit can have a successful optimizer. Singularity
    # is tracked separately and must remain eligible for the fallback pass.
    messages <- messages[!grepl("boundary (singular) fit", messages, fixed = TRUE)]
    optimizer_codes <- unlist(fit@optinfo$conv$opt)
    optimizer_ok <- all(optimizer_codes == 0)
    converged <- optimizer_ok && length(messages) == 0
    if (!optimizer_ok) {
        messages <- c(messages, paste("Optimizer convergence code:",
                                     paste(optimizer_codes, collapse = ", ")))
    }


    list(
        converged = converged,

        convergence_message =
            if (converged) {

                NA_character_

            } else {

                paste(
                    messages,
                    collapse = "; "
                )
            }
    )
}




# ============================================================
# 5. Retry optimizers
#
# Handles convergence only.
# Does NOT change random effects.
# ============================================================

failure_message <- function(history) {
    values <- unlist(history, recursive = TRUE, use.names = TRUE)
    messages <- as.character(values[grepl("(^|\\.)(message|error)$", names(values))])
    messages <- unique(messages[!is.na(messages) & nzchar(messages)])
    if (!length(messages)) return("Model fitting failed without a diagnostic message.")
    paste(messages, collapse = "; ")
}

fit_with_optimizer_retries <- function(
    formula,
    data,
    REML = FALSE
) {

    validate_fixed_effects_rank(lme4::nobars(formula), data)

    optimizers <- c(
        "nloptwrap",
        "bobyqa",
        "Nelder_Mead"
    )


    history <- list()

    last_fit <- NULL
    last_optimizer <- NA_character_
    last_message <- NA_character_


    for (optimizer in optimizers) {

        error_message <- NULL


        fit <- tryCatch(

            lmer(
                formula,
                data = data,
                REML = REML,

                control =
                    lmerControl(
                        optimizer =
                            optimizer
                    )
            ),

            error = function(e) {

                error_message <<-
                    conditionMessage(e)

                NULL
            }
        )


        if (is.null(fit)) {

            last_optimizer <- optimizer
            last_message <- error_message
            history[[optimizer]] <- list(
                optimizer = optimizer,
                fitted = FALSE,
                converged = FALSE,
                message = error_message
            )

            next
        }


        conv <- check_lmer_convergence(
            fit
        )


        history[[optimizer]] <- list(
            optimizer = optimizer,
            fitted = TRUE,
            converged = conv$converged,
            message =
                conv$convergence_message
        )


        last_fit <- fit
        last_optimizer <- optimizer
        last_message <-
            conv$convergence_message


        if (conv$converged) {

            return(
                list(
                    fit = fit,
                    converged = TRUE,
                    optimizer = optimizer,
                    convergence_message =
                        conv$convergence_message,
                    history = history
                )
            )
        }
    }


    list(
        fit = last_fit,
        converged = FALSE,
        optimizer = last_optimizer,
        convergence_message =
            last_message,
        history = history
    )
}


# ============================================================
# 6. Build random-effect term
# ============================================================

build_random_term <- function(
    subject_id,
    time,
    random_slope
) {

    if (random_slope) {

        return(
            paste0(
                "(1 + ",
                time,
                " | ",
                subject_id,
                ")"
            )
        )
    }


    paste0(
        "(1 | ",
        subject_id,
        ")"
    )
}



# ============================================================
# Run one longitudinal feature
# ============================================================

run_longitudinal_feature_impl <- function(
    data,
    outcome,
    predictor,
    predictor_type,
    time,
    subject_id,
    analysis_goal = c(
        "overall_predictor_association",
        "time_effect",
        "trajectory_difference"
    ),
    covariates = character(0),
    ordered_levels = NULL,
    random_slope = TRUE
) {

    analysis_goal <- match.arg(
        analysis_goal
    )


    # --------------------------------------------------------
    # Numeric outcome + numeric time
    # --------------------------------------------------------

    data[[outcome]] <- as.numeric(
        data[[outcome]]
    )

    data[[time]] <- as.numeric(
        data[[time]]
    )


    # --------------------------------------------------------
    # Prepare predictor
    #
    # strategy:
    # coefficient
    # ordered_trend
    # overall_and_pairwise
    # --------------------------------------------------------

    prepared <- prepare_primary_predictor(
        data = data,
        predictor = predictor,
        predictor_type = predictor_type,
        ordered_levels = ordered_levels
    )

    dat <- prepared$data

    model_predictor <-
        prepared$model_predictor

    strategy <-
        prepared$strategy


    # --------------------------------------------------------
    # Fixed random structure for this pass
    # --------------------------------------------------------

    random_term <- build_random_term(
        subject_id = subject_id,
        time = time,
        random_slope = random_slope
    )


    # ========================================================
    # CASE 1
    # Numeric / binary categorical / ordered trend
    # ========================================================

    if (
        strategy %in%
        c(
            "coefficient",
            "ordered_trend"
        )
    ) {

        # ----------------------------------------------------
        # Overall association
        # ----------------------------------------------------

        if (
            analysis_goal %in%
            c(
                "overall_predictor_association",
                "time_effect"
            )
        ) {

            # ------------------------------------------------
            # Same main-effects model for:
            # - overall predictor association
            # - time effect
            # ------------------------------------------------

            rhs <- build_rhs(
                terms = c(
                    model_predictor,
                    time
                ),
                covariates = covariates
            )

        } else {

            # ------------------------------------------------
            # Trajectory difference
            #
            # predictor * time
            # ------------------------------------------------

            interaction <- paste0(
                model_predictor,
                " * ",
                time
            )

            rhs <- build_rhs(
                terms = interaction,
                covariates = covariates
            )
        }


        formula <- as.formula(
            paste(
                outcome,
                "~",
                rhs,
                "+",
                random_term
            )
        )


        fit_result <- fit_with_optimizer_retries(
            formula = formula,
            data = dat,
            REML = FALSE
        )

        fit <- fit_result$fit


        # ----------------------------------------------------
        # Failed model
        # ----------------------------------------------------

        if (
            is.null(fit) ||
            !fit_result$converged
        ) {

            return(
                list(

                    overall = data.frame(
                        feature = outcome,
                        predictor = predictor,
                        predictor_type = predictor_type,
                        analysis_goal = analysis_goal,
                        test_type = strategy,
                        term = NA_character_,
                        estimate = NA_real_,
                        std_error = NA_real_,
                        p_value = NA_real_,
                        n_obs = if (
                            is.null(fit)
                        ) {
                            NA_real_
                        } else {
                            nobs(fit)
                        },
                        converged = FALSE,
                        optimizer = fit_result$optimizer,
                        singular = NA,
                        random_slope = random_slope,
                        stringsAsFactors = FALSE
                    ),

                    pairwise = NULL,

                    diagnostics =
                        fit_result$history
                )
            )
        }


        singular <- isSingular(
            fit,
            tol = 1e-4
        )


        # ----------------------------------------------------
        # Extract appropriate statistic
        # ----------------------------------------------------

        if (
            analysis_goal ==
            "overall_predictor_association"
        ) {

            stat <- extract_primary_coefficient(
                fit = fit,
                predictor = model_predictor,
                strategy = strategy
            )

        } else if (
            analysis_goal ==
            "time_effect"
        ) {

            stat <- extract_time_coefficient(
                fit = fit,
                time = time
            )

        } else {

            stat <- extract_interaction_coefficient(
                fit = fit,
                model_predictor = model_predictor,
                time = time
            )
        }


        return(
            list(

                overall = data.frame(
                    feature = outcome,
                    predictor = predictor,
                    predictor_type = predictor_type,
                    analysis_goal = analysis_goal,
                    test_type = strategy,
                    term = stat$term,
                    estimate = stat$estimate,
                    std_error = stat$std_error,
                    p_value = stat$p_value,
                    n_obs = nobs(fit),
                    converged = TRUE,
                    optimizer = fit_result$optimizer,
                    singular = singular,
                    random_slope = random_slope,
                    stringsAsFactors = FALSE
                ),

                pairwise = NULL,

                diagnostics =
                    fit_result$history
            )
        )
    }


    # ========================================================
    # CASE 2
    # Multi-category unordered predictor
    # ========================================================

    if (
        strategy ==
        "overall_and_pairwise"
    ) {

        # ====================================================
        # TIME EFFECT
        #
        # The multi-category predictor is included only as an
        # adjustment variable. No overall predictor test and
        # no pairwise comparisons are needed.
        # ====================================================

        if (
            analysis_goal ==
            "time_effect"
        ) {

            rhs <- build_rhs(
                terms = c(
                    model_predictor,
                    time
                ),
                covariates = covariates
            )

            formula <- as.formula(
                paste(
                    outcome,
                    "~",
                    rhs,
                    "+",
                    random_term
                )
            )

            fit_result <- fit_with_optimizer_retries(
                formula = formula,
                data = dat,
                REML = FALSE
            )

            fit <- fit_result$fit

            if (
                is.null(fit) ||
                !fit_result$converged
            ) {

                return(
                    list(

                        overall = data.frame(
                            feature = outcome,
                            predictor = predictor,
                            predictor_type = predictor_type,
                            analysis_goal = analysis_goal,
                            test_type = "time_effect",
                            term = time,
                            estimate = NA_real_,
                            std_error = NA_real_,
                            p_value = NA_real_,
                            n_obs = if (
                                is.null(fit)
                            ) {
                                NA_real_
                            } else {
                                nobs(fit)
                            },
                            converged = FALSE,
                            optimizer = fit_result$optimizer,
                            singular = NA,
                            random_slope = random_slope,
                            stringsAsFactors = FALSE
                        ),

                        pairwise = NULL,

                        diagnostics =
                            fit_result$history
                    )
                )
            }

            singular <- isSingular(
                fit,
                tol = 1e-4
            )

            stat <- extract_time_coefficient(
                fit = fit,
                time = time
            )

            return(
                list(

                    overall = data.frame(
                        feature = outcome,
                        predictor = predictor,
                        predictor_type = predictor_type,
                        analysis_goal = analysis_goal,
                        test_type = "time_effect",
                        term = stat$term,
                        estimate = stat$estimate,
                        std_error = stat$std_error,
                        p_value = stat$p_value,
                        n_obs = nobs(fit),
                        converged = TRUE,
                        optimizer = fit_result$optimizer,
                        singular = singular,
                        random_slope = random_slope,
                        stringsAsFactors = FALSE
                    ),

                    pairwise = NULL,

                    diagnostics =
                        fit_result$history
                )
            )
        }


        # ====================================================
        # Overall association / trajectory difference
        # ====================================================

        if (
            analysis_goal ==
            "overall_predictor_association"
        ) {

            full_rhs <- build_rhs(
                terms = c(
                    model_predictor,
                    time
                ),
                covariates = covariates
            )


            reduced_rhs <- build_rhs(
                terms = time,
                covariates = covariates
            )


        } else {

            # =================================================
            # Trajectory difference
            #
            # Full:
            # predictor * time
            #
            # Reduced:
            # predictor + time
            # =================================================

            full_rhs <- build_rhs(
                terms = paste0(
                    model_predictor,
                    " * ",
                    time
                ),
                covariates = covariates
            )


            reduced_rhs <- build_rhs(
                terms = c(
                    model_predictor,
                    time
                ),
                covariates = covariates
            )
        }


        full_formula <- as.formula(
            paste(
                outcome,
                "~",
                full_rhs,
                "+",
                random_term
            )
        )


        reduced_formula <- as.formula(
            paste(
                outcome,
                "~",
                reduced_rhs,
                "+",
                random_term
            )
        )


        full_result <- fit_with_optimizer_retries(
            formula = full_formula,
            data = dat,
            REML = FALSE
        )


        reduced_result <- fit_with_optimizer_retries(
            formula = reduced_formula,
            data = dat,
            REML = FALSE
        )


        fit1 <- full_result$fit
        fit0 <- reduced_result$fit


        # ----------------------------------------------------
        # Require both models to converge
        # ----------------------------------------------------

        if (
            is.null(fit0) ||
            is.null(fit1) ||
            !full_result$converged ||
            !reduced_result$converged
        ) {

            return(
                list(

                    overall = data.frame(
                        feature = outcome,
                        predictor = predictor,
                        predictor_type = predictor_type,
                        analysis_goal = analysis_goal,
                        test_type = if (
                            analysis_goal ==
                            "overall_predictor_association"
                        ) {
                            "overall_predictor_association"
                        } else {
                            "trajectory_difference"
                        },
                        term = if (
                		analysis_goal ==
                		"overall_predictor_association"
            			) {predictor} else {
                		paste0(
                    		predictor,
                    		":",
                    		time
                		)
            			},
                        estimate = NA_real_,
                        std_error = NA_real_,
                        p_value = NA_real_,
                        n_obs = NA_real_,
                        converged = FALSE,
                        optimizer = full_result$optimizer,
                        singular = NA,
                        random_slope = random_slope,
                        stringsAsFactors = FALSE
                    ),

                    pairwise = NULL,

                    diagnostics = list(
                        full = full_result$history,
                        reduced = reduced_result$history
                    )
                )
            )
        }


        # ----------------------------------------------------
        # Overall LRT
        # ----------------------------------------------------

        test <- anova(
            fit0,
            fit1,
            test = "LRT"
        )

        p <- test$`Pr(>Chisq)`[2]


        overall_result <- data.frame(
            feature = outcome,
            predictor = predictor,
            predictor_type = predictor_type,
            analysis_goal = analysis_goal,

            test_type = if (
                analysis_goal ==
                "overall_predictor_association"
            ) {
                "overall_predictor_association"
            } else {
                "trajectory_difference"
            },

            term = if (
                analysis_goal ==
                "overall_predictor_association"
            ) {
                predictor
            } else {
                paste0(
                    predictor,
                    ":",
                    time
                )
            },

            estimate = NA_real_,
            std_error = NA_real_,
            p_value = p,
            n_obs = nobs(fit1),
            converged = TRUE,
            optimizer = full_result$optimizer,

            singular = isSingular(
                fit1,
                tol = 1e-4
            ),

            random_slope = random_slope,

            stringsAsFactors = FALSE
        )


        # ----------------------------------------------------
        # Pairwise only for overall-association analysis
        #
        # For trajectory analysis, simple pairwise group
        # comparisons do not directly represent the
        # interaction, so leave them out for MVP.
        # ----------------------------------------------------

        pairwise_result <- NULL


        if (
            analysis_goal ==
            "overall_predictor_association"
        ) {

            emm <- emmeans(
                fit1,
                specs = predictor
            )


            pairwise_result <- as.data.frame(
                pairs(
                    emm,
                    adjust = "none"
                )
            )


            pairwise_result$feature <- outcome
            pairwise_result$predictor <- predictor
            pairwise_result$test_type <- "pairwise"
        }


        return(
            list(

                overall =
                    overall_result,

                pairwise =
                    pairwise_result,

                diagnostics =
                    list(
                        full = full_result$history,
                        reduced = reduced_result$history
                    )
            )
        )
    }


    stop(
        paste(
            "Unknown predictor strategy:",
            strategy
        )
    )
}


run_longitudinal_pass <- function(
    metadata,
    data_matrix,
    predictor,
    predictor_type,
    time,
    subject_id,
    analysis_goal,
    covariates = character(0),
    ordered_levels = NULL,
    random_slope = TRUE
) {

    aligned <- align_sample_tables(metadata, data_matrix)
    metadata <- aligned$metadata
    data_matrix <- aligned$data_matrix

    # ========================================================
    # Containers
    # ========================================================

    overall_results <- list()
    pairwise_results <- list()
    diagnostics <- list()


    # ========================================================
    # Outcomes / features
    #
    # Each column of data_matrix is treated as one
    # numeric outcome.
    # ========================================================

    features <- setdiff(
        colnames(
            data_matrix
        ),
        "SampleID"
    )


    # ========================================================
    # Loop through every outcome
    # ========================================================

    for (feature in features) {

        message(
            "Running: ",
            feature
        )


        # ----------------------------------------------------
        # Add current outcome to metadata.
        #
        # Sample tables were aligned before entering the feature loop.
        # ----------------------------------------------------

        dat <- metadata

        outcome_column <- make.unique(c(names(dat), "analysis_outcome"))[ncol(dat) + 1L]

        dat[[outcome_column]] <- (
            data_matrix[[feature]]
        )


        # ====================================================
        # Fit one longitudinal model
        #
        # run_longitudinal_feature() handles:
        #
        # - overall association
        # - trajectory difference
        #
        # - numeric predictor
        # - binary categorical predictor
        # - ordered categorical trend
        # - multi-category categorical predictor
        #
        # - optimizer retries
        # - convergence reporting
        # - singularity reporting
        #
        # It DOES NOT change random effects feature-by-feature.
        # ====================================================

        result <- tryCatch(

            {

                run_longitudinal_feature(

                    data = dat,

                    outcome = outcome_column,

                    predictor = predictor,

                    predictor_type =
                        predictor_type,

                    time = time,

                    subject_id =
                        subject_id,

                    analysis_goal =
                        analysis_goal,

                    covariates =
                        covariates,

                    ordered_levels =
                        ordered_levels,

                    random_slope =
                        random_slope
                )
            },


            # =================================================
            # If one feature fails completely, do not stop
            # the entire omics analysis.
            #
            # Instead return an NA result for that feature.
            # =================================================

            error = function(e) {

                list(

                    overall =
                        data.frame(

                            feature =
                                feature,

                            predictor =
                                predictor,

                            predictor_type =
                                predictor_type,

                            analysis_goal =
                                analysis_goal,

                            test_type =
                                NA_character_,

                            term =
                                NA_character_,

                            estimate =
                                NA_real_,

                            std_error =
                                NA_real_,

                            p_value =
                                NA_real_,

                            n_obs =
                                NA_real_,

                            converged =
                                FALSE,

                            optimizer =
                                NA_character_,

                            singular =
                                NA,

                            random_slope =
                                random_slope,

                            error =
                                conditionMessage(e),

                            stringsAsFactors =
                                FALSE
                        ),


                    pairwise =
                        NULL,


                    diagnostics =
                        list(error = conditionMessage(e))
                )
            }
        )


        # ====================================================
        # Store main result
        #
        # Depending on predictor / analysis type, this could be:
        #
        # numeric:
        #   predictor coefficient + p
        #
        # binary categorical:
        #   group coefficient + p
        #
        # ordered categorical:
        #   ordered trend coefficient + p
        #
        # multi-category:
        #   overall association / interaction p
        # ====================================================

        result$overall$feature <- feature
        if (!is.null(result$pairwise)) result$pairwise$feature <- feature

        overall_results[[feature]] <- (
            result$overall
        )


        # ====================================================
        # Store pairwise comparisons
        #
        # Pairwise results are only generated when appropriate,
        # such as a multi-category unordered predictor.
        #
        # At the feature level we use:
        #
        # pairs(
        #     emm,
        #     adjust = "none"
        # )
        #
        # so these p-values are RAW pairwise p-values.
        #
        # We deliberately do NOT apply Tukey adjustment here,
        # because the omics-wide multiplicity correction is
        # performed after all features are combined.
        # ====================================================

        if (
            !is.null(
                result$pairwise
            )
        ) {

            pairwise_results[[feature]] <- (
                result$pairwise
            )
        }


        # ====================================================
        # Store model diagnostics
        #
        # This can contain:
        #
        # - optimizers attempted
        # - convergence messages
        # - optimizer that succeeded
        #
        # Diagnostics are kept separately from the primary
        # results table.
        # ====================================================

        diagnostics[[feature]] <- (
            result$diagnostics
        )
    }


    # ========================================================
    # Combine primary / overall results
    # ========================================================

    overall_results <- bind_feature_results(overall_results)


    rownames(
        overall_results
    ) <- NULL


    # ========================================================
    # Omics-wide FDR for the PRIMARY test
    #
    # There is one primary p-value per outcome.
    #
    # Examples:
    #
    # Numeric predictor:
    #   beta_predictor p-value
    #
    # Binary predictor:
    #   group coefficient p-value
    #
    # Ordered predictor:
    #   ordered-trend p-value
    #
    # Multi-category predictor:
    #   overall predictor p-value
    #
    # Trajectory analysis:
    #   interaction p-value
    #
    # BH correction is therefore performed across features.
    # ========================================================

    overall_results$FDR <- p.adjust(
        overall_results$p_value,
        method = "BH"
    )

    # ========================================================
    # Convergence summary across features
    # ========================================================

    n_features <- nrow(
    overall_results)

    n_converged <- sum(
    overall_results$converged == TRUE,
    na.rm = TRUE)

    n_failed <- (
    n_features - n_converged)

    convergence_fraction <- if (
    n_features > 0) {

    n_converged / n_features

    } else {

    NA_real_}


    # ========================================================
    # Combine pairwise comparisons
    # ========================================================

    if (
        length(
            pairwise_results
        ) > 0
    ) {

        pairwise_results <- do.call(
            rbind,
            pairwise_results
        )


        rownames(
            pairwise_results
        ) <- NULL


        # ====================================================
        # Pairwise multiplicity correction
        #
        # IMPORTANT:
        #
        # At the feature level, pairwise comparisons were run:
        #
        # emm <- emmeans(
        #     fit1,
        #     specs = predictor
        # )
        #
        # pairwise_result <- as.data.frame(
        #     pairs(
        #         emm,
        #         adjust = "none"
        #     )
        # )
        #
        # Therefore pairwise_results$p.value contains the
        # unadjusted pairwise p-values.
        #
        #
        # We now apply BH FDR across omics features,
        # SEPARATELY for each contrast.
        #
        #
        # Example:
        #
        # Predictor levels:
        # A / B / C
        #
        # Pairwise contrasts:
        #
        # A - B
        # A - C
        # B - C
        #
        #
        # For A - B:
        #
        # Gene1 p
        # Gene2 p
        # Gene3 p
        # ...
        # Gene5000 p
        #
        #        ↓
        #
        # BH FDR across genes
        #
        #
        # Then separately:
        #
        # A - C
        #
        # Gene1 p
        # Gene2 p
        # ...
        #
        #        ↓
        #
        # another BH correction.
        #
        #
        # This preserves the interpretation that each
        # biological contrast represents a separate
        # hypothesis family across the omics features.
        # ====================================================

        pairwise_results$FDR <- ave(

            pairwise_results$p.value,

            pairwise_results$contrast,

            FUN = function(p) {

                p.adjust(
                    p,
                    method = "BH"
                )
            }
        )


    } else {

        pairwise_results <- NULL
    }


    # ========================================================
    # Return one complete longitudinal pass
    #
    # This function does NOT decide whether to switch from
    # random slope to random intercept.
    #
    # It simply runs every feature using the same requested
    # random-effects structure and reports singularity.
    #
    # The higher-level run_longitudinal_omics() function can
    # inspect singularity across all features and decide
    # whether a second pass should be performed.
    # ========================================================

    return(
        list(

            overall =
                overall_results,

            pairwise =
                pairwise_results,

            diagnostics =
                diagnostics,

	    n_features =
            n_features,

            n_converged =
            n_converged,

            n_failed =
            n_failed,

            convergence_fraction =
            convergence_fraction,


            random_slope =
                random_slope,

            analysis_goal =
                analysis_goal
        )
    )
}


run_longitudinal_omics <- function(
    metadata,
    data_matrix,
    predictor,
    predictor_type,
    time,
    subject_id,
    analysis_goal,
    covariates = character(0),
    ordered_levels = NULL,
    random_slope = TRUE,
    singularity_threshold = 0.20
) {

    # ========================================================
    # PASS 1
    #
    # Run all features using the random-effects structure
    # requested by the user.
    #
    # Example:
    #
    # random_slope = TRUE
    #
    # outcome ~ fixed effects +
    #           (1 + time | subject)
    #
    # Every feature gets the SAME random-effects structure.
    # Optimizers may differ by feature if retries are needed.
    # ========================================================

    first_pass <- run_longitudinal_pass(

        metadata = metadata,

        data_matrix = data_matrix,

        predictor = predictor,

        predictor_type = predictor_type,

        time = time,

        subject_id = subject_id,

        analysis_goal = analysis_goal,

        covariates = covariates,

        ordered_levels = ordered_levels,

        random_slope = random_slope
    )


    first_results <- (
        first_pass$overall
    )


    # ========================================================
    # CASE 1
    #
    # User requested random intercept only.
    #
    # No singularity-based random-slope decision is needed.
    # The first pass is already the final analysis.
    # ========================================================

    if (!random_slope) {

        return(
            list(

                results =
                    first_results,

                pairwise =
                    first_pass$pairwise,

                diagnostics =
                    first_pass$diagnostics,

                analysis_goal =
                    analysis_goal,

                random_slope_requested =
                    FALSE,

                final_random_slope =
                    FALSE,

                rerun =
                    FALSE,

	        # --------------------------------------------
                # Convergence summary
                # --------------------------------------------

                n_features =
                    first_pass$n_features,

                n_converged =
                    first_pass$n_converged,

                n_failed =
                    first_pass$n_failed,

                convergence_fraction =
                    first_pass$convergence_fraction,

                singular_fraction = converged_singular_fraction(first_results)
            )
        )
    }


    # ========================================================
    # Evaluate singularity
    #
    # IMPORTANT:
    #
    # Only evaluate singularity among models that actually
    # converged.
    #
    # A non-converged model should not be used to decide
    # whether the random-slope structure is supported.
    # ========================================================

    valid <- (
        first_results$converged == TRUE &
        !is.na(
            first_results$singular
        )
    )


    n_converged_for_singularity <- sum(
        valid
    )


    # ========================================================
    # No usable converged models
    #
    # We cannot make an omics-wide singularity decision.
    # Return the first pass and flag the issue.
    # ========================================================

    if (
        n_converged_for_singularity == 0
    ) {

        return(
            list(

                results =
                    first_results,

                pairwise =
                    first_pass$pairwise,

                diagnostics =
                    first_pass$diagnostics,

                analysis_goal =
                    analysis_goal,

                random_slope_requested =
                    TRUE,

                final_random_slope =
                    TRUE,

                rerun =
                    FALSE,
	       
	        # --------------------------------------------
                # Convergence summary
                # --------------------------------------------

                n_features =
                    first_pass$n_features,

                n_converged =
                    first_pass$n_converged,

                n_failed =
                    first_pass$n_failed,

                convergence_fraction =
                    first_pass$convergence_fraction,


                singular_fraction =
                    NA_real_,

                warning =
                    paste(
                        "No converged models were available",
                        "for evaluating random-slope singularity."
                    )
            )
        )
    }


    # ========================================================
    # Fraction of converged models that are singular
    # ========================================================

    singular_fraction <- mean(
        first_results$singular[
            valid
        ]
    )


    # ========================================================
    # CASE 2
    #
    # Singularity is below configured threshold.
    #
    # Keep the random-slope analysis.
    #
    # Individual singular features remain flagged in the
    # result table, but we do NOT change their models
    # feature-by-feature.
    # ========================================================

    if (
        singular_fraction <
        singularity_threshold
    ) {

        return(
            list(

                results =
                    first_results,

                pairwise =
                    first_pass$pairwise,

                diagnostics =
                    first_pass$diagnostics,

                analysis_goal =
                    analysis_goal,

                random_slope_requested =
                    TRUE,

                final_random_slope =
                    TRUE,

                rerun =
                    FALSE,

	    	# --------------------------------------------
                # Convergence summary
                # --------------------------------------------

                n_features =
                    first_pass$n_features,

                n_converged =
                    first_pass$n_converged,

                n_failed =
                    first_pass$n_failed,

                convergence_fraction =
                    first_pass$convergence_fraction,

                singular_fraction =
                    singular_fraction,

                singularity_threshold =
                    singularity_threshold
            )
        )
    }


    # ========================================================
    # CASE 3
    #
    # Random-slope singularity is widespread.
    #
    # Rather than changing only the affected features,
    # rerun ALL features using random intercept only.
    #
    # This keeps the final omics analysis statistically
    # consistent across outcomes.
    # ========================================================

    second_pass <- run_longitudinal_pass(

        metadata = metadata,

        data_matrix = data_matrix,

        predictor = predictor,

        predictor_type = predictor_type,

        time = time,

        subject_id = subject_id,

        analysis_goal = analysis_goal,

        covariates = covariates,

        ordered_levels = ordered_levels,

        random_slope = FALSE
    )


    # ========================================================
    # Return second pass as FINAL analysis.
    #
    # Keep the first-pass results for diagnostics/audit trail.
    # ========================================================

    return(
        list(

            results =
                second_pass$overall,

            pairwise =
                second_pass$pairwise,

            diagnostics =
                second_pass$diagnostics,

            analysis_goal =
                analysis_goal,

            random_slope_requested =
                TRUE,

            final_random_slope =
                FALSE,

            rerun =
                TRUE,

	    # ------------------------------------------------
            # FINAL convergence summary
            # ------------------------------------------------

            n_features =
                second_pass$n_features,

            n_converged =
                second_pass$n_converged,

            n_failed =
                second_pass$n_failed,

            convergence_fraction =
                second_pass$convergence_fraction,

            singular_fraction = converged_singular_fraction(second_pass$overall),

            first_pass_singular_fraction = singular_fraction,

            singularity_threshold =
                singularity_threshold,

            rerun_reason =
                paste(
                    "Random-slope singularity exceeded",
                    "the configured omics-wide threshold."
                ),

            first_pass_results =
                first_results,

            first_pass_pairwise =
                first_pass$pairwise,

            first_pass_diagnostics =
                first_pass$diagnostics
        )
    )
}
