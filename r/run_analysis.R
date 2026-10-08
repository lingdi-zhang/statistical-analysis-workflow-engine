# ============================================================
# run_analysis.R
#
# Entry point called by Python.
#
# Responsibilities:
#   1. Read metadata
#   2. Read outcome / feature matrix
#   3. Read analysis request JSON
#   4. Load model_runner.R
#   5. Run cross-sectional or longitudinal analysis
#   6. Write results for Python
# ============================================================


library(jsonlite)


# ============================================================
# 1. Load modeling functions
# ============================================================

args_full <- commandArgs(trailingOnly = FALSE)

script_arg <- grep(
    "^--file=",
    args_full,
    value = TRUE
)

script_path <- sub(
    "^--file=",
    "",
    script_arg
)

script_dir <- dirname(
    normalizePath(script_path)
)

model_runner_path <- file.path(
    script_dir,
    "model_runner.R"
)

prepare_model_variables_path <- file.path(
    script_dir,
    "prepare_model_variables.R"
)

if (!file.exists(model_runner_path)) {
    stop(
        paste(
            "model_runner.R not found at:",
            model_runner_path
        )
    )
}

if (!file.exists(prepare_model_variables_path)) {
    stop(
        paste(
            "prepare_model_variables.R not found at:",
            prepare_model_variables_path
        )
    )
}

source(model_runner_path)
source(prepare_model_variables_path)

# ============================================================
# 2. Read command-line arguments
#
# Python will call:
#
# Rscript R/run_analysis.R \
#     metadata.csv \
#     data_matrix.csv \
#     request.json \
#     overall_results.csv \
#     pairwise_results.csv \
#     summary.json
# ============================================================

args <- commandArgs(
    trailingOnly = TRUE
)


if (length(args) != 6) {

    stop(
        paste(
            "Expected 6 arguments:",
            "metadata_path,",
            "data_matrix_path,",
            "request_path,",
            "overall_output,",
            "pairwise_output,",
            "summary_output"
        )
    )
}


metadata_path <- args[1]

data_matrix_path <- args[2]

request_path <- args[3]

overall_output <- args[4]

pairwise_output <- args[5]

summary_output <- args[6]


# ============================================================
# 3. Read input data
# ============================================================

request <- fromJSON(request_path, simplifyVector = TRUE)

# Python supplies a collision-checked marker; retain compatibility with direct
# command-line CSV inputs, whose blank cells represent missing values.
missing_values <- request$transfer_na_marker
if (is.null(missing_values)) missing_values <- ""

metadata <- read.csv(
    metadata_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    na.strings = missing_values
)

data_matrix <- read.csv(
    data_matrix_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    na.strings = missing_values
)


# ============================================================
# 4. Extract common request parameters
# ============================================================

if (!all(vapply(list(metadata, data_matrix), function(table) {
    "SampleID" %in% names(table)
}, logical(1)))) stop("Both CSV inputs must contain SampleID.")

analysis_type <- request$analysis_type


# ------------------------------------------------------------
# Primary predictor
#
# Python request stores this as:
#
# "primary_predictors": ["FA"]
# ------------------------------------------------------------

predictor <- request$primary_predictors[[1]]


# ------------------------------------------------------------
# Covariates
# ------------------------------------------------------------

covariates <- request$covariates


if (is.null(covariates)) {

    covariates <- character(0)
}


# ------------------------------------------------------------
# Predictor type
#
# Example:
#
# "variable_types": {
#     "FA": "categorical",
#     "Age": "numeric"
# }
# ------------------------------------------------------------

predictor_type <- request$variable_types[[predictor]]


if (is.null(predictor_type)) {

    stop(
        paste(
            "Predictor type was not specified for",
            predictor
        )
    )
}


# ============================================================
# 5. Ordered levels
#
# Only needed if predictor is ordered categorical.
# ============================================================

ordered_levels <- NULL


if (
    predictor_type ==
    "ordered_categorical"
) {

    ordered_levels <- (
        request$ordered_levels[[predictor]]
    )


    if (is.null(ordered_levels)) {

        stop(
            paste(
                "Ordered predictor",
                predictor,
                "does not have ordered_levels."
            )
        )
    }
}

# ============================================================
# 6. prepare variables 
# ============================================================
if (analysis_type == "longitudinal" && !is.null(request$time)) {
    request$variable_types[[request$time]] <- "numeric"
}
metadata <- prepare_model_variables(
    dat = metadata,
    variable_types = request$variable_types,
    ordered_levels = request$ordered_levels
)

# ============================================================
# 6. CROSS-SECTIONAL ANALYSIS
# ============================================================

if (
    analysis_type ==
    "cross_sectional"
) {

    message(
        "Running cross-sectional analysis."
    )


    # --------------------------------------------------------
    # Run analysis
    # --------------------------------------------------------

    results <- run_cross_sectional_omics(

        metadata =
            metadata,

        data_matrix =
            data_matrix,

        predictor =
            predictor,

        predictor_type =
            predictor_type,

        covariates =
            covariates,

        ordered_levels =
            ordered_levels
    )


    # --------------------------------------------------------
    # Normalize main results
    # --------------------------------------------------------

    final_results <- (
        results$overall
    )


    final_pairwise <- (
        results$pairwise
    )


    # --------------------------------------------------------
    # Summary returned to Python
    # --------------------------------------------------------

    analysis_summary <- list(

        analysis_type =
            "cross_sectional",

        predictor =
            predictor,

        predictor_type =
            predictor_type,

        n_features = nrow(final_results),
        n_successful = sum(is.finite(final_results$p_value)),
        n_estimate_only = sum(final_results$status == "estimate_only", na.rm = TRUE),
        n_failed = sum(!is.finite(final_results$p_value) &
                       final_results$status != "estimate_only", na.rm = TRUE)
    )



# ============================================================
# 7. LONGITUDINAL ANALYSIS
# ============================================================

} else if (
    analysis_type ==
    "longitudinal"
) {

    message(
        "Running longitudinal analysis."
    )


    # --------------------------------------------------------
    # Required longitudinal parameters
    # --------------------------------------------------------

    time <- request$time

    subject_id <- request$subject_id

    analysis_goal <- (
        request$analysis_goal
    )


    if (is.null(time)) {

        stop(
            "Time variable was not specified."
        )
    }


    if (is.null(subject_id)) {

        stop(
            "Subject ID was not specified."
        )
    }


    if (is.null(analysis_goal)) {

        stop(
            "Longitudinal analysis goal was not specified."
        )
    }


    # --------------------------------------------------------
    # Random slope
    #
    # Default to TRUE if Python request does not provide it.
    # --------------------------------------------------------

    random_slope <- TRUE


    if (
        !is.null(
            request$random_slope
        )
    ) {

        random_slope <- (
            request$random_slope
        )
    }


    # --------------------------------------------------------
    # Singularity threshold
    #
    # Default for MVP = 0.20
    #
    # This is configurable rather than hard-coded into
    # model_runner.R.
    # --------------------------------------------------------

    singularity_threshold <- 0.20


    if (
        !is.null(
            request$singularity_threshold
        )
    ) {

        singularity_threshold <- (
            request$singularity_threshold
        )
    }


    # --------------------------------------------------------
    # Run longitudinal omics model
    # --------------------------------------------------------

    results <- run_longitudinal_omics(

        metadata =
            metadata,

        data_matrix =
            data_matrix,

        predictor =
            predictor,

        predictor_type =
            predictor_type,

        time =
            time,

        subject_id =
            subject_id,

        analysis_goal =
            analysis_goal,

        covariates =
            covariates,

        ordered_levels =
            ordered_levels,

        random_slope =
            random_slope,

        singularity_threshold =
            singularity_threshold
    )


    # --------------------------------------------------------
    # Normalize main results
    #
    # run_longitudinal_omics() returns:
    #
    # results$results
    # --------------------------------------------------------

    final_results <- (
        results$results
    )


    final_pairwise <- (
        results$pairwise
    )


    # --------------------------------------------------------
    # Build summary for Python / Streamlit
    # --------------------------------------------------------

    analysis_summary <- list(

        analysis_type =
            "longitudinal",

        analysis_goal =
            analysis_goal,

        predictor =
            predictor,

        predictor_type =
            predictor_type,

        time =
            time,

        subject_id =
            subject_id,

        n_features =
            results$n_features,

        n_converged =
            results$n_converged,

        n_failed =
            results$n_failed,

        convergence_fraction =
            results$convergence_fraction,

        random_slope_requested =
            results$random_slope_requested,

        final_random_slope =
            results$final_random_slope,

        rerun =
            results$rerun
    )


    # --------------------------------------------------------
    # Add singularity information if available
    # --------------------------------------------------------

    if (
        !is.null(
            results$singular_fraction
        )
    ) {

        analysis_summary$singular_fraction <- (
            results$singular_fraction
        )
    }


    if (
        !is.null(
            results$first_pass_singular_fraction
        )
    ) {

        analysis_summary$
            first_pass_singular_fraction <- (
                results$
                    first_pass_singular_fraction
            )
    }


    if (
        !is.null(
            results$rerun_reason
        )
    ) {

        analysis_summary$rerun_reason <- (
            results$rerun_reason
        )
    }


    if (
        !is.null(
            results$warning
        )
    ) {

        analysis_summary$warning <- (
            results$warning
        )
    }



# ============================================================
# 8. Unknown analysis type
# ============================================================

} else {

    stop(
        paste(
            "Unknown analysis type:",
            analysis_type
        )
    )
}


# ============================================================
# 9. Write primary / overall results
# ============================================================

write.csv(

    final_results,

    overall_output,

    row.names = FALSE,
    na = missing_values[[1]]
)


# ============================================================
# 10. Write pairwise results
#
# Pairwise results only exist when relevant.
#
# If none exist, create an empty file so Python knows
# the analysis completed successfully but there are no
# pairwise results.
# ============================================================

if (
    !is.null(
        final_pairwise
    ) &&
    nrow(final_pairwise) > 0
) {

    write.csv(

        final_pairwise,

        pairwise_output,

        row.names = FALSE,
        na = missing_values[[1]]
    )

} else {

    file.create(
        pairwise_output
    )
}


# ============================================================
# 11. Write analysis summary
# ============================================================

write_json(

    analysis_summary,

    summary_output,

    pretty = TRUE,

    auto_unbox = TRUE,

    na = "null"
)


# Preserve final and first-pass optimizer histories independently of the summary.
diagnostics_output <- file.path(dirname(summary_output), "diagnostics.json")
write_json(list(final = results$diagnostics,
                first_pass = results$first_pass_diagnostics),
           diagnostics_output, pretty = TRUE, auto_unbox = TRUE, na = "null")

# ============================================================
# 12. Console completion message
# ============================================================

message(
    "Analysis completed successfully."
)
