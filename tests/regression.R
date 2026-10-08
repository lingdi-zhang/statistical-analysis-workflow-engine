# Run from the repository root: Rscript tests/regression.R
# Base-R checks run even when the optional mixed-model packages are absent.
for (expression in parse("r/model_runner.R")) {
    if (is.call(expression) && identical(expression[[1]], as.name("<-"))) {
        eval(expression)
    }
}

# Intercept-only summaries count singularity among converged fits.
local({
    original_pass <- run_longitudinal_pass
    on.exit(assign("run_longitudinal_pass", original_pass, envir = .GlobalEnv))
    assign("run_longitudinal_pass", function(...) list(
        overall = data.frame(converged = c(TRUE, TRUE, FALSE),
                             singular = c(TRUE, FALSE, NA)),
        n_features = 3L, n_converged = 2L, n_failed = 1L,
        convergence_fraction = 2 / 3), envir = .GlobalEnv)
    result <- run_longitudinal_omics(NULL, NULL, "x", "numeric", "time", "id",
                                      "time_effect", random_slope = FALSE)
    stopifnot(identical(result$singular_fraction, 0.5))
})

set.seed(42)
metadata <- data.frame(SampleID = paste0("S", 1:60),
                       x = rnorm(60), age = rnorm(60))
outcome <- 2 * metadata$x + metadata$age + rnorm(60)
features <- data.frame(SampleID = metadata$SampleID, age = outcome)
actual <- run_cross_sectional_omics(metadata, features, "x", "numeric", "age")$overall
expected <- summary(lm(outcome ~ x + age, data = metadata))$coefficients["x", ]
stopifnot(actual$feature == "age", actual$term == "x",
          isTRUE(all.equal(actual$estimate, unname(expected["Estimate"]))))

metadata[["dose (mg)"]] <- metadata$x
metadata[["gene-1"]] <- outcome
result <- run_cross_sectional_feature(metadata, "gene-1", "dose (mg)", "numeric", "age")
stopifnot(result$overall$feature == "gene-1", result$overall$term == "dose (mg)",
          isTRUE(all.equal(result$overall$estimate, actual$estimate)))

metadata$group <- rep(c("a", "b", "c"), 20)
metadata$group[c(1, 4)] <- NA
prepared <- prepare_feature_model(metadata, "gene-1", "group", "categorical", "age")
m <- prepared$mapping
full <- lm(reformulate(unname(m[c("group", "age")]), unname(m["gene-1"])), prepared$data)
reduced <- lm(reformulate(unname(m["age"]), unname(m["gene-1"])), prepared$data)
stopifnot(nobs(full) == 58, nobs(reduced) == 58,
          is.finite(anova(reduced, full)$`Pr(>F)`[2]))

lost_category <- metadata
lost_category[["gene-1"]][!is.na(lost_category$group) & lost_category$group == "c"] <- NA
failure <- tryCatch(run_cross_sectional_feature(lost_category, "gene-1", "group", "categorical"),
                    error = function(e) conditionMessage(e))
stopifnot(is.character(failure), grepl("No usable observations", failure, fixed = TRUE))
visits <- data.frame(y = 1:6, x = 1:6, time = factor(rep(c(0, 12), 3)),
                     subject = rep(1:3, each = 2))
prepared_visits <- prepare_feature_model(visits, "y", "x", "numeric", character(0),
                                         "time", "subject")
stopifnot(identical(prepared_visits$data[[prepared_visits$mapping[["time"]]]],
                    rep(c(0, 12), 3)))

# Rank is checked per feature, after missing outcomes select its model rows.
rank_metadata <- metadata
rank_metadata$covariate <- rank_metadata$x
rank_metadata$covariate[31:60] <- rnorm(30)
rank_features <- data.frame(good = outcome, bad = c(outcome[1:30], rep(NA_real_, 30)))
rank_results <- run_cross_sectional_omics(rank_metadata, rank_features, "x", "numeric",
                                         "covariate")$overall
stopifnot(is.finite(rank_results$p_value[rank_results$feature == "good"]),
          is.na(rank_results$p_value[rank_results$feature == "bad"]),
          grepl("rank deficient", rank_results$error[rank_results$feature == "bad"]))
confounded <- data.frame(y = outcome, group = rep(c("a", "b", "c"), 20))
confounded$batch <- factor(confounded$group)
failure <- tryCatch(run_cross_sectional_feature(confounded, "y", "group", "categorical", "batch"),
                    error = function(e) conditionMessage(e))
stopifnot(is.character(failure), grepl("rank deficient", failure))

# Factor labels containing punctuation must not be mistaken for interactions.
set.seed(81)
for (label in c("treated:high", "treated high", "treated(time)-high")) {
    punctuated <- data.frame(y = rnorm(40),
                             group = rep(c("control", label), 20),
                             time = rep(0:3, each = 10))
    simple <- punctuated
    simple$group <- rep(c("control", "treated"), 20)
    actual <- run_cross_sectional_feature(punctuated, "y", "group", "categorical")$overall
    expected <- run_cross_sectional_feature(simple, "y", "group", "categorical")$overall
    stopifnot(isTRUE(all.equal(actual$estimate, expected$estimate)),
              isTRUE(all.equal(actual$p_value, expected$p_value)))
    fit <- lm(y ~ group * time, punctuated)
    baseline <- lm(y ~ group * time, simple)
    for (variables in list("group", c("group", "time"), c("time", "group"))) {
        term <- coefficient_for_variables(fit, variables)
        reference <- coefficient_for_variables(baseline, variables)
        stopifnot(isTRUE(all.equal(unname(coef(fit)[term]),
                                  unname(coef(baseline)[reference]))))
    }
}

# Internal-looking names and category labels survive label restoration.
set.seed(115)
collision_data <- data.frame(y = rnorm(40), modelvar000003X = rnorm(40), age = rnorm(40))
collision_result <- run_cross_sectional_feature(collision_data, "y", "modelvar000003X",
                                               "numeric", "age")$overall
stopifnot(collision_result$predictor == "modelvar000003X",
          collision_result$term == "modelvar000003X")
collision_data$group <- rep(c("control", "modelvarX000003X:treated"), 20)
prepared <- prepare_feature_model(collision_data, "y", "group", "categorical", "age")
stopifnot(!any(vapply(prepared$mapping, function(token) {
    any(grepl(token, c(names(collision_data), collision_data$group), fixed = TRUE))
}, logical(1))))
group_result <- run_cross_sectional_feature(collision_data, "y", "group", "categorical", "age")$overall
stopifnot(group_result$predictor == "group",
          group_result$term == "groupmodelvarX000003X:treated")
# Restoration is simultaneous, including interaction terms and missing labels.
manual_mapping <- c(modelvar000003X = "modelvar000002X", age = "modelvar000003X")
restored <- restore_model_labels(list(overall = data.frame(
    feature = c("modelvar000002X", NA), predictor = c("modelvar000002X", NA),
    term = c("modelvar000002X:modelvar000003X", NA))), manual_mapping)$overall
stopifnot(restored$feature[1] == "modelvar000003X",
          restored$predictor[1] == "modelvar000003X",
          restored$term[1] == "modelvar000003X:age", is.na(restored$term[2]))

# Failure messages preserve errors across all optimizer attempts and nested fits.
history <- list(full = list(nloptwrap = list(message = "optimizer failed")),
                reduced = list(bobyqa = list(message = "model did not converge")))
stopifnot(grepl("optimizer failed", failure_message(history), fixed = TRUE),
          grepl("model did not converge", failure_message(history), fixed = TRUE))

# Saturated regressions retain estimates while withholding unsupported inference.
small <- data.frame(y = c(2, 3), x = c(0, 1))
for (type in c("numeric", "categorical", "ordered_categorical")) {
    result <- run_cross_sectional_feature(small, "y", "x", type,
        ordered_levels = if (type == "ordered_categorical") c("0", "1") else NULL)$overall
    stopifnot(result$status == "estimate_only", isTRUE(all.equal(result$estimate, 1)),
              is.na(result$std_error), is.na(result$p_value),
              grepl("insufficient residual degrees", result$warning, fixed = TRUE))
}
saturated_metadata <- data.frame(x = c(0, 1, 2, 3, 4, 5),
                                 covariate = c(1, 0, 1, 0, 1, 0))
saturated_features <- data.frame(good = c(2, 4, 3, 7, 5, 8),
                                 limited = c(2, 4, 3, NA, NA, NA))
saturated <- run_cross_sectional_omics(saturated_metadata, saturated_features,
                                      "x", "numeric", "covariate")$overall
stopifnot(is.finite(saturated$p_value[saturated$feature == "good"]),
          saturated$status[saturated$feature == "limited"] == "estimate_only",
          is.finite(saturated$estimate[saturated$feature == "limited"]),
          is.na(saturated$FDR[saturated$feature == "limited"]))

# Unsupported slopes are rejected before fitting, with an actionable message.
two_visits <- data.frame(y = 1:8, x = c(1, 3, 4, 2, 6, 5, 8, 7),
                         time = rep(0:1, 4), subject = rep(1:4, each = 2))
unsupported <- tryCatch(run_longitudinal_feature(two_visits, "y", "x", "numeric",
    "time", "subject", "time_effect", random_slope = TRUE),
    error = function(e) conditionMessage(e))
stopifnot(is.character(unsupported),
          grepl("8 usable observations for 4 subjects", unsupported, fixed = TRUE),
          grepl("Select random intercept only", unsupported, fixed = TRUE))

# Sample IDs, rather than row order, determine the outcome/predictor pairing.
set.seed(94)
alignment_metadata <- data.frame(SampleID = sprintf("%03d", 1:20), x = rnorm(20))
alignment_features <- data.frame(SampleID = alignment_metadata$SampleID,
                                 gene = 2 * alignment_metadata$x + rnorm(20, sd = 0.1))
reference <- run_cross_sectional_omics(alignment_metadata, alignment_features,
                                       "x", "numeric")$overall
shuffled <- run_cross_sectional_omics(alignment_metadata, alignment_features[20:1, ],
                                      "x", "numeric")$overall
stopifnot(isTRUE(all.equal(reference, shuffled)))
aligned <- suppressWarnings(align_sample_tables(alignment_metadata,
                                                 alignment_features[20:2, ]))
stopifnot(identical(aligned$metadata$SampleID, aligned$data_matrix$SampleID),
          nrow(aligned$metadata) == 19L)
for (bad_ids in list(c(NA, alignment_features$SampleID[-1]),
                     c("", alignment_features$SampleID[-1]),
                     rep("001", 20), paste0("other", 1:20))) {
    bad_features <- alignment_features
    bad_features$SampleID <- bad_ids
    failure <- tryCatch(align_sample_tables(alignment_metadata, bad_features),
                        error = function(e) conditionMessage(e))
    stopifnot(is.character(failure))
}

# Fallback exposes final and first-pass singularity separately.
local({
    original_pass <- run_longitudinal_pass
    on.exit(assign("run_longitudinal_pass", original_pass, envir = .GlobalEnv))
    assign("run_longitudinal_pass", function(..., random_slope) list(
        overall = data.frame(converged = c(TRUE, TRUE, FALSE),
                             singular = c(random_slope, random_slope, NA)),
        n_features = 3L, n_converged = 2L, n_failed = 1L,
        convergence_fraction = 2 / 3), envir = .GlobalEnv)
    result <- run_longitudinal_omics(NULL, NULL, "x", "numeric", "time", "id",
                                      "time_effect", random_slope = TRUE)
    stopifnot(isTRUE(result$rerun), identical(result$singular_fraction, 0),
              identical(result$first_pass_singular_fraction, 1))
    stopifnot(is.na(converged_singular_fraction(data.frame(
        converged = FALSE, singular = NA))))
})

# Strict numeric conversion rejects malformed data while retaining true missing values.
source("r/prepare_model_variables.R")
for (bad in c("invalid", "NaN", "Inf", "-Inf")) {
    numeric_metadata <- data.frame(x = c("1", bad, NA_character_),
                                    age = c("20", "21", "22"))
    for (variable in c("x", "age")) {
        numeric_metadata[[variable]][2] <- bad
        failure <- tryCatch(prepare_model_variables(numeric_metadata,
            setNames(list("numeric"), variable)), error = function(e) conditionMessage(e))
        stopifnot(is.character(failure), grepl(variable, failure, fixed = TRUE))
    }
    bad_outcomes <- alignment_features
    bad_outcomes$gene <- as.character(bad_outcomes$gene)
    bad_outcomes$gene[1] <- bad
    bad_outcomes$good <- alignment_features$gene
    checked <- run_cross_sectional_omics(alignment_metadata, bad_outcomes, "x", "numeric")$overall
    stopifnot(checked$status[checked$feature == "gene"] == "failed",
              nzchar(checked$error[checked$feature == "gene"]),
              checked$status[checked$feature == "good"] == "success")
}
converted <- prepare_model_variables(data.frame(x = c("1", NA, "3")), list(x = "numeric"))
stopifnot(identical(converted$x, c(1, NA_real_, 3)))
factor_values <- prepare_model_variables(data.frame(x = factor(c("10", "20"))), list(x = "numeric"))
stopifnot(identical(factor_values$x, c(10, 20)))
missing_outcome <- alignment_features
missing_outcome$gene[1] <- NA_real_
checked <- run_cross_sectional_omics(alignment_metadata, missing_outcome, "x", "numeric")$overall
stopifnot(checked$status == "success", checked$n_obs == 19)

# Exercise the entry point's real CSV-reader expressions without optional packages.
local({
    directory <- tempfile("csv-label-regression-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE))
    metadata_path <- file.path(directory, "metadata.csv")
    data_matrix_path <- file.path(directory, "features.csv")
    writeLines(c("SampleID,group,x", "NA,NA,1", "null,null,2", "S3,,"), metadata_path)
    writeLines(c("SampleID,gene", "NA,1", "null,2", "S3,"), data_matrix_path)
    request <- list()
    for (expression in parse("r/run_analysis.R")) {
        if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
            as.character(expression[[2]]) %in% c("missing_values", "metadata", "data_matrix")) {
            eval(expression)
        } else if (is.call(expression) && identical(expression[[1]], as.name("if")) &&
                   grepl("is.null(missing_values)", paste(deparse(expression), collapse = ""), fixed = TRUE)) {
            eval(expression)
        }
    }
    stopifnot(identical(metadata$SampleID, c("NA", "null", "S3")),
              identical(metadata$group, c("NA", "null", NA_character_)),
              is.na(metadata$x[3]), is.na(data_matrix$gene[3]),
              identical(data_matrix$SampleID, metadata$SampleID))
})

packages <- c("lme4", "lmerTest", "emmeans", "jsonlite")
if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
    source("r/model_runner.R")
    saturated_groups <- data.frame(y = c(1, 2, 4), group = c("a", "b", "c"))
    saturated_group_result <- run_cross_sectional_feature(saturated_groups, "y", "group", "categorical")
    stopifnot(saturated_group_result$overall$status == "estimate_only",
              is.na(saturated_group_result$overall$p_value),
              all(is.finite(saturated_group_result$pairwise$estimate)),
              all(is.na(saturated_group_result$pairwise$SE)),
              all(is.na(saturated_group_result$pairwise$p.value)))

    # A feature with only one usable subject fails, while other features continue.
    set.seed(102)
    failure_metadata <- data.frame(subject = rep(1:8, each = 4),
                                   time = rep(0:3, 8), x = rnorm(32))
    failure_features <- data.frame(good = rnorm(32),
                                   bad = c(rnorm(4), rep(NA_real_, 28)))
    failure_pass <- run_longitudinal_pass(failure_metadata, failure_features,
        "x", "numeric", "time", "subject", "time_effect", random_slope = FALSE)
    stopifnot(is.finite(failure_pass$overall$p_value[failure_pass$overall$feature == "good"]),
              nzchar(failure_pass$overall$error[failure_pass$overall$feature == "bad"]),
              length(failure_pass$diagnostics$bad) > 0)

    # A successful singular fit must count toward random-slope recovery.
    set.seed(4)
    recovery_data <- data.frame(subject = rep(1:40, each = 3),
                                time = rep(0:2, 40), x = rep(rnorm(40), each = 3))
    recovery_data$y <- rep(rnorm(40, sd = 3), each = 3) + recovery_data$x +
        0.3 * recovery_data$time + rnorm(120)
    singular_fit <- lmer(y ~ x + time + (1 + time | subject), recovery_data, REML = FALSE)
    stopifnot(isSingular(singular_fit), check_lmer_convergence(singular_fit)$converged)
    bad_fit <- singular_fit
    bad_fit@optinfo$conv$opt <- 1L
    stopifnot(!check_lmer_convergence(bad_fit)$converged)
    bad_fit <- singular_fit
    bad_fit@optinfo$conv$lme4$messages <- "Model failed to converge with max|grad|"
    stopifnot(!check_lmer_convergence(bad_fit)$converged)
    recovered <- run_longitudinal_omics(recovery_data, data.frame(gene = recovery_data$y),
                                        "x", "numeric", "time", "subject", "time_effect")
    stopifnot(recovered$rerun, !recovered$final_random_slope,
              recovered$n_converged == 1, is.finite(recovered$results$p_value))
    recovery_data$covariate <- recovery_data$x
    recovery_data$covariate[61:120] <- rnorm(60)
    rank_mixed <- run_longitudinal_omics(
        recovery_data, data.frame(good = recovery_data$y,
                                  bad = c(recovery_data$y[1:60], rep(NA_real_, 60))),
        "x", "numeric", "time", "subject", "time_effect",
        covariates = "covariate", random_slope = FALSE)
    stopifnot(rank_mixed$n_failed == 1, rank_mixed$n_converged == 1,
              grepl("rank deficient", rank_mixed$results$error[rank_mixed$results$feature == "bad"]),
              grepl("rank deficient", rank_mixed$diagnostics$bad$error))
    categorical <- run_cross_sectional_feature(metadata, "gene-1", "group", "categorical", "age")
    stopifnot(categorical$overall$n_obs == 58, is.finite(categorical$overall$p_value))
    metadata$subject <- rep(1:20, each = 3)
    metadata$time <- rep(0:2, 20)
    # Keep group independent of visit; the original repeating a/b/c fixture
    # would make group perfectly confounded with time.
    metadata$group <- rep(rep(c("a", "b", "c"), each = 6), length.out = 60)
    metadata$group[c(1, 4)] <- NA
    metadata[["gene-1"]] <- rep(rnorm(20, sd = 3), each = 3) +
        0.5 * metadata$time + rnorm(60)
    longitudinal <- run_longitudinal_feature(metadata, "gene-1", "group", "categorical",
                                             "time", "subject", random_slope = FALSE)
    stopifnot(longitudinal$overall$n_obs == 58, is.finite(longitudinal$overall$p_value))

    # Exercise the CLI and its CSV/JSON output contract with generated sample data.
    directory <- tempfile("analysis-regression-")
    dir.create(directory)
    request <- file.path(directory, "request.json")
    set.seed(210)
    demo_metadata <- data.frame(SampleID = sprintf("S%03d", 1:60), phenotype = rnorm(60))
    demo_features <- data.frame(SampleID = demo_metadata$SampleID,
                                gene1 = 2 * demo_metadata$phenotype + rnorm(60),
                                gene2 = -demo_metadata$phenotype + rnorm(60))
    metadata_input <- file.path(directory, "metadata.csv")
    features_input <- file.path(directory, "features.csv")
    write.csv(demo_metadata, metadata_input, row.names = FALSE)
    write.csv(demo_features, features_input, row.names = FALSE)
    predictor <- "phenotype"
    predictor_type <- if (is.numeric(demo_metadata[[predictor]])) "numeric" else "categorical"
    jsonlite::write_json(list(analysis_type = "cross_sectional",
                             primary_predictors = list(predictor), covariates = list(),
                             variable_types = setNames(list(predictor_type), predictor)), request,
                         auto_unbox = TRUE)
    outputs <- file.path(directory, c("overall.csv", "pairwise.csv", "summary.json"))
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      shQuote(c("r/run_analysis.R", metadata_input,
                                features_input, request, outputs)))
    demo_results <- read.csv(outputs[1])
    stopifnot(status == 0, all(file.exists(outputs)), nrow(demo_results) > 0,
              all(is.finite(demo_results$p_value)))
    unlink(directory, recursive = TRUE)
} else {
    message("SKIP: full R integration checks require lme4, lmerTest, emmeans, and jsonlite.")
}
message("Regression checks passed.")
