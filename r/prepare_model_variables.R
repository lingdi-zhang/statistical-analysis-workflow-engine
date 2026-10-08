prepare_model_variables <- function(
    dat,
    variable_types,
    ordered_levels = list()
) {

    for (variable in names(variable_types)) {

        variable_type <- variable_types[[variable]]

        if (variable_type == "numeric") {

            values <- dat[[variable]]
            if (is.factor(values)) values <- as.character(values)
            converted <- suppressWarnings(as.numeric(values))
            if (any(!is.na(values) & is.na(converted))) {
                stop(paste0("Numeric variable '", variable, "' contains nonnumeric non-missing values."))
            }
            if (any(is.infinite(converted))) {
                stop(paste0("Numeric variable '", variable, "' contains infinite values."))
            }
            dat[[variable]] <- converted

        } else if (variable_type == "categorical") {

            dat[[variable]] <- factor(
                dat[[variable]]
            )

        } else if (
            variable_type == "ordered_categorical"
        ) {

            dat[[variable]] <- factor(
                dat[[variable]],
                levels = ordered_levels[[variable]],
                ordered = TRUE
            )

        } else {

            stop(
                paste0(
                    "Unsupported variable type for ",
                    variable,
                    ": ",
                    variable_type
                )
            )
        }
    }

    return(dat)
}
