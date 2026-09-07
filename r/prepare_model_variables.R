prepare_model_variables <- function(
    dat,
    variable_types,
    ordered_levels = list()
) {

    for (variable in names(variable_types)) {

        variable_type <- variable_types[[variable]]

        if (variable_type == "numeric") {

            dat[[variable]] <- as.numeric(
                dat[[variable]]
            )

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
