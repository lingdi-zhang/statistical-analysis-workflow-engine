# ==========================================================
# Streamlit request builder
#
# Output request structures:
#
# Cross-sectional:
# {
#     "analysis_type": "cross_sectional",
#     "primary_predictors": [predictor],
#     "covariates": covariates,
#     "variable_types": variable_types,
#     "ordered_levels": ordered_levels,
#     "subset": {
#         "enabled": use_subset == "Yes",
#         "variable": subset_variable,
#         "value": subset_value
#     }
# }
#
# Longitudinal:
# {
#     "analysis_type": "longitudinal",
#     "analysis_goal": analysis_goal,
#     "subject_id": subject_id,
#     "time": time_variable,
#     "primary_predictors": [predictor],
#     "covariates": covariates,
#     "variable_types": variable_types,
#     "ordered_levels": ordered_levels
# }
# ==========================================================

import streamlit as st


# ==========================================================
# Helpers
# ==========================================================

def safe_unique_values(series):
    """
    Return unique non-missing values as strings.
    """
    return (
        series
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )


def get_metadata_columns(metadata):
    """
    Return metadata column names as a list.
    """
    if metadata is None:
        return []

    return metadata.columns.tolist()


def build_variable_type_ui(
    metadata,
    variables,
    prefix=""
):
    """
    Dynamically build datatype controls for selected variables.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Uploaded metadata.

    variables : list
        Variables for which datatype controls should be shown.

    prefix : str
        Prefix used to ensure unique Streamlit widget keys.

    Returns
    -------
    variable_types : dict
        Example:
        {
            "FA": "categorical",
            "Age": "numeric",
            "Race": "auto"
        }

    ordered_levels : dict
        Example:
        {
            "Severity": [
                "mild",
                "moderate",
                "severe"
            ]
        }
    """

    variable_types = {}
    ordered_levels = {}

    variables = list(
        dict.fromkeys(
            variable
            for variable in variables
            if variable is not None
            and variable != ""
        )
    )

    if not variables:
        st.caption(
            "Select a primary predictor and/or covariates "
            "to define their variable types."
        )

        return (
            variable_types,
            ordered_levels
        )

    for variable in variables:

        col1, col2 = st.columns(
            [2, 3]
        )

        with col1:
            st.write(
                f"**{variable}**"
            )

        with col2:
            type_choice = st.selectbox(
                "Type",
                options=[
                    "Auto",
                    "Numeric",
                    "Categorical",
                    "Ordered categorical"
                ],
                key=f"{prefix}_type_{variable}",
                label_visibility="collapsed"
            )

        resolved_choice = (
            type_choice
            .lower()
            .replace(" ", "_")
        )

        variable_types[
            variable
        ] = resolved_choice

        if type_choice == "Ordered categorical":

            if variable in metadata.columns:

                levels = safe_unique_values(
                    metadata[
                        variable
                    ]
                )

                if levels:

                    st.caption(
                        f"Observed values for {variable}: "
                        + ", ".join(levels)
                    )

                    ordered_levels[
                        variable
                    ] = st.multiselect(
                        f"Order for {variable}",
                        options=levels,
                        default=levels,
                        help=(
                            "Arrange the selected levels "
                            "from lowest to highest."
                        ),
                        key=(
                            f"{prefix}_order_"
                            f"{variable}"
                        )
                    )

                else:
                    st.caption(
                        f"No non-missing values "
                        f"were found for "
                        f"'{variable}'."
                    )

                    ordered_levels[
                        variable
                    ] = []

            else:
                st.caption(
                    f"Column '{variable}' "
                    "was not found in metadata."
                )

    return (
        variable_types,
        ordered_levels
    )


# ==========================================================
# Cross-sectional request builder
# ==========================================================

def build_cross_sectional_request(
    metadata
):
    """
    Render the cross-sectional request UI.

    Returns
    -------
    tuple
        Current request and whether "Inspect Data →" was clicked.
    """

    metadata_columns = (
        get_metadata_columns(
            metadata
        )
    )

    st.header(
        "3. Cross-sectional Analysis"
    )

    st.info(
        "Cross-sectional analysis can use the full dataset "
        "or a user-defined subset."
    )

    # ======================================================
    # Variables
    # ======================================================

    st.subheader(
        "Define Variables"
    )

    predictor = st.selectbox(
        "Primary predictor",
        options=metadata_columns,
        index=None,
        placeholder=(
            "Select the primary predictor"
        ),
        key="cross_predictor"
    )

    covariate_options = [
        column
        for column in metadata_columns
        if column != predictor
    ]

    covariates = st.multiselect(
        "Covariates",
        options=covariate_options,
        placeholder=(
            "Select optional covariates"
        ),
        key="cross_covariates"
    )

    # ======================================================
    # Optional subset
    # ======================================================

    st.subheader(
        "Optional Subset"
    )

    use_subset = st.radio(
        "Subset the uploaded dataset?",
        options=[
            "No",
            "Yes"
        ],
        horizontal=True,
        key="cross_subset_choice"
    )

    subset_variable = None
    subset_value = None

    if use_subset == "Yes":

        subset_options = [
            column
            for column in metadata_columns
            if column != predictor
        ]

        subset_variable = (
            st.selectbox(
                "Subset column",
                options=subset_options,
                index=None,
                placeholder=(
                    "Select a column"
                ),
                key=(
                    "cross_subset_variable"
                )
            )
        )

        if (
            subset_variable
            is not None
            and subset_variable
            in metadata.columns
        ):

            possible_values = (
                metadata[
                    subset_variable
                ]
                .dropna()
                .unique()
                .tolist()
            )

            if possible_values:

                subset_value = (
                    st.selectbox(
                        "Value to keep",
                        options=possible_values,
                        key=(
                            "cross_subset_value"
                        )
                    )
                )

            else:
                st.caption(
                    "No non-missing values "
                    f"were found in "
                    f"'{subset_variable}'."
                )

    # ======================================================
    # Variable types
    # ======================================================

    st.subheader(
        "Variable Types"
    )

    st.write(
        "Specify the type of each selected variable, "
        "or leave it as Auto to detect the type "
        "during data inspection."
    )

    variables_to_define = (
        [predictor]
        + covariates
    )

    (
        variable_types,
        ordered_levels
    ) = build_variable_type_ui(
        metadata=metadata,
        variables=variables_to_define,
        prefix="cross"
    )

    st.divider()

    # ======================================================
    # Build request
    # ======================================================

    inspect_clicked = st.button(
        "Inspect Data →",
        type="primary",
        key="cross_inspect"
    )

    request = {
        "analysis_type":
            "cross_sectional",

        "primary_predictors":
            (
                [predictor]
                if predictor
                is not None
                else []
            ),

        "covariates":
            covariates,

        "variable_types":
            variable_types,

        "ordered_levels":
            ordered_levels,

        "subset": {
            "enabled":
                use_subset
                == "Yes",

            "variable":
                subset_variable,

            "value":
                subset_value
        }
    }

    return request, inspect_clicked




# ==========================================================
# Longitudinal request builder
# ==========================================================

def build_longitudinal_request(
    metadata
):
    """
    Render the longitudinal request UI.

    Returns
    -------
    tuple
        Current request and whether "Inspect Data →" was clicked.
    """

    metadata_columns = (
        get_metadata_columns(
            metadata
        )
    )

    st.header(
        "3. Longitudinal Analysis"
    )

    st.info(
        "Longitudinal analysis requires a subject identifier "
        "and a time variable."
    )

    # ======================================================
    # Longitudinal structure
    # ======================================================

    st.subheader(
        "Define Longitudinal Structure"
    )

    subject_id = st.selectbox(
        "Subject ID",
        options=metadata_columns,
        index=None,
        placeholder=(
            "Select the subject ID column"
        ),
        key="long_subject"
    )

    time_options = [
        column
        for column in metadata_columns
        if column != subject_id
    ]

    time_variable = st.selectbox(
        "Time variable",
        options=time_options,
        index=None,
        placeholder=(
            "Select the time variable"
        ),
        key="long_time"
    )

    # ======================================================
    # Predictor / covariates
    # ======================================================

    st.subheader(
        "Define Variables"
    )

    predictor_options = [
        column
        for column in metadata_columns
        if column
        not in {
            subject_id,
            time_variable
        }
    ]

    predictor = st.selectbox(
        "Primary predictor",
        options=predictor_options,
        index=None,
        placeholder=(
            "Select the primary predictor"
        ),
        key="long_predictor"
    )

    covariate_options = [
        column
        for column in metadata_columns
        if column
        not in {
            subject_id,
            time_variable,
            predictor
        }
    ]

    covariates = st.multiselect(
        "Covariates",
        options=covariate_options,
        placeholder=(
            "Select optional covariates"
        ),
        key="long_covariates"
    )

    # ======================================================
    # Analysis goal
    # ======================================================

    st.subheader(
        "Analysis Goal"
    )

    analysis_goal_label = (
        st.selectbox(
            "What do you want to test?",
            options=[
                (
                    "Do groups change differently over time? "
                    "- trajectory difference"
                ),
                (
                    "Is the predictor associated with the feature "
                    "overall? - overall predictor association"
                ),
                (
                    "Does the feature change over time? "
                    "- time effect"
                )
            ],
            key="long_goal"
        )
    )

    goal_mapping = {
        (
            "Do groups change differently over time? "
            "- trajectory difference"
        ):
            "trajectory_difference",

        (
            "Is the predictor associated with the feature "
            "overall? - overall predictor association"
        ):
            "overall_predictor_association",

        (
            "Does the feature change over time? "
            "- time effect"
        ):
            "time_effect"
    }

    analysis_goal = (
        goal_mapping[
            analysis_goal_label
        ]
    )

    # ======================================================
    # Variable types
    # ======================================================

    st.subheader(
        "Variable Types"
    )

    st.write(
        "Specify the type of each selected variable, "
        "or leave it as Auto to detect the type "
        "during data inspection."
    )

    variables_to_define = (
        [predictor]
        + covariates
    )

    (
        variable_types,
        ordered_levels
    ) = build_variable_type_ui(
        metadata=metadata,
        variables=variables_to_define,
        prefix="long"
    )

    st.subheader("Random Effects")
    random_effects = st.radio(
        "Subject-specific effects",
        options=["Random intercept only", "Random intercept + time slope"],
        index=1,
        key="long_random_effects",
        help="Random intercepts allow subject baselines to differ. Time slopes also allow subject trajectories to differ and require more usable observations."
    )
    random_slope = random_effects == "Random intercept + time slope"

    st.divider()

    # ======================================================
    # Build request
    # ======================================================

    inspect_clicked = st.button(
        "Inspect Data →",
        type="primary",
        key="long_inspect"
    )

    request = {
        "analysis_type":
            "longitudinal",

        "random_slope": random_slope,

        "analysis_goal":
            analysis_goal,

        "subject_id":
            subject_id,

        "time":
            time_variable,

        "primary_predictors":
            (
                [predictor]
                if predictor
                is not None
                else []
            ),

        "covariates":
            covariates,

        "variable_types":
            variable_types,

        "ordered_levels":
            ordered_levels
    }

    return request, inspect_clicked




# ==========================================================
# Main request-builder entry point
# ==========================================================

def render_request_builder(
    metadata
):
    """
    Render analysis-type selection and the appropriate
    request UI.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Uploaded sample-level metadata.

    Returns
    -------
    tuple
        Current analysis request and whether inspection was requested.

    Notes
    -----
    This module does NOT:

    - validate requested columns,
    - infer variable types,
    - inspect missingness,
    - inspect longitudinal structure,
    - match metadata to the feature matrix,
    - fit models,
    - display DataInspector results.

    Those responsibilities belong downstream.
    """

    st.header(
        "2. Analysis Type"
    )

    analysis_type = st.radio(
        "Choose analysis type",
        options=[
            "Cross-sectional",
            "Longitudinal"
        ],
        horizontal=True,
        key="analysis_type"
    )

    st.divider()

    if (
        analysis_type
        == "Cross-sectional"
    ):

        return (
            build_cross_sectional_request(
                metadata=metadata
            )
        )

    return (
        build_longitudinal_request(
            metadata=metadata
        )
    )
