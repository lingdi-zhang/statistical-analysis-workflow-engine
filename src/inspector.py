import pandas as pd


class DataInspector:

    def __init__(
        self,
        metadata,
        data_matrix,
        request
    ):
        self.metadata = (
            metadata.copy()
            if metadata is not None
            else None
        )

        self.data_matrix = (
            data_matrix.copy()
            if data_matrix is not None
            else None
        )

        self.request = (
            request.copy()
            if request is not None
            else {}
        )

        self.warnings = []
        self.errors = []

    # =========================================================
    # Public entry point
    # =========================================================

    def inspect(self):

        report = {
            "status": "success",
            "analysis_type":
                self.request.get(
                    "analysis_type"
                ),
            "variables": {},
            "longitudinal": None,
            "sample_matching": None,
            "warnings": [],
            "errors": []
        }

        # -----------------------------------------------------
        # 1. Validate request and metadata
        # -----------------------------------------------------

        self._check_metadata()

        if self.errors:
            return self._finish_report(
                report
            )

        # -----------------------------------------------------
        # 2. Match metadata and data_matrix by SampleID
        # -----------------------------------------------------

        if self.data_matrix is not None:

            report["sample_matching"] = (
                self._match_samples()
            )

            if self.errors:
                return self._finish_report(
                    report
                )

        # -----------------------------------------------------
        # 3. Optional cross-sectional subset
        # -----------------------------------------------------

        self._apply_cross_sectional_subset()

        if self.errors:
            return self._finish_report(
                report
            )

        # -----------------------------------------------------
        # 4. Variable inspection
        # -----------------------------------------------------

        variables = (
            self._get_variables_to_inspect()
        )

        for variable in variables:

            report["variables"][
                variable
            ] = (
                self._inspect_variable(
                    variable
                )
            )

        # -----------------------------------------------------
        # 5. Update request with resolved variable types
        # -----------------------------------------------------

        self._update_variable_types(
            report["variables"]
        )

        # -----------------------------------------------------
        # 6. Ordered categorical validation
        # -----------------------------------------------------

        self._validate_ordered_levels(
            report["variables"]
        )

        if self.errors:
            return self._finish_report(report)

        self._validate_feature_model_inputs()

        if self.errors:
            return self._finish_report(report)

        # -----------------------------------------------------
        # 7. Longitudinal-specific inspection
        # -----------------------------------------------------

        if (
            self.request.get(
                "analysis_type"
            )
            == "longitudinal"
        ):

            report["longitudinal"] = (
                self._inspect_longitudinal_structure()
            )

        # -----------------------------------------------------
        # Finish
        # -----------------------------------------------------

        return self._finish_report(
            report
        )

    # =========================================================
    # Request / metadata validation
    # =========================================================

    def _check_metadata(self):

        # -----------------------------------------------------
        # Metadata itself
        # -----------------------------------------------------

        if self.metadata is None:

            self.errors.append(
                "Metadata was not provided."
            )

            return

        if self.metadata.empty:

            self.errors.append(
                "Metadata is empty."
            )

            return

        # -----------------------------------------------------
        # Require SampleID in first column
        # -----------------------------------------------------

        if (
            self.metadata.columns[0]
            != "SampleID"
        ):

            self.errors.append(
                "The first column of metadata "
                "must be 'SampleID'."
            )

        # -----------------------------------------------------
        # Analysis type
        # -----------------------------------------------------

        analysis_type = (
            self.request.get(
                "analysis_type"
            )
        )

        if analysis_type not in [
            "cross_sectional",
            "longitudinal"
        ]:

            self.errors.append(
                "Analysis type must be either "
                "'cross_sectional' or 'longitudinal'."
            )

            return

        # -----------------------------------------------------
        # Primary predictor
        # -----------------------------------------------------

        predictors = (
            self.request.get(
                "primary_predictors",
                []
            )
        )

        predictors = [
            x
            for x in predictors
            if x is not None
            and x != ""
        ]

        if len(predictors) == 0:

            self.errors.append(
                "Primary predictor was not provided."
            )

        # -----------------------------------------------------
        # Longitudinal required fields
        # -----------------------------------------------------

        if analysis_type == "longitudinal":

            subject_id = (
                self.request.get(
                    "subject_id"
                )
            )

            time_variable = (
                self.request.get(
                    "time"
                )
            )

            if not subject_id:

                self.errors.append(
                    "Subject ID was not provided."
                )

            if not time_variable:

                self.errors.append(
                    "Time variable was not provided."
                )

            analysis_goal = (
                self.request.get(
                    "analysis_goal"
                )
            )

            valid_goals = [
                "trajectory_difference",
                "overall_predictor_association",
                "time_effect"
            ]

            if analysis_goal not in valid_goals:

                self.errors.append(
                    "Longitudinal analysis goal "
                    "is missing or invalid."
                )

        # -----------------------------------------------------
        # Cross-sectional subset request
        # -----------------------------------------------------

        if analysis_type == "cross_sectional":

            subset = (
                self.request.get(
                    "subset",
                    {}
                )
            )

            if subset.get(
                "enabled",
                False
            ):

                if not subset.get(
                    "variable"
                ):

                    self.errors.append(
                        "Subset variable was not provided."
                    )

                if (
                    subset.get(
                        "value"
                    )
                    is None
                ):

                    self.errors.append(
                        "Subset value was not provided."
                    )

        if self.errors:
            return

        # -----------------------------------------------------
        # Required metadata columns
        # -----------------------------------------------------

        required = []

        required += predictors

        required += (
            self.request.get(
                "covariates",
                []
            )
        )

        if analysis_type == "longitudinal":

            required += [
                self.request.get(
                    "subject_id"
                ),
                self.request.get(
                    "time"
                )
            ]

        if analysis_type == "cross_sectional":

            subset = (
                self.request.get(
                    "subset",
                    {}
                )
            )

            if subset.get(
                "enabled",
                False
            ):

                required.append(
                    subset.get(
                        "variable"
                    )
                )

        required = list(
            dict.fromkeys(
                x
                for x in required
                if x is not None
                and x != ""
            )
        )

        missing = [
            x
            for x in required
            if x not in self.metadata.columns
        ]

        if missing:

            self.errors.append(
                "Missing metadata columns: "
                + ", ".join(
                    missing
                )
            )

    # =========================================================
    # Sample matching
    # =========================================================

    def _match_samples(self):

        """
        Match metadata and data_matrix using SampleID.

        Assumptions
        -----------
        - First column of metadata is SampleID.
        - First column of data_matrix is SampleID.
        - Only shared samples are retained.
        - data_matrix is reordered to match metadata.
        """

        # -----------------------------------------------------
        # Require SampleID as first column of data_matrix
        # -----------------------------------------------------

        if (
            self.data_matrix.columns[0]
            != "SampleID"
        ):

            self.errors.append(
                "The first column of the feature matrix "
                "must be 'SampleID'."
            )

            return None

        # -----------------------------------------------------
        # Missing SampleID
        # -----------------------------------------------------

        if (
            self.metadata[
                "SampleID"
            ]
            .isna()
            .any()
        ):

            self.errors.append(
                "Metadata contains missing SampleID values."
            )

        if (
            self.data_matrix[
                "SampleID"
            ]
            .isna()
            .any()
        ):

            self.errors.append(
                "Feature matrix contains missing SampleID values."
            )

        if self.errors:
            return None

        # -----------------------------------------------------
        # Convert SampleID to string
        # -----------------------------------------------------

        self.metadata[
            "SampleID"
        ] = (
            self.metadata[
                "SampleID"
            ]
            .astype(str)
        )

        self.data_matrix[
            "SampleID"
        ] = (
            self.data_matrix[
                "SampleID"
            ]
            .astype(str)
        )

        # -----------------------------------------------------
        # Duplicate SampleID
        # -----------------------------------------------------

        metadata_duplicates = int(
            self.metadata[
                "SampleID"
            ]
            .duplicated()
            .sum()
        )

        data_duplicates = int(
            self.data_matrix[
                "SampleID"
            ]
            .duplicated()
            .sum()
        )

        if metadata_duplicates > 0:

            self.errors.append(
                f"Metadata contains "
                f"{metadata_duplicates} duplicated "
                "SampleID values."
            )

        if data_duplicates > 0:

            self.errors.append(
                f"Feature matrix contains "
                f"{data_duplicates} duplicated "
                "SampleID values."
            )

        if self.errors:
            return None

        # -----------------------------------------------------
        # Shared samples
        # -----------------------------------------------------

        metadata_samples = set(
            self.metadata[
                "SampleID"
            ]
        )

        data_samples = set(
            self.data_matrix[
                "SampleID"
            ]
        )

        shared_samples = (
            metadata_samples
            & data_samples
        )

        if len(shared_samples) == 0:

            self.errors.append(
                "Metadata and feature matrix "
                "have no shared SampleID values."
            )

            return None

        # -----------------------------------------------------
        # Warnings if unmatched samples are removed
        # -----------------------------------------------------

        if len(
            metadata_samples
            - data_samples
        ) > 0:

            self.warnings.append(
                "Some metadata samples were not found "
                "in the feature matrix and were removed."
            )

        if len(
            data_samples
            - metadata_samples
        ) > 0:

            self.warnings.append(
                "Some feature-matrix samples were not found "
                "in metadata and were removed."
            )

        # -----------------------------------------------------
        # Keep shared samples in metadata
        # Preserve metadata order
        # -----------------------------------------------------

        self.metadata = (
            self.metadata[
                self.metadata[
                    "SampleID"
                ].isin(
                    shared_samples
                )
            ]
            .copy()
            .reset_index(
                drop=True
            )
        )

        # -----------------------------------------------------
        # Keep shared samples in data_matrix
        # Align to metadata order
        # -----------------------------------------------------

        self.data_matrix = (
            self.data_matrix[
                self.data_matrix[
                    "SampleID"
                ].isin(
                    shared_samples
                )
            ]
            .set_index(
                "SampleID"
            )
            .loc[
                self.metadata[
                    "SampleID"
                ]
            ]
            .reset_index()
        )

        return {
            "n_shared_samples":
                len(
                    shared_samples
                )
        }

    # =========================================================
    # Cross-sectional subsetting
    # =========================================================

    def _apply_cross_sectional_subset(self):

        if (
            self.request.get(
                "analysis_type"
            )
            != "cross_sectional"
        ):

            return

        subset = (
            self.request.get(
                "subset",
                {}
            )
        )

        if not subset.get(
            "enabled",
            False
        ):

            return

        variable = (
            subset.get(
                "variable"
            )
        )

        value = (
            subset.get(
                "value"
            )
        )

        if variable not in self.metadata.columns:

            self.errors.append(
                f"Subset variable '{variable}' "
                "was not found."
            )

            return

        # -----------------------------------------------------
        # SampleIDs to retain
        # -----------------------------------------------------

        keep_samples = (
            self.metadata.loc[
                self.metadata[
                    variable
                ] == value,
                "SampleID"
            ]
            .tolist()
        )

        # -----------------------------------------------------
        # Subset metadata
        # -----------------------------------------------------

        self.metadata = (
            self.metadata[
                self.metadata[
                    "SampleID"
                ].isin(
                    keep_samples
                )
            ]
            .copy()
            .reset_index(
                drop=True
            )
        )

        # -----------------------------------------------------
        # Subset data_matrix
        # -----------------------------------------------------

        if self.data_matrix is not None:

            self.data_matrix = (
                self.data_matrix[
                    self.data_matrix[
                        "SampleID"
                    ].isin(
                        keep_samples
                    )
                ]
                .set_index(
                    "SampleID"
                )
                .loc[
                    self.metadata[
                        "SampleID"
                    ]
                ]
                .reset_index()
            )

        after = len(
            self.metadata
        )

        if after == 0:

            self.errors.append(
                f"Subsetting {variable} = {value} "
                "returned zero samples."
            )

        elif after < 10:

            self.warnings.append(
                f"Only {after} samples remain after "
                f"subsetting {variable} = {value}."
            )

    # =========================================================
    # Which variables should be inspected
    # =========================================================

    def _get_variables_to_inspect(self):

        variables = []

        variables += (
            self.request.get(
                "primary_predictors",
                []
            )
        )

        variables += (
            self.request.get(
                "covariates",
                []
            )
        )

        if (
            self.request.get(
                "analysis_type"
            )
            == "longitudinal"
        ):

            variables.append(
                self.request.get(
                    "time"
                )
            )

        variables = list(
            dict.fromkeys(
                x
                for x in variables
                if x is not None
                and x != ""
            )
        )

        return variables

    # =========================================================
    # Resolve requested / inferred type
    # =========================================================

    def _resolve_variable_type(
        self,
        variable
    ):

        # Time is measured numerically even when only two visits are observed.
        if (self.request.get("analysis_type") == "longitudinal"
                and variable == self.request.get("time")):
            return "numeric"

        requested_types = (
            self.request.get(
                "variable_types",
                {}
            )
        )

        requested = (
            requested_types.get(
                variable,
                "auto"
            )
        )

        if requested != "auto":

            return requested

        return self._infer_variable_type(
            self.metadata[
                variable
            ]
        )

    # =========================================================
    # Infer type
    # =========================================================

    def _infer_variable_type(
        self,
        series
    ):

        x = (
            series
            .dropna()
        )

        if len(x) == 0:

            return "unknown"

        if pd.api.types.is_bool_dtype(
            x
        ):

            return "categorical"

        if (pd.api.types.is_numeric_dtype(x)
                or pd.to_numeric(x, errors="coerce").notna().all()):

            n_unique = (
                x.nunique()
            )

            if n_unique == 2:

                return "categorical"

            return "numeric"

        return "categorical"

    # =========================================================
    # Inspect one variable
    # =========================================================

    def _inspect_variable(
        self,
        variable
    ):

        series = (
            self.metadata[
                variable
            ]
        )

        resolved_type = (
            self._resolve_variable_type(
                variable
            )
        )

        # Auto detection can classify a two-valued numeric column as categorical;
        # that must not turn numeric infinity into an acceptable category.
        requested_type = self.request.get("variable_types", {}).get(variable, "auto")
        if resolved_type != "numeric" and (
                pd.api.types.is_numeric_dtype(series) or requested_type == "auto"):
            numeric_values = pd.to_numeric(series, errors="coerce")
            infinite_n = int(numeric_values.isin([float("inf"), float("-inf")]).sum())
            if infinite_n:
                self.errors.append(f"{variable} contains {infinite_n} infinite numeric values.")

        missing_n = int(
            series
            .isna()
            .sum()
        )

        missing_fraction = float(
            series
            .isna()
            .mean()
        )

        unique_n = int(
            series
            .nunique(
                dropna=True
            )
        )

        out = {
            "type":
                resolved_type,

            "n_unique":
                unique_n,

            "missing_n":
                missing_n,

            "missing_fraction":
                missing_fraction
        }

        # -----------------------------------------------------
        # Unknown
        # -----------------------------------------------------

        if resolved_type == "unknown":

            self.errors.append(
                f"{variable} contains no observed "
                "non-missing values."
            )

            return out

        # -----------------------------------------------------
        # Numeric
        # -----------------------------------------------------

        if resolved_type == "numeric":

            numeric = (
                pd.to_numeric(
                    series,
                    errors="coerce"
                )
            )

            infinite_values = numeric.isin([float("inf"), float("-inf")])
            infinite_n = int(infinite_values.sum())
            if infinite_n:
                self.errors.append(
                    f"{variable} contains {infinite_n} infinite numeric values."
                )
                # Calculate descriptive statistics on finite values only.
                numeric = numeric.mask(infinite_values)

            out.update({
                "mean":
                    self._safe_float(
                        numeric.mean()
                    ),

                "sd":
                    self._safe_float(
                        numeric.std()
                    ),

                "min":
                    self._safe_float(
                        numeric.min()
                    ),

                "median":
                    self._safe_float(
                        numeric.median()
                    ),

                "max":
                    self._safe_float(
                        numeric.max()
                    )
            })

            invalid_numeric = int(
                (
                    series.notna()
                    & numeric.isna()
                    & ~infinite_values
                )
                .sum()
            )

            if invalid_numeric > 0:

                self.errors.append(
                    f"{variable} was defined as numeric "
                    f"but {invalid_numeric} non-missing values "
                    "could not be converted to numeric."
                )

            self.metadata[variable] = numeric

            if unique_n < 2:

                self.errors.append(
                    f"{variable} has fewer than "
                    "2 observed values."
                )

        # -----------------------------------------------------
        # Categorical / ordered categorical
        # -----------------------------------------------------

        elif resolved_type in [
            "categorical",
            "ordered_categorical"
        ]:

            counts = (
                series
                .astype("string")
                .value_counts(
                    dropna=True
                )
            )

            out[
                "levels"
            ] = (
                counts
                .index
                .tolist()
            )

            out[
                "counts"
            ] = {
                str(k):
                    int(v)

                for k, v
                in counts.items()
            }

            if unique_n < 2:

                self.errors.append(
                    f"{variable} has fewer than "
                    "2 observed levels."
                )

            if unique_n == 2:

                props = (
                    counts
                    / counts.sum()
                )

                min_prop = float(
                    props.min()
                )

                out[
                    "smallest_group_fraction"
                ] = min_prop

                if min_prop < 0.10:

                    self.warnings.append(
                        f"{variable} is highly imbalanced; "
                        f"the smallest group contains "
                        f"{min_prop:.1%} of samples."
                    )

        # -----------------------------------------------------
        # Unsupported type
        # -----------------------------------------------------

        else:

            self.errors.append(
                f"Unsupported variable type "
                f"'{resolved_type}' for {variable}."
            )

        # -----------------------------------------------------
        # Missingness warning
        # -----------------------------------------------------

        if missing_fraction >= 0.20:

            self.warnings.append(
                f"{variable} has "
                f"{missing_fraction:.1%} missing values."
            )

        return out

    # =========================================================
    # Update request with resolved variable types
    # =========================================================

    def _update_variable_types(
        self,
        variable_reports
    ):

        resolved_types = {}

        for variable, info in (
            variable_reports.items()
        ):

            resolved_types[
                variable
            ] = (
                info.get(
                    "type"
                )
            )

        self.request[
            "variable_types"
        ] = (
            resolved_types
        )

    # =========================================================
    # Ordered level validation
    # =========================================================

    def _validate_ordered_levels(
        self,
        variable_reports
    ):

        requested_orders = (
            self.request.get(
                "ordered_levels",
                {}
            )
        )

        for variable, info in (
            variable_reports.items()
        ):

            if (
                info.get(
                    "type"
                )
                != "ordered_categorical"
            ):

                continue

            requested_order = (
                requested_orders.get(
                    variable
                )
            )

            if not requested_order:

                self.errors.append(
                    f"{variable} was defined as ordered "
                    "categorical but no level order was provided."
                )

                continue

            observed = set(
                str(x)
                for x in (
                    self.metadata[
                        variable
                    ]
                    .dropna()
                    .unique()
                )
            )

            requested = set(
                str(x)
                for x in requested_order
            )

            missing_from_order = (
                observed
                - requested
            )

            extra_in_order = (
                requested
                - observed
            )

            if missing_from_order:

                self.errors.append(
                    f"Ordered levels for {variable} "
                    "are missing observed values: "
                    + ", ".join(
                        sorted(
                            missing_from_order
                        )
                    )
                )

            if extra_in_order:

                self.warnings.append(
                    f"Ordered levels for {variable} "
                    "contain values not observed in the data: "
                    + ", ".join(
                        sorted(
                            extra_in_order
                        )
                    )
                )

            if (
                len(
                    requested_order
                )
                != len(
                    set(
                        str(x)
                        for x in requested_order
                    )
                )
            ):

                self.errors.append(
                    f"Ordered levels for {variable} "
                    "contain duplicate values."
                )

            info[
                "ordered_levels"
            ] = (
                requested_order
            )

    # =========================================================
    # Longitudinal inspection
    # =========================================================

    def _inspect_longitudinal_structure(
        self
    ):

        subject = (
            self.request[
                "subject_id"
            ]
        )

        time = (
            self.request[
                "time"
            ]
        )

        n_subjects = int(
            self.metadata[
                subject
            ]
            .nunique(
                dropna=True
            )
        )

        counts = (
            self.metadata
            .dropna(
                subset=[
                    subject
                ]
            )
            .groupby(
                subject
            )
            .size()
        )

        duplicate_rows = int(
            self.metadata
            .duplicated(
                subset=[
                    subject,
                    time
                ]
            )
            .sum()
        )

        timepoints = (
            self.metadata[
                time
            ]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )

        if len(counts) == 0:

            median_obs = None
            min_obs = None
            max_obs = None

        else:

            median_obs = (
                self._safe_float(
                    counts.median()
                )
            )

            min_obs = int(
                counts.min()
            )

            max_obs = int(
                counts.max()
            )

        report = {
            "n_subjects":
                n_subjects,

            "n_timepoints":
                len(
                    timepoints
                ),

            "timepoints":
                timepoints,

            "median_observations_per_subject":
                median_obs,

            "min_observations_per_subject":
                min_obs,

            "max_observations_per_subject":
                max_obs,

            "duplicate_subject_time_rows":
                duplicate_rows
        }

        if n_subjects < 2:

            self.errors.append(
                "Longitudinal analysis requires "
                "at least two subjects."
            )

        if len(timepoints) < 2:

            self.errors.append(
                "Longitudinal analysis requires "
                "at least two observed time points."
            )

        if duplicate_rows > 0:

            self.warnings.append(
                f"{duplicate_rows} duplicated "
                f"{subject} + {time} combinations were found."
            )

        return report

    # =========================================================
    # Return finalized inputs for downstream analysis
    # =========================================================

    def get_analysis_inputs(
        self
    ):

        return {
            "metadata":
                self.metadata,

            "data_matrix":
                self.data_matrix,

            "request":
                self.request
        }

    def _validate_feature_model_inputs(self):
        """Reject collisions and category loss before any models are run."""
        if self.data_matrix is None:
            return

        features = [name for name in self.data_matrix.columns if name != "SampleID"]
        for feature in features:
            if feature in self.metadata.columns:
                self.errors.append(
                    f"Feature '{feature}' conflicts with metadata column '{feature}'."
                )

        predictor = self.request["primary_predictors"][0]
        types = self.request["variable_types"]
        categorical = types.get(predictor) in ("categorical", "ordered_categorical")

        if not features:
            self.errors.append("Feature matrix contains no outcome columns.")
            return

        required = list(dict.fromkeys(
            [predictor] + self.request.get("covariates", [])
            + ([self.request["time"], self.request["subject_id"]]
               if self.request["analysis_type"] == "longitudinal" else [])
        ))
        model_metadata = self.metadata[required].copy()
        for variable in required:
            if types.get(variable) == "numeric":
                model_metadata[variable] = pd.to_numeric(
                    model_metadata[variable], errors="coerce"
                )

        # Alignment is by SampleID, independent of DataFrame row indexes.
        model_metadata.index = self.metadata["SampleID"]
        expected = set(model_metadata[predictor].dropna().astype(str))
        metadata_complete = model_metadata.notna().all(axis=1)
        outcomes = self.data_matrix.set_index("SampleID")
        for feature in features:
            raw = outcomes[feature]
            outcome = pd.to_numeric(raw, errors="coerce")
            invalid = int((raw.notna() & outcome.isna()).sum())
            infinite = outcome.isin([float("inf"), float("-inf")])
            if invalid:
                self.errors.append(
                    f"Feature '{feature}' contains {invalid} nonnumeric non-missing values."
                )
            if infinite.any():
                self.errors.append(f"Feature '{feature}' contains infinite values.")
            outcome = outcome.mask(infinite)
            usable = metadata_complete & outcome.reindex(model_metadata.index).notna()
            usable_outcome = outcome.reindex(model_metadata.index)[usable]
            if usable_outcome.empty:
                self.errors.append(
                    f"Feature '{feature}' has no usable observations after excluding missing model values."
                )
            elif usable_outcome.nunique() < 2:
                self.errors.append(
                    f"Feature '{feature}' is constant after excluding missing model values."
                )
            if (self.request["analysis_type"] == "longitudinal"
                    and self.request.get("random_slope", True)):
                n_subjects = model_metadata.loc[usable, self.request["subject_id"]].nunique()
                n_observations = int(usable.sum())
                if n_subjects and n_observations <= 2 * n_subjects:
                    self.warnings.append(
                        f"Feature '{feature}': random intercept + time slope is unsupported "
                        f"with {n_observations} usable observations for {n_subjects} subjects "
                        f"({2 * n_subjects} random effects). Select 'Random intercept only' "
                        "and inspect again. The selected model will not be changed automatically."
                    )
            if not categorical:
                continue
            observed = set(model_metadata.loc[usable, predictor].astype(str))
            missing = sorted(expected - observed)
            if missing:
                self.errors.append(
                    f"Feature '{feature}' has no usable observations for predictor "
                    f"'{predictor}' categories: {', '.join(missing)} after excluding "
                    "missing model values."
                )

    # =========================================================
    # Helpers
    # =========================================================

    @staticmethod
    def _safe_float(
        value
    ):

        if pd.isna(
            value
        ):

            return None

        return float(
            value
        )

    # =========================================================
    # Finish report
    # =========================================================

    def _finish_report(
        self,
        report
    ):

        report[
            "warnings"
        ] = (
            self.warnings
        )

        report[
            "errors"
        ] = (
            self.errors
        )

        if self.errors:

            report[
                "status"
            ] = "failed"

        elif self.warnings:

            report[
                "status"
            ] = (
                "passed_with_warnings"
            )

        else:

            report[
                "status"
            ] = "passed"

        return report
