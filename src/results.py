from io import BytesIO

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st


# ============================================================
# Plot download helpers
# ============================================================

def safe_file_name(value):
    value = str(value)
    safe = "".join(ch if ch.isalnum() or ch in ["-", "_", "."] else "_" for ch in value)
    return safe.strip("_") or "plot"


def add_pdf_download(fig, file_name, key):
    pdf_buffer = BytesIO()
    fig.savefig(pdf_buffer, format="pdf", bbox_inches="tight")
    pdf_buffer.seek(0)
    st.download_button(
        label="Download plot as PDF",
        data=pdf_buffer.getvalue(),
        file_name=file_name,
        mime="application/pdf",
        key=key
    )


# ============================================================
# Feature Explorer color palette
#
# Intentionally separate from volcano-plot colors.
# ============================================================

def get_feature_group_colors(n_groups):
    """
    Return distinct colors for Feature Explorer groups.

    This palette is intentionally used only for individual
    feature plots so it remains visually separate from the
    volcano plots.
    """

    # Purple, green, pink, gray, olive, brown.
    # Deliberately avoids blue and orange.
    palette = [
        "#9467BD",
        "#2CA02C",
        "#E377C2",
        "#7F7F7F",
        "#8C9A3C",
        "#8C564B"
    ]

    return [
        palette[i % len(palette)]
        for i in range(n_groups)
    ]


# ============================================================
# Main results display
# ============================================================

def display_analysis_results(
    analysis_output,
    request,
    metadata=None,
    data_matrix=None
):

    results = analysis_output["overall"]

    pairwise = analysis_output.get(
        "pairwise"
    )

    summary = analysis_output.get(
        "summary",
        {}
    )

    n_features = len(
        results
    )

    predictor = (
        request["primary_predictors"][0]
        if request.get("primary_predictors")
        else None
    )

    predictor_type = (
        request
        .get("variable_types", {})
        .get(predictor, "auto")
    )

    analysis_type = (
        request["analysis_type"]
    )

    analysis_goal = (
        request.get(
            "analysis_goal"
        )
    )

    # ========================================================
    # 1. Model summary
    # ========================================================

    st.header(
        "Analysis Results"
    )

    display_model_summary(
        summary=summary,
        request=request
    )

    # ========================================================
    # 2. Main results
    # ========================================================

    if n_features < 5:

        st.subheader(
            "Results"
        )

        display_results_table(
            results
        )

        if (
            pairwise is not None
            and
            len(pairwise) > 0
        ):

            st.subheader(
                "Pairwise Comparisons"
            )

            display_pairwise_table(
                pairwise
            )

    else:

        # ====================================================
        # Cross-sectional
        # ====================================================

        if (
            analysis_type
            == "cross_sectional"
        ):

            multi_category = (
                predictor_type == "categorical"
                and
                "estimate" in results.columns
                and
                results["estimate"].isna().all()
            )

            if multi_category:

                st.subheader(
                    "Pairwise Associations"
                )

                if (
                    pairwise is not None
                    and
                    len(pairwise) > 0
                ):

                    plot_pairwise_volcano(
                        pairwise
                    )

                else:

                    st.info(
                        "No pairwise results are available."
                    )

            else:

                st.subheader(
                    "Feature Associations"
                )

                plot_volcano(
                    results=results,
                    title="Feature Associations"
                )

        # ====================================================
        # Longitudinal
        # ====================================================

        elif (
            analysis_type
            == "longitudinal"
        ):

            if (
                analysis_goal
                == "overall_predictor_association"
            ):

                multi_category = (
                    predictor_type == "categorical"
                    and
                    "estimate" in results.columns
                    and
                    results["estimate"].isna().all()
                )

                if multi_category:

                    st.subheader(
                        "Pairwise Associations"
                    )

                    if (
                        pairwise is not None
                        and
                        len(pairwise) > 0
                    ):

                        plot_pairwise_volcano(
                            pairwise
                        )

                    else:

                        st.info(
                            "No pairwise results are available."
                        )

                else:

                    st.subheader(
                        "Feature Associations"
                    )

                    plot_volcano(
                        results=results,
                        title="Overall Predictor Association"
                    )

            elif (
                analysis_goal
                == "time_effect"
            ):

                st.subheader(
                    "Time Effects"
                )

                plot_volcano(
                    results=results,
                    title="Longitudinal Time Effect"
                )

            elif (
                analysis_goal
                == "trajectory_difference"
            ):

                st.subheader(
                    "Trajectory Differences"
                )

                has_effect_estimate = (
                    "estimate" in results.columns
                    and
                    results["estimate"].notna().any()
                )

                if has_effect_estimate:

                    plot_volcano(
                        results=results,
                        title="Predictor × Time Interaction"
                    )

                else:

                    st.info(
                        "The primary test is an overall "
                        "predictor × time interaction test. "
                        "Because the predictor has multiple "
                        "unordered categories, there is no "
                        "single effect estimate for a "
                        "volcano plot."
                    )

                    display_results_table(
                        results
                    )

            else:

                st.warning(
                    "Unknown longitudinal analysis goal."
                )

    # ========================================================
    # 3. Feature Explorer
    # ========================================================

    if (
        metadata is not None
        and
        data_matrix is not None
    ):

        display_feature_explorer(
            results=results,
            pairwise=pairwise,
            metadata=metadata,
            data_matrix=data_matrix,
            request=request
        )

    # ========================================================
    # 4. Full results tables
    # ========================================================

    if n_features >= 5:

        with st.expander(
            "View full results table"
        ):

            display_results_table(
                results
            )

        if (
            pairwise is not None
            and
            len(pairwise) > 0
        ):

            with st.expander(
                "View pairwise results table"
            ):

                display_pairwise_table(
                    pairwise
                )


# ============================================================
# Model summary
# ============================================================

def display_model_summary(
    summary,
    request
):

    analysis_type = (
        request["analysis_type"]
    )

    if (
        analysis_type
        == "cross_sectional"
    ):

        st.subheader(
            "Model Summary"
        )

        st.write(
            "**Model used:** Linear regression"
        )

        return

    if (
        analysis_type
        == "longitudinal"
    ):

        st.subheader(
            "Model Summary"
        )

        analysis_goal = (
            summary.get(
                "analysis_goal"
            )
            or
            request.get(
                "analysis_goal"
            )
        )

        if (
            analysis_goal
            == "overall_predictor_association"
        ):

            st.write(
                "**Primary test:** "
                "Overall predictor association"
            )

        elif (
            analysis_goal
            == "time_effect"
        ):

            st.write(
                "**Primary test:** "
                "Longitudinal time effect"
            )

        elif (
            analysis_goal
            == "trajectory_difference"
        ):

            st.write(
                "**Primary test:** "
                "Predictor × time interaction"
            )

        final_random_slope = (
            summary.get(
                "final_random_slope"
            )
        )

        if final_random_slope is True:

            model_used = (
                "Linear mixed model with "
                "random intercept + random slope"
            )

        elif final_random_slope is False:

            model_used = (
                "Linear mixed model with "
                "random intercept only"
            )

        else:

            model_used = (
                "Linear mixed model"
            )

        st.write(
            f"**Model used:** {model_used}"
        )

        convergence_fraction = (
            summary.get(
                "convergence_fraction"
            )
        )

        if (
            convergence_fraction
            is not None
        ):

            st.write(
                "**Convergence fraction:** "
                f"{convergence_fraction:.1%}"
            )

        singular_fraction = (
            summary.get(
                "singular_fraction"
            )
        )

        if (
            singular_fraction
            is None
        ):

            singular_fraction = (
                summary.get(
                    "first_pass_singular_fraction"
                )
            )

        if (
            singular_fraction
            is not None
        ):

            st.write(
                "**Random-slope singularity fraction:** "
                f"{singular_fraction:.1%}"
            )

        if summary.get(
            "rerun",
            False
        ):

            st.warning(
                "The initial random-slope analysis showed "
                "widespread singularity. All features were "
                "rerun using a common random-intercept model."
            )

        n_converged = (
            summary.get(
                "n_converged"
            )
        )

        n_failed = (
            summary.get(
                "n_failed"
            )
        )

        if (
            n_converged is not None
            and
            n_failed is not None
        ):

            st.write(
                f"**Features converged:** "
                f"{n_converged}"
            )

            st.write(
                f"**Features failed:** "
                f"{n_failed}"
            )


# ============================================================
# Main results table
# ============================================================

def display_results_table(
    results
):

    preferred_columns = [
        "feature",
        "term",
        "estimate",
        "std_error",
        "p_value",
        "FDR",
        "n_obs",
        "converged",
        "singular",
        "optimizer"
    ]

    columns = [
        column
        for column in preferred_columns
        if column in results.columns
    ]

    extra_columns = [
        column
        for column in results.columns
        if column not in columns
    ]

    display_df = results[
        columns + extra_columns
    ].copy()

    if (
        "FDR"
        in display_df.columns
    ):

        display_df = (
            display_df.sort_values(
                "FDR",
                na_position="last"
            )
        )

    elif (
        "p_value"
        in display_df.columns
    ):

        display_df = (
            display_df.sort_values(
                "p_value",
                na_position="last"
            )
        )

    st.dataframe(
        display_df,
        use_container_width=True
    )


# ============================================================
# Pairwise results table
# ============================================================

def display_pairwise_table(
    pairwise
):

    preferred_columns = [
        "feature",
        "contrast",
        "estimate",
        "SE",
        "p.value",
        "FDR"
    ]

    columns = [
        column
        for column in preferred_columns
        if column in pairwise.columns
    ]

    extra_columns = [
        column
        for column in pairwise.columns
        if column not in columns
    ]

    display_df = pairwise[
        columns + extra_columns
    ].copy()

    if (
        "contrast" in display_df.columns
        and
        "FDR" in display_df.columns
    ):

        display_df = (
            display_df.sort_values(
                [
                    "contrast",
                    "FDR"
                ],
                na_position="last"
            )
        )

    st.dataframe(
        display_df,
        use_container_width=True
    )


# ============================================================
# Label top significant features
# ============================================================

def label_top_features(
    ax,
    plot_data,
    top_n=10,
    rank_by="FDR",
    fdr_cutoff=0.05
):

    if (
        len(plot_data)
        == 0
    ):

        return

    significant_data = (
        plot_data[
            plot_data["FDR"]
            < fdr_cutoff
        ]
        .copy()
    )

    if (
        len(significant_data)
        == 0
    ):

        return

    n_to_label = min(
        top_n,
        len(significant_data)
    )

    if (
        rank_by
        == "FDR"
    ):

        top = (
            significant_data
            .sort_values(
                "FDR",
                ascending=True
            )
            .head(
                n_to_label
            )
        )

    elif (
        rank_by
        == "effect"
    ):

        top = (
            significant_data
            .assign(
                abs_effect=(
                    significant_data[
                        "estimate"
                    ].abs()
                )
            )
            .sort_values(
                "abs_effect",
                ascending=False
            )
            .head(
                n_to_label
            )
        )

    else:

        raise ValueError(
            "rank_by must be "
            "'FDR' or 'effect'."
        )

    for _, row in top.iterrows():

        ax.annotate(
            str(
                row["feature"]
            ),
            (
                row["estimate"],
                row[
                    "minus_log10_FDR"
                ]
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8
        )


# ============================================================
# Volcano plot
# ============================================================

def plot_volcano(
    results,
    fdr_cutoff=0.05,
    title="Feature Associations"
):
    required = {"estimate", "FDR", "feature"}
    if not required.issubset(results.columns):
        st.info("The result table does not contain the columns required for a volcano plot.")
        return

    plot_data = results.copy()
    plot_data["estimate"] = pd.to_numeric(plot_data["estimate"], errors="coerce")
    plot_data["FDR"] = pd.to_numeric(plot_data["FDR"], errors="coerce")
    plot_data = plot_data[plot_data["estimate"].notna() & plot_data["FDR"].notna()].copy()
    if len(plot_data) == 0:
        st.info("No features have both an estimate and an FDR value.")
        return

    minimum_fdr = np.nextafter(0, 1)
    plot_data["FDR_plot"] = plot_data["FDR"].clip(lower=minimum_fdr)
    plot_data["minus_log10_FDR"] = -np.log10(plot_data["FDR_plot"])
    plot_data["significant"] = plot_data["FDR"] < fdr_cutoff

    rank_option = st.radio(
        "Label significant features by",
        ["Smallest FDR", "Largest absolute estimate"],
        horizontal=True,
        key=f"volcano_rank_{title}"
    )
    rank_by = "FDR" if rank_option == "Smallest FDR" else "effect"

    fig, ax = plt.subplots(figsize=(8, 6))
    nonsig = plot_data[~plot_data["significant"]]
    sig = plot_data[plot_data["significant"]]
    ax.scatter(nonsig["estimate"], nonsig["minus_log10_FDR"], alpha=0.5, label="FDR ≥ 0.05")
    ax.scatter(sig["estimate"], sig["minus_log10_FDR"], alpha=0.8, label="FDR < 0.05")
    ax.axhline(-np.log10(fdr_cutoff), linestyle="--")
    ax.axvline(0, linestyle="--")
    label_top_features(ax=ax, plot_data=plot_data, top_n=10, rank_by=rank_by, fdr_cutoff=fdr_cutoff)
    ax.set_xlabel("Estimate")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    st.pyplot(fig)
    add_pdf_download(
        fig=fig,
        file_name=f"{safe_file_name(title)}_volcano.pdf",
        key=f"download_volcano_{safe_file_name(title)}"
    )
    plt.close(fig)
    st.caption(
        f"{len(sig)} of {len(plot_data)} features have FDR < {fdr_cutoff}. "
        f"At most 10 significant features are labeled."
    )


# ============================================================
# Pairwise volcano
# ============================================================

def plot_pairwise_volcano(
    pairwise,
    fdr_cutoff=0.05
):

    if (
        "contrast"
        not in pairwise.columns
    ):

        st.info(
            "Pairwise contrast information "
            "is not available."
        )

        return

    contrasts = (
        pairwise[
            "contrast"
        ]
        .dropna()
        .unique()
        .tolist()
    )

    if (
        len(contrasts)
        == 0
    ):

        st.info(
            "No pairwise contrasts are available."
        )

        return

    selected_contrast = (
        st.selectbox(
            "Select pairwise comparison",
            contrasts,
            key="pairwise_contrast"
        )
    )

    plot_data = (
        pairwise[
            pairwise[
                "contrast"
            ]
            == selected_contrast
        ]
        .copy()
    )

    plot_data[
        "estimate"
    ] = pd.to_numeric(
        plot_data[
            "estimate"
        ],
        errors="coerce"
    )

    plot_data[
        "FDR"
    ] = pd.to_numeric(
        plot_data[
            "FDR"
        ],
        errors="coerce"
    )

    plot_data = (
        plot_data[
            plot_data[
                "estimate"
            ].notna()
            &
            plot_data[
                "FDR"
            ].notna()
        ]
        .copy()
    )

    if (
        len(plot_data)
        == 0
    ):

        st.info(
            "No valid results are available "
            "for this pairwise comparison."
        )

        return

    minimum_fdr = (
        np.nextafter(
            0,
            1
        )
    )

    plot_data[
        "FDR_plot"
    ] = (
        plot_data[
            "FDR"
        ]
        .clip(
            lower=minimum_fdr
        )
    )

    plot_data[
        "minus_log10_FDR"
    ] = (
        -np.log10(
            plot_data[
                "FDR_plot"
            ]
        )
    )

    plot_data[
        "significant"
    ] = (
        plot_data[
            "FDR"
        ]
        < fdr_cutoff
    )

    rank_option = (
        st.radio(
            "Label significant features by",
            [
                "Smallest FDR",
                "Largest effect size"
            ],
            horizontal=True,
            key=(
                f"pairwise_rank_"
                f"{selected_contrast}"
            )
        )
    )

    rank_by = (
        "FDR"
        if (
            rank_option
            == "Smallest FDR"
        )
        else
        "effect"
    )

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    nonsig = (
        plot_data[
            ~plot_data[
                "significant"
            ]
        ]
    )

    sig = (
        plot_data[
            plot_data[
                "significant"
            ]
        ]
    )

    ax.scatter(
        nonsig[
            "estimate"
        ],
        nonsig[
            "minus_log10_FDR"
        ],
        alpha=0.5,
        label="FDR ≥ 0.05"
    )

    ax.scatter(
        sig[
            "estimate"
        ],
        sig[
            "minus_log10_FDR"
        ],
        alpha=0.8,
        label="FDR < 0.05"
    )

    ax.axhline(
        -np.log10(
            fdr_cutoff
        ),
        linestyle="--"
    )

    ax.axvline(
        0,
        linestyle="--"
    )

    label_top_features(
        ax=ax,
        plot_data=plot_data,
        top_n=10,
        rank_by=rank_by,
        fdr_cutoff=fdr_cutoff
    )

    ax.set_xlabel(
        f"Effect estimate: "
        f"{selected_contrast}"
    )

    ax.set_ylabel(
        "-log10(FDR)"
    )

    ax.set_title(
        f"Pairwise Association: "
        f"{selected_contrast}"
    )

    ax.legend()

    fig.tight_layout()

    st.pyplot(
        fig
    )

    add_pdf_download(
        fig=fig,
        file_name=f"{safe_file_name(selected_contrast)}_pairwise_volcano.pdf",
        key=f"download_pairwise_volcano_{safe_file_name(selected_contrast)}"
    )

    plt.close(
        fig
    )

    st.caption(
        f"{len(sig)} of "
        f"{len(plot_data)} features "
        f"have FDR < {fdr_cutoff} "
        f"for {selected_contrast}."
    )


# ============================================================
# Feature Explorer
# ============================================================

def display_feature_explorer(
    results,
    pairwise,
    metadata,
    data_matrix,
    request
):

    if (
        "feature"
        not in results.columns
    ):

        return

    available_features = (
        results[
            "feature"
        ]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    available_features = [
        feature
        for feature in available_features
        if feature
        in data_matrix.columns
    ]

    if (
        len(available_features)
        == 0
    ):

        return

    st.divider()

    st.subheader(
        "Feature Explorer"
    )

    st.write(
        "Select a feature to inspect its model statistics "
        "and the underlying data pattern."
    )

    if (
        "FDR"
        in results.columns
    ):

        feature_order = (
            results[
                [
                    "feature",
                    "FDR"
                ]
            ]
            .copy()
        )

        feature_order[
            "FDR"
        ] = pd.to_numeric(
            feature_order[
                "FDR"
            ],
            errors="coerce"
        )

        feature_order = (
            feature_order
            .sort_values(
                "FDR",
                na_position="last"
            )
        )

        ordered = (
            feature_order[
                "feature"
            ]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )

        available_features = [
            feature
            for feature in ordered
            if feature
            in available_features
        ]

    selected_feature = (
        st.selectbox(
            "Select feature",
            available_features,
            key="feature_explorer_feature"
        )
    )

    display_feature_statistics(
        results=results,
        feature=selected_feature
    )

    plot_feature_data(
        metadata=metadata,
        data_matrix=data_matrix,
        feature=selected_feature,
        request=request
    )


# ============================================================
# Feature model statistics
# ============================================================

def display_feature_statistics(
    results,
    feature
):

    feature_results = (
        results[
            results[
                "feature"
            ].astype(str)
            == str(feature)
        ]
    )

    if (
        len(feature_results)
        == 0
    ):

        return

    row = (
        feature_results
        .iloc[0]
    )

    metric_columns = (
        st.columns(3)
    )

    estimate = (
        row.get(
            "estimate"
        )
    )

    if (
        pd.notna(
            estimate
        )
    ):

        metric_columns[
            0
        ].metric(
            "Effect estimate",
            f"{float(estimate):.4g}"
        )

    else:

        metric_columns[
            0
        ].metric(
            "Effect estimate",
            "N/A"
        )

    p_value = (
        row.get(
            "p_value"
        )
    )

    if (
        pd.notna(
            p_value
        )
    ):

        metric_columns[
            1
        ].metric(
            "P-value",
            f"{float(p_value):.3g}"
        )

    else:

        metric_columns[
            1
        ].metric(
            "P-value",
            "N/A"
        )

    fdr = (
        row.get(
            "FDR"
        )
    )

    if (
        pd.notna(
            fdr
        )
    ):

        metric_columns[
            2
        ].metric(
            "FDR",
            f"{float(fdr):.3g}"
        )

    else:

        metric_columns[
            2
        ].metric(
            "FDR",
            "N/A"
        )


# ============================================================
# Prepare feature-level raw data
# ============================================================

def prepare_feature_data(
    metadata,
    data_matrix,
    feature
):

    if (
        feature
        not in data_matrix.columns
    ):

        return None

    if (
        "SampleID"
        in metadata.columns
        and
        "SampleID"
        in data_matrix.columns
    ):

        sample_id = (
            "SampleID"
        )

    else:

        first_column = (
            data_matrix
            .columns[0]
        )

        if (
            first_column
            not in metadata.columns
        ):

            return None

        sample_id = (
            first_column
        )

    feature_data = (
        data_matrix[
            [
                sample_id,
                feature
            ]
        ]
        .copy()
    )

    plot_data = (
        metadata.merge(
            feature_data,
            on=sample_id,
            how="inner"
        )
    )

    plot_data[
        feature
    ] = pd.to_numeric(
        plot_data[
            feature
        ],
        errors="coerce"
    )

    plot_data = (
        plot_data[
            plot_data[
                feature
            ].notna()
        ]
        .copy()
    )

    return plot_data


# ============================================================
# Choose correct feature-level graph
# ============================================================

def plot_feature_data(
    metadata,
    data_matrix,
    feature,
    request
):

    plot_data = (
        prepare_feature_data(
            metadata=metadata,
            data_matrix=data_matrix,
            feature=feature
        )
    )

    if (
        plot_data is None
        or
        len(plot_data) == 0
    ):

        st.info(
            "No raw data are available "
            "for this feature."
        )

        return

    predictor = (
        request[
            "primary_predictors"
        ][0]
    )

    predictor_type = (
        request
        .get(
            "variable_types",
            {}
        )
        .get(
            predictor,
            "auto"
        )
    )

    analysis_type = (
        request[
            "analysis_type"
        ]
    )

    analysis_goal = (
        request.get(
            "analysis_goal"
        )
    )

    if (
        analysis_type
        == "cross_sectional"
    ):

        if (
            predictor
            not in plot_data.columns
        ):

            st.info(
                "The predictor is not available "
                "for feature plotting."
            )

            return

        if (
            predictor_type
            == "numeric"
        ):

            plot_cross_sectional_numeric(
                plot_data=plot_data,
                feature=feature,
                predictor=predictor
            )

        else:

            plot_cross_sectional_categorical(
                plot_data=plot_data,
                feature=feature,
                predictor=predictor
            )

        return

    if (
        analysis_type
        == "longitudinal"
    ):

        time_variable = (
            request.get(
                "time"
            )
        )

        if (
            time_variable
            not in plot_data.columns
        ):

            st.info(
                "The time variable is not available "
                "for feature plotting."
            )

            return

        if (
            analysis_goal
            == "time_effect"
        ):

            plot_longitudinal_time_effect(
                plot_data=plot_data,
                feature=feature,
                time_variable=time_variable
            )

            return

        if (
            predictor_type
            == "numeric"
        ):

            plot_longitudinal_numeric_predictor(
                plot_data=plot_data,
                feature=feature,
                predictor=predictor,
                time_variable=time_variable
            )

            return

        plot_longitudinal_categorical(
            plot_data=plot_data,
            feature=feature,
            predictor=predictor,
            time_variable=time_variable
        )


# ============================================================
# Cross-sectional categorical plot
# ============================================================

def plot_cross_sectional_categorical(
    plot_data,
    feature,
    predictor
):

    data = (
        plot_data[
            [
                predictor,
                feature
            ]
        ]
        .dropna()
        .copy()
    )

    groups = (
        data[predictor]
        .astype(str)
        .unique()
        .tolist()
    )

    if len(groups) == 0:
        return

    values = [
        data.loc[
            data[predictor].astype(str) == group,
            feature
        ].values
        for group in groups
    ]

    group_colors = get_feature_group_colors(
        len(groups)
    )

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    box = ax.boxplot(
        values,
        tick_labels=groups,
        patch_artist=True,
        medianprops={
            "color": "black",
            "linewidth": 2
        }
    )

    for patch, color in zip(
        box["boxes"],
        group_colors
    ):
        patch.set_facecolor(
            color
        )

        patch.set_alpha(
            0.85
        )

    ax.set_xlabel(
        predictor
    )

    ax.set_ylabel(
        feature
    )

    ax.set_title(
        f"{feature} by {predictor}"
    )

    fig.tight_layout()

    st.pyplot(
        fig
    )

    add_pdf_download(
        fig=fig,
        file_name=(
            f"{safe_file_name(feature)}_"
            f"by_{safe_file_name(predictor)}.pdf"
        ),
        key=(
            f"download_cross_categorical_"
            f"{safe_file_name(feature)}_"
            f"{safe_file_name(predictor)}"
        )
    )

    plt.close(
        fig
    )


# ============================================================
# Cross-sectional numeric plot
# ============================================================

def plot_cross_sectional_numeric(
    plot_data,
    feature,
    predictor
):

    data = (
        plot_data[
            [
                predictor,
                feature
            ]
        ]
        .copy()
    )

    data[
        predictor
    ] = pd.to_numeric(
        data[
            predictor
        ],
        errors="coerce"
    )

    data = (
        data.dropna()
    )

    if (
        len(data)
        < 2
    ):

        return

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    ax.scatter(
        data[
            predictor
        ],
        data[
            feature
        ],
        alpha=0.5
    )

    if (
        data[
            predictor
        ].nunique()
        > 1
    ):

        coefficients = (
            np.polyfit(
                data[
                    predictor
                ],
                data[
                    feature
                ],
                1
            )
        )

        x_line = np.linspace(
            data[
                predictor
            ].min(),
            data[
                predictor
            ].max(),
            100
        )

        y_line = (
            coefficients[
                0
            ]
            * x_line
            +
            coefficients[
                1
            ]
        )

        ax.plot(
            x_line,
            y_line
        )

    ax.set_xlabel(
        predictor
    )

    ax.set_ylabel(
        feature
    )

    ax.set_title(
        f"{feature} vs {predictor}"
    )

    fig.tight_layout()

    st.pyplot(
        fig
    )

    add_pdf_download(
        fig=fig,
        file_name=f"{safe_file_name(feature)}_vs_{safe_file_name(predictor)}.pdf",
        key=f"download_cross_numeric_{safe_file_name(feature)}_{safe_file_name(predictor)}"
    )

    plt.close(
        fig
    )


# ============================================================
# Longitudinal categorical predictor
# ============================================================

def plot_longitudinal_categorical(
    plot_data,
    feature,
    predictor,
    time_variable
):

    data = (
        plot_data[
            [predictor, time_variable, feature]
        ]
        .dropna()
        .copy()
    )

    if len(data) == 0:
        return

    predictor_levels = (
        data[predictor]
        .astype(str)
        .unique()
        .tolist()
    )

    time_levels = sort_time_values(
        data[time_variable]
    )

    group_colors = get_feature_group_colors(
        len(predictor_levels)
    )

    fig, ax = plt.subplots(figsize=(9, 6))
    rng = np.random.default_rng(123)

    for group_index, group in enumerate(predictor_levels):

        group_data = data[
            data[predictor].astype(str) == group
        ]

        color = group_colors[group_index]
        means = []
        sems = []
        positions = []

        for time_index, time_value in enumerate(time_levels):

            values = (
                group_data.loc[
                    group_data[time_variable] == time_value,
                    feature
                ]
                .dropna()
                .values
            )

            if len(values) == 0:
                continue

            jitter = rng.normal(
                0,
                0.045,
                size=len(values)
            )

            ax.scatter(
                np.full(len(values), time_index) + jitter,
                values,
                alpha=0.35,
                s=24,
                color=color
            )

            means.append(np.mean(values))
            sems.append(
                np.std(values, ddof=1) / np.sqrt(len(values))
                if len(values) > 1
                else 0.0
            )
            positions.append(time_index)

        if len(positions) == 0:
            continue

        positions = np.asarray(positions, dtype=float)
        means = np.asarray(means, dtype=float)
        sems = np.asarray(sems, dtype=float)

        ax.plot(
            positions,
            means,
            marker="o",
            linewidth=2,
            color=color,
            label=str(group)
        )

        ax.fill_between(
            positions,
            means - 1.96 * sems,
            means + 1.96 * sems,
            alpha=0.15,
            color=color
        )

    ax.set_xticks(range(len(time_levels)))
    ax.set_xticklabels([str(value) for value in time_levels])
    ax.set_xlabel(time_variable)
    ax.set_ylabel(feature)
    ax.set_title(f"{feature}: longitudinal pattern")
    ax.legend(title=predictor)

    fig.tight_layout()
    st.pyplot(fig)

    add_pdf_download(
        fig=fig,
        file_name=(
            f"{safe_file_name(feature)}_"
            f"longitudinal_by_{safe_file_name(predictor)}.pdf"
        ),
        key=(
            f"download_longitudinal_categorical_"
            f"{safe_file_name(feature)}_"
            f"{safe_file_name(predictor)}"
        )
    )

    plt.close(fig)


# ============================================================
# Longitudinal time effect
# ============================================================

def plot_longitudinal_time_effect(
    plot_data,
    feature,
    time_variable
):

    data = (
        plot_data[
            [time_variable, feature]
        ]
        .dropna()
        .copy()
    )

    if len(data) == 0:
        return

    time_levels = sort_time_values(
        data[time_variable]
    )

    color = get_feature_group_colors(1)[0]

    fig, ax = plt.subplots(figsize=(8, 5))
    rng = np.random.default_rng(123)

    means = []
    sems = []
    positions = []

    for time_index, time_value in enumerate(time_levels):

        values = (
            data.loc[
                data[time_variable] == time_value,
                feature
            ]
            .dropna()
            .values
        )

        if len(values) == 0:
            continue

        jitter = rng.normal(
            0,
            0.045,
            size=len(values)
        )

        ax.scatter(
            np.full(len(values), time_index) + jitter,
            values,
            alpha=0.35,
            s=24,
            color=color
        )

        means.append(np.mean(values))
        sems.append(
            np.std(values, ddof=1) / np.sqrt(len(values))
            if len(values) > 1
            else 0.0
        )
        positions.append(time_index)

    if len(positions) == 0:
        plt.close(fig)
        return

    positions = np.asarray(positions, dtype=float)
    means = np.asarray(means, dtype=float)
    sems = np.asarray(sems, dtype=float)

    ax.plot(
        positions,
        means,
        marker="o",
        linewidth=2,
        color=color
    )

    ax.fill_between(
        positions,
        means - 1.96 * sems,
        means + 1.96 * sems,
        alpha=0.15,
        color=color
    )

    ax.set_xticks(range(len(time_levels)))
    ax.set_xticklabels([str(value) for value in time_levels])
    ax.set_xlabel(time_variable)
    ax.set_ylabel(feature)
    ax.set_title(f"{feature}: longitudinal time pattern")

    fig.tight_layout()
    st.pyplot(fig)

    add_pdf_download(
        fig=fig,
        file_name=f"{safe_file_name(feature)}_time_effect.pdf",
        key=f"download_time_effect_{safe_file_name(feature)}"
    )

    plt.close(fig)


# ============================================================
# Longitudinal numeric predictor
# ============================================================

def plot_longitudinal_numeric_predictor(
    plot_data,
    feature,
    predictor,
    time_variable
):

    data = (
        plot_data[
            [
                predictor,
                time_variable,
                feature
            ]
        ]
        .copy()
    )

    data[
        predictor
    ] = pd.to_numeric(
        data[
            predictor
        ],
        errors="coerce"
    )

    data = (
        data.dropna()
    )

    if (
        len(data)
        == 0
    ):

        return

    time_levels = (
        sort_time_values(
            data[
                time_variable
            ]
        )
    )

    time_map = {
        value: index
        for index, value
        in enumerate(
            time_levels
        )
    }

    x = (
        data[
            time_variable
        ]
        .map(
            time_map
        )
        .astype(float)
    )

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    scatter = ax.scatter(
        x,
        data[
            feature
        ],
        c=data[
            predictor
        ],
        alpha=0.5
    )

    colorbar = (
        fig.colorbar(
            scatter,
            ax=ax
        )
    )

    colorbar.set_label(
        predictor
    )

    ax.set_xticks(
        range(
            len(
                time_levels
            )
        )
    )

    ax.set_xticklabels(
        [
            str(value)
            for value
            in time_levels
        ]
    )

    ax.set_xlabel(
        time_variable
    )

    ax.set_ylabel(
        feature
    )

    ax.set_title(
        f"{feature}: longitudinal pattern"
    )

    fig.tight_layout()

    st.pyplot(
        fig
    )

    add_pdf_download(
        fig=fig,
        file_name=f"{safe_file_name(feature)}_longitudinal_{safe_file_name(predictor)}.pdf",
        key=f"download_longitudinal_numeric_{safe_file_name(feature)}_{safe_file_name(predictor)}"
    )

    plt.close(
        fig
    )


# ============================================================
# Sort longitudinal time values
# ============================================================

def sort_time_values(
    series
):

    values = (
        series
        .dropna()
        .unique()
        .tolist()
    )

    try:

        return sorted(
            values,
            key=float
        )

    except (
        TypeError,
        ValueError
    ):

        return sorted(
            values,
            key=lambda x: str(x)
        )
