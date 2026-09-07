import streamlit as st
import pandas as pd

from request_parse_streamlit import render_request_builder
from inspector import DataInspector
from r_runner import run_r_analysis
from results import display_analysis_results


# ============================================================
# Page setup
# ============================================================

st.set_page_config(
    page_title="Cross Sectional and Longitudinal Analysis Agent",
    layout="wide"
)

st.markdown(
    """
    <style>

    /* ---------- Main text ---------- */
    .stApp {
        font-size: 20px;
    }

    .stMarkdown p {
        font-size: 20px !important;
        line-height: 1.6 !important;
    }

    /* ---------- Headings ---------- */
    h1 {
        font-size: 42px !important;
    }

    h2 {
        font-size: 32px !important;
    }

    h3 {
        font-size: 25px !important;
    }

    /* ---------- Widget labels ---------- */
    [data-testid="stWidgetLabel"] p {
        font-size: 19px !important;
        font-weight: 500 !important;
    }

    /* ---------- Selectbox / multiselect ---------- */
    div[data-baseweb="select"] * {
        font-size: 18px !important;
    }

    /* ---------- Radio buttons ---------- */
    [data-testid="stRadio"] label p {
        font-size: 18px !important;
    }

    /* ---------- Checkboxes ---------- */
    [data-testid="stCheckbox"] label p {
        font-size: 18px !important;
    }

    /* ---------- Buttons ---------- */
    .stButton button p {
        font-size: 19px !important;
    }

    .stDownloadButton button p {
        font-size: 18px !important;
    }

    /* ---------- File uploader ---------- */
    [data-testid="stFileUploader"] {
        font-size: 18px !important;
    }

    [data-testid="stFileUploader"] * {
        font-size: 18px !important;
    }

    /* ---------- Metrics ---------- */
    [data-testid="stMetricLabel"] p {
        font-size: 18px !important;
    }

    [data-testid="stMetricValue"] {
        font-size: 30px !important;
    }

    /* ---------- Alerts / info / warning ---------- */
    [data-testid="stAlert"] p {
        font-size: 18px !important;
    }

    /* ---------- Expanders ---------- */
    [data-testid="stExpander"] summary p {
        font-size: 18px !important;
    }

    /* ---------- Captions ---------- */
    [data-testid="stCaptionContainer"] p {
        font-size: 16px !important;
    }

    </style>
    """,
    unsafe_allow_html=True
)

st.title("Statistical Analysis Agent")

st.write(
    "Upload your metadata and feature matrix, "
    "define the analysis, inspect the data, "
    "and run the statistical model."
)


# ============================================================
# Session-state helpers
# ============================================================

def clear_analysis_results():
    """
    Clear model results while keeping the current inspection.
    """
    for key in [
        "analysis_output",
        "analysis_metadata",
        "analysis_data_matrix",
        "analysis_request"
    ]:
        st.session_state.pop(
            key,
            None
        )


def clear_inspection_and_results():
    """
    Clear both inspection and model results.
    """
    for key in [
        "inspection_report",
        "inspection_metadata",
        "inspection_data_matrix",
        "inspection_request",
        "analysis_output",
        "analysis_metadata",
        "analysis_data_matrix",
        "analysis_request"
    ]:
        st.session_state.pop(
            key,
            None
        )


def display_inspection_report(
    report
):
    """
    Display a previously completed DataInspector report.
    """

    st.header(
        "4. Data Inspection"
    )

    if (
        report["status"]
        == "failed"
    ):

        st.error(
            "Data inspection failed."
        )

    elif (
        report["status"]
        == "passed_with_warnings"
    ):

        st.warning(
            "Data inspection passed with warnings."
        )

    else:

        st.success(
            "Data inspection passed."
        )

    # --------------------------------------------------------
    # Errors
    # --------------------------------------------------------

    for error in report.get(
        "errors",
        []
    ):

        st.error(
            error
        )

    # --------------------------------------------------------
    # Warnings
    # --------------------------------------------------------

    for warning in report.get(
        "warnings",
        []
    ):

        st.warning(
            warning
        )

    # --------------------------------------------------------
    # Sample matching
    # --------------------------------------------------------

    sample_matching = (
        report.get(
            "sample_matching"
        )
    )

    if (
        sample_matching
        is not None
    ):

        st.subheader(
            "Sample Matching"
        )

        st.write(
            f"Shared samples: "
            f"{sample_matching['n_shared_samples']}"
        )

    # --------------------------------------------------------
    # Variable inspection
    # --------------------------------------------------------

    variables = (
        report.get(
            "variables"
        )
    )

    if variables:

        st.subheader(
            "Variable Inspection"
        )

        st.json(
            variables
        )

    # --------------------------------------------------------
    # Longitudinal structure
    # --------------------------------------------------------

    longitudinal = (
        report.get(
            "longitudinal"
        )
    )

    if (
        longitudinal
        is not None
    ):

        st.subheader(
            "Longitudinal Structure"
        )

        st.json(
            longitudinal
        )


# ============================================================
# 1. Upload data
# ============================================================

st.header(
    "1. Upload Data"
)

col1, col2 = st.columns(
    2
)

with col1:

    metadata_file = st.file_uploader(
        "Upload metadata",
        type=[
            "csv"
        ],
        help=(
            "CSV file with SampleID as the first column."
        )
    )

with col2:

    data_file = st.file_uploader(
        "Upload feature matrix",
        type=[
            "csv"
        ],
        help=(
            "CSV file with SampleID as the first column. "
            "All remaining columns should be features."
        )
    )


# ============================================================
# Continue after both files are uploaded
# ============================================================

if (
    metadata_file is not None
    and
    data_file is not None
):

    metadata = pd.read_csv(
        metadata_file
    )

    data_matrix = pd.read_csv(
        data_file
    )

    # --------------------------------------------------------
    # Detect changed uploads
    #
    # If the uploaded data change, previous inspection/model
    # results should no longer be considered valid.
    # --------------------------------------------------------

    current_upload_signature = (
        metadata_file.name,
        metadata_file.size,
        data_file.name,
        data_file.size
    )

    previous_upload_signature = (
        st.session_state.get(
            "upload_signature"
        )
    )

    if (
        previous_upload_signature
        is not None
        and
        previous_upload_signature
        != current_upload_signature
    ):

        clear_inspection_and_results()

    st.session_state[
        "upload_signature"
    ] = current_upload_signature

    # --------------------------------------------------------
    # Metadata preview
    # --------------------------------------------------------

    st.success(
        f"Metadata loaded: "
        f"{metadata.shape[0]} rows × "
        f"{metadata.shape[1]} columns"
    )

    with st.expander(
        "Preview metadata"
    ):

        st.dataframe(
            metadata.head(
                10
            ),
            use_container_width=True
        )

    # --------------------------------------------------------
    # Feature matrix preview
    # --------------------------------------------------------

    st.success(
        f"Feature matrix loaded: "
        f"{data_matrix.shape[0]} rows × "
        f"{data_matrix.shape[1]} columns"
    )

    with st.expander(
        "Preview feature matrix"
    ):

        st.dataframe(
            data_matrix.head(
                10
            ),
            use_container_width=True
        )

    st.divider()

    # ========================================================
    # 2-3. Build analysis request
    # ========================================================

    request = render_request_builder(
        metadata
    )

    # ========================================================
    # 4. Run inspection only when user clicks Inspect Data
    # ========================================================

    if (
        request is not None
    ):

        # A newly submitted analysis request invalidates
        # previous model results.
        clear_analysis_results()

        inspector = DataInspector(
            metadata=metadata,
            data_matrix=data_matrix,
            request=request
        )

        report = (
            inspector.inspect()
        )

        # Store the inspection report so it survives
        # Streamlit reruns caused by later UI interactions.
        st.session_state[
            "inspection_report"
        ] = report

        if (
            report["status"]
            != "failed"
        ):

            analysis_inputs = (
                inspector.get_analysis_inputs()
            )

            st.session_state[
                "inspection_metadata"
            ] = analysis_inputs[
                "metadata"
            ]

            st.session_state[
                "inspection_data_matrix"
            ] = analysis_inputs[
                "data_matrix"
            ]

            st.session_state[
                "inspection_request"
            ] = analysis_inputs[
                "request"
            ]

        else:

            for key in [
                "inspection_metadata",
                "inspection_data_matrix",
                "inspection_request"
            ]:

                st.session_state.pop(
                    key,
                    None
                )

    # ========================================================
    # 5. Display stored inspection
    # ========================================================

    inspection_report = (
        st.session_state.get(
            "inspection_report"
        )
    )

    if (
        inspection_report
        is not None
    ):

        display_inspection_report(
            inspection_report
        )

        # ====================================================
        # 6. Run analysis only when explicitly requested
        # ====================================================

        if (
            inspection_report[
                "status"
            ]
            != "failed"
            and
            "inspection_request"
            in st.session_state
        ):

            st.header(
                "5. Run Analysis"
            )

            st.write(
                "The data have already been inspected. "
                "Run the statistical models when ready."
            )

            if st.button(
                "Run Analysis",
                type="primary",
                key="run_analysis"
            ):

                with st.spinner(
                    "Running statistical models..."
                ):

                    analysis_output = (
                        run_r_analysis(
                            metadata=(
                                st.session_state[
                                    "inspection_metadata"
                                ]
                            ),
                            data_matrix=(
                                st.session_state[
                                    "inspection_data_matrix"
                                ]
                            ),
                            request=(
                                st.session_state[
                                    "inspection_request"
                                ]
                            )
                        )
                    )

                # --------------------------------------------
                # Persist completed analysis
                # --------------------------------------------

                st.session_state[
                    "analysis_output"
                ] = analysis_output

                st.session_state[
                    "analysis_metadata"
                ] = st.session_state[
                    "inspection_metadata"
                ]

                st.session_state[
                    "analysis_data_matrix"
                ] = st.session_state[
                    "inspection_data_matrix"
                ]

                st.session_state[
                    "analysis_request"
                ] = st.session_state[
                    "inspection_request"
                ]

    # ========================================================
    # 7. Display stored results
    #
    # Changing volcano labels, pairwise contrasts, or Feature
    # Explorer selections now reruns only the Streamlit display.
    # It does NOT rerun DataInspector or the R analysis.
    # ========================================================

    if (
        "analysis_output"
        in st.session_state
    ):

        st.divider()

        display_analysis_results(
            analysis_output=(
                st.session_state[
                    "analysis_output"
                ]
            ),
            request=(
                st.session_state[
                    "analysis_request"
                ]
            ),
            metadata=(
                st.session_state[
                    "analysis_metadata"
                ]
            ),
            data_matrix=(
                st.session_state[
                    "analysis_data_matrix"
                ]
            )
        )

else:

    st.info(
        "Upload both the metadata and feature matrix "
        "CSV files to begin."
    )
