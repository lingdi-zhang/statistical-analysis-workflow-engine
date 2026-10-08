# Omics Analysis Workbench

**Cross-sectional and longitudinal statistical modeling.**

A laptop-friendly Streamlit application that automates high-throughput cross-sectional and longitudinal omics modeling while keeping statistical specifications, diagnostics, and model-recovery decisions explicit.

The application provides a Streamlit interface for defining an analysis, validating metadata and feature matrices, fitting statistical models in R, handling common mixed-model failures, and generating interpretable results and downloadable figures.

---

## Why I Built This

Longitudinal omics analyses often require repeatedly writing custom code to:

- Align metadata and feature matrices
- Specify statistical models
- Handle categorical and continuous predictors
- Fit hundreds or thousands of feature-level models
- Diagnose convergence and singularity problems
- Correct for multiple testing
- Generate publication-ready visualizations

This project automates repetitive implementation and model-fitting tasks while keeping statistical specifications and model-recovery decisions explicit and visible to the user.

---

## Statistical Guardrails

The workflow automates model execution rather than statistical decision-making.

- Analysis type is explicitly selected by the user.
- Predictor type determines the valid statistical test.
- Model specifications remain visible to the user.
- Failed and singular models are reported rather than silently discarded.
- Random-effects fallback follows a predefined rule: if ≥20% of converged feature-level random-slope models are singular, the analysis is refit using random intercepts without random slopes.
- Multiple-testing correction is applied across feature-level tests.

---

## Key Capabilities

### Cross-Sectional Analysis

Supports associations between molecular features and:

- Continuous predictors
- Binary categorical predictors
- Multi-category predictors
- Ordered categorical predictors

Depending on the predictor type, the workflow performs the appropriate coefficient, overall, trend, or pairwise tests.

### Longitudinal Analysis

Uses linear mixed-effects models to evaluate:

- **Overall association** — whether a predictor is associated with the feature across repeated measurements
- **Time effect** — whether the feature changes over time
- **Trajectory difference** — whether longitudinal changes differ between predictor groups

The longitudinal UI lets you explicitly select random intercepts only or random
intercepts plus time slopes. The initial selection remains intercepts plus slopes.
Inspection warns per feature when usable observations do not exceed twice the
number of subjects; select random intercepts only and inspect again for such
features, including datasets with exactly two visits per subject. Unsupported
slope models are reported as feature-level failures without changing the
selection. The existing singularity fallback applies to converged slope models.

Models can include:

- Random subject intercepts
- Random slopes for time
- Covariates
- Predictor × time interactions

Example longitudinal model specifications:

#### Trajectory Difference

Tests the `Phenotype × Time` interaction:

```text
Feature ~ Phenotype * Time + covariates + (1 + Time | SubjectID)
```

#### Overall Association

Tests the phenotype effect while adjusting for time and specified covariates:

```text
Feature ~ Phenotype + Time + covariates + (1 + Time | SubjectID)
```

#### Time Effect

Tests the time effect while adjusting for phenotype and specified covariates:

```text
Feature ~ Phenotype + Time + covariates + (1 + Time | SubjectID)
```

### Robust Model Fitting

High-throughput mixed models frequently encounter numerical problems. The workflow therefore includes:

- Optimizer retries
- Convergence tracking
- Singular-fit detection
- Automatic fallback from random intercept + slope to random intercept models when singularity is widespread
- Per-feature failure handling so one failed model does not terminate the entire analysis

After model-level diagnostics determine the random-effects structure, the final model specification is applied consistently across features to preserve comparability.

### Multiple-Testing Correction

P-values are corrected using the **Benjamini-Hochberg false discovery rate (FDR)**.
Overall tests are corrected across features. Pairwise tests are corrected across
features separately for each contrast.

### Model Diagnostics

For high-throughput analyses, the application reports:

- Number of features analyzed
- Convergence fraction
- Model failure count
- Singularity fraction for the final models (with the initial fraction reported separately after fallback)
- Final random-effects specification

This makes model-fitting problems visible rather than silently ignoring them.

### Results and Visualization

Results are automatically summarized using visualizations appropriate for the selected analysis.

Outputs include:

- Volcano plots for high-dimensional association results
- Feature-level plots for selected signals
- Longitudinal trajectory plots
- Compact result tables for small feature sets
- Downloadable PDF figures

---

## Workflow

```text
┌─────────────────────┐
│     Upload Data     │
│ Metadata + Features │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Analysis Request   │
│ Predictor / Time /  │
│     Covariates      │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│    Data Inspector   │
│  Alignment + QC     │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ Statistical Engine  │
│  R + Mixed Models   │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│    Results + QC     │
│ FDR + Diagnostics   │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│    Visualization    │
│  Tables + Figures   │
└─────────────────────┘
```

---

## Technology Stack

### Frontend and Orchestration

- Python
- Streamlit
- pandas
- NumPy
- Matplotlib

### Statistical Modeling

- R
- `lme4`
- `lmerTest`
- `emmeans`

Python handles user interaction, data validation, analysis configuration, and visualization, while R performs the statistical modeling.

---

## Input Data

The application expects two CSV files. Each file must start with its column
headers; leading blank lines (including whitespace-only lines) are rejected.
Use empty cells for missing data. UI uploads containing literal `NA` or `null`
cell values are rejected with an error asking you to replace them with empty
cells if they mean missing data, or rename them if they are intentional labels.

### 1. Metadata

One row per sample containing variables such as:

```text
SampleID,SubjectID,Visit,Disease,Age,Sex
S001,P001,0,Control,35,F
S002,P001,2,Control,35,F
S003,P002,0,Disease,42,M
S004,P002,2,Disease,42,M
```

### 2. Feature Matrix

One row per sample and one column per molecular feature:

```text
SampleID,Gene1,Gene2,Gene3
S001,4.32,7.81,2.14
S002,4.91,7.22,2.67
S003,6.14,5.93,3.82
S004,6.83,5.41,4.12
```

Samples are aligned using `SampleID` before analysis, including direct R calls
and the CSV/JSON command-line workflow. Missing or duplicate IDs are rejected;
unmatched samples are removed with a warning. The command-line workflow requires
`SampleID` in both files. Low-level R calls without IDs use row order and require
equal row counts.

---

## Project Structure

```text
statistical-analysis-workflow-engine/
│
├── src/
│   ├── app.py
│   ├── data_loading.py
│   ├── upload_signature.py
│   ├── inspector.py
│   ├── request_parse_streamlit.py
│   ├── r_runner.py
│   └── results.py
│
├── r/
│   ├── run_analysis.R
│   ├── model_runner.R
│   └── prepare_model_variables.R
│
├── tests/
│   ├── regression.R
│   ├── test_data_loading.py
│   ├── test_inspector.py
│   ├── test_model_summary.py
│   ├── test_plot_time.py
│   ├── test_r_transfer.py
│   └── test_upload_signature.py
│
├── environment.yml
└── README.md
```

---

## Running Locally

### 1. Clone the Repository

```bash
git clone https://github.com/lingdi-zhang/statistical-analysis-workflow-engine.git
cd statistical-analysis-workflow-engine
```

### 2. Create and Activate the Environment

```bash
micromamba create -f environment.yml
micromamba activate longitudinal-omics-agent
```

### 3. Launch the Application

From the repository root:

```bash
streamlit run src/app.py
```

Open the local Streamlit URL shown in the terminal.

### Regression Checks

From the repository root, run:

```bash
python -m unittest discover -s tests
Rscript tests/regression.R
```


---

## Design Philosophy

The goal of this project is **not to replace statistical judgment with a black-box model**.

Instead, the application separates analysis into explicit stages:

```text
Analysis Specification
      ↓
Data Validation
      ↓
Statistical Modeling
      ↓
Model Diagnostics
      ↓
Multiple-Testing Correction
      ↓
Visualization
```

This separation keeps statistical decisions transparent while automating repetitive implementation and model-fitting tasks.

---

## Current Scope

The current implementation focuses on:

- Numeric molecular features
- Cross-sectional regression
- Longitudinal linear mixed-effects models
- Continuous and categorical predictors
- Omics-scale multiple testing
- Automated model diagnostics
- Automated visualization

---

## License

This project is intended for research and educational use.
