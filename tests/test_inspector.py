import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from inspector import DataInspector


class InspectorTests(unittest.TestCase):
    def setUp(self):
        self.metadata = pd.DataFrame({
            "SampleID": [f"S{i}" for i in range(6)],
            "group": ["a", "b", "c", "a", "b", "c"],
            "age": [20, 21, 22, 23, 24, 25],
            "subject": [1, 1, 2, 2, 3, 3],
            "time": [0, 12, 0, 12, 0, 12],
        })
        self.features = pd.DataFrame({"SampleID": self.metadata.SampleID,
                                      "gene": [1., 2., 3., 4., 5., 6.]})
        self.request = {"analysis_type": "cross_sectional",
                        "primary_predictors": ["group"], "covariates": ["age"],
                        "variable_types": {"group": "categorical", "age": "numeric"}}

    def inspect(self):
        inspector = DataInspector(self.metadata, self.features, self.request)
        return inspector, inspector.inspect()

    def test_auto_text_infinity_is_rejected_for_predictors_and_covariates(self):
        for variable in ("group", "age"):
            for infinity in ("inf", "-inf"):
                with self.subTest(variable=variable, infinity=infinity):
                    self.metadata[variable] = pd.Series(["1", infinity] * 3, dtype="string")
                    self.request["variable_types"][variable] = "auto"
                    report = self.inspect()[1]
                    self.assertEqual(report["status"], "failed")
                    self.assertTrue(any(f"{variable} contains 3 infinite" in error
                                        for error in report["errors"]))

    def test_explicit_categorical_infinity_label_is_preserved(self):
        self.metadata["group"] = pd.Series(["control", "inf"] * 3, dtype="string")
        inspector, report = self.inspect()
        self.assertNotEqual(report["status"], "failed")
        self.assertIn("inf", report["variables"]["group"]["levels"])

    def test_auto_binary_numeric_text_remains_categorical(self):
        self.metadata["group"] = pd.Series(["0", "1"] * 3, dtype="string")
        self.request["variable_types"]["group"] = "auto"
        _, report = self.inspect()
        self.assertNotEqual(report["status"], "failed")
        self.assertEqual(report["variables"]["group"]["type"], "categorical")

    def test_two_visit_random_effects_warning_and_explicit_choice(self):
        self.request.update(analysis_type="longitudinal", time="time",
                            subject_id="subject", analysis_goal="time_effect",
                            random_slope=True)
        inspector, report = self.inspect()
        self.assertEqual(report["status"], "passed_with_warnings")
        self.assertTrue(any("6 usable observations for 3 subjects" in w for w in report["warnings"]))
        self.assertTrue(inspector.get_analysis_inputs()["request"]["random_slope"])
        self.request["random_slope"] = False
        inspector, report = self.inspect()
        self.assertNotEqual(report["status"], "failed")
        self.assertFalse(any("random intercept + time slope is unsupported" in w for w in report["warnings"]))
        self.assertFalse(inspector.get_analysis_inputs()["request"]["random_slope"])

    def test_slope_warning_uses_each_features_complete_cases(self):
        self.metadata = pd.DataFrame({
            "SampleID": [f"S{i}" for i in range(9)],
            "subject": [1, 1, 1, 2, 2, 2, 3, 3, 3],
            "time": [0, 1, 2] * 3,
            "x": [1, 2, 4, 3, 5, 8, 6, 9, 7],
        })
        self.features = pd.DataFrame({"SampleID": self.metadata.SampleID,
                                      "good": list(range(9)),
                                      "limited": [1., 2., None, 3., 4., None, 5., 6., None]})
        self.request = {"analysis_type": "longitudinal", "primary_predictors": ["x"],
                        "time": "time", "subject_id": "subject", "analysis_goal": "time_effect",
                        "covariates": [], "variable_types": {"x": "numeric"}, "random_slope": True}
        report = self.inspect()[1]
        slope_warnings = [w for w in report["warnings"] if "time slope is unsupported" in w]
        self.assertEqual(len(slope_warnings), 1)
        self.assertIn("Feature 'limited'", slope_warnings[0])

    def test_valid_and_reordered_samples(self):
        self.features = self.features.iloc[::-1]
        self.assertNotEqual(self.inspect()[1]["status"], "failed")

    def test_infinite_numeric_metadata_stops_inspection(self):
        for variable in ("age", "time"):
            for value in (float("inf"), float("-inf")):
                with self.subTest(variable=variable, value=value):
                    self.setUp()
                    if variable == "time":
                        self.request.update(analysis_type="longitudinal", time="time",
                                            subject_id="subject", analysis_goal="time_effect")
                    else:
                        self.request.update(primary_predictors=["age"], covariates=[])
                    self.metadata[variable] = self.metadata[variable].astype(float)
                    self.metadata.loc[0, variable] = value
                    report = self.inspect()[1]
                    self.assertEqual(report["status"], "failed")
                    self.assertTrue(any(f"{variable} contains 1 infinite" in e for e in report["errors"]))
        self.setUp()
        self.metadata["age"] = self.metadata["age"].astype(float)
        self.metadata.loc[0, "age"] = float("inf")
        self.assertEqual(self.inspect()[1]["status"], "failed")

    def test_auto_detection_does_not_hide_infinity(self):
        self.metadata["age"] = [1., float("inf")] * 3
        self.request["variable_types"]["age"] = "auto"
        self.assertEqual(self.inspect()[1]["status"], "failed")

    def test_numeric_predictor_rejects_invalid_outcomes(self):
        self.request.update(primary_predictors=["age"], covariates=[])
        for values, expected in [(["invalid"] * 6, "nonnumeric"),
                                 ([float("inf")] * 6, "infinite"),
                                 ([float("nan")] * 6, "no usable observations"),
                                 ([2.] * 6, "constant")]:
            with self.subTest(expected=expected):
                self.features["gene"] = values
                report = self.inspect()[1]
                self.assertEqual(report["status"], "failed")
                self.assertTrue(any(expected in e for e in report["errors"]))

    def test_constant_outcome_after_complete_case_filtering(self):
        self.features["gene"] = [1., 1., 1., 2., 3., 4.]
        self.metadata.loc[[3, 4, 5], "age"] = float("nan")
        report = self.inspect()[1]
        self.assertEqual(report["status"], "failed")
        self.assertTrue(any("constant" in e for e in report["errors"]))

    def test_empty_feature_matrix(self):
        self.features = self.features[["SampleID"]]
        self.assertEqual(self.inspect()[1]["status"], "failed")

    def test_collision_with_unselected_metadata_column(self):
        self.features.rename(columns={"gene": "time"}, inplace=True)
        report = self.inspect()[1]
        self.assertEqual(report["status"], "failed")
        self.assertTrue(any("conflicts with metadata" in e for e in report["errors"]))

    def test_missing_outcomes_remove_category(self):
        self.features.loc[[2, 5], "gene"] = float("nan")
        report = self.inspect()[1]
        self.assertEqual(report["status"], "failed")
        self.assertTrue(any("categories: c" in e for e in report["errors"]))

    def test_missing_covariate_removes_category(self):
        self.metadata.loc[[2, 5], "age"] = float("nan")
        self.assertEqual(self.inspect()[1]["status"], "failed")

    def test_invalid_feature_values_remove_category(self):
        self.features["gene"] = self.features["gene"].astype(object)
        self.features.loc[[2, 5], "gene"] = "invalid"
        self.assertEqual(self.inspect()[1]["status"], "failed")

    def test_expected_categories_follow_subset(self):
        self.metadata["subset"] = ["keep", "keep", "drop"] * 2
        self.request["subset"] = {"enabled": True, "variable": "subset", "value": "keep"}
        self.assertNotEqual(self.inspect()[1]["status"], "failed")

    def test_ordered_predictor_category_loss(self):
        self.request["variable_types"]["group"] = "ordered_categorical"
        self.request["ordered_levels"] = {"group": ["a", "b", "c"]}
        self.features.loc[[2, 5], "gene"] = float("nan")
        self.assertEqual(self.inspect()[1]["status"], "failed")

    def test_two_visit_time_preserves_units(self):
        self.request.update(analysis_type="longitudinal", time="time",
                            subject_id="subject", analysis_goal="time_effect")
        inspector, report = self.inspect()
        self.assertNotEqual(report["status"], "failed")
        self.assertEqual(inspector.request["variable_types"]["time"], "numeric")
        self.assertEqual(inspector.metadata.time.tolist(), [0, 12, 0, 12, 0, 12])

    def test_text_time_is_rejected(self):
        self.request.update(analysis_type="longitudinal", time="time",
                            subject_id="subject", analysis_goal="time_effect")
        self.metadata["time"] = ["baseline", "followup"] * 3
        self.assertEqual(self.inspect()[1]["status"], "failed")

    def test_missing_subject_removes_category(self):
        self.request.update(analysis_type="longitudinal", time="time",
                            subject_id="subject", analysis_goal="time_effect")
        self.metadata.loc[[2, 5], "subject"] = float("nan")
        self.assertEqual(self.inspect()[1]["status"], "failed")


if __name__ == "__main__":
    unittest.main()
