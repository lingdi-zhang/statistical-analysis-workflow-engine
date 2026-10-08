import os
import json
import shutil
import sys
import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from inspector import DataInspector
from r_runner import run_r_analysis


class ResultLabelTests(unittest.TestCase):
    def test_numeric_missing_markers_are_not_truncated(self):
        metadata = pd.DataFrame({'SampleID': ['001', '002', '003'],
                                 'x': [1.0, float('nan'), 3.0]})
        features = pd.DataFrame({'SampleID': metadata.SampleID,
                                 'gene': [float('nan'), 2.0, 3.0]})
        originals = [table.copy(deep=True) for table in (metadata, features)]
        request = {'variable_types': {'x': 'numeric'}}

        def fake_r(command, **kwargs):
            marker = json.loads(Path(command[4]).read_text())['transfer_na_marker']
            self.assertGreater(len(marker), 32)
            for path, variable, row in ((command[2], 'x', 1), (command[3], 'gene', 0)):
                transferred = pd.read_csv(path, dtype='string', keep_default_na=False)
                self.assertEqual(transferred.loc[row, variable], marker)
                self.assertEqual(transferred.SampleID.tolist(), ['001', '002', '003'])
            Path(command[5]).write_text('feature,p_value\ngene,0.01\n')
            Path(command[6]).write_text('')
            Path(command[7]).write_text('{}')
            return subprocess.CompletedProcess(command, 0, '', '')

        with patch('r_runner.subprocess.run', side_effect=fake_r):
            run_r_analysis(metadata, features, request)
        for actual, original in zip((metadata, features), originals):
            pd.testing.assert_frame_equal(actual, original)
        self.assertNotIn('transfer_na_marker', request)

    def test_literal_na_labels_missing_statistics_and_diagnostics(self):
        def fake_r(command, **kwargs):
            marker = json.loads(Path(command[4]).read_text())["transfer_na_marker"]
            Path(command[5]).write_text(
                f'feature,predictor,term,estimate,p_value,error\nNA,null,NA,1.5,0.01,{marker}\n'
                f'001,null,{marker},{marker},{marker},optimizer failed\n')
            Path(command[6]).write_text(
                f'feature,predictor,contrast,estimate,p.value\nnull,NA,null,{marker},{marker}\n')
            Path(command[7]).write_text('{"warning":"Some models failed"}')
            Path(command[7]).with_name("diagnostics.json").write_text(
                '{"final":{"001":{"nloptwrap":{"message":"optimizer failed"}}},"first_pass":null}')
            return subprocess.CompletedProcess(command, 0, '', '')

        with patch('r_runner.subprocess.run', side_effect=fake_r):
            output = run_r_analysis(pd.DataFrame({'SampleID': ['a']}),
                                    pd.DataFrame({'NA': [1], 'null': [2]}), {})
        self.assertEqual(output['overall'].feature.tolist(), ['NA', '001'])
        self.assertEqual(output['overall'].iloc[0]['predictor'], 'null')
        self.assertEqual(output['overall'].iloc[0]['term'], 'NA')
        self.assertTrue(pd.isna(output['overall'].iloc[1]['term']))
        self.assertTrue(pd.isna(output['overall'].iloc[1]['estimate']))
        self.assertTrue(pd.isna(output['overall'].iloc[1]['p_value']))
        self.assertTrue(pd.isna(output['overall'].iloc[0]['error']))
        self.assertEqual(output['pairwise'].iloc[0]['feature'], 'null')
        self.assertEqual(output['pairwise'].iloc[0]['contrast'], 'null')
        self.assertTrue(pd.isna(output['pairwise'].iloc[0]['p.value']))
        self.assertEqual(output['diagnostics']['final']['001']['nloptwrap']['message'],
                         'optimizer failed')

    def test_missing_metadata_survives_r_reading(self):
        executable = shutil.which("Rscript")
        if not executable:
            self.skipTest("Rscript is required for the CSV transfer check")
        real_run = subprocess.run
        metadata = pd.DataFrame({
            "SampleID": ["a", "b", "c", "d", "e"],
            "group": ["A", None, "B", "A", "NA"],
            "covariate": ["X", "X", None, "Y", ""],
            "subject": ["001", "001", "1", None, "1"],
        })
        request = {"variable_types": {"group": "categorical", "covariate": "categorical"}}
        original = metadata.copy(deep=True)

        def fake_r(command, **kwargs):
            transferred = json.loads(Path(command[4]).read_text())
            marker = transferred["transfer_na_marker"]
            # Exercise base R's actual CSV parser without optional model packages.
            expression = (
                'd <- read.csv(' + json.dumps(command[2]) +
                ', colClasses="character", na.strings=' + json.dumps(marker) + '); '
                'stopifnot(is.na(d$group[2]), is.na(d$covariate[3]), '
                'is.na(d$subject[4]), sum(complete.cases(d)) == 2, '
                'd$group[5] == "NA", d$covariate[5] == "", d$subject[1] == "001")'
            )
            checked = real_run([executable, "-e", expression], capture_output=True, text=True)
            self.assertEqual(checked.returncode, 0, checked.stderr)
            Path(command[5]).write_text('feature,p_value\ngene,0.01\n')
            Path(command[6]).write_text('')
            Path(command[7]).write_text('{}')
            return subprocess.CompletedProcess(command, 0, '', '')

        with patch('r_runner.subprocess.run', side_effect=fake_r):
            run_r_analysis(metadata, pd.DataFrame({"gene": [1, 2, 3, 4, 5]}), request)
        pd.testing.assert_frame_equal(metadata, original)
        self.assertNotIn("transfer_na_marker", request)

    def test_result_identifiers_survive_transfer(self):
        def fake_r(command, **kwargs):
            Path(command[5]).write_text(
                'feature,predictor,term,estimate,p_value\n001,002,003,1.5,0.01\n')
            Path(command[6]).write_text(
                'feature,predictor,contrast,estimate,p.value\n001,002,004,2.5,0.02\n')
            Path(command[7]).write_text('{}')
            return subprocess.CompletedProcess(command, 0, '', '')

        with patch('r_runner.subprocess.run', side_effect=fake_r):
            output = run_r_analysis(pd.DataFrame({'SampleID': ['a']}),
                                    pd.DataFrame({'SampleID': ['a'], '001': [1]}), {})
        self.assertEqual(output['overall'].iloc[0]['feature'], '001')
        self.assertEqual(output['overall'].iloc[0]['predictor'], '002')
        self.assertEqual(output['overall'].iloc[0]['term'], '003')
        self.assertEqual(output['pairwise'].iloc[0]['contrast'], '004')
        self.assertEqual(output['pairwise'].iloc[0]['feature'], '001')
        self.assertEqual(output['overall'].iloc[0]['estimate'], 1.5)
        self.assertEqual(output['pairwise'].iloc[0]['p.value'], 0.02)


class RTransferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")
        executable = shutil.which("Rscript", path=path)
        if not executable:
            raise unittest.SkipTest("Rscript is required for integration checks")
        check = subprocess.run([executable, "-e",
            'packages <- c("lme4", "lmerTest", "emmeans", "jsonlite"); quit(status=if (all(vapply(packages, requireNamespace, logical(1), quietly=TRUE))) 0 else 1)'],
            capture_output=True, text=True)
        if check.returncode:
            raise unittest.SkipTest("R integration checks require lme4, lmerTest, emmeans, and jsonlite")

    def test_numeric_missing_values_use_complete_observations(self):
        import numpy as np
        rng = np.random.default_rng(103)
        metadata = pd.DataFrame({'SampleID': [f'S{i:03}' for i in range(24)],
                                 'x': rng.normal(size=24),
                                 'age': rng.normal(size=24)})
        features = pd.DataFrame({'SampleID': metadata.SampleID,
                                 'gene': 2 * metadata.x + rng.normal(size=24)})
        metadata.loc[0, 'age'] = np.nan
        features.loc[1, 'gene'] = np.nan
        originals = [table.copy(deep=True) for table in (metadata, features)]
        inspector = DataInspector(metadata, features, {
            'analysis_type': 'cross_sectional', 'primary_predictors': ['x'],
            'covariates': ['age'], 'variable_types': {'x': 'numeric', 'age': 'numeric'}})
        self.assertNotEqual(inspector.inspect()['status'], 'failed')
        path = str(Path(sys.executable).parent) + os.pathsep + os.environ.get('PATH', '')
        with patch.dict(os.environ, {'PATH': path}):
            output = run_r_analysis(**inspector.get_analysis_inputs())
        row = output['overall'].iloc[0]
        self.assertEqual(row['status'], 'success')
        self.assertEqual(row['n_obs'], 22)
        self.assertTrue(pd.notna(row['p_value']))
        for actual, original in zip((metadata, features), originals):
            pd.testing.assert_frame_equal(actual, original)

    def test_two_visit_choice_runs_intercepts_and_reports_unsupported_slopes(self):
        import numpy as np
        rng = np.random.default_rng(38)
        metadata = pd.DataFrame({"SampleID": [f"S{i}" for i in range(60)],
                                 "subject": np.repeat(np.arange(30), 2),
                                 "time": np.tile([0., 1.], 30),
                                 "x": np.repeat(rng.normal(size=30), 2)})
        features = pd.DataFrame({"SampleID": metadata.SampleID,
                                 "gene": np.repeat(rng.normal(scale=3, size=30), 2)
                                         + metadata.x + 0.5 * metadata.time + rng.normal(size=60)})
        request = {"analysis_type": "longitudinal", "analysis_goal": "time_effect",
                   "primary_predictors": ["x"], "covariates": [],
                   "subject_id": "subject", "time": "time",
                   "variable_types": {"x": "numeric"}, "random_slope": False}
        path = str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")
        for random_slope in (False, True):
            request["random_slope"] = random_slope
            inspector = DataInspector(metadata, features, request)
            report = inspector.inspect()
            self.assertNotEqual(report["status"], "failed")
            with patch.dict(os.environ, {"PATH": path}):
                output = run_r_analysis(**inspector.get_analysis_inputs())
            self.assertEqual(output["summary"]["final_random_slope"], random_slope)
            self.assertFalse(output["summary"]["rerun"])
            if random_slope:
                self.assertEqual(output["summary"]["n_failed"], 1)
                self.assertIn("Select random intercept only", output["overall"].iloc[0]["error"])
            else:
                self.assertEqual(output["summary"]["n_converged"], 1)
                self.assertTrue(pd.notna(output["overall"].iloc[0]["p_value"]))

    def test_rank_failure_is_per_feature_after_inspection(self):
        path = str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")
        metadata = pd.DataFrame({"SampleID": [f"S{i}" for i in range(18)],
                                 "x": list(range(18)),
                                 "adjustment": list(range(9)) + list(range(18, 9, -1))})
        outcomes = pd.DataFrame({"SampleID": metadata.SampleID,
                                 "good": [1., 3., 4., 2., 4., 3.] * 3,
                                 "bad": [1., 3., 4.] * 3 + [float("nan")] * 9})
        inspector = DataInspector(metadata, outcomes, {
            "analysis_type": "cross_sectional", "primary_predictors": ["x"],
            "covariates": ["adjustment"],
            "variable_types": {"x": "numeric", "adjustment": "numeric"},
        })
        self.assertNotEqual(inspector.inspect()["status"], "failed")
        with patch.dict(os.environ, {"PATH": path}):
            output = run_r_analysis(**inspector.get_analysis_inputs())
        rows = output["overall"].set_index("feature")
        self.assertTrue(pd.notna(rows.loc["good", "p_value"]))
        self.assertTrue(pd.isna(rows.loc["bad", "p_value"]))
        self.assertIn("rank deficient", rows.loc["bad", "error"])
        self.assertEqual(output["summary"]["n_failed"], 1)

    def test_numeric_ordered_labels_survive_csv_round_trip(self):
        executable_directory = str(Path(sys.executable).parent)
        path = executable_directory + os.pathsep + os.environ.get("PATH", "")
        if not shutil.which("Rscript", path=path):
            self.skipTest("Rscript is required for the integration check")
        metadata = pd.DataFrame({"SampleID": [f"S{i}" for i in range(18)],
                                 "severity": [1., 2., 3.] * 6})
        outcomes = pd.DataFrame({"SampleID": metadata.SampleID,
                                 "gene": [1., 3., 4., 2., 4., 3.] * 3})
        request = {"analysis_type": "cross_sectional",
                   "primary_predictors": ["severity"], "covariates": [],
                   "variable_types": {"severity": "ordered_categorical"},
                   "ordered_levels": {"severity": ["1.0", "2.0", "3.0"]}}
        inspector = DataInspector(metadata, outcomes, request)
        self.assertNotEqual(inspector.inspect()["status"], "failed")
        original = metadata.copy(deep=True)
        with patch.dict(os.environ, {"PATH": path}):
            output = run_r_analysis(**inspector.get_analysis_inputs())
        row = output["overall"].iloc[0]
        self.assertEqual(row["n_obs"], 18)
        self.assertEqual(row["term"], "severity_trend")
        self.assertTrue(pd.notna(row["p_value"]))
        self.assertAlmostEqual(row["estimate"], 1.)
        pd.testing.assert_frame_equal(metadata, original)


if __name__ == "__main__":
    unittest.main()
