"""Verify summary output without requiring a running Streamlit app."""
import ast
from pathlib import Path
import unittest
from unittest.mock import MagicMock


class ModelSummaryTests(unittest.TestCase):
    def render(self, **summary):
        source = Path(__file__).resolve().parents[1] / 'src' / 'results.py'
        tree = ast.parse(source.read_text())
        function = next(node for node in tree.body
                        if isinstance(node, ast.FunctionDef)
                        and node.name == 'display_model_summary')
        st = MagicMock()
        namespace = {'st': st}
        exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), 'exec'),
             namespace)
        namespace['display_model_summary'](summary, {'analysis_type': 'longitudinal'})
        return st, [call.args[0] for call in st.write.call_args_list]

    def test_fallback_shows_initial_and_final_fractions(self):
        st, output = self.render(rerun=True, first_pass_singular_fraction=0.4,
                                 singular_fraction=0.05, final_random_slope=False)
        self.assertIn('**Initial random-slope singularity:** 40.0%', output)
        self.assertIn('**Final random-intercept singularity:** 5.0%', output)
        self.assertIn('fallback threshold', st.warning.call_args.args[0])

    def test_no_fallback_shows_only_selected_model_fraction(self):
        for slope in (True, False):
            with self.subTest(slope=slope):
                st, output = self.render(rerun=False, singular_fraction=0,
                                         final_random_slope=slope)
                self.assertIn('**Singularity fraction:** 0.0%', output)
                self.assertFalse(any('Initial random-slope singularity' in line for line in output))
                st.warning.assert_not_called()

    def test_missing_final_fraction_does_not_display_initial_as_final(self):
        _, output = self.render(rerun=True, first_pass_singular_fraction=0.4,
                                singular_fraction=None)
        self.assertIn('**Final random-intercept singularity:** Unavailable', output)
        self.assertNotIn('**Final random-intercept singularity:** 40.0%', output)


if __name__ == '__main__':
    unittest.main()
