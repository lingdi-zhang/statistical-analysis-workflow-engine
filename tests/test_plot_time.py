"""Check rendered coordinates without requiring Streamlit or a GUI backend."""
import ast
from pathlib import Path
import unittest
from unittest.mock import MagicMock

import numpy as np
import pandas as pd


class LongitudinalTimeTests(unittest.TestCase):
    def test_unequal_visit_intervals_are_preserved(self):
        names = {'plot_longitudinal_categorical', 'plot_longitudinal_time_effect',
                 'plot_longitudinal_numeric_predictor', 'sort_time_values'}
        source = Path(__file__).resolve().parents[1] / 'src' / 'results.py'
        tree = ast.parse(source.read_text())
        functions = ast.Module(body=[node for node in tree.body
                                    if isinstance(node, ast.FunctionDef) and node.name in names],
                               type_ignores=[])
        data = pd.DataFrame({'time': [0, 1, 12], 'gene': [1., 2., 3.],
                             'group': ['control'] * 3, 'dose': [2., 3., 4.]})
        for name in sorted(names - {'sort_time_values'}):
            with self.subTest(plot=name):
                fig, ax = MagicMock(), MagicMock()
                plt = MagicMock()
                plt.subplots.return_value = (fig, ax)
                namespace = {'np': np, 'pd': pd, 'plt': plt, 'st': MagicMock(),
                             'get_feature_group_colors': lambda n: ['blue'] * n,
                             'add_pdf_download': MagicMock(), 'safe_file_name': str}
                exec(compile(functions, str(source), 'exec'), namespace)
                kwargs = dict(plot_data=data, feature='gene', time_variable='time')
                if name.endswith('categorical'):
                    kwargs['predictor'] = 'group'
                elif name.endswith('numeric_predictor'):
                    kwargs['predictor'] = 'dose'
                namespace[name](**kwargs)
                np.testing.assert_allclose(ax.set_xticks.call_args.args[0], [0, 1, 12])
                if name.endswith('numeric_predictor'):
                    np.testing.assert_allclose(ax.scatter.call_args.args[0], [0, 1, 12])
                else:
                    np.testing.assert_allclose(ax.plot.call_args.args[0], [0, 1, 12])


if __name__ == '__main__':
    unittest.main()
