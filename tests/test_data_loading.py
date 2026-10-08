import sys
import unittest
from tempfile import TemporaryDirectory
from io import StringIO, BytesIO
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from data_loading import read_sample_csv, read_metadata_csv
from inspector import DataInspector


class DataLoadingTests(unittest.TestCase):
    def test_leading_blank_lines_are_rejected(self):
        for loader in (read_sample_csv, read_metadata_csv):
            for prefix in ("\n", "\r\n", " \t\n", "\ufeff\n"):
                for header in ("SampleID,gene", "SampleID,gene,gene"):
                    content = prefix + header + "\na,1,2\n"
                    with self.subTest(loader=loader.__name__, prefix=prefix, header=header):
                        for source in (StringIO(content), BytesIO(content.encode())):
                            with self.assertRaisesRegex(ValueError, "first line must contain"):
                                loader(source)
                            self.assertEqual(source.tell(), 0)
                        with TemporaryDirectory() as directory:
                            path = Path(directory) / "input.csv"
                            path.write_bytes(content.encode())
                            with self.assertRaisesRegex(ValueError, "first line must contain"):
                                loader(path)

    def test_bom_header_without_leading_blanks_is_valid(self):
        for loader in (read_sample_csv, read_metadata_csv):
            table = loader(BytesIO("\ufeffSampleID,gene\na,1\n".encode()))
            self.assertEqual(table.columns.tolist(), ["SampleID", "gene"])

    def test_duplicate_headers_are_rejected_in_both_tables(self):
        for loader in (read_sample_csv, read_metadata_csv):
            for header in ('SampleID,gene,gene', 'SampleID,SampleID,x',
                           'SampleID,"gene,one","gene,one"'):
                with self.subTest(loader=loader.__name__, header=header):
                    source = BytesIO((header + "\na,1,2\n").encode())
                    with self.assertRaisesRegex(ValueError, "Duplicate CSV column names"):
                        loader(source)
                    self.assertEqual(source.tell(), 0)

    def test_distinct_and_multiline_headers_remain_valid(self):
        source = BytesIO(b'SampleID,gene,gene.1,"gene\nname"\na,1,2,3\n')
        table = read_sample_csv(source)
        self.assertEqual(table.columns.tolist(), ['SampleID', 'gene', 'gene.1', 'gene\nname'])
        self.assertEqual(table.iloc[0]['gene.1'], 2)

    def test_sample_ids_retain_leading_zeros(self):
        table = read_sample_csv(StringIO("SampleID,x\n001,1\n1,2\n002,3\n"))
        self.assertEqual(table.SampleID.tolist(), ["001", "1", "002"])
        self.assertEqual(table.x.tolist(), [1, 2, 3])

    def test_literal_na_and_null_cells_block_both_uploads(self):
        for loader in (read_metadata_csv, read_sample_csv):
            for column in ("SampleID", "value"):
                for value in ("NA", "null", " NA ", " null "):
                    with self.subTest(loader=loader.__name__, column=column, value=value):
                        row = f"{value},1" if column == "SampleID" else f"a,{value}"
                        with self.assertRaises(ValueError) as error:
                            loader(StringIO("SampleID,value\n" + row + "\n"))
                        self.assertIn(column, str(error.exception))
                        self.assertIn("Replace NA and null with empty cells", str(error.exception))

    def test_empty_cells_remain_missing_and_other_labels_remain_valid(self):
        metadata = read_metadata_csv(StringIO(
            "SampleID,group,SubjectID\na,NAD,null_subject\nb,,\n"))
        self.assertEqual(metadata.group.iloc[0], "NAD")
        self.assertEqual(metadata.SubjectID.iloc[0], "null_subject")
        self.assertTrue(metadata.group.isna().iloc[1])
        self.assertTrue(metadata.SubjectID.isna().iloc[1])
        features = read_sample_csv(StringIO("SampleID,gene\na,1\nb,\n"))
        self.assertTrue(features.gene.isna().iloc[1])

    def test_subject_labels_and_numeric_inspection(self):
        metadata = read_metadata_csv(StringIO(
            "SampleID,SubjectID,time,x\na,001,0,1\nb,001,1,2\n"
            "c,1,0,3\nd,1,1,4\n"))
        features = read_sample_csv(StringIO("SampleID,gene\na,1\nb,3\nc,2\nd,4\n"))
        inspector = DataInspector(metadata, features, {
            "analysis_type": "longitudinal", "primary_predictors": ["x"],
            "subject_id": "SubjectID", "time": "time", "covariates": [],
            "analysis_goal": "time_effect", "variable_types": {"x": "auto"},
        })
        report = inspector.inspect()
        self.assertNotEqual(report["status"], "failed")
        self.assertEqual(report["longitudinal"]["n_subjects"], 2)
        inputs = inspector.get_analysis_inputs()
        self.assertEqual(inputs["metadata"].SubjectID.tolist(), ["001", "001", "1", "1"])
        self.assertEqual(inputs["request"]["variable_types"]["x"], "numeric")
        self.assertEqual(inputs["metadata"].x.tolist(), [1, 2, 3, 4])

    def test_distinct_ids_do_not_match(self):
        metadata = read_sample_csv(StringIO("SampleID,x\n001,1\n002,2\n003,3\n"))
        features = read_sample_csv(StringIO("SampleID,gene\n1,2\n2,3\n3,4\n"))
        report = DataInspector(metadata, features, {
            "analysis_type": "cross_sectional", "primary_predictors": ["x"],
            "covariates": [], "variable_types": {"x": "numeric"},
        }).inspect()
        self.assertEqual(report["status"], "failed")
        self.assertTrue(any("no shared SampleID" in e for e in report["errors"]))


if __name__ == "__main__":
    unittest.main()
