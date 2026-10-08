import csv
import os
from collections import Counter

import pandas as pd


def _validate_header(source):
    """Check original names before pandas can rename duplicate columns."""
    if isinstance(source, (str, os.PathLike)):
        with open(source, encoding="utf-8-sig", newline="") as handle:
            header = next(csv.reader(handle), [])
    else:
        position = source.tell()
        try:
            def lines():
                first = True
                while True:
                    line = source.readline()
                    if not line:
                        break
                    if isinstance(line, bytes):
                        line = line.decode("utf-8")
                    if first:
                        line = line.lstrip("\ufeff")
                        first = False
                    yield line
            header = next(csv.reader(lines()), [])
        finally:
            source.seek(position)
    if not header or (len(header) == 1 and not header[0].strip()):
        raise ValueError("The first line must contain the CSV column headers; leading blank lines are not allowed.")
    duplicates = [name for name, count in Counter(header).items() if count > 1]
    if duplicates:
        raise ValueError("Duplicate CSV column names: " + ", ".join(repr(name) for name in duplicates))


def _validate_missing_labels(table, table_name):
    """Require users to resolve ambiguous missing-value text before modeling."""
    affected = []
    for column in table.columns:
        values = table[column].astype("string").str.strip()
        count = int(values.isin(["NA", "null"]).sum())
        if count:
            affected.append(f"{column!r} ({count} cells)")
    if affected:
        raise ValueError(
            f"{table_name} contains literal 'NA' or 'null' values in columns: "
            + ", ".join(affected)
            + ". Replace NA and null with empty cells if they should be treated "
              "as missing data. If they are intentional labels or identifiers, "
              "rename them before uploading."
        )
    return table


def read_sample_csv(source):
    """Load a sample table while preserving the exact SampleID text."""
    _validate_header(source)
    table = pd.read_csv(source, dtype={"SampleID": "string"},
                        keep_default_na=False, na_values=[""])
    return _validate_missing_labels(table, "Feature matrix")


def read_metadata_csv(source):
    """Preserve metadata labels; inspection converts numeric model variables."""
    _validate_header(source)
    table = pd.read_csv(source, dtype="string", keep_default_na=False, na_values=[""])
    return _validate_missing_labels(table, "Metadata")
