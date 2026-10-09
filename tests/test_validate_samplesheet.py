import unittest
from pathlib import Path
from scripts.validate_samplesheet import validate_samplesheet
ROOT=Path(__file__).resolve().parents[1]
class SamplesheetTests(unittest.TestCase):
    def test_valid_fixture_passes(self):
        self.assertEqual(validate_samplesheet(ROOT/"tests/fixtures/samplesheet_valid.csv"),[])
    def test_missing_columns_reported(self):
        self.assertTrue(any("Missing required column" in e for e in validate_samplesheet(ROOT/"tests/fixtures/samplesheet_invalid.csv")))
    def test_duplicate_sample_reported(self):
        src=ROOT/"tests/fixtures/samplesheet_valid.csv"; tmp=ROOT/"tests/fixtures/.duplicate.tmp.csv"
        try:
            tmp.write_text(src.read_text().replace("P02_T,","P01_T,"))
            self.assertTrue(any("duplicate sample ID" in e for e in validate_samplesheet(tmp)))
        finally: tmp.unlink(missing_ok=True)
if __name__=="__main__": unittest.main()
