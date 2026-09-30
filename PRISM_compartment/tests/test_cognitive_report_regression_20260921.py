"""The negative label contains 'dementia'; it must be checked first."""
import importlib.util
from pathlib import Path
import sys
import unittest

p = Path(__file__).parents[1] / "src/kmlee_bam/preprocessing/data/make_pathology_aware_donor_split.py"
spec = importlib.util.spec_from_file_location("cognition_split_module", p)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


class CognitionReportTest(unittest.TestCase):
    def test_negative_before_positive_substring(self):
        for label, expected in [("No dementia", 0), ("Dementia", 1), ("Reference", 0), ("MCI", 1), (None, None), ("Unknown", None)]:
            with self.subTest(label=label):
                self.assertEqual(module.cognitive_impaired(label), expected)

    def test_is_report_only(self):
        import pandas as pd
        df = pd.DataFrame({"Cognitive status": ["No dementia", "Dementia"]}, index=["a", "b"])
        cfg = module.PathologyAwareSplitConfigV4(donor_pathology_table_path='unused_for_unit_test')
        labels, primary, weak, report = module.build_labels(df, cfg)
        name = "axis_Cognitive_impaired_report"
        self.assertEqual(list(labels[name]), [0, 1])
        self.assertIn(name, report)
        self.assertNotIn(name, primary + weak)
