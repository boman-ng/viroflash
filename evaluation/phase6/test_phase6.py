#!/usr/bin/env python3
import unittest

from phase6_common import ReportHtmlParser, sample_id_from_path
from phase6_scorer import percentile, wilson


class Phase6Tests(unittest.TestCase):
    def test_wilson_zero_and_complete_counts_are_bounded(self):
        zero = wilson(0, 10)
        complete = wilson(10, 10)
        self.assertEqual(zero["rate"], 0)
        self.assertEqual(complete["rate"], 1)
        self.assertGreaterEqual(zero["lower"], 0)
        self.assertLessEqual(complete["upper"], 1)

    def test_percentile_uses_linear_interpolation(self):
        self.assertEqual(percentile([1, 2, 3], 0.5), 2)
        self.assertEqual(percentile([0, 10], 0.95), 9.5)

    def test_html_parser_separates_run_and_target_fields(self):
        parser = ReportHtmlParser()
        parser.feed('<section id="run-integrity"><table><tr data-field="sample_id"><td>a&amp;b</td></tr></table></section><section class="target-signal" data-group="g"><tr data-field="target_group_id"><td>g</td></tr></section>')
        self.assertEqual(parser.run, {"sample_id": "a&b"})
        self.assertEqual(parser.targets, [{"data_group": "g", "target_group_id": "g"}])

    def test_sample_id_matches_binary_suffix_order(self):
        self.assertEqual(sample_id_from_path("sample_R1.fastq.gz"), "sample")
        self.assertEqual(sample_id_from_path("SRR1.fastq.gz"), "SRR1")


if __name__ == "__main__":
    unittest.main()
