#!/usr/bin/env python3
import unittest

from phase6_common import ReportHtmlParser, sample_id_from_path
from phase6_scorer import adjudicate, percentile, performance_hard_gates, wilson


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

    def test_non_evaluable_manifest_status_overrides_label_classification(self):
        run = {
            "expectation_kind": "EXPECTED_GROUP",
            "expected_group_key": "K02718.1",
            "dataset_id": "external-59",
            "evaluability_status": "LABEL_SEQUENCE_AMBIGUITY",
        }
        target = {
            "target_group_id": "group", "representative_id": "K02718.1",
            "member_ids": "K02718.1", "evidence_status": "REFERENCE_SIGNAL_OBSERVED",
            "attribution_status": "RESOLVED_TO_REFERENCE_GROUP",
            "supporting_selected_fragments": "2", "selected_fragment_denominator": "10",
            "attributed_fragment_fraction": "0.2", "interval_lower": "0.1",
            "interval_upper": "0.3", "covered_bases": "20", "coverage_fraction": "0.1",
            "occupied_windows": "2", "host_confounded_fragments": "0",
            "cross_group_ambiguous_fragments": "0", "integration_status": "NONE",
            "split_events": "0", "discordant_fragments": "0",
        }
        decision = adjudicate(run, {"targets": [target]}, {})
        self.assertEqual(decision["classification"], "NOT_EVALUABLE")
        self.assertTrue(decision["expected_observed"])
        self.assertEqual(len(decision["evidence_summary"]), 1)

    def test_performance_false_conditions_surface_in_hard_gates(self):
        performance = {
            "datasets": {"internal-68": {"new_samples_over_historical_max": ["sample"]}},
            "cohorts": {"internal": {"systematic_wall_regression": True}},
            "internal_index_build_memory_comparison": {"within_historical_max": True},
        }
        gates = performance_hard_gates(performance)
        self.assertFalse(gates["sample_peak_rss_within_corresponding_historical_max"])
        self.assertFalse(gates["no_unaccepted_systematic_wall_regression"])
        self.assertFalse(gates["performance_converged"])


if __name__ == "__main__":
    unittest.main()
