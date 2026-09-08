#!/usr/bin/env python3
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from phase6_common import (
    RETROSPECTIVE_PROVENANCE_LIMIT, SCORER_CORRECTIONS, ReportHtmlParser,
    sample_id_from_path, validate_report_invariants, verify_manifest_identity,
    verify_scoring_provenance,
)
from phase6_scorer import (
    adjudicate, fasta_labels, percentile, performance_hard_gates,
    verified_internal_labels, wilson,
)


class Phase6Tests(unittest.TestCase):
    @staticmethod
    def target(representative, evidence_status, attribution_status="RESOLVED_TO_REFERENCE_GROUP"):
        return {
            "target_group_id": representative, "representative_id": representative,
            "member_ids": representative, "evidence_status": evidence_status,
            "attribution_status": attribution_status,
            "supporting_selected_fragments": "2", "selected_fragment_denominator": "10",
            "attributed_fragment_fraction": "0.2", "interval_lower": "0.1",
            "interval_upper": "0.3", "covered_bases": "20", "coverage_fraction": "0.1",
            "occupied_windows": "2", "host_confounded_fragments": "0",
            "cross_group_ambiguous_fragments": "0", "integration_status": "NONE",
            "split_events": "0", "discordant_fragments": "0",
        }

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

    def test_manifest_digest_and_frozen_identity_are_both_enforced(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text(json.dumps({"frozen_evaluation_digest": "frozen-a"}))
            prepared = {
                "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "frozen_evaluation_digest": "frozen-a",
            }
            verify_manifest_identity(path, prepared)
            path.write_text(json.dumps({"frozen_evaluation_digest": "frozen-b"}))
            with self.assertRaisesRegex(ValueError, "manifest digest changed"):
                verify_manifest_identity(path, prepared)
            prepared["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError, "frozen evaluation identity changed"):
                verify_manifest_identity(path, prepared)

    def test_scoring_provenance_rejects_input_drift(self):
        with tempfile.TemporaryDirectory() as directory:
            inputs = {"scorer": {"sha256": "a" * 64}}
            ledger = {
                "schema_id": "viroflash.phase6.scoring-provenance.v1",
                "contract_id": "contract",
                "mode": "RETROSPECTIVE_PRE_RESCORING",
                "created_at": "now",
                "inputs": inputs,
                "prior_evaluation_only_corrections": list(SCORER_CORRECTIONS),
                "provenance_limitation": RETROSPECTIVE_PROVENANCE_LIMIT,
            }
            Path(directory, "scoring-provenance.json").write_text(json.dumps(ledger))
            with patch("phase6_common.scoring_input_records", return_value=inputs):
                verify_scoring_provenance(directory, {}, {"contract_id": "contract"})
            with patch("phase6_common.scoring_input_records", return_value={"scorer": {"sha256": "b" * 64}}):
                with self.assertRaisesRegex(ValueError, "scoring input digest ledger changed"):
                    verify_scoring_provenance(directory, {}, {"contract_id": "contract"})

    def test_taxon_matching_starts_at_description_and_observes_boundary(self):
        terms = {
            "EBV": ["Epstein-Barr virus", "Human gammaherpesvirus 4"],
            "HBV": ["Hepatitis B virus"],
            "HPV16": ["Human papillomavirus type 16"],
            "HPV18": ["Human papillomavirus type 18", "Human papillomavirus 18"],
        }
        fasta = """>hbv |Hepatitis B virus isolate human\nA\n>heron |UNVERIFIED: Heron hepatitis B virus isolate 32\nA\n>ebv |Human gammaherpesvirus 4, complete genome\nA\n>hpv16 |Human papillomavirus type 16 isolate x\nA\n>hpv161 |Human papillomavirus type 161 isolate x\nA\n>hpv18 |Human papillomavirus 18 isolate x\nA\n"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "target.fa"
            path.write_text(fasta)
            self.assertEqual(
                fasta_labels(path, terms),
                {"hbv": "HBV", "ebv": "EBV", "hpv16": "HPV16", "hpv18": "HPV18"},
            )

    def test_internal_labels_require_frozen_fasta_bytes_and_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "target.fa"
            path.write_text(">hbv |Hepatitis B virus isolate human\nA\n")
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            record = {"path": str(path), "compressed_bytes": path.stat().st_size, "sha256": digest}
            manifest = {"runs": [{"dataset_id": "internal-68", "target_reference": record}]}
            config = {"internal_label_header_terms": {"HBV": ["Hepatitis B virus"]}}
            self.assertEqual(verified_internal_labels(manifest, config), {"hbv": "HBV"})
            path.write_text(path.read_text() + "A\n")
            with self.assertRaisesRegex(ValueError, "bytes differ"):
                verified_internal_labels(manifest, config)

    def test_report_invariants_bind_index_and_public_arithmetic(self):
        run = {
            "input_fragments": "1000", "selected_fragments": "1000",
            "prescreen_passed_fragments": "10", "aligned_fragments": "8",
            "unassigned_fragments": "3", "target_family_size": "1",
            "selection_probability": "1", "minimum_relevant_fraction": "0.001",
            "familywise_miss_probability": "0.05", "interval_level": "0.95",
            "analysis_status": "CONFORMANT_COMPLETE", "reason_codes": "",
        }
        target = {
            "target_group_id": "group", "representative_id": "representative",
            "member_ids": "representative;member", "representative_length": "100",
            "supporting_selected_fragments": "2", "selected_fragment_denominator": "1000",
            "attributed_fragment_fraction": "0.002", "interval_lower": "0.001",
            "interval_upper": "0.003", "interval_level": "0.95",
            "interval_method": "EQUAL_TAILED_EXACT_HYPERGEOMETRIC_INVERSION",
            "estimated_input_supporting_fragments": "2", "covered_bases": "20",
            "coverage_fraction": "0.2", "evidence_status": "REFERENCE_SIGNAL_OBSERVED",
            "attribution_status": "RESOLVED_TO_REFERENCE_GROUP",
            "occupied_windows": "1", "host_confounded_fragments": "0",
            "cross_group_ambiguous_fragments": "0", "split_events": "0",
            "discordant_fragments": "0",
        }
        index = {
            "target_family_size": 1,
            "groups": {"group": {"ordinal": 0, "representative_id": "representative", "member_ids": ["representative", "member"], "representative_length": 100}},
        }
        validate_report_invariants(run, [target], index)
        cases = (
            ("target membership", target, "target_group_id", "other", "absent from index ledger"),
            ("fraction arithmetic", target, "attributed_fragment_fraction", "0.2", "fraction arithmetic mismatch"),
            ("coverage arithmetic", target, "coverage_fraction", "0.3", "coverage fraction arithmetic mismatch"),
            ("status coherence", run, "analysis_status", "CONFORMANT_WITH_LIMITATIONS", "status and run limitations"),
        )
        for name, source, field, value, message in cases:
            with self.subTest(name=name):
                changed_run = copy.deepcopy(run)
                changed_target = copy.deepcopy(target)
                destination = changed_run if source is run else changed_target
                destination[field] = value
                with self.assertRaisesRegex(ValueError, message):
                    validate_report_invariants(changed_run, [changed_target], index)

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

    def test_any_observed_matching_reference_group_is_concordant(self):
        run = {
            "expectation_kind": "EXPECTED_GROUP", "expected_group_key": "EBV",
            "dataset_id": "internal-68", "evaluability_status": "EVALUABLE",
        }
        targets = [
            self.target("ebv-indeterminate", "INDETERMINATE", "INDETERMINATE"),
            self.target("ebv-observed", "REFERENCE_SIGNAL_OBSERVED"),
        ]
        decision = adjudicate(
            run, {"targets": targets},
            {"ebv-indeterminate": "EBV", "ebv-observed": "EBV"},
        )
        self.assertEqual(decision["classification"], "LABEL_CONCORDANT_SIGNAL")
        self.assertTrue(decision["expected_observed"])
        self.assertEqual(decision["expected_signals"], targets)
        self.assertEqual(decision["observed_expected_signals"], [targets[1]])
        self.assertFalse(decision["wrong_group_resolved"])

    def test_mock_signal_is_unexpected_but_not_wrong_group(self):
        run = {
            "expectation_kind": "MOCK", "expected_group_key": None,
            "dataset_id": "external-59", "evaluability_status": "EVALUABLE",
        }
        decision = adjudicate(
            run, {"targets": [self.target("unexpected", "REFERENCE_SIGNAL_OBSERVED")]}, {},
        )
        self.assertEqual(decision["classification"], "LABEL_DISCORDANT_WITH_SEQUENCE_EVIDENCE")
        self.assertEqual(decision["observed_signal_count"], 1)
        self.assertEqual(decision["label_scope_observed_signal_count"], 1)
        self.assertFalse(decision["wrong_group_resolved"])

    def test_internal_mock_ignores_off_scope_signal_for_label_adjudication(self):
        run = {
            "expectation_kind": "MOCK", "expected_group_key": None,
            "dataset_id": "internal-68", "evaluability_status": "EVALUABLE",
        }
        off_scope = self.target("off-panel-virus", "REFERENCE_SIGNAL_OBSERVED")
        decision = adjudicate(run, {"targets": [off_scope]}, {"ebv": "EBV"})
        self.assertEqual(decision["classification"], "LABEL_CONCORDANT_NO_SIGNAL")
        self.assertEqual(decision["observed_signal_count"], 1)
        self.assertEqual(decision["label_scope_observed_signal_count"], 0)
        self.assertFalse(decision["wrong_group_resolved"])
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
