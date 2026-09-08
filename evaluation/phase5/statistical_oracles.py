#!/usr/bin/env python3
import json
import subprocess


TESTS = [
    "sampling_design::tests::phase5_production_probability_meets_exact_rational_miss_budget",
    "sampling_design::tests::phase5_bernoulli_selection_conditioned_on_realized_n_is_uniform",
    "evidence::tests::phase5_production_endpoints_bias_and_family_coverage_match_exact_enumeration",
]


def main():
    results = []
    for test in TESTS:
        command = ["cargo", "test", test, "--locked"]
        completed = subprocess.run(command, capture_output=True, text=True)
        results.append(
            {
                "test": test,
                "returncode": completed.returncode,
                "status": "COMPLETE" if completed.returncode == 0 else "FAILED",
            }
        )
        if completed.returncode:
            raise SystemExit(
                json.dumps(
                    {
                        "results": results,
                        "stdout": completed.stdout,
                        "stderr": completed.stderr,
                    },
                    indent=2,
                )
            )
    print(json.dumps({"production_connected_exact_oracles": results}, indent=2))


if __name__ == "__main__":
    main()
