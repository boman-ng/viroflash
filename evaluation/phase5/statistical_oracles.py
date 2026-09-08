#!/usr/bin/env python3
import argparse
import itertools
import json
import math
import random
from decimal import Decimal, getcontext
from fractions import Fraction
from functools import cache
from pathlib import Path


getcontext().prec = 80


def load_profile(path):
    with path.open(encoding="utf-8") as handle:
        profile = json.load(handle, parse_float=Decimal)
    parameters = profile["parameters"]
    delta = parameters["minimum_relevant_fraction"]["value"]
    beta = parameters["familywise_miss_probability"]["value"]
    alpha = parameters["familywise_interval_error"]["value"]
    assert delta == Decimal("0.00001")
    assert beta == Decimal("0.05")
    assert alpha == Decimal("0.05")
    return delta, beta, alpha


def enumerate_sampling_power(delta, beta):
    cases = []
    for population, family_size in [(80_000, 1), (250_000, 4), (1_000_000, 20)]:
        relevant = math.ceil(Decimal(population) * delta)
        group_miss = beta / family_size
        probability = Decimal(1) - (group_miss.ln() / relevant).exp()
        while (Decimal(1) - probability) ** relevant > group_miss:
            probability = probability.next_plus()
        detection = Decimal(0)
        for pattern in itertools.product((False, True), repeat=relevant):
            selected = sum(pattern)
            weight = probability**selected * (Decimal(1) - probability) ** (
                relevant - selected
            )
            if selected:
                detection += weight
        required = Decimal(1) - group_miss
        assert detection >= required
        cases.append(
            {
                "population": population,
                "target_family_size": family_size,
                "minimum_relevant_fragments": relevant,
                "selection_probability": str(probability),
                "enumerated_detection_probability": str(detection),
                "required_detection_probability": str(required),
            }
        )
    return cases


def enumerate_estimator_bias():
    cases = []
    population = 8
    for sample in (1, 3, 5, 8):
        for total_successes in range(population + 1):
            values = []
            labels = [index < total_successes for index in range(population)]
            for selected in itertools.combinations(range(population), sample):
                observed = sum(labels[index] for index in selected)
                values.append(Fraction(observed, sample))
            mean = sum(values, Fraction()) / len(values)
            truth = Fraction(total_successes, population)
            assert mean == truth
            cases.append(
                {
                    "population": population,
                    "sample": sample,
                    "total_successes": total_successes,
                    "enumerated_mean": str(mean),
                    "truth": str(truth),
                    "bias": "0",
                }
            )
    return cases


def hypergeometric_probability(population, total_successes, sample, observed):
    failures = population - total_successes
    if observed < 0 or observed > total_successes or sample - observed > failures:
        return Fraction()
    return Fraction(
        math.comb(total_successes, observed)
        * math.comb(failures, sample - observed),
        math.comb(population, sample),
    )


@cache
def exact_interval(population, sample, observed, alpha):
    accepted = []
    for total_successes in range(observed, population - (sample - observed) + 1):
        lower_tail = sum(
            (
                hypergeometric_probability(
                    population, total_successes, sample, candidate
                )
                for candidate in range(observed + 1)
            ),
            Fraction(),
        )
        upper_tail = sum(
            (
                hypergeometric_probability(
                    population, total_successes, sample, candidate
                )
                for candidate in range(observed, sample + 1)
            ),
            Fraction(),
        )
        if lower_tail >= alpha / 2 and upper_tail >= alpha / 2:
            accepted.append(total_successes)
    assert accepted
    return min(accepted), max(accepted)


def category_compositions(total, categories):
    if categories == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for remainder in category_compositions(total - first, categories - 1):
            yield (first,) + remainder


def enumerate_simultaneous_coverage(alpha):
    population = 12
    sample = 6
    group_alpha = Fraction(alpha) / 2
    minimum = Fraction(1)
    minimum_composition = None
    checked = 0
    for counts in category_compositions(population, 4):
        masks = []
        for mask, count in enumerate(counts):
            masks.extend([mask] * count)
        group_totals = (
            sum(mask & 1 != 0 for mask in masks),
            sum(mask & 2 != 0 for mask in masks),
        )
        covered = 0
        selections = 0
        for selected in itertools.combinations(range(population), sample):
            observed = (
                sum(masks[index] & 1 != 0 for index in selected),
                sum(masks[index] & 2 != 0 for index in selected),
            )
            intervals = tuple(
                exact_interval(population, sample, value, group_alpha)
                for value in observed
            )
            covered += all(
                lower <= truth <= upper
                for truth, (lower, upper) in zip(group_totals, intervals)
            )
            selections += 1
        coverage = Fraction(covered, selections)
        assert coverage >= Fraction(1) - Fraction(alpha)
        if coverage < minimum:
            minimum = coverage
            minimum_composition = counts
        checked += 1
    return {
        "population": population,
        "sample": sample,
        "groups": 2,
        "population_compositions_checked": checked,
        "minimum_simultaneous_coverage": str(minimum),
        "required_simultaneous_coverage": str(Fraction(1) - Fraction(alpha)),
        "minimum_coverage_category_counts_00_01_10_11": minimum_composition,
    }


def deterministic_simulation(alpha):
    def run():
        generator = random.Random(0x5649524F464C4153)
        population = 100
        sample = 20
        trials = 20_000
        labels = [(index < 7, index % 11 == 0) for index in range(population)]
        estimates = [Fraction(), Fraction()]
        covered = 0
        group_alpha = Fraction(alpha) / 2
        truths = [sum(label[group] for label in labels) for group in range(2)]
        for _ in range(trials):
            selected = generator.sample(range(population), sample)
            observed = [
                sum(labels[index][group] for index in selected) for group in range(2)
            ]
            for group in range(2):
                estimates[group] += Fraction(observed[group], sample)
            intervals = [
                exact_interval(population, sample, value, group_alpha)
                for value in observed
            ]
            covered += all(
                lower <= truth <= upper
                for truth, (lower, upper) in zip(truths, intervals)
            )
        return {
            "seed": "0x5649524f464c4153",
            "population": population,
            "sample": sample,
            "trials": trials,
            "observed_bias": [
                str(estimates[group] / trials - Fraction(truths[group], population))
                for group in range(2)
            ],
            "observed_simultaneous_coverage": str(Fraction(covered, trials)),
            "acceptance_role": "measurement_only_no_tolerance",
        }

    first = run()
    assert first == run()
    return first


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--profile",
        type=Path,
        default=Path("evaluation/phase0/analysis-profile.json"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    delta, beta, alpha = load_profile(args.profile)
    result = {
        "profile": {"delta": str(delta), "beta": str(beta), "alpha": str(alpha)},
        "sampling_power_enumeration": enumerate_sampling_power(delta, beta),
        "estimator_bias_enumeration": enumerate_estimator_bias(),
        "simultaneous_interval_coverage_enumeration": enumerate_simultaneous_coverage(
            Fraction(alpha)
        ),
        "deterministic_simulation": deterministic_simulation(Fraction(alpha)),
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
