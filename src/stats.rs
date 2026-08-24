//! Statistical models and multiple-testing adjustment for candidate evidence.
//!
//! The primary path uses a one-sided exact conditional two-Poisson rate test with reference length
//! as exposure, retaining log probabilities through floating-point underflow. Every fixed discovery
//! hypothesis enters a Benjamini-Hochberg running-min adjustment. Negative-binomial and pi0 helpers
//! remain diagnostic and do not silently replace the primary model.

const FPMIN: f64 = 1e-300;
const ITMAX: usize = 200;
const EPS: f64 = 3e-14;

/// log Gamma(x) using the Lanczos approximation from Numerical Recipes `gammaln`.
pub fn ln_gamma(x: f64) -> f64 {
    const COF: [f64; 6] = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5,
    ];
    let tmp = x + 5.5;
    let tmp = tmp - (x + 0.5) * tmp.ln();
    let mut ser = 1.000000000190015;
    for (i, c) in COF.iter().enumerate() {
        ser += c / (x + i as f64 + 1.0);
    }
    -tmp + (2.5066282746310005 * ser / x).ln()
}

fn gamma_p_series(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..ITMAX {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * EPS {
            break;
        }
    }
    sum * (-x + a * x.ln() - ln_gamma(a)).exp()
}

fn gamma_q_cf(a: f64, x: f64) -> f64 {
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / FPMIN;
    let mut d = 1.0 / b.max(FPMIN);
    let mut h = d;
    for i in 1..=ITMAX {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = b + an / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            break;
        }
    }
    (-x + a * x.ln() - ln_gamma(a)).exp() * h
}

pub fn regularized_gamma_p(a: f64, x: f64) -> f64 {
    if x < 0.0 || a <= 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        gamma_p_series(a, x)
    } else {
        1.0 - gamma_q_cf(a, x)
    }
}

pub fn poisson_upper_tail(c: u64, lambda: f64) -> f64 {
    if lambda <= 0.0 {
        return if c > 0 { 0.0 } else { 1.0 };
    }
    if c == 0 {
        return 1.0;
    }
    regularized_gamma_p(c as f64, lambda)
}

pub fn digamma(x: f64) -> f64 {
    if x < 6.0 {
        return digamma(x + 1.0) - 1.0 / x;
    }
    let r = 1.0 / x;
    let r2 = r * r;
    x.ln()
        - 0.5 * r
        - r2 * (1.0 / 12.0
            - r2 * (1.0 / 120.0
                - r2 * (1.0 / 252.0
                    - r2 * (1.0 / 240.0 - r2 * (1.0 / 132.0 - r2 * 691.0 / 32760.0)))))
}

pub fn trigamma(x: f64) -> f64 {
    if x < 6.0 {
        return trigamma(x + 1.0) + 1.0 / (x * x);
    }
    let r = 1.0 / x;
    let r2 = r * r;
    r + 0.5 * r2
        + r2 * r
            * (1.0 / 6.0
                - r2 * (1.0 / 30.0 - r2 * (1.0 / 42.0 - r2 * (1.0 / 30.0 - r2 * 5.0 / 66.0))))
}

pub fn nb_alpha_mom(mean: f64, var: f64) -> f64 {
    ((var - mean) / (mean * mean)).max(0.0)
}

/// `ell` is proportional to `sum(logGamma(v+r)) - n*logGamma(r) + n*r*log(r/(r+mu)) - sum(v)*log(r+mu)`.
/// `s2 = sum(trigamma(v+r)) + n*mu/(r*(r+mu)) - n*trigamma(r)`.
pub fn nb_mle_score_and_hessian(values: &[f64], mean: f64, r: f64) -> (f64, f64) {
    let n = values.len() as f64;
    let mut s1 = 0.0;
    let mut s2 = 0.0;
    for &v in values {
        s1 += digamma(v + r);
        s2 += trigamma(v + r);
    }
    s1 += n * (r / (r + mean)).ln() - n * digamma(r);
    s2 += n * mean / (r * (r + mean)) - n * trigamma(r);
    (s1, s2)
}

pub fn nb_alpha_mle(values: &[f64], init_alpha: f64) -> Option<f64> {
    let n = values.len();
    if n == 0 {
        return None;
    }
    let mean = values.iter().sum::<f64>() / n as f64;
    if mean <= 0.0 {
        return None;
    }
    let mut r = (1.0 / init_alpha.max(1e-6)).max(1e-4);
    for _ in 0..10 {
        let (s1, s2) = nb_mle_score_and_hessian(values, mean, r);
        if s2.abs() < 1e-12 {
            return None;
        }
        let step = s1 / s2;
        if step.abs() < 1e-6 * r {
            return Some(1.0 / r);
        }
        r -= step;
        if r < 1e-4 {
            return None;
        }
    }
    None
}

fn betacf(a: f64, b: f64, x: f64) -> f64 {
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=ITMAX {
        let m2 = 2.0 * m as f64;
        let mut aa = m as f64 * (b - m as f64) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m as f64) * (qab + m as f64) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 3e-7 {
            break;
        }
    }
    h
}

pub fn regularized_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * betacf(a, b, x) / a
    } else {
        1.0 - bt * betacf(b, a, 1.0 - x) / b
    }
}

pub fn nb_upper_tail(c: u64, mean: f64, alpha: f64) -> f64 {
    if mean <= 0.0 {
        return if c > 0 { 0.0 } else { 1.0 };
    }
    if alpha <= 1e-10 {
        return poisson_upper_tail(c, mean);
    }
    if c == 0 {
        return 1.0;
    }
    let r = 1.0 / alpha;
    let p = mean * alpha / (1.0 + mean * alpha);
    regularized_beta(p, c as f64, r)
}

pub fn trimmed_mean_20(values: &[f64]) -> Option<f64> {
    let n = values.len();
    if n == 0 {
        return None;
    }
    let mut sorted: Vec<f64> = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let trim = (n as f64 * 0.1).floor() as usize;
    let kept = &sorted[trim..n - trim];
    Some(kept.iter().sum::<f64>() / kept.len() as f64)
}

pub fn eb_shrink_layer(layer_mean: f64, layer_n: usize, global_mean: f64, kappa: f64) -> f64 {
    (kappa * global_mean + layer_n as f64 * layer_mean) / (kappa + layer_n as f64)
}

pub fn expected_hits(contig_len: u64, n_total_reads: u64) -> f64 {
    contig_len as f64 * (n_total_reads as f64 / 1e7)
}

/// Probability represented in ordinary and log space, with explicit underflow status.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TailProbability {
    pub probability: f64,
    pub ln_probability: f64,
    pub underflow: bool,
}

fn log_add_exp(left: f64, right: f64) -> f64 {
    if left == f64::NEG_INFINITY {
        return right;
    }
    if right == f64::NEG_INFINITY {
        return left;
    }
    let high = left.max(right);
    let low = left.min(right);
    high + (low - high).exp().ln_1p()
}

/// H0: `target_count / target_exposure == background_count / background_exposure`.
pub fn exact_poisson_rate_upper_tail(
    target_count: u64,
    target_exposure: f64,
    background_count: u64,
    background_exposure: f64,
) -> Result<TailProbability, String> {
    if !target_exposure.is_finite() || target_exposure <= 0.0 {
        return Err("target exposure must be finite and positive".to_string());
    }
    if !background_exposure.is_finite() || background_exposure <= 0.0 {
        return Err("background exposure must be finite and positive".to_string());
    }
    let total_count = target_count
        .checked_add(background_count)
        .ok_or_else(|| "Poisson conditional-test total-count overflow".to_string())?;
    if target_count == 0 || total_count == 0 {
        return Ok(TailProbability {
            probability: 1.0,
            ln_probability: 0.0,
            underflow: false,
        });
    }

    let target_probability = target_exposure / (target_exposure + background_exposure);
    if !(0.0..1.0).contains(&target_probability) {
        return Err("Poisson conditional-test exposure ratio is degenerate".to_string());
    }

    let direct = regularized_beta(
        target_probability,
        target_count as f64,
        background_count as f64 + 1.0,
    );
    let ln_probability = if direct.is_finite() && direct >= f64::MIN_POSITIVE {
        direct.ln()
    } else {
        let n = total_count;
        let k = target_count;
        let ln_p = target_probability.ln();
        let ln_one_minus_p = (-target_probability).ln_1p();
        let mut x = k;
        let mut ln_term =
            ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0)
                + k as f64 * ln_p
                + (n - k) as f64 * ln_one_minus_p;
        let mut ln_sum = ln_term;
        while x < n {
            ln_term += ((n - x) as f64).ln() - ((x + 1) as f64).ln() + ln_p - ln_one_minus_p;
            ln_sum = log_add_exp(ln_sum, ln_term);
            x += 1;
        }
        ln_sum.min(0.0)
    };
    if !ln_probability.is_finite() {
        return Err("Poisson conditional test produced a non-finite log-p".to_string());
    }
    let probability = ln_probability.exp();
    Ok(TailProbability {
        probability,
        ln_probability,
        underflow: probability < f64::MIN_POSITIVE && ln_probability < 0.0,
    })
}

pub fn benjamini_hochberg_from_log(ln_p_values: &[f64]) -> Result<Vec<TailProbability>, String> {
    if ln_p_values.is_empty() {
        return Ok(Vec::new());
    }
    for &ln_p in ln_p_values {
        if !ln_p.is_finite() || ln_p > 0.0 {
            return Err(format!("Invalid BH input log-p: {ln_p}"));
        }
    }
    let m = ln_p_values.len();
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&left, &right| ln_p_values[left].total_cmp(&ln_p_values[right]));
    let mut ln_q_values = vec![0.0; m];
    let mut running_min = 0.0f64;
    for (rank, &index) in order.iter().enumerate().rev() {
        let adjusted = ln_p_values[index] + (m as f64).ln() - ((rank + 1) as f64).ln();
        running_min = running_min.min(adjusted.min(0.0));
        ln_q_values[index] = running_min;
    }
    Ok(ln_q_values
        .into_iter()
        .map(|ln_probability| {
            let probability = ln_probability.exp();
            TailProbability {
                probability,
                ln_probability,
                underflow: probability < f64::MIN_POSITIVE && ln_probability < 0.0,
            }
        })
        .collect())
}

fn monotone_storey_q(p_values: &[f64], pi0: f64) -> Vec<f64> {
    let m = p_values.len();
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| {
        p_values[a]
            .partial_cmp(&p_values[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut q = vec![1.0f64; m];
    let mut min_q = 1.0f64;
    for (rank, &i) in order.iter().enumerate().rev() {
        let val = pi0 * (m as f64) * p_values[i] / (rank as f64 + 1.0);
        min_q = min_q.min(val.min(1.0));
        q[i] = min_q;
    }
    q
}

/// When `survival(k >= 1) = 0`, return `k* = 1` and `lambda = 0`.
fn poisson_support_point_le(r: f64, tau: f64) -> f64 {
    if r.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
        return 0.0;
    }
    let mut hi: u64 = 1;
    while poisson_upper_tail(hi, r) > tau {
        hi *= 2;
    }
    let mut lo: u64 = 0;
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        if poisson_upper_tail(mid, r) > tau {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    poisson_upper_tail(hi, r)
}

/// Estimate pi0 using the CDH 2018 discrete adjustment.
pub fn cdh2018_pi0(p_values: &[f64], rates: &[f64]) -> f64 {
    let m = p_values.len();
    debug_assert_eq!(m, rates.len());
    if m == 0 {
        return 1.0;
    }
    let n = 100usize;
    let tau_1 = 0.25f64;
    let tau_n = 0.5f64;
    let step = (tau_n - tau_1) / (n as f64 - 1.0);
    let mut total = 0.0f64;
    for j in 0..n {
        let tau = tau_1 + step * j as f64;
        let mut est = 1.0 / ((1.0 - tau) * m as f64);
        for t in 0..m {
            let lam = poisson_support_point_le(rates[t], tau);
            if p_values[t] > lam {
                est += 1.0 / ((1.0 - lam) * m as f64);
            }
        }
        total += est.min(1.0);
    }
    (total / n as f64).clamp(0.0, 1.0)
}

pub fn discrete_q_values(p_values: &[f64], rates: &[f64]) -> (f64, Vec<f64>) {
    let pi0 = cdh2018_pi0(p_values, rates);
    let q = monotone_storey_q(p_values, pi0);
    (pi0, q)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

    #[test]
    fn digamma_known_values() {
        assert!((digamma(1.0) + EULER_GAMMA).abs() < 1e-10);
        assert!((digamma(2.0) - (1.0 - EULER_GAMMA)).abs() < 1e-10);
        assert!((digamma(10.0) - 2.251_752_589_066_721).abs() < 1e-9);
    }

    #[test]
    fn trigamma_known_values() {
        let pi2_6 = std::f64::consts::PI * std::f64::consts::PI / 6.0;
        assert!((trigamma(1.0) - pi2_6).abs() < 1e-9);
        assert!((trigamma(10.0) - 0.105_166_335_681_685).abs() < 1e-9);
    }

    #[test]
    fn poisson_upper_tail_matches_partial_sum() {
        // P(X≥5 | λ=2) = 1 − e^{−2}(1+2+2+4/3+2/3)
        let expect = 1.0 - (-2.0f64).exp() * (1.0 + 2.0 + 2.0 + 4.0 / 3.0 + 2.0 / 3.0);
        let got = poisson_upper_tail(5, 2.0);
        assert!((got - expect).abs() < 1e-12, "got {got} expect {expect}");
        assert_eq!(poisson_upper_tail(0, 2.0), 1.0);
        assert_eq!(poisson_upper_tail(3, 0.0), 0.0);
    }

    #[test]
    fn exact_poisson_rate_test_handles_zero_background_without_zero_p() {
        let one = exact_poisson_rate_upper_tail(1, 1.0, 0, 9.0).unwrap();
        assert!((one.probability - 0.1).abs() < 1e-12, "{one:?}");
        assert!((one.ln_probability - 0.1f64.ln()).abs() < 1e-12);
        assert!(!one.underflow);

        let two = exact_poisson_rate_upper_tail(2, 1.0, 0, 9.0).unwrap();
        assert!((two.probability - 0.01).abs() < 1e-12, "{two:?}");

        // n=3, k=2: C(3,2)*0.1^2*0.9 + 0.1^3 = 0.028.
        let mixed = exact_poisson_rate_upper_tail(2, 1.0, 1, 9.0).unwrap();
        assert!((mixed.probability - 0.028).abs() < 1e-12, "{mixed:?}");
    }

    #[test]
    fn exact_poisson_rate_test_preserves_log_tail_after_underflow() {
        let tail = exact_poisson_rate_upper_tail(1_000, 1.0, 0, 999.0).unwrap();
        assert_eq!(tail.probability, 0.0);
        assert!(tail.underflow);
        assert!((tail.ln_probability - 1_000.0 * 0.001f64.ln()).abs() < 1e-8);
    }

    #[test]
    fn exact_poisson_rate_test_recomputes_subnormal_tail_in_log_space() {
        let tail = exact_poisson_rate_upper_tail(1_050, 1.0, 0, 1.0).unwrap();
        assert!(tail.probability > 0.0 && tail.probability < f64::MIN_POSITIVE);
        assert!(tail.underflow);
        assert!((tail.ln_probability - 1_050.0 * 0.5f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn bh_adjustment_uses_log_probabilities() {
        let ln_p = [0.001f64.ln(), 0.01f64.ln(), 0.2f64.ln(), 1.0f64.ln()];
        let q = benjamini_hochberg_from_log(&ln_p).unwrap();
        let expected = [0.004, 0.02, 0.266_666_666_666_666_66, 1.0];
        for (got, want) in q.iter().zip(expected) {
            assert!((got.probability - want).abs() < 1e-12, "{got:?} != {want}");
        }

        let tiny = benjamini_hochberg_from_log(&[-1_000.0, -900.0]).unwrap();
        assert!(tiny[0].underflow && tiny[1].underflow);
        assert!(tiny[0].ln_probability < tiny[1].ln_probability);
    }

    #[test]
    fn nb_upper_tail_approaches_poisson_when_alpha_tiny() {
        let nb = nb_upper_tail(5, 2.0, 1e-12);
        let po = poisson_upper_tail(5, 2.0);
        assert!((nb - po).abs() < 1e-12);
    }

    #[test]
    fn nb_upper_tail_overdispersed_is_larger() {
        let po = poisson_upper_tail(20, 10.0);
        let nb = nb_upper_tail(20, 10.0, 0.2);
        assert!(
            nb > po,
            "overdispersion should increase the upper tail: nb={nb} po={po}"
        );
    }

    #[test]
    fn regularized_beta_symmetry() {
        assert!((regularized_beta(0.5, 2.0, 2.0) - 0.5).abs() < 1e-12);
        assert!((regularized_beta(0.0, 2.0, 3.0)).abs() < 1e-15);
        assert!((regularized_beta(1.0, 2.0, 3.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn trimmed_mean_20_and_eb_shrink_and_expected_hits() {
        let v = [0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 100.0];
        let m = trimmed_mean_20(&v).unwrap();
        assert!(
            (m - (1.0 + 2.0 + 3.0 + 4.0 + 5.0 + 6.0 + 7.0 + 0.0) / 8.0).abs() < 1e-12,
            "trimmed mean {m}"
        );

        let sorted: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let m = trimmed_mean_20(&sorted).unwrap();
        assert!((m - 4.5).abs() < 1e-12, "{m}");

        assert!((trimmed_mean_20(&[1.0, 2.0, 3.0]).unwrap() - 2.0).abs() < 1e-12);
        assert!(trimmed_mean_20(&[]).is_none());

        let v2 = [1.0, 1.0, 1.0, 1.0, 100.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let m2 = trimmed_mean_20(&v2).unwrap();
        assert!((m2 - 1.0).abs() < 1e-12, "outlier was not trimmed: {m2}");

        assert_eq!(eb_shrink_layer(0.5, 0, 1.0, 2.0), 1.0);
        let s = eb_shrink_layer(0.5, 8, 1.0, 2.0);
        assert!((s - (2.0 * 1.0 + 8.0 * 0.5) / 10.0).abs() < 1e-12, "{s}");
        let s = eb_shrink_layer(0.5, 1_000_000, 1.0, 2.0);
        assert!((s - 0.5).abs() < 1e-6, "{s}");

        assert!((expected_hits(300_000, 10_000_000) - 300_000.0).abs() < 1e-9);
        assert!((expected_hits(300_000, 20_000_000) - 600_000.0).abs() < 1e-9);
        assert_eq!(expected_hits(0, 10_000_000), 0.0);
    }

    fn surv_sum(k: u64, r: f64) -> f64 {
        if k == 0 {
            return 1.0;
        }
        let mut s = (-r).exp();
        let mut term = s;
        for i in 1..k {
            term *= r / i as f64;
            s += term;
        }
        1.0 - s
    }

    #[test]
    fn poisson_support_point_le_matches_brute_force() {
        for r in [0.0, 0.5, 2.0, 5.0, 20.0, 300.0] {
            for tau in [0.25, 0.3, 0.3233235838169365, 0.4, 0.5] {
                let got = poisson_support_point_le(r, tau);
                if r <= 0.0 {
                    assert_eq!(got, 0.0);
                    continue;
                }
                let mut k = 1u64;
                loop {
                    let s = poisson_upper_tail(k, r);
                    if s <= tau {
                        assert!(
                            (got - s).abs() < 1e-12,
                            "r={r} τ={tau} support point mismatch: {got} vs {s}"
                        );
                        assert!(
                            poisson_upper_tail(k - 1, r) > tau,
                            "k* is not the first: {k}"
                        );
                        break;
                    }
                    k += 1;
                    assert!(k < 100_000, "k* was not found");
                }
                if r <= 20.0 {
                    let mut k2 = 1u64;
                    while surv_sum(k2, r) > tau {
                        k2 += 1;
                    }
                    let want = surv_sum(k2, r);
                    assert!(
                        (got - want).abs() < 1e-9,
                        "r={r} τ={tau}: gamma path {got} vs summation {want}"
                    );
                }
            }
        }
    }

    #[test]
    fn cdh2018_pi0_matches_reference_values() {
        let r = 2.0;

        let p_a: Vec<f64> = (1..=16)
            .map(|k| poisson_upper_tail(k, r))
            .chain([1e-12; 4])
            .collect();
        let pi0 = cdh2018_pi0(&p_a, &[r; 20]);
        assert!((pi0 - 0.2370636320139451).abs() < 1e-9, "case A: {pi0}");
        assert!((0.0..=1.0).contains(&pi0));

        let pi0 = cdh2018_pi0(&[1e-12; 20], &[r; 20]);
        assert!((pi0 - 0.08111565855756667).abs() < 1e-9, "case B: {pi0}");

        assert_eq!(cdh2018_pi0(&[1.0; 5], &[0.0; 5]), 1.0);

        let p_f: Vec<f64> = (0..80)
            .map(|c| poisson_upper_tail((c % 4) as u64, r))
            .chain([1e-12; 20])
            .collect();
        let pi0 = cdh2018_pi0(&p_f, &vec![r; 100]);
        assert!((pi0 - 0.916910180402218).abs() < 1e-9, "case F: {pi0}");
        assert!(
            pi0 >= 0.8,
            "must be conservative for a plausible null sample: {pi0}"
        );

        assert_eq!(cdh2018_pi0(&[], &[]), 1.0);
    }

    #[test]
    fn cdh2018_pi0_all_null_plausible_is_conservative() {
        let r = 2.0;
        let p: Vec<f64> = (0..30)
            .map(|c| poisson_upper_tail((c % 4) as u64, r))
            .collect();
        let pi0 = cdh2018_pi0(&p, &vec![r; 30]);
        assert!(
            (pi0 - 1.0).abs() < 1e-12,
            "all-null estimate should be approximately 1: {pi0}"
        );
    }

    #[test]
    fn discrete_q_values_monotone_and_bh_formula() {
        let r = 2.0;
        let p: Vec<f64> = (0..80)
            .map(|c| poisson_upper_tail((c % 4) as u64, r))
            .chain([1e-12; 20])
            .collect();
        let (pi0, q) = discrete_q_values(&p, &vec![r; 100]);
        assert!((pi0 - 0.916910180402218).abs() < 1e-9, "π̂0^G: {pi0}");
        assert!(q.iter().all(|&v| (0.0..=1.0).contains(&v)));
        for i in 0..p.len() {
            for j in 0..p.len() {
                if p[i] < p[j] {
                    assert!(q[i] <= q[j], "q is not monotone in p: p[{i}] p[{j}]");
                }
            }
        }
        let min_q = q.iter().copied().fold(f64::INFINITY, f64::min);

        assert!(
            (min_q - pi0 * 100.0 * 1e-12 / 20.0).abs() < 1e-20,
            "rank-1 q formula: {min_q}"
        );
    }

    #[test]
    fn nb_upper_tail_large_c_is_finite_and_monotone() {
        let p1 = nb_upper_tail(200_000, 100_000.0, 0.1);
        let p2 = nb_upper_tail(200_001, 100_000.0, 0.1);
        assert!(p1.is_finite() && p1 > 0.0 && p1 <= 1.0, "p1={p1}");
        assert!(p2 <= p1, "upper tail should decrease with c: {p2} > {p1}");
    }

    #[test]
    fn nb_mle_score_and_hessian_finite_difference() {
        let values = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let r = 3.0;
        let h = 1e-4;
        let (s1_lo, _) = nb_mle_score_and_hessian(&values, mean, r - h);
        let (s1_hi, _) = nb_mle_score_and_hessian(&values, mean, r + h);
        let (_, s2) = nb_mle_score_and_hessian(&values, mean, r);
        let fd = (s1_hi - s1_lo) / (2.0 * h);
        assert!(
            (fd - s2).abs() < 1e-4 + 1e-6 * s2.abs(),
            "finite difference {fd} does not match Hessian {s2}"
        );
    }

    #[test]
    fn nb_alpha_mle_stationary_and_mom_baseline() {
        assert_eq!(nb_alpha_mom(10.0, 10.0), 0.0);
        assert!((nb_alpha_mom(10.0, 30.0) - 0.2).abs() < 1e-12);

        let values: Vec<f64> = (0..12).map(|i| i as f64 % 5.0).collect();
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let var = values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1.0);
        let alpha =
            nb_alpha_mle(&values, nb_alpha_mom(mean, var).max(1e-4)).expect("MLE should converge");
        let r = 1.0 / alpha;
        let (s1, _) = nb_mle_score_and_hessian(&values, mean, r);
        assert!(
            s1.abs() < 1e-4,
            "score at stationary point should be approximately 0: alpha={alpha} s1={s1}"
        );
    }

    #[test]
    fn nb_alpha_mle_degenerate_returns_none() {
        assert!(nb_alpha_mle(&[], 0.2).is_none());
        assert!(nb_alpha_mle(&[0.0, 0.0, 0.0], 0.2).is_none());
    }
}
