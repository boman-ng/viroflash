//! 统计计算：两 Poisson 率精确条件检验与 log-space BH（主路径），以及
//! Poisson/NB/CDH 离散方法（诊断和已覆盖的数学原语）。
//!
//! - 正则化 gamma（级数 + Lentz 连分数）：Numerical Recipes 标准实现；
//! - NB 上尾 = I_p(c, r)（r=1/α），α 估计 MoM 优先（Robinson & Smyth 2008, Biostatistics）；
//! - 经验 decoy p 的 +1 修正：North/Curtis/Sham 2002, AJHG；
//! - BH running-min：Benjamini & Hochberg 1995, JRSS-B。

const FPMIN: f64 = 1e-300;
const ITMAX: usize = 200;
const EPS: f64 = 3e-14;

/// log Γ(x)（Lanczos, Numerical Recipes gammaln）。
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

/// 下正则化 gamma P(a,x)（级数展开，x < a+1 时收敛快）。
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

/// 上正则化 gamma Q(a,x) = Γ(a,x)/Γ(a)（Lentz 连分数）。
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

/// 下正则化 gamma P(a,x) = γ(a,x)/Γ(a)；Poisson 上尾 P(X ≥ c | λ) = P(c, λ)。
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

/// Poisson 上尾：P(X ≥ c | λ)。
pub fn poisson_upper_tail(c: u64, lambda: f64) -> f64 {
    if lambda <= 0.0 {
        return if c > 0 { 0.0 } else { 1.0 };
    }
    if c == 0 {
        return 1.0;
    }
    regularized_gamma_p(c as f64, lambda)
}

/// digamma（x<6 递归上移 + 渐近展开）。
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

/// trigamma（x<6 递归上移 + 渐近展开）。
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

/// NB overdispersion 矩估计：α = max((s²−μ)/μ², 0)。
pub fn nb_alpha_mom(mean: f64, var: f64) -> f64 {
    ((var - mean) / (mean * mean)).max(0.0)
}

/// NB α 的 Newton-Raphson MLE 得分与信息（均值固定，r = 1/α）。
/// 模型（Robinson & Smyth 2008 固定均值 profile 似然）：
/// ℓ ∝ Σ logΓ(v+r) − n·logΓ(r) + n·r·log(r/(r+μ)) − Σv·log(r+μ)。
/// d/dr 中 n·μ/(r+μ)（r·log(r/(r+μ)) 项）与 −Σv/(r+μ) = −n·μ/(r+μ) 相消，
/// 故 s1 = Σ digamma(v+r) + n·log(r/(r+μ)) − n·digamma(r)；
/// 再对 r 求导（末两项 n/(r+μ)² 对消）：
/// s2 = Σ trigamma(v+r) + n·μ/(r(r+μ)) − n·trigamma(r)。
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

/// NB α 的 Newton-Raphson MLE（均值固定，迭代 r=1/α，≤10 次，初值来自 init_alpha）。
/// 收敛判据 |step| < 1e-6·r；不收敛或退化返回 None（调用方回退 MoM）。
/// 依据 Robinson & Smyth 2008；小样本 MLE 不稳，调用方须优先使用 MoM。
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
            return None; // 越界：不收敛
        }
    }
    None // 迭代耗尽未收敛
}

/// 正则化 beta 的连分数（Lentz）。
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

/// 正则化 beta I_x(a,b)。
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

/// NB 上尾 P(X≥c | mean μ, overdispersion α)：I_p(c, r)，r=1/α，p=μα/(1+μα)。
/// 恒用正则化 beta，无近似分支。
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

/// 20% 截尾均值：两端各裁 ⌊0.1·n⌋ 个，n < 10 时退化普通均值；
/// 小层由 EB 收缩兜底。空输入返回 None。
pub fn trimmed_mean_20(values: &[f64]) -> Option<f64> {
    let n = values.len();
    if n == 0 {
        return None;
    }
    let mut sorted: Vec<f64> = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let trim = (n as f64 * 0.1).floor() as usize; // 每端
    let kept = &sorted[trim..n - trim];
    Some(kept.iter().sum::<f64>() / kept.len() as f64)
}

/// 经验贝叶斯收缩：λ̂_ℓ* = (κ·λ̂_glob + n_ℓ·λ̂_ℓ) / (κ + n_ℓ)。
/// n_ℓ=0 → λ̂_glob（无诱饵层退化）；n_ℓ→∞ → λ̂_ℓ。
/// κ 为伪计数权重（假想全局层贡献 κ 个"伪诱饵"）。
pub fn eb_shrink_layer(layer_mean: f64, layer_n: usize, global_mean: f64, kappa: f64) -> f64 {
    (kappa * global_mean + layer_n as f64 * layer_mean) / (kappa + layer_n as f64)
}

/// 暴露量：E_t = L_t · (N_total / 10^7)。
/// λ̂（每碱基期望 reads）与 E_t（碱基×深度系数）相乘 = 目标 t 的 Poisson 期望，
/// 使 λ̂ 不随样本 reads 总量漂移、跨样本可比。
pub fn expected_hits(contig_len: u64, n_total_reads: u64) -> f64 {
    contig_len as f64 * (n_total_reads as f64 / 1e7)
}

/// 可同时保留普通概率和对数概率的尾概率。
///
/// `underflow` 表示普通概率进入 f64 subnormal/zero 区间；此时
/// `ln_probability` 由 log-space 递推保留主要数值信息，不表示数学概率为 0。
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

/// 两个独立 Poisson 率的单侧精确条件检验。
///
/// H0：`target_count / target_exposure == background_count / background_exposure`；
/// H1：target 率更高。对总事件数条件化后，target 计数服从
/// `Binomial(total, target_exposure / total_exposure)`，与 R `poisson.test` 的两样本
/// 精确检验同构。该检验能在 background_count=0 时保留有限样本不确定性，不把
/// 观测零背景率当作已知的精确零。
pub fn exact_poisson_rate_upper_tail(
    target_count: u64,
    target_exposure: f64,
    background_count: u64,
    background_exposure: f64,
) -> Result<TailProbability, String> {
    if !target_exposure.is_finite() || target_exposure <= 0.0 {
        return Err("target exposure 必须为有限正数".to_string());
    }
    if !background_exposure.is_finite() || background_exposure <= 0.0 {
        return Err("background exposure 必须为有限正数".to_string());
    }
    let total_count = target_count
        .checked_add(background_count)
        .ok_or_else(|| "Poisson 条件检验总计数溢出".to_string())?;
    if target_count == 0 || total_count == 0 {
        return Ok(TailProbability {
            probability: 1.0,
            ln_probability: 0.0,
            underflow: false,
        });
    }

    let target_probability = target_exposure / (target_exposure + background_exposure);
    if !(0.0..1.0).contains(&target_probability) {
        return Err("Poisson 条件检验 exposure 比例退化".to_string());
    }

    // P[X >= target_count] = I_p(target_count, background_count + 1)。通常直接走
    // regularized beta；极小尾部下溢时再以 log-PMF 递推求 log-sum-exp。
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
        return Err("Poisson 条件检验产生非有限 log-p".to_string());
    }
    let probability = ln_probability.exp();
    Ok(TailProbability {
        probability,
        ln_probability,
        underflow: probability < f64::MIN_POSITIVE && ln_probability < 0.0,
    })
}

/// 固定检验族的 Benjamini-Hochberg 调整；全程以 ln(p) 排序和单调化，避免极小
/// p 在进入多重校正前先下溢为 0。
pub fn benjamini_hochberg_from_log(ln_p_values: &[f64]) -> Result<Vec<TailProbability>, String> {
    if ln_p_values.is_empty() {
        return Ok(Vec::new());
    }
    for &ln_p in ln_p_values {
        if !ln_p.is_finite() || ln_p > 0.0 {
            return Err(format!("BH 输入 log-p 非法: {ln_p}"));
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

/// BH 型单调 q：排序后 q_(i) = min(π̂0·m·p_(i)/i, q_(i+1))（Storey 2002 §3）。
/// π̂0=1 时退化为 BH q。离散 p 主路径由 `discrete_q_values` 提供 π̂0（CDH2018）。
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
        // rank+1 == m 时 m/(rank+1)=1 → q_(m)=π̂0·p_(m)
        let val = pi0 * (m as f64) * p_values[i] / (rank as f64 + 1.0);
        min_q = min_q.min(val.min(1.0));
        q[i] = min_q;
    }
    q
}

/// max{λ ∈ S_t : λ ≤ τ}：S_t = {P(Pois(r) ≥ k) : k ≥ 0}，survival 随 k 严格递减
/// （r>0），故 λ = survival(k*)，k* = min{k≥1 : survival(k) ≤ τ}
/// （survival(0)=1 > τ 恒成立，τ∈(0,1)）。倍增加二分定位 k*，
/// 每 (t,τ) 约 2·log₂(k*) 次尾部求值。r ≤ 0 或非有限视作 rate 0：
/// survival(k≥1)=0 → k*=1、λ=0。
fn poisson_support_point_le(r: f64, tau: f64) -> f64 {
    if r.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
        return 0.0;
    }
    let mut hi: u64 = 1;
    while poisson_upper_tail(hi, r) > tau {
        hi *= 2; // survival(k)→0 单调递减，必然终止
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

/// CDH2018 离散感知 π0 估计（Chen–Doerge–Heyse 2018, DOI 10.1002/bimj.201700157，
/// Algorithm 1 + §3.2 引导值；与作者 R 包 fdrDiscreteNull `GenEstProp` 逐项核对）。
/// p_t = Poisson 上尾 P(Pois(r_t) ≥ n_t)（r_t = λ̂_ℓ*·E_t），支持 S_t 如上。
/// - q_t = inf S_t = 0（有限速率下 survival(k) > 0 恒正、下确界 0）→ γ = 0，
///   Dirac 集 C = {t : q_t = 1} 恒空，|C|/m 项省略；
/// - 引导值：τ_1 = γ + 0.5·(0.5−γ) = 0.25，τ_n = 0.5，n = 100 等距
///   （步长 (τ_n−τ_1)/99，`seq(0.25, 0.5, length.out=100)` 语义）；
/// - 试估计 β(τ_j) = 1/((1−τ_j)m) + (1/m)·Σ_t 1{p_t > λ_tj}/(1−λ_tj)，
///   逐项截断 ≤ 1（原文 "Truncate β(τ_j) at 1"）；π̂0^G = mean_j β(τ_j)。
///
/// Theorem 1：E[β(τ_j)] = π0 + 非负偏（保守）；Theorem 2：aBH = BH@α/π̂0^G 保守。
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

/// 离散感知 Storey q：π̂0 = `cdh2018_pi0`（离散 p 支持集感知），
/// q 单调化用 `monotone_storey_q`（Storey 2002 §3；aBH 语义：BH@α/π̂0^G，Theorem 2 保守）。
/// `rates` = 各测试 Poisson 期望（λ̂_ℓ*·E_t），长度须与 p_values 一致。
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
        // 条件于总计数：exposure 1:9 时 X~Binomial(n, 0.1)。
        let one = exact_poisson_rate_upper_tail(1, 1.0, 0, 9.0).unwrap();
        assert!((one.probability - 0.1).abs() < 1e-12, "{one:?}");
        assert!((one.ln_probability - 0.1f64.ln()).abs() < 1e-12);
        assert!(!one.underflow);

        let two = exact_poisson_rate_upper_tail(2, 1.0, 0, 9.0).unwrap();
        assert!((two.probability - 0.01).abs() < 1e-12, "{two:?}");

        // n=3、k=2：C(3,2)·0.1²·0.9 + 0.1³ = 0.028。
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
        // 0.5^1050 是正的 subnormal f64；普通值精度受限，但 log-p 应保持解析值。
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
        assert!(nb > po, "overdispersion 应抬升上尾: nb={nb} po={po}");
    }

    #[test]
    fn regularized_beta_symmetry() {
        assert!((regularized_beta(0.5, 2.0, 2.0) - 0.5).abs() < 1e-12);
        assert!((regularized_beta(0.0, 2.0, 3.0)).abs() < 1e-15);
        assert!((regularized_beta(1.0, 2.0, 3.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn trimmed_mean_20_and_eb_shrink_and_expected_hits() {
        // 20% 截尾：n=10 两端各裁 1 → 中 8 项均值
        let v = [0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 100.0];
        let m = trimmed_mean_20(&v).unwrap();
        assert!(
            (m - (1.0 + 2.0 + 3.0 + 4.0 + 5.0 + 6.0 + 7.0 + 0.0) / 8.0).abs() < 1e-12,
            "截尾均值 {m}"
        );
        // 已排序输入同样正确（排序无破坏）
        let sorted: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let m = trimmed_mean_20(&sorted).unwrap();
        assert!((m - 4.5).abs() < 1e-12, "{m}");
        // n<10 退化普通均值；空 → None
        assert!((trimmed_mean_20(&[1.0, 2.0, 3.0]).unwrap() - 2.0).abs() < 1e-12);
        assert!(trimmed_mean_20(&[]).is_none());
        // 抗异常值：截尾后不受 100 影响（20% 截尾的关键性质）
        let v2 = [1.0, 1.0, 1.0, 1.0, 100.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let m2 = trimmed_mean_20(&v2).unwrap();
        assert!((m2 - 1.0).abs() < 1e-12, "异常值未被裁掉: {m2}");

        // EB 收缩：κ=2，n=0 → 全局；n=8 → 加权；n→∞ 逼近层均值
        assert_eq!(eb_shrink_layer(0.5, 0, 1.0, 2.0), 1.0);
        let s = eb_shrink_layer(0.5, 8, 1.0, 2.0);
        assert!((s - (2.0 * 1.0 + 8.0 * 0.5) / 10.0).abs() < 1e-12, "{s}");
        let s = eb_shrink_layer(0.5, 1_000_000, 1.0, 2.0);
        assert!((s - 0.5).abs() < 1e-6, "{s}");

        // 暴露量：L=3e5、N=1e7 → E=3e5；N 翻倍 E 翻倍
        assert!((expected_hits(300_000, 10_000_000) - 300_000.0).abs() < 1e-9);
        assert!((expected_hits(300_000, 20_000_000) - 600_000.0).abs() < 1e-9);
        assert_eq!(expected_hits(0, 10_000_000), 0.0);
    }

    /// 测试参考实现：Poisson 上尾的逐项求和（与 gamma 路径独立；k 小）。
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
        // λ 必须是支持点：λ = survival(k*)，k* = 首个 survival ≤ τ 的计数。
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
                            "r={r} τ={tau} 支持点不符: {got} vs {s}"
                        );
                        assert!(poisson_upper_tail(k - 1, r) > tau, "k* 非首个: {k}");
                        break;
                    }
                    k += 1;
                    assert!(k < 100_000, "k* 未找到");
                }
                if r <= 20.0 {
                    // 交叉验证：独立求和路径（r 大时求和会溢出，仅小 r）
                    let mut k2 = 1u64;
                    while surv_sum(k2, r) > tau {
                        k2 += 1;
                    }
                    let want = surv_sum(k2, r);
                    assert!(
                        (got - want).abs() < 1e-9,
                        "r={r} τ={tau}: gamma 路径 {got} vs 求和 {want}"
                    );
                }
            }
        }
    }

    #[test]
    fn cdh2018_pi0_matches_reference_values() {
        // 参考值由独立 Python 实现逐项计算（算法镜像，栅格 seq(0.25,0.5,len=100)）。
        // 注意：p 值必须与 λ 同走 gamma 路径（poisson_upper_tail）——生产路径同源，
        // n_t == k* 时 p == λ 精确相等（严格 > 判 0）；混用求和路径会因 ulp 差翻转边界。
        let r = 2.0;
        // Case A: m=20，16 个 Pois(2) 上尾 p（计数 1..16）+ 4 备择 p=1e-12。
        let p_a: Vec<f64> = (1..=16)
            .map(|k| poisson_upper_tail(k, r))
            .chain([1e-12; 4])
            .collect();
        let pi0 = cdh2018_pi0(&p_a, &[r; 20]);
        assert!((pi0 - 0.2370636320139451).abs() < 1e-9, "case A: {pi0}");
        assert!((0.0..=1.0).contains(&pi0));
        // Case B: 全备择 → 技术项主导的保守小值。
        let pi0 = cdh2018_pi0(&[1e-12; 20], &[r; 20]);
        assert!((pi0 - 0.08111565855756667).abs() < 1e-9, "case B: {pi0}");
        // Case D: rate 0 退化（p 恒 1）→ 全 null → 1.0。
        assert_eq!(cdh2018_pi0(&[1.0; 5], &[0.0; 5]), 1.0);
        // Case F: m=100，80 null（计数循环 0..3，null 可信）+ 20 备择。
        // π0_true = 0.8；样本值须保守（≥ π0_true）。
        let p_f: Vec<f64> = (0..80)
            .map(|c| poisson_upper_tail((c % 4) as u64, r))
            .chain([1e-12; 20])
            .collect();
        let pi0 = cdh2018_pi0(&p_f, &vec![r; 100]);
        assert!((pi0 - 0.916910180402218).abs() < 1e-9, "case F: {pi0}");
        assert!(pi0 >= 0.8, "对可信 null 样本须保守: {pi0}");
        // 空输入
        assert_eq!(cdh2018_pi0(&[], &[]), 1.0);
    }

    #[test]
    fn cdh2018_pi0_all_null_plausible_is_conservative() {
        // 全 null、计数 null 可信（Pois(2) 支撑 0..3）：逐 β 截断 → π̂0^G = 1.0。
        let r = 2.0;
        let p: Vec<f64> = (0..30)
            .map(|c| poisson_upper_tail((c % 4) as u64, r))
            .collect();
        let pi0 = cdh2018_pi0(&p, &vec![r; 30]);
        assert!((pi0 - 1.0).abs() < 1e-12, "全 null 应≈1: {pi0}");
    }

    #[test]
    fn discrete_q_values_monotone_and_bh_formula() {
        // Case F 数据：q 单调、界内；最小 p 的 q = π̂0·m·p_min（rank-1 BH 公式）。
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
                    assert!(q[i] <= q[j], "q 不随 p 单调: p[{i}] p[{j}]");
                }
            }
        }
        let min_q = q.iter().copied().fold(f64::INFINITY, f64::min);
        // 20 个 p=1e-12 并列：BH 公式按并列组最大秩 20 计算。
        assert!(
            (min_q - pi0 * 100.0 * 1e-12 / 20.0).abs() < 1e-20,
            "rank-1 q 公式: {min_q}"
        );
    }

    #[test]
    fn nb_upper_tail_large_c_is_finite_and_monotone() {
        // 回归：大 c 走正则化 beta 路径（无近似分支），须有限且随 c 递减
        let p1 = nb_upper_tail(200_000, 100_000.0, 0.1);
        let p2 = nb_upper_tail(200_001, 100_000.0, 0.1);
        assert!(p1.is_finite() && p1 > 0.0 && p1 <= 1.0, "p1={p1}");
        assert!(p2 <= p1, "上尾应随 c 递减: {p2} > {p1}");
    }

    #[test]
    fn nb_mle_score_and_hessian_finite_difference() {
        // Hessian 项 n·μ/(r(r+μ)) 的回归防线：有限差分 (s1(r+h)−s1(r−h))/2h ≈ s2(r)
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
            "有限差分 {fd} 与 Hessian {s2} 不符"
        );
    }

    #[test]
    fn nb_alpha_mle_stationary_and_mom_baseline() {
        assert_eq!(nb_alpha_mom(10.0, 10.0), 0.0);
        assert!((nb_alpha_mom(10.0, 30.0) - 0.2).abs() < 1e-12);
        // 中等过离散：MLE 收敛到驻点（s1(1/α) ≈ 0）
        let values: Vec<f64> = (0..12).map(|i| i as f64 % 5.0).collect();
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let var = values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1.0);
        let alpha = nb_alpha_mle(&values, nb_alpha_mom(mean, var).max(1e-4)).expect("MLE 应收敛");
        let r = 1.0 / alpha;
        let (s1, _) = nb_mle_score_and_hessian(&values, mean, r);
        assert!(s1.abs() < 1e-4, "驻点处得分应≈0: alpha={alpha} s1={s1}");
    }

    #[test]
    fn nb_alpha_mle_degenerate_returns_none() {
        assert!(nb_alpha_mle(&[], 0.2).is_none());
        assert!(nb_alpha_mle(&[0.0, 0.0, 0.0], 0.2).is_none());
    }
}
