# Exact enclosure for one finite entropy game

```sh
python certificates/finite_entropy/verify.py --enable
python -m unittest discover -s certificates/finite_entropy -p 'test_*.py' -v
```

Python 3.12 standard library only. The verifier prints its result to stdout;
it does not run an optimizer or rewrite scientific files. The API also requires
`certify(case, enabled=True)`. Validation remains active under `python -O`.

For **lambda = 3/5, E = 8, c = 1/20, levels i = 0,...,320**, the exact finite
game value is enclosed in **[4.9836437909, 4.9836437910] bits**. The rational
parameters specify the mathematical target; this is not a claim about exact
IEEE interpretations of the historical program's `.6` and `.05` literals.

## Problem and certificate

Let x_i = i/8, q = 1-lambda and m = lambda p + q t. Maximize in probability
vectors p with sum x_i p_i <= 1, minimize in probability vectors t with
sum x_i t_i <= c, and use

    Phi(p,t) = lambda H(p) + H(m) - q H(t) - h2(lambda).

Entropy is in bits, with 0 log 0 = 0. The certificate contains strictly positive
rational p,t with exact normalization and exact means 1,c. Their denominator
is 10^130. The objective is concave in p and convex in t: the diagonal second
derivatives, multiplied by ln(2), are respectively

    -lambda/p_i - lambda^2/m_i,
    q/t_i - q^2/m_i = lambda*q*p_i/(t_i*m_i).

Continuity extends the tangent inequalities to the boundary. At the certified
pair define

    g_i = -lambda (ln p_i + ln m_i + 2)/ln(2),
    h_i = q (ln t_i - ln m_i)/ln(2).

For rational slopes B >= 0 and b <= 0, choose A and a so that all coordinates
satisfy g_i <= A+B*x_i and h_i >= a+b*x_i. Then the max-min value V obeys

    Phi + a + b*c - sum h_i*t_i <= V
        <= Phi + A + B - sum g_i*p_i.

Indeed, the left side bounds min_t Phi(p_cert,t), and the right side bounds
max_p Phi(p,t_cert). Weak minimax inequalities suffice; no minimax equality
or optimal multiplier assumption is needed. The verifier obtains A from the
maximum of upper gradient bounds minus B*x_i, and a from the minimum of lower
gradient bounds minus b*x_i. It uses conservative endpoints for Phi and both
weighted gradient sums. These are 642 coordinate inequalities, not a floating
LP optimality assertion.

## Rigorous logarithms

For rational z>0, exact power-of-two scaling gives z=2^k*r with 1<=r<2.
With u=(r-1)/(r+1), use

    ln r = 2 sum_{j=0}^{79} u^(2j+1)/(2j+1) + R,
    0 <= R <= 9 / (4*161*3^161).

The same bound permits r=2 to enclose ln(2). Integer fixed-point arithmetic
at scale 10^60 rounds every lower operation down and every upper operation
up. The geometric remainder is included in the upper endpoint. Negative k
reverses the endpoints of k*ln(2). Remaining operations use exact Fraction
arithmetic, including division by the strictly positive interval for ln(2).
The displayed decimal endpoints are rounded outwards to 10^-10 bits.

No platform logarithm, Decimal approximation, floating-point solver or stored
decimal result is trusted by this verifier. The tests include a second direct
Fraction series with its own remainder, parameter/mean/mass/sign mutations,
and rejection under optimized Python. These software controls are not Lean
formalization or independent mathematical review.

## Provenance and limits

The source is the first archived row in `finite-game-results50.json` at
commit `a916576da3c770601732bfe5252353989bcb844f` (E=8,c=.05,dimension321).
Its alpha,beta,delta parameters were used to propose probabilities using
160-digit Decimal arithmetic; no optimization was repeated. For i>=2 the
normalized proposed probabilities were rounded down to multiples of 10^-130.
Coordinates 1 and 0 were then solved exactly from mean and normalization:

    p_1 = 8 - sum_{i>=2} i*p_i; p_0 = 1-p_1-sum_{i>=2}p_i,
    t_1 = 2/5 - sum_{i>=2} i*t_i; t_0 = 1-t_1-sum_{i>=2}t_i.

The verifier checks positivity and both equalities afresh. Thus this is a
**new auxiliary feasible pair**, not an assertion that historical floating
arrays were exactly feasible. The complete certificate has SHA256
`fde590462bae8b16d085ebe966156acc223b6fc40c8a0a0eaed45e9c28c9bff5`.
Its validity depends on the stored rational vectors and the checks, not on
reproducing the Decimal proposal or trusting its optimizer metadata.

This is a standard tangent-duality certificate for one finite game. It does
not establish an infinite-cutoff limit, quantum capacity, tensorization,
uniqueness, methodological novelty, a manuscript assessment, or a new paper.
Historical algorithms, outputs and scientific manifests remain unchanged.

## CI scope

The exact certificate has its own dependency-free Windows/Linux workflow.
The three historical diagnostic workflows retain their commands and manual
dispatch. Their automatic triggers now select their scientific inputs,
runner scripts and dependency files, so certificate-only and documentation
changes do not rerun the 28 optimizations or the Holevo/tensorization campaigns.
Changes to a diagnostic workflow itself can be checked through manual dispatch;
they do not automatically rerun those scientific campaigns.
