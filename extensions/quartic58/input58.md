# Root candidate58: quartic separation of assisted and single-letter losses

2026-10-24. UNREVIEWED corollary candidate of54+56, conditional on their valid fixed-ratio limits. No grade, priority or physical claim for the interpolating functional.

Let a in[0,1], F_a=a lambda H(p)+H(lambda p+qt)-qH(t)-h(lambda). It is concave in p, convex in t. At a0 it is the single-letter reduction; at a1 the assisted reduction. Other a are only a mathematical interpolation, not a claimed capacity formula. Baseline (a+1)lambda g_H(E). Consider continuum mean constraints1 and c=1-epsilon, counting nu_d, g=G_d^(-1) exp(-dx).

Conjectured local smooth saddle, uniformly a in[0,1], no atoms near c1. Let v=d(1-x)g, Q2=d^2 x^2-2d(d+1)x+d(d+1). Orthogonality: int Q2 g dnu=0, int x Q2 g dnu=0, int Q2^2 g dnu=2d(d+1).

Write p_a=g+epsilon^2 u_a+O(epsilon^3), t_a=g+epsilon v+epsilon^2 w_a+O(epsilon^3) in weighted moments. Project second-order KKT log equations onto orthogonal complement of span{1,x} in L2(g dnu):
input: (a+lambda)u_a+q w_a-(q^2/2)Q2 g=0;
reset: w_a-u_a=((1+q)/2)Q2 g.
Therefore u_a=-q Q2 g/[2(a+1)], w_a=[(1+q)/2-q/(2(a+1))]Q2 g.
Check a0: p=(gamma_(lambda+qc)-q gamma_c)/lambda -> u0=-q Q2g/2, w0=Q2g/2 (thermal reset). At a1: u1=-q Q2g/4 and w1=(2+q)Q2g/4, so t_EA-gamma_c=(q/4)epsilon^2 Q2 g+O(epsilon^3). Distinct adversarial shape starts at second order despite identical quadratic loss contact.

At fixed constraints, Gibbs entropy deficit g_entropy-H(p_a)=D(p_a||g)=epsilon^4 q^2 d(d+1)/[4(a+1)^2]+o(epsilon^4), nats. Global saddle envelope gives partial_a D_a=lambda D(p_a||g). Integrate a0to1; int_0^1 (a+1)^(-2) da=1/2. Proposed operational contrast:

D_EA(c)-D_C1(c) ~ [lambda q^2 d(d+1)/(8 ln2)](1-c)^4.

Both curves already have h2(lambda)-lambda q d(1-c)^2/(2ln2)+O((1-c)^3); the new claim says cubic coefficients coincide too and first difference is quartic positive. All limits E->infinity first, then c->1, parameters lambda,d fixed. Regularly varying spectra only if57 closes separately. Do not claim a finite-energy quartic law or uniform joint scaling.

## Analytic obligations / falsifiers

Construct local F_a saddle with fixed mass/means and uniform positive-domain moment IFT for a in[0,1]. Input stationary profile for k=a/(a+1): p proportional exp(-alpha x)(1-z exp(-beta x))^(1/(a+1)), t proportional exp(-(alpha+beta)x)(1-z exp(-beta x))^(-k); normalize with appropriate lambda/q prefactors. At c1, beta0, zq, alpha d. Roots/parameter conventions must be derived carefully rather than reusing54 jacobian with wrong factors. Projection removes only mass and mean components; verify equations/signs and coefficient conventions epsilon^2 versus second derivatives. Establish enough weighted regularity and uniform entropy expansion (unbounded x tails matter). Envelope is for saddles, not illicit inf/sup exchange; can instead use global tangent inequalities to bound a increments. Failure of local uniform control invalidates the equivalent even if formal Taylor algebra is right. Numerical fit alone is not proof.
