# Root-authored candidate 56: product-state classical versus assisted reset games

2026-10-24. UNREVIEWED proof candidate, not a paper or priority claim. Root derivation; no grade. Preserve manuscript53 unchanged. Possible companion to spectral candidate54, preferably combined if validated.

## Scope and falsifiers

H_j=j^(1/d), d>0; lambda in (0,1), q=1-lambda. Reset R_tau(rho)=lambda rho+q tau. Input energy E, reset budget M. Study single-letter Holevo capacity C1 with product-state classical encoding, NOT regularized classical capacity. Infinite ensembles/finite approximations, uniform reduction, saddle signs, boundary positivity and joint parameter limits require audit. No multiplicativity-p2 to Shannon shortcut. Spectral exponent d is not physical dimension.

## Uniform reduction candidate

For a pure ensemble with average rho, output entropy of each letter equals h(lambda)+q S(tau)-S(Z|B). Binary discrimination/Fano gives 0<=S(Z|B)_psi<=h(q <psi|tau|psi>). Jensen bounds its ensemble average by h(q Tr rho tau). Let J(rho,tau)=S(lambda rho+q tau)-q S(tau)-h(lambda). Thus chi=J+average S(Z|B), and sup J<=C1. Flag upper bound chi<=lambda S(rho). A Gibbs trial gives sup J>=lambda g_d(E)-h(lambda), so near-optimal inputs have uniformly bounded entropy deficit, hence ||rho||infinity=O(1/log E). Jensen correction then o(1), uniformly over each finite reset budget (including M depending on E). Need justify pure finite ensembles approximate sup J under energy constraint. J is convex in tau via weighted Jensen-Shannon divergence. Enumeration phase twirls U_N|j>=exp(2pi i j/N)|j> commute with H; dephasing and finite-branch continuity should give, without swapping inf/sup,

W(E,M)=inf_t sup_p [H(lambda p+q t)-q H(t)-h(lambda)],

and inf_tau C1=W+o(1). Entropies in nats below; divide by ln2 for bits.

## Thermal branch

Counting measure nu_d(dx)=d x^(d-1)dx; mean input normalized to1, reset mean c. mu=lambda+q c. m=Gibbs(mean mu), t=Gibbs(mean c), p=(m-q t)/lambda. Positivity iff (c/mu)^d>=q, equivalently c>=cstar=lambda q^(1/d)/(1-q^(1+1/d)). For cstar<c<1 the finite-E Gibbs construction is eventually positive, since t_j/m_j decreases with H_j and minimum p positivity occurs at j0. Input derivative -lambda(log m+1) affine in energy; reset derivative q log(t/m) affine with negative slope. Entropy/KL tangent inequalities should certify global saddle. At cstar use monotone/continuity squeeze, not unproved finite positivity. For c>=1 plateau h(lambda).

Deficit relative lambda g_d(E): D_C(c)=[h(lambda)-d ln(lambda+q c)+q d ln c]/ln2 on [cstar,1].

## Atom branch, all d>0

For 0<c<cstar let v=(q c/(lambda+q c))^(1/(d+1)), r=(1-v^d)/(d(1-v^(d+1))), s=v r, alpha=1/r, beta=1/s-1/r, A=lambda/[Gamma(d+1)(r^d-s^d)]. Relative to nu_d,

m_cont=A exp(-alpha x),
t_cont=(A/q) exp(-(alpha+beta)x),
f=(A/lambda)[exp(-alpha x)-exp(-(alpha+beta)x)].

Reset atom w=(q-v^d)/[q(1-v^d)]>0. Input mass and mean1; reset continuous mass lambda v^d/[q(1-v^d)], mean c=lambda v^(d+1)/[q(1-v^(d+1))]. Reset derivative q log(t/m)=-q ln q-q beta x matches atom intercept -q ln q. Input gradient affine on x>0.

Discrete certificate j>=1, x_j=H_j/L:
u_j=exp(-alpha x_j)-exp(-(alpha+beta)x_j), Z=sum u_j;
p_j=u_j/Z, p0=0;
t_j=(lambda/q) exp(-(alpha+beta)x_j)/Z, t0=1-sum_{j>=1} t_j;
m_j=lambda exp(-alpha x_j)/Z, m0=q t0.
For fixed c<cstar, t0->w>0 and Z~L^d integral u. Input ground gradient inequality requires q t0>=lambda/Z, true eventually. Reset gradient intercept agrees including ground. Ordinary Riemann sums, no singular t continuum at zero. L=(1-eta)E then eta->0 and tangent upper/lower squeezes should give exact limit.

Continuum value lambda ln(lambda/A)+lambda alpha-q beta c, with lambda alpha-q beta c=lambda d. Deficit:
D_C(c)=lambda[d ln(1-v^(d+1))-(d+1)ln(1-v^d)]/ln2.

## Consequences to verify

Small c: D_C(c)~(d+1)lambda^(1/(d+1))q^(d/(d+1))c^(d/(d+1))/ln2.
Near1: D_C(c)=h2(lambda)-lambda q d(1-c)^2/(2ln2)+O((1-c)^3).
At d1: cstar=q/(1+q), r+s=1, D_C=-lambda ln(1-2s)/ln2; s=[-q c+sqrt(q^2 c^2+lambda q c)]/lambda.
EA candidate54 has nonthermal optimum for every fixed0<c<1; this C1 candidate has thermal saddle above positive cstar. C1 atomic low-c branch exists for every d>0, whereas candidate54's constructed EA atom exists only for d>1/2. This is a contrast between constructions; uniqueness/all weak limits not claimed.
EA/C1 small-c coefficient ratio [2 C_d/Gamma(d)]^(1/(d+1)), C_d=integral y^(d-1)(1-sqrt(1-exp(-y)))dy, lies between1 and2^(1/(d+1)).

## Bibliographic uncertainty

GLR2005 quant-ph0408103v3 SectionIV.B Eq34 already defines the channel. Baitenov2024 DOI10.1134/S0081543824010048 studies one-shot shifted-depolarizing capacity; only abstract obtained so far, overlap unresolved. Search robust contamination/Jensen-Shannon games and quantum one-shot capacities before claiming priority. General unassisted regularized capacity is not addressed; entangled codewords invalidate naive tensorization of pure-letter Fano bound.
