# Root analytic closure proposed for58, separate immutable addendum

UNREVIEWED, supplements but does not overwrite58. Put k=a/(a+1), e=1/(a+1), z=exp(-delta). Profile kernels relative nu_d:
u=exp(-alpha x)(1-z exp(-beta x))^e,
v=exp(-alpha x)z exp(-beta x)(1-z exp(-beta x))^(-k),
Z=int u, p=u/Z, t=lambda v/(qZ), m=lambda exp(-alpha x)(1-z exp(-beta x))^(-k)/Z.
Then p^a m=lambda exp(-(a+1)alpha x)/Z^(a+1), t/m=exp(-delta-beta x)/q. These identities supply global concave/convex saddle inequalities at given moments; no inf/sup swap. For finite positive energies all log kernels grow at most linearly. Continuum entropy integrals finite.

Moment map T=(ln(lambda/q)+ln(V/U),U1/U,V1/V), target(0,1,1-epsilon). At alpha=d,beta0,delta=-lnq its Jacobian rows are
[0,-1/lambda,-1/lambda],
[-1/d,e q/(lambda d),0],
[-1/d,-(lambda+k q)/(lambda d),0].
Determinant -1/(lambda^2 d^2), independent a. Thus alpha'=e qd,beta'=lambda d,delta'=-lambda d; p'=0,t'=d(1-x)g. The scalar a is an auxiliary smooth parameter on a neighborhood of[0,1]. We may use its neighborhood only locally; global concavity and capacity interpretation needed only a in[0,1].

For alpha>=d/2,delta>=(-lnq)/2,beta>=0, all fixed-order derivatives in alpha,beta,delta,a of u,v and entropy integrands are dominated by a polynomial in x times x^(d-1) exp(-d x/2), uniformly a in[0,1]. Powers of log(1-z exp(-beta x)) are uniformly bounded since z bounded below1. Logs p/g only add bounded constants and linear x; normalizers bounded away0 on compact parameter sets. A one-sided C^r moment map can be extended to beta<0 by its boundary Taylor polynomial as in54; no physical density at beta<0. The uniform nonzero determinant, compactness in a and beta'>0 give a uniform local physical saddle for positive epsilon. Need articulate uniqueness of overlapping IFT solutions, which agree in common neighborhoods by the IFT; only a finite subcover is required.

Mass/mean constraints imply int u_a= int x u_a=0 and similarly for w_a, using the coefficient convention p=g+epsilon^2 u_a+O(epsilon^3), t=g+epsilon v+epsilon^2 w_a+O(epsilon^3). Project log KKT equations as58. The polynomial v^2/g, after subtracting its span{1,x} projection, is Q2g. All projected second-order components are fixed by the invertible two-by-two system with determinant a+1. Hence u_a=-q Q2g/[2(a+1)], w_a=u_a+(1+q)Q2g/2.

The relative entropy D(p_a||g) is C^5 in epsilon uniformly compact a by these dominated derivatives (one-sided physical domain). Its zeroth through third Taylor coefficients vanish because p'=0 and mass constraints. Fourth coefficient is one-half int u_a^2/g =q^2 d(d+1)/[4(a+1)^2]. This avoids assuming a pointwise uniform expansion of p_a/g on unbounded x. Can use C4 with uniform Taylor remainder instead of C5.

At the stationary constrained saddle, derivatives of p,t with respect to a have zero mass and means, so their gradient terms vanish, giving partial_a V_a=lambda H(p_a). Consequently partial_a D_a=lambda[g_entropy-H(p_a)]=lambda D(p_a||g). A derivative proof can be replaced by two saddle inequalities for nearby a; differentiability holds for this local IFT branch. Integrating the uniform entropy expansion gives the58 quartic coefficient. Exact continuum value, useful as a numerical check:
V_a=(a+1)lambda lnZ+(a+1)lambda alpha-q beta c-q delta;
D_a=(a+1)lambda(d+lnG_d)-V_a (nats).

For E->infinity correspondence for a interior, if needed for the mathematical envelope, the same discrete normalization and moment Riemann sums as54 apply locally beta>=0,delta bounded away0; no singular kernel or atom. Gamma baseline and exact entropy identity are the displayed formulas with Z_L. The slack squeeze retains both constraints. The interpolation at a interior remains solely a mathematical game.
