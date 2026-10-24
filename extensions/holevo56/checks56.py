"""Independent finite diagnostics for candidate56. These are not a proof or solver.
All logarithms are natural. Deterministic seed; no modification of inputs.
"""
from pathlib import Path
import hashlib, json, math, platform
import numpy as np
import scipy
from scipy.special import xlogy, gammaln, gamma
from scipy.integrate import quad
from scipy.optimize import brentq
import sympy as sp

ROOT = Path(__file__).resolve().parent
RNG = np.random.default_rng(2026102456)

def entropy(p):
    p = np.asarray(p)
    return float(-np.sum(xlogy(p, p)))

def h(lam):
    return entropy([lam, 1-lam])

def params(d, lam, c):
    q = 1-lam
    v = (q*c/(lam+q*c))**(1/(d+1))
    r = (1-v**d)/(d*(1-v**(d+1)))
    s = v*r
    a, b = 1/r, 1/s-1/r
    A = lam/(gamma(d+1)*(r**d-s**d))
    w = (q-v**d)/(q*(1-v**d))
    return v,r,s,a,b,A,w

def dc_atom(d, lam, c):
    v = params(d,lam,c)[0]
    return lam*(d*np.log1p(-v**(d+1))-(d+1)*np.log1p(-v**d))

def dc_thermal(d,lam,c):
    return h(lam)-d*math.log(lam+(1-lam)*c)+(1-lam)*d*math.log(c)

def matrix_entropy(mat):
    ev = np.linalg.eigvalsh(mat)
    return entropy(np.maximum(ev,0))

out = {"purpose":"finite diagnostics only; no exact-arithmetic or formal verification", "seed":2026102456,
       "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__},
       "input_sha256":hashlib.sha256((ROOT/'input56.md').read_bytes()).hexdigest()}

# Noncommuting finite quantum ensembles: pure-letter identity, Fano/Jensen,
# flag upper bound. No numerical optimization of a quantum capacity.
qr=[]
for n in [2,3,7,15]:
  for lam in [.07,.4,.83]:
    q=1-lam
    for rep in range(4):
      z=RNG.normal(size=(n,n))+1j*RNG.normal(size=(n,n))
      tau=z@z.conj().T; tau/=np.trace(tau)
      k=n+3
      z=RNG.normal(size=(n,k))+1j*RNG.normal(size=(n,k))
      z/=np.linalg.norm(z,axis=0)
      weights=RNG.dirichlet(np.ones(k))
      letters=[np.outer(z[:,i],z[:,i].conj()) for i in range(k)]
      rho=sum(weights[i]*letters[i] for i in range(k))
      st=matrix_entropy(tau)
      conditional=np.array([h(lam)+q*st-matrix_entropy(lam*x+q*tau) for x in letters])
      indiv=np.array([h(q*float(np.trace(x@tau).real)) for x in letters])
      mix=lam*rho+q*tau
      J=matrix_entropy(mix)-q*st-h(lam)
      chi=matrix_entropy(mix)-sum(weights[i]*matrix_entropy(lam*letters[i]+q*tau) for i in range(k))
      avg=float(weights@conditional)
      overlap=float(np.trace(rho@tau).real)
      qr.append({"n":n,"lambda":lam,"rep":rep,"identity_error":abs(chi-J-avg),
         "minimum_conditional":float(min(conditional)),"max_individual_fano_violation":float(max(conditional-indiv)),
         "jensen_violation":avg-h(q*overlap),"flag_violation":chi-lam*matrix_entropy(rho)})
out['quantum_cases']=qr

# Symbolic checks with arbitrary v and lambda, exact positive integer d.
v,l=sp.symbols('v l',positive=True)
symbolic=[]
for dd in range(1,7):
    r=(1-v**dd)/(dd*(1-v**(dd+1)))
    s=v*r
    cc=l*v**(dd+1)/((1-l)*(1-v**(dd+1)))
    alpha,beta=1/r,1/s-1/r
    symbolic.append({'d':dd,'input_mean_minus_one':str(sp.factor(dd*r*(1-v**(dd+1))/(1-v**dd)-1)),
        'value_energy_identity':str(sp.factor(l*alpha-(1-l)*beta*cc-l*dd)),
        'reset_mass_plus_atom_minus_one':str(sp.factor(l*v**dd/((1-l)*(1-v**dd))+((1-l)-v**dd)/((1-l)*(1-v**dd))-1))})
out['symbolic_integer_d']=symbolic

# Continuum identities are additionally integrated directly. General proof in report.
continuum=[]
atoms=[]
for d in [.25,.5,1.,2.,3.]:
  for lam in [.2,.6,.9]:
    q=1-lam
    cs=lam*q**(1/d)/(1-q**(1+1/d))
    for frac in [.5,.95]:
      c=frac*cs
      v,r,s,a,b,A,w=params(d,lam,c)
      # Integrate after x=y^(1/d), making nu_d(dx)=dy. Exponential decay.
      fun=lambda y: (A/lam)*math.exp(-a*y**(1/d))*(-math.expm1(-b*y**(1/d)))
      tcon=lambda y: (A/q)*math.exp(-(a+b)*y**(1/d))
      vals=[quad(fun,0,np.inf,epsabs=2e-10)[0],quad(lambda y:y**(1/d)*fun(y),0,np.inf,epsabs=2e-10)[0],
            quad(tcon,0,np.inf,epsabs=2e-10)[0],quad(lambda y:y**(1/d)*tcon(y),0,np.inf,epsabs=2e-10)[0]]
      continuum.append({'d':d,'lambda':lam,'c':c,'cstar':cs,'fraction':frac,'v':v,'atom':w,
          'quadrature_input_mass_error':vals[0]-1,'quadrature_input_mean_error':vals[1]-1,
          'quadrature_reset_mass_error':vals[2]+w-1,'quadrature_reset_mean_error':vals[3]-c,
          'energy_identity_error':lam*a-q*b*c-lam*d,
          'join_error_nats':dc_atom(d,lam,cs)-dc_thermal(d,lam,cs)})
      for N in [64,256]:
        L=N**(1/d)
        Jmax=max(100,int(math.ceil(N*(42/a)**d)))
        x=(np.arange(1,Jmax+1,dtype=float)/N)**(1/d)
        ea=np.exp(-a*x); eb=np.exp(-(a+b)*x)
        u=ea*(-np.expm1(-b*x)); Z=float(np.sum(u))
        p=np.r_[0,u/Z]
        tc=(lam/q)*eb/Z; t=np.r_[1-float(sum(tc)),tc]
        m=lam*p+q*t
        ep=float(np.dot(p[1:],x)); et=float(np.dot(t[1:],x))
        F=entropy(m)-q*entropy(t)-h(lam)
        identity=lam*math.log(Z)+lam*a*ep-q*b*et
        K=lam*math.log(lam/A)+lam*d
        atoms.append({'d':d,'lambda':lam,'fraction':frac,'c':c,'N':N,'L':L,'cutoff_j':Jmax,
           't0':float(t[0]),'t0_limit':w,'input_mean_over_L':ep,'reset_mean_over_L':et,
           'ground_KKT_margin':q*float(t[0])-lam/Z,'identity_error':abs(F-identity),
           'value_limit_error_nats':F-lam*d*math.log(L)-K,
           'finite_truncation_note':'cutoff alpha*x >=42, not an infinite-tail interval certificate'})
out['continuum_quadrature']=continuum
out['atomic_finite_sums']=atoms

# Exact finite-support Gibbs saddles, with all competitors subjected to same
# finite support and energy budgets. Global inequalities are proved in report.
thermal=[]
tangents=[]
for d in [.5,1.,2.]:
  for lam in [.2,.6,.9]:
    q=1-lam
    cs=lam*q**(1/d)/(1-q**(1+1/d)); c=(1+cs)/2; mu=lam+q*c
    N=128; E=N**(1/d); jj=max(128,int(math.ceil(N*(40/d)**d)))
    H=np.arange(jj+1,dtype=float)**(1/d)
    def gibbs(mean):
      def state(bb):
        a=np.exp(-bb*H); return a/sum(a)
      bb=brentq(lambda z:float(np.dot(H,state(z)))-mean,1e-10,100,xtol=1e-15,rtol=1e-14)
      return bb,state(bb)
    bt,t=gibbs(c*E); bm,m=gibbs(mu*E); p=(m-q*t)/lam
    F=entropy(m)-q*entropy(t)-h(lam)
    gp=-lam*(np.log(m)+1)
    gt=q*np.log(t/m)
    thermal.append({'d':d,'lambda':lam,'c':c,'E':E,'cutoff_j':jj,
       'ground_positivity_ratio':float(q*t[0]/m[0]),'p_min':float(min(p)),
       'input_mass_error':float(sum(p)-1),'input_energy_error':float(np.dot(H,p)-E),
       'reset_energy_error':float(np.dot(H,t)-c*E),
       'input_affine_gradient_error':float(max(abs(gp-(gp[0]+lam*bm*H)))),
       'reset_affine_gradient_error':float(max(abs(gt-(gt[0]-q*(bt-bm)*H))))})
    # Random competitors use only the first 100 levels, scaled toward ground to meet budget.
    for rep in range(10):
      pc=np.zeros(jj+1); tc=np.zeros(jj+1)
      pc[:100]=RNG.dirichlet(np.ones(100)); tc[:100]=RNG.dirichlet(np.ones(100))
      for arr,budget in [(pc,E),(tc,c*E)]:
        energy=float(np.dot(H,arr))
        z=min(1,budget/energy); arr*=z; arr[0]+=1-z
      Fpc=entropy(lam*pc+q*t)-q*entropy(t)-h(lam)
      Ftc=entropy(lam*p+q*tc)-q*entropy(tc)-h(lam)
      tangents.append({'d':d,'lambda':lam,'rep':rep,'maximizer_violation':Fpc-F,'minimizer_violation':F-Ftc})
out['thermal_finite_support']=thermal
out['saddle_competitors']=tangents

# Endpoint ratios are diagnostic only: E is not present in these continuum functions.
endpoints=[]
for d in [.25,.5,1.,2.,3.]:
  for lam in [.2,.6,.9]:
    q=1-lam; cs=lam*q**(1/d)/(1-q**(1+1/d))
    cc=cs*1e-9
    small=(d+1)*lam**(1/(d+1))*q**(d/(d+1))*cc**(d/(d+1))
    eps=1e-4
    C=quad(lambda y:y**(d-1)*math.exp(-y)/(1+math.sqrt(-math.expm1(-y))),0,np.inf,epsabs=2e-10)[0]
    endpoints.append({'d':d,'lambda':lam,'small_c_ratio':dc_atom(d,lam,cc)/small,
       'near_one_quadratic_ratio':(h(lam)-dc_thermal(d,lam,1-eps))/(lam*q*d*eps*eps/2),
       'conditional_EA_coefficient_ratio':(2*C/gamma(d))**(1/(d+1)),
       'strict_ratio_lower':1.,'strict_ratio_upper':2**(1/(d+1))})
out['endpoint_function_checks']=endpoints
out['summary']={
 'quantum_case_count':len(qr),'max_quantum_identity_error':max(z['identity_error'] for z in qr),
 'max_fano_violation':max(z['max_individual_fano_violation'] for z in qr),
 'max_jensen_violation':max(z['jensen_violation'] for z in qr),'max_flag_violation':max(z['flag_violation'] for z in qr),
 'continuum_case_count':len(continuum),'atomic_case_count':len(atoms),
 'min_atomic_ground_margin':min(z['ground_KKT_margin'] for z in atoms),
 'max_atomic_identity_error':max(z['identity_error'] for z in atoms),
 'thermal_case_count':len(thermal),'finite_competitor_count':len(tangents),
 'max_thermal_ground_ratio':max(z['ground_positivity_ratio'] for z in thermal),
 'max_saddle_violation':max(max(z['maximizer_violation'],z['minimizer_violation']) for z in tangents),
 'note':'Finite numerical successes cannot prove theorem, positivity eventually, tails, or global optimality.'}
(ROOT/'checks56-results.json').write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(out['summary'],indent=2))
