"""Finite diagnostics of the written RV extension; not a limit proof.
All probabilities below are per eigenstate. Multiplicities weight every sum.
"""
import pathlib,json,math,warnings
import numpy as np
import scipy
from scipy.integrate import quad,IntegrationWarning
from scipy.optimize import brentq
from scipy.special import gamma,xlogy
ROOT=pathlib.Path(__file__).resolve().parent
lam=.6;q=1-lam;W=[]
def integrate(fun,a,b,**kw):
 with warnings.catch_warnings(record=True) as ws:
  warnings.simplefilter('always',IntegrationWarning)
  val,err=quad(fun,a,b,epsabs=2e-10,epsrel=2e-10,limit=300,**kw)
 W.extend(str(w.message) for w in ws)
 return val
def U(k,s,delta):
 def ker(y):return math.exp(-s*y-delta-y)/(1+math.sqrt(-math.expm1(-delta-y)))
 co=integrate(ker,0,1,weight='alg',wvar=(k-1,0))+integrate(lambda y:y**(k-1)*ker(y),1,np.inf)
 return gamma(k)*s**(-k)-co
def V(k,s,delta):
 if delta==0:
  def reg(y):return 1 if y==0 else math.exp(-(s+1)*y)/math.sqrt(-math.expm1(-y)/y)
  a=integrate(reg,0,1,weight='alg',wvar=(k-1.5,0))
 else:
  a=integrate(lambda y:math.exp(-(s+1)*y-delta)/math.sqrt(-math.expm1(-delta-y)),0,1,weight='alg',wvar=(k-1,0))
 return a+integrate(lambda y:y**(k-1)*math.exp(-(s+1)*y-delta)/math.sqrt(-math.expm1(-delta-y)),1,np.inf)
def profile(d,atomic):
 delta=0. if atomic else .2
 s=.08 if atomic else brentq(lambda ss:lam/q*V(d,ss,delta)/U(d,ss,delta)-1,1e-5,100)
 ud=U(d,s,delta);ud1=U(d+1,s,delta);beta=ud1/ud;alpha=s*beta
 I=d*ud/beta**d;c=lam/q*V(d+1,s,delta)/ud1
 w=1-lam/q*V(d,s,delta)/ud
 phi=2*lam*math.log(I)+2*lam*alpha-q*beta*c-q*delta
 assert w>=-1e-10
 return dict(d=d,delta=delta,s=s,alpha=alpha,beta=beta,I=I,c=c,w=w,phi=phi)
def counting(x,d,A,mode,r):
 n=np.floor(x);positive=np.where(n>0,A*np.maximum(n,1)**d*(np.log(np.maximum(n,1))+1 if mode=='log' else 1),0)
 return r+np.floor(positive)
models=[(.25,False,'log',3.,7),(.5,False,'log',3.,7),(1.,True,'log',1.,1),(1.,True,'log',1.,7),(1.,True,'power',3.,7),(1.,False,'log',1.,7)]
rows=[];profiles=[]
for d,atomic,mode,A,r in models:
 pr=profile(d,atomic);profiles.append(pr|{'mode':mode,'A':A,'r':r});alpha=pr['alpha'];beta=pr['beta']
 for L in [64,256,1024,4096]:
  cutoff=max(45/alpha,45/d);levels=np.arange(int(math.ceil(cutoff*L))+1,dtype=float)
  counts=counting(levels,d,A,mode,r);mult=np.diff(counts,prepend=0)
  mask=mult>0;levels=levels[mask];mult=mult[mask];x=levels/L
  NL=float(counting(L,d,A,mode,r))
  def terms(logdelta):
   delta=math.exp(logdelta);h=delta+beta*x;sq=np.sqrt(-np.expm1(-h));u=np.exp(-alpha*x)*sq;v=np.exp(-alpha*x-h)/sq
   return delta,u,v
  def balance(logdelta):
   delta,u,v=terms(logdelta);return math.log(lam/q)+math.log(mult@v)-math.log(mult@u)
  ld=brentq(balance,-100,math.log(3.),xtol=1e-12);delta,u,v=terms(ld)
  Z=mult@u;p=u/Z;t=lam/q*v/Z;m=lam*p+q*t
  mp=mult@(x*p);mt=mult@(x*t)
  hh=-(lam*math.log(lam)+q*math.log(q))
  ent=lambda a:-float(mult@xlogy(a,a))
  F=lam*ent(p)+ent(m)-q*ent(t)-hh
  exact=2*lam*math.log(Z)+2*lam*alpha*mp-q*beta*mt-q*delta
  identity=max(abs(F-exact),abs(mult@p-1),abs(mult@t-1),abs(balance(ld)),float(np.max(np.abs((p*m)/(lam/Z**2*np.exp(-2*alpha*x))-1))),float(np.max(np.abs((t/m)/(np.exp(-delta-beta*x)/q)-1))))
  assert identity<1e-9,identity
  ww=np.exp(-d*x);ZG=mult@ww;EG=float(mult@(levels*ww)/ZG);SG=d*EG/L+math.log(ZG)
  G=gamma(d+1)/d**d
  rec={'d':d,'atomic':atomic,'mode':mode,'A':A,'ground_degeneracy':r,'L':L,'N_L':NL,'cutoff_scaled':cutoff,
       'positive_distinct_levels_in_truncation':len(x)-1,'delta_L':delta,'normalizer_ratio_to_limit':float(Z/NL/pr['I']),
       'input_scaled_mean':float(mp),'reset_scaled_mean':float(mt),'target_reset_mean':pr['c'],
       'total_ground_reset_mass':float(r*t[0]),'per_ground_reset_mass':float(t[0]),'target_ground_mass':pr['w'],
       'centered_value_residual':float(F-2*lam*math.log(NL)-pr['phi']),
       'ground_functional_contribution':float(lam*(-r*xlogy(p[0],p[0]))-r*xlogy(m[0],m[0])+q*r*xlogy(t[0],t[0])),
       'ground_contribution_limit':float(-q*pr['w']*math.log(q)),
       'gibbs_partition_ratio':float(ZG/NL/G),'gibbs_scaled_energy':EG/L,
       'gibbs_entropy_constant_residual':SG-math.log(float(counting(EG,d,A,mode,r)))-d-math.log(G),
       'max_exact_identity_residual':identity}
  if atomic:
   rec['singular_positive_level_sums']=[{'a':a,'sum':float(mult[(x>0)&(x<=a)]@(x[(x>0)&(x<=a)]**-.5)/NL),'limit':d/(d-.5)*a**(d-.5)} for a in [.01,.05,.1]]
  rows.append(rec)

# Exact falsifier of any use of x^-1/2 UI at d=1/2: H_j=j^2,
# L=n^2 gives sqrt(L)/(floor(sqrt(L))+1)*sum_{j<=sqrt(aL)}1/j.
critical=[]
for n in [100,1000,10000]:
 a=.01;k=int(math.sqrt(a)*n)
 critical.append({'n':n,'L':n*n,'a':a,'singular_sum':n/(n+1)*float(np.sum(1/np.arange(1,k+1)))})
assert all(critical[i+1]['singular_sum']>critical[i]['singular_sum'] for i in range(2))

# A slowly varying oscillatory factor need not converge. This falsifies
# replacing ln N(E) by d ln E plus a fixed additive constant.
slow_osc=[]
for t in [math.pi/2,3*math.pi/2,5*math.pi/2,7*math.pi/2]:
 logE=math.exp(t)
 slow_osc.append({'ln_ln_E':t,'ln_ell_E':math.sin(t),'ln_ell_2E_minus_ln_ell_E':math.sin(math.log(logE+math.log(2)))-math.sin(t)})

out={'audit_id':'regular57','scope':'Finite diagnostics only. Infinite tails and convergence are proved analytically in report57. Fixed cutoff is not a proof.',
 'versions':{'numpy':np.__version__,'scipy':scipy.__version__},'profiles':profiles,'rows':rows,
 'critical_d_half_wrong_UI_falsifier':critical,'slow_oscillation_wrong_centering_falsifier':slow_osc,
 'quadrature_warnings':W,'summary':{'models':len(models),'cases':len(rows),'maximum_exact_identity_residual':max(r['max_exact_identity_residual'] for r in rows),'warnings':len(W),'all_assertions_passed':True}}
(ROOT/'results/checks57.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out['summary'],indent=2))
print(json.dumps([{k:r[k] for k in ['d','mode','ground_degeneracy','L','normalizer_ratio_to_limit','centered_value_residual','total_ground_reset_mass','target_ground_mass','gibbs_entropy_constant_residual']} for r in rows if r['L']==4096],indent=2))
