"""Finite, high precision diagnostics of the exact continuum profiles in58.
Not a proof, not interval arithmetic, and no physical meaning for intermediate a.
"""
from pathlib import Path
import json, hashlib, platform, time
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parent
mp.mp.dps=65
START=time.time()
NTERMS=480

la,dd,aa=sp.symbols('lambda d a',positive=True)
qq=1-la; ee=1/(aa+1); kk=aa/(aa+1)
J=sp.Matrix([[0,-1/la,-1/la],[-1/dd,ee*qq/(la*dd),0],[-1/dd,-(la+kk*qq)/(la*dd),0]])
target=sp.Matrix([0,0,-1])
direction=J.inv()*target
y=sp.symbols('y')
Q=y*y-2*(dd+1)*y+dd*(dd+1)
def gamma_average(poly):
    return sp.factor(sum(co*sp.rf(dd,pow_[0]) for pow_,co in sp.Poly(sp.expand(poly),y).terms()))
u,w=sp.symbols('u w')
solution=sp.solve([(aa+la)*u+qq*w-qq**2/2,w-u-(1+qq)/2],[u,w])
out={'purpose':'Finite diagnostics only; analytic proof is in report58.md.',
 'input_sha256':hashlib.sha256((ROOT/'input58.md').read_bytes()).hexdigest(),
 'addendum_sha256':hashlib.sha256((ROOT/'input58-addendum.md').read_bytes()).hexdigest(),
 'runtime':{'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__},
 'precision_decimal_digits':mp.mp.dps,'series_terms':NTERMS,
 'symbolic':{'jacobian_determinant':str(sp.factor(J.det())),
    'parameter_first_derivatives':[str(sp.factor(z)) for z in direction],
    'Q_mass':str(gamma_average(Q)),'xQ_mean':str(gamma_average(y*Q)/dd),
    'Q_squared_norm':str(gamma_average(Q*Q)),
    'second_order_u_coefficient':str(sp.factor(solution[u])),
    'second_order_w_coefficient':str(sp.factor(solution[w])),
    'integrated_entropy_factor':str(sp.integrate(1/(aa+1)**2,(aa,0,1)))},
 'numerical_limitations':['No interval certification of truncated series tails.','Finite grid of lambda,d,a,epsilon; no global parameter validation.','Exact profiles evaluated, not fitted to determine theorem coefficient.']}

def solve_profile(d,lam,a,eps):
    d,lam,a,eps=map(mp.mpf,map(str,[d,lam,a,eps]))
    q=1-lam; e=1/(a+1); k=a/(a+1)
    cu=[mp.mpf(1)]; cv=[mp.mpf(0),mp.mpf(1)]
    for n in range(1,NTERMS): cu.append(cu[-1]*(n-1-e)/n)
    for n in range(2,NTERMS): cv.append(cv[-1]*(k+n-2)/(n-1))
    def moments(alpha,beta,delta):
        z=mp.exp(-delta); zpow=mp.mpf(1)
        U=U1=V=V1=mp.mpf(0)
        for n in range(NTERMS):
            den=alpha+n*beta
            common=zpow/den**d
            uc=cu[n]*common; vc=cv[n]*common
            U+=uc; U1+=uc/den; V+=vc; V1+=vc/den
            zpow*=z
        return mp.gamma(d+1)*U,d*U1/U,mp.gamma(d+1)*V,d*V1/V
    def eq(alpha,beta,delta):
        U,meanp,V,meant=moments(alpha,beta,delta)
        return (lam*V/(q*U)-1,meanp-1,meant-(1-eps))
    init=(d+e*q*d*eps,lam*d*eps,-mp.log(q)-lam*d*eps)
    pars=mp.findroot(eq,init,tol=mp.mpf('1e-40'),maxsteps=30)
    alpha,beta,delta=pars
    U,meanp,V,meant=moments(*pars)
    G=mp.gamma(d+1)/d**d
    value=(a+1)*lam*mp.log(U)+(a+1)*lam*alpha-q*beta*(1-eps)-q*delta
    deficit=(a+1)*lam*(d+mp.log(G))-value
    row={'d':str(d),'lambda':str(lam),'a':str(a),'epsilon':str(eps),
       'alpha':mp.nstr(alpha,45),'beta':mp.nstr(beta,45),'delta':mp.nstr(delta,45),
       'normalization_error':mp.nstr(lam*V/(q*U)-1,8),'input_mean_error':mp.nstr(meanp-1,8),
       'reset_mean_error':mp.nstr(meant-(1-eps),8),
       'deficit_nats':mp.nstr(deficit,50),'positive_domain':bool(alpha>0 and beta>0 and delta>0)}
    if a==0:
        mu=lam+q*(1-eps)
        expected=-lam*mp.log(lam)-q*mp.log(q)-d*mp.log(mu)+q*d*mp.log(1-eps)
        row['C1_closed_formula_error']=mp.nstr(deficit-expected,12)
    # A direct entropy integral at intermediate a, independent of value subtraction.
    if a==mp.mpf('.5'):
        def integrand(x):
            b=1-mp.exp(-delta-beta*x)
            p=mp.exp(-alpha*x)*b**e/U
            logratio=mp.log(G/U)+(d-alpha)*x+e*mp.log(b)
            return d*x**(d-1)*p*logratio
        rel=mp.quad(integrand,[0,1,mp.inf])
        coeff=q*q*d*(d+1)/(4*(a+1)**2)
        row['direct_relative_entropy']=mp.nstr(rel,45)
        row['entropy_quartic_ratio']=mp.nstr(rel/(coeff*eps**4),30)
    return row,deficit

cases=[]; contrasts=[]
for d in [.5,1.,2.]:
  for lam in [.3,.7]:
    for eps in [.1,.05,.025]:
      ds={}
      for a in [0.,.5,1.]:
        row,val=solve_profile(d,lam,a,eps)
        cases.append(row); ds[a]=val
      dmp,lmp,emp=map(mp.mpf,map(str,[d,lam,eps])); qmp=1-lmp
      coef=lmp*qmp*qmp*dmp*(dmp+1)/8
      contrasts.append({'d':d,'lambda':lam,'epsilon':eps,
        'difference_nats':mp.nstr(ds[1.]-ds[0.],50),'predicted_quartic_coefficient_nats':mp.nstr(coef,35),
        'ratio_to_predicted_quartic':mp.nstr((ds[1.]-ds[0.])/(coef*emp**4),30)})
      print('completed',d,lam,eps,flush=True)
out['profiles']=cases; out['contrasts']=contrasts
out['summary']={'profiles_count':len(cases),'contrasts_count':len(contrasts),
 'all_positive_domain':all(z['positive_domain'] for z in cases),
 'max_moment_residual':mp.nstr(max(abs(mp.mpf(z[k])) for z in cases for k in ['normalization_error','input_mean_error','reset_mean_error']),12),
 'max_C1_closed_formula_error':mp.nstr(max(abs(mp.mpf(z['C1_closed_formula_error'])) for z in cases if 'C1_closed_formula_error' in z),12),
 'max_quartic_ratio_at_epsilon_0025':mp.nstr(max(mp.mpf(z['ratio_to_predicted_quartic']) for z in contrasts if z['epsilon']==.025),20),
 'min_quartic_ratio_at_epsilon_0025':mp.nstr(min(mp.mpf(z['ratio_to_predicted_quartic']) for z in contrasts if z['epsilon']==.025),20),
 'elapsed_seconds':time.time()-START}
(ROOT/'checks58-results.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(json.dumps(out['summary'],indent=2),flush=True)
