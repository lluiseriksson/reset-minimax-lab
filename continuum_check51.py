"""Independent integral evaluations and larger finite checks; not a limit proof."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import math,json,hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares,brentq
from scipy.special import digamma,beta as beta_fn
import finite_game50 as game
ROOT=Path(__file__).parent
LAM=game.LAM;Q=game.Q;LN2=math.log(2)
original=json.loads((ROOT/'finite-game-results50.json').read_text())
continuum=[];large=[]
def moments(theta):
    alpha,beta,delta=np.exp(theta)
    def f(x):return math.exp(-alpha*x)*math.sqrt(-math.expm1(-delta-beta*x))
    def t(x):return math.exp(-(alpha+beta)*x-delta)/math.sqrt(-math.expm1(-delta-beta*x))
    ip=quad(f,0,np.inf,epsabs=1e-11,epsrel=1e-11)[0]
    it=quad(t,0,np.inf,epsabs=1e-11,epsrel=1e-11)[0]
    mp=quad(lambda x:x*f(x),0,np.inf,epsabs=1e-11,epsrel=1e-11)[0]/ip
    mt=quad(lambda x:x*t(x),0,np.inf,epsabs=1e-11,epsrel=1e-11)[0]/it
    return ip,it,mp,mt
for c in (.2,.5,.8,.99):
    initial=next(x for x in original['rows'] if x['E']==512 and x['c']==c)
    start=np.log([initial['alpha'],initial['beta'],initial['delta']])
    def residual(theta):
        ip,it,mp,mt=moments(theta)
        return [math.log(LAM/Q*it/ip),mp-1,mt-c]
    result=least_squares(residual,start,xtol=1e-12,gtol=1e-12,ftol=1e-12,max_nfev=200)
    alpha,beta,delta=np.exp(result.x);ip,it,mp,mt=moments(result.x)
    phi=2*LAM*math.log(ip)+2*LAM*alpha-Q*beta*c-Q*delta
    continuum.append({'branch':'diffuse','c':c,'alpha':float(alpha),'beta':float(beta),'delta':float(delta),'residual_max':float(np.max(np.abs(residual(result.x)))),'candidate_gap_bits':(2*LAM-phi)/LN2,'solver_success':bool(result.success)})
sstar=Q/(2*LAM)
A=lambda s:float(digamma(s+1.5)-digamma(s))
c_of_s=lambda s:2*LAM/Q*(s-1/A(s))
cstar=c_of_s(sstar)
for c in (.05,.1,cstar):
    s=brentq(lambda u:c_of_s(u)-c,1e-8,sstar+1e-8,xtol=1e-13)
    b=A(s);w=max(0.,1-2*LAM*s/Q)
    continuum.append({'branch':'atom_or_boundary','c':c,'s':s,'beta':b,'candidate_reset_atom':w,'candidate_gap_bits':2*LAM*math.log2(b/beta_fn(s,1.5)),'mean_residual':abs(c_of_s(s)-c)})
for c in (.05,.1,.2):
    previous=next(x for x in original['rows'] if x['E']==512 and x['c']==c)
    guess=np.log([previous['alpha'],previous['beta'],previous['delta']])
    for E in (2048,8192):
        row,guess=game.solve(E,c,guess);large.append(row)
        print(json.dumps({'E':E,'c':c,'checks_pass':row['passed_numerical_checks'],'reset_p0':row['reset_mass_at_zero'],'gap':row['oscillator_benchmark_minus_game_bits']}),flush=True)
packet={'scope':'Numerical integral and finite-grid comparison only; continuum proof51 separately under audit. No infinite-dimensional capacity certificate or score.','input_finite_results_sha256':hashlib.sha256((ROOT/'finite-game-results50.json').read_bytes()).hexdigest(),'lambda':LAM,'candidate_boundary_ratio':cstar,'continuum_integral_evaluations':continuum,'larger_finite_grid_checks':large,'all_large_checks_passed':all(x['passed_numerical_checks'] for x in large)}
out=ROOT/'continuum-results51.json';out.write_text(json.dumps(packet,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'file':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'all_large_pass':packet['all_large_checks_passed']}),flush=True)
