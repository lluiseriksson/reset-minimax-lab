"""Original finite entropy-game diagnostic; no infinite-dimensional capacity certificate."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,math,time,hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,linprog
from scipy.special import logsumexp,digamma

LAM=.6;Q=1-LAM;LN2=math.log(2)
def entropy(p):
    p=p[p>0]
    return float(-np.dot(p,np.log(p))/LN2)
def g(E):return (math.log1p(E)+E*math.log1p(1/E))/LN2
def binary(x):return -(x*math.log(x)+(1-x)*math.log1p(-x))/LN2
def distributions(parameters,x):
    alpha,beta,delta=np.exp(parameters)
    v=delta+beta*x
    log_one_minus_z=np.log(-np.expm1(-v))
    logp=-alpha*x+.5*log_one_minus_z
    logt=-alpha*x-v-.5*log_one_minus_z
    zp=logsumexp(logp);zt=logsumexp(logt)
    p=np.exp(logp-zp);t=np.exp(logt-zt)
    return p,t,math.log(LAM/Q)+zt-zp
def solve(E,c,guess=None,xmax=40):
    x=np.arange(int(xmax*E)+1,dtype=float)/E
    def residual(parameters):
        p,t,norm=distributions(parameters,x)
        return np.array([norm,np.dot(p,x)-1,np.dot(t,x)-c])
    if guess is None:guess=np.log([.8,max(.02,1/c-1),max(1e-6,c*c)])
    attempts=[]
    starts=[guess,np.log([.8,2,.001]),np.log([1,1,.001]),np.log([.5,.5,.05])]
    for initial in starts:
        result=least_squares(residual,initial,bounds=([-10,-15,-60],[10,12,5]),xtol=2e-13,ftol=2e-13,gtol=2e-13,max_nfev=2500)
        attempt_residual=float(np.max(np.abs(residual(result.x))))
        attempts.append({'initial_parameters':initial.tolist(),'solution_parameters':result.x.tolist(),'solver_success':bool(result.success),'residual_max':attempt_residual,'nfev':int(result.nfev)})
        if result.success and attempt_residual<1e-8:break
    p,t,norm=distributions(result.x,x);m=LAM*p+Q*t
    alpha,beta,delta=np.exp(result.x)
    residual_max=float(np.max(np.abs(residual(result.x))))
    # Compute log probabilities directly, including tails below floating underflow.
    v=delta+beta*x;lo=np.log(-np.expm1(-v))
    lp=-alpha*x+.5*lo;lp-=logsumexp(lp)
    lt=-alpha*x-v-.5*lo;lt-=logsumexp(lt)
    lm=np.logaddexp(math.log(LAM)+lp,math.log(Q)+lt)
    input_kkt=lp+lm+2*alpha*x
    reset_kkt=lt-lm+beta*x
    kkt_shape=max(float(np.ptp(input_kkt)),float(np.ptp(reset_kkt)))
    gp=-LAM*(lp+lm+2)/LN2;gt=Q*(lt-lm)/LN2
    input_lp=linprog(-gp,A_ub=[x],b_ub=[1],A_eq=[np.ones(len(x))],b_eq=[1],bounds=(0,None),method='highs')
    reset_lp=linprog(gt,A_ub=[x],b_ub=[c],A_eq=[np.ones(len(x))],b_eq=[1],bounds=(0,None),method='highs')
    input_gap=float(-input_lp.fun-np.dot(gp,p)) if input_lp.success else None
    reset_gap=float(np.dot(gt,t)-reset_lp.fun) if reset_lp.success else None
    gap_pass=(input_gap is not None and reset_gap is not None and abs(input_gap)<1e-7 and abs(reset_gap)<1e-7)
    phi=LAM*entropy(p)+entropy(m)-Q*entropy(t)-binary(LAM)
    # This purification-test correction applies to diagonal input p and reset t.
    e=Q*float(np.dot(p*p,t))
    row={'E':E,'c':c,'xmax':xmax,'dimension':len(x),'solver_success':bool(result.success),'nfev':int(result.nfev),'attempts':attempts,'constraint_residual_max':residual_max,'kkt_shape_residual':kkt_shape,'input_linear_best_response_gap':input_gap,'reset_linear_best_response_gap':reset_gap,'alpha':float(alpha),'beta':float(beta),'delta':float(delta),'input_mass_at_zero':float(p[0]),'reset_mass_at_zero':float(t[0]),'input_scaled_mean':float(np.dot(p,x)),'reset_scaled_mean':float(np.dot(t,x)),'finite_entropy_game_value_bits':phi,'oscillator_benchmark_minus_game_bits':2*LAM*g(E)-phi,'pointwise_quantum_correction_upper_bits':binary(e) if e>0 else 0.,'pointwise_purification_error':e,'last_level_reset_probability':float(t[-1]),'last_level_input_probability':float(p[-1]),'passed_numerical_checks':bool(result.success and residual_max<1e-8 and kkt_shape<1e-8 and gap_pass)}
    return row,result.x
def main():
    started=time.time();rows=[]
    for c in (.05,.1,.12823752863469068,.2,.5,.8,.99):
        guess=None
        for E in (8,32,128,512):
            row,guess=solve(E,c,guess);rows.append(row)
            print(json.dumps({'E':E,'c':c,'pass':row['passed_numerical_checks'],'residual':row['constraint_residual_max'],'reset_p0':row['reset_mass_at_zero'],'gap':row['oscillator_benchmark_minus_game_bits']}),flush=True)
    s=Q/(2*LAM)
    cstar=float((digamma(s+1.5)-digamma(s+1))/(digamma(s+1.5)-digamma(s)))
    packet={'status':'finite_numerical_diagnostic_not_a_capacity_or_limit_theorem','lambda':LAM,'candidate_continuum_boundary_ratio':cstar,'scope':'Finite probability simplex game: maximize in p and minimize in t, both constrained by mean energy and normalized. Positive stationary points satisfy formal saddle KKT equations. Floating residuals are not interval certificates; cutoff and continuum/quantum minimax limits remain unproved.','not_claimed':['global quantum capacity curve','optimal infinite-dimensional reset','physical phase transition','novelty','AIRR score'],'rows':rows,'all_numerical_checks_passed':all(r['passed_numerical_checks'] for r in rows),'elapsed_seconds':time.time()-started}
    out=Path(__file__).with_name('finite-game-results50.json');out.write_text(json.dumps(packet,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'result_file':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'all_pass':packet['all_numerical_checks_passed']}),flush=True)
if __name__=='__main__':main()
