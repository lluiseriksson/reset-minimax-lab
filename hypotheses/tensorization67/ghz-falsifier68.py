"""Larger GHZ blocks: exact diagonal-plus-rank-one reduction, finite tests only."""
from pathlib import Path
import hashlib,itertools,json,math,time,importlib.util
import numpy as np
from scipy.optimize import minimize
C=Path(__file__).parent
spec=importlib.util.spec_from_file_location('f66',C/'entropy-falsifier66.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
rng=np.random.default_rng(2026102468)

def fast(weights,tau,lam,n,d,perms):
    # A partial trace over any GHZ site kills all off-diagonal terms. Only the
    # no-replacement branch retains coherence, of rank one on d code words.
    coords=np.array(perms).T
    diag=np.zeros([d]*n)
    diag[tuple(coords.T)]=weights
    K=lam*np.eye(d)+(1-lam)*np.array(tau)[:,None]*np.ones((1,d))
    for i in range(n):
        diag=np.moveaxis(np.tensordot(K,diag,axes=(1,i)),0,i)
    support=np.ravel_multi_index(tuple(coords.T),[d]*n)
    flat=diag.ravel();base=flat[support]-lam**n*weights
    block=np.diag(base)+lam**n*np.outer(np.sqrt(weights),np.sqrt(weights))
    ev=np.linalg.eigvalsh(block)
    remaining=np.ones(d**n,bool);remaining[support]=False
    spectrum=np.concatenate((flat[remaining],ev))
    assert spectrum.min()>-2e-12 and abs(spectrum.sum()-1)<2e-11
    positive=spectrum[spectrum>1e-15];lhs=float(-positive@np.log(positive))
    t=[float(np.dot(weights,np.array(tau)[perms[i]])) for i in range(n)]
    st=float(-sum(s*math.log(s) for s in tau if s>0))
    rhs=n*((1-lam)*st+m.h(lam))-sum(m.h((1-lam)*v) for v in t)
    return dict(n=n,d=d,lam=float(lam),weights=list(map(float,weights)),tau=list(map(float,tau)),
                permutations=[list(map(int,p)) for p in perms],overlaps=t,lhs=lhs,rhs=rhs,margin=lhs-rhs)

start=time.time();cross=[]
for n,d in [(2,2),(2,3),(3,2),(3,3)]:
    for j in range(8):
        w=rng.dirichlet(np.ones(d));t=rng.dirichlet(np.ones(d));lam=rng.uniform(.01,.99)
        p=[rng.permutation(d) for _ in range(n)];r=fast(w,t,lam,n,d,p)
        x=np.zeros(d**n,complex)
        for k in range(d):x[np.ravel_multi_index(tuple(p[i][k] for i in range(n)),[d]*n)]=math.sqrt(w[k])
        exact=m.evaluate(x,t,lam,n,d,check=True)
        err=abs(exact['lhs']-r['lhs']);assert err<2e-11
        assert abs(exact['rhs']-r['rhs'])<2e-11
        cross.append(dict(n=n,d=d,error=err))
rows=[];optimized=[]
for n,d in [(4,3),(5,3),(6,2),(4,4),(8,2),(5,4)]:
    for j in range(480):
        w=rng.dirichlet(np.ones(d)*[.03,.3,3][j%3])
        tau=rng.dirichlet(np.ones(d)*[.03,.3,3][(j//3)%3])
        if j%7==0:tau=np.eye(d)[j%d]
        perms=[rng.permutation(d) for _ in range(n)]
        lam=float(rng.choice([.001,.01,.05,.2,.5,.8,.95,.99,.999]))
        rows.append(fast(w,tau,lam,n,d,perms))
    print('GHZ block',n,d,'minimum',min(r['margin'] for r in rows if r['n']==n and r['d']==d),flush=True)
    for j in range(6):
        perms=[rng.permutation(d) for _ in range(n)]
        def unpack(z):
            a=np.exp(z[:d]-max(z[:d]));b=np.exp(z[d:2*d]-max(z[d:2*d]))
            return a/a.sum(),b/b.sum(),float(z[-1])
        def objective(z):return fast(*unpack(z),n,d,perms)['margin']
        z=rng.normal(size=2*d+1);z[-1]=rng.uniform(.01,.99)
        opt=minimize(objective,z,method='L-BFGS-B',bounds=[(-12,12)]*(2*d)+[(.001,.999)],
                     options={'maxiter':200,'ftol':1e-12,'gtol':2e-7,'maxls':30})
        r=fast(*unpack(opt.x),n,d,perms)
        r.update(success=bool(opt.success),message=str(opt.message),nfev=int(opt.nfev),nit=int(opt.nit))
        optimized.append(r)
out=dict(scope='Finite GHZ-family falsifier search up to8 uses; no proof, capacity or score claim',
 seed=2026102468,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 crosschecks=cross,structured=rows,optimized=optimized,elapsed_seconds=time.time()-start)
out['summary']=dict(crosschecks=len(cross),cases=len(rows),optimization_runs=len(optimized),
 minimum_margin=min(r['margin'] for r in rows+optimized),negative_cases=sum(r['margin'] < -2e-10 for r in rows+optimized),
 optimizer_nonconvergences=sum(not r['success'] for r in optimized))
(C/'ghz-falsifier68-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(out['summary']),flush=True)
