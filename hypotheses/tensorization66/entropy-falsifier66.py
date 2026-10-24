"""Structured and optimized falsifier search. Numerical evidence, never a proof."""
from pathlib import Path
import itertools, json, math, time, hashlib
import numpy as np
from scipy.optimize import minimize

ROOT=Path(__file__).parent
RNG=np.random.default_rng(2026102466)

def ent(a):
    w=np.linalg.eigvalsh((a+a.conj().T)/2)
    if w.min() < -2e-10 or abs(w.sum()-1)>2e-9: raise ValueError('invalid density')
    w=w[w>1e-15]
    return float(-w@np.log(w))

def h(x):
    x=float(np.clip(x,0,1))
    return 0. if x==0 or x==1 else -x*math.log(x)-(1-x)*math.log1p(-x)

def partial(rho, keep, n, d):
    # Direct index contraction, independent from the earlier sequential trace code.
    labels=list(range(n)); right=[i if i not in keep else n+i for i in labels]
    t=np.einsum(rho.reshape([d]*(2*n)),labels+right,
                list(keep)+[n+i for i in keep],optimize=True)
    return t.reshape(d**len(keep),d**len(keep))

def branch_output(rho,tau,lam,n,d):
    out=np.zeros_like(rho,dtype=complex)
    for mask in itertools.product((0,1),repeat=n):
        keep=[i for i in range(n) if mask[i]]
        gone=[i for i in range(n) if not mask[i]]
        block=partial(rho,keep,n,d)
        order=keep+gone
        for _ in gone: block=np.kron(block,tau)
        p=[order.index(i) for i in range(n)]
        block=block.reshape([d]*(2*n)).transpose(p+[n+i for i in p]).reshape(d**n,d**n)
        out+=lam**len(keep)*(1-lam)**len(gone)*block
    return out

def sequential(rho,tau,lam,n,d):
    out=rho.copy()
    for i in range(n):
        rest=[j for j in range(n) if j!=i]
        block=np.kron(tau,partial(out,rest,n,d))
        order=[i]+rest; p=[order.index(j) for j in range(n)]
        rep=block.reshape([d]*(2*n)).transpose(p+[n+j for j in p]).reshape(d**n,d**n)
        out=lam*out+(1-lam)*rep
    return out

def evaluate(x,spec,lam,n,d,check=False):
    x=x/np.linalg.norm(x); tau=np.diag(spec); rho=np.outer(x,x.conj())
    out=branch_output(rho,tau,lam,n,d)
    if check:
        assert np.max(np.abs(out-sequential(rho,tau,lam,n,d)))<2e-12
    overlaps=[float(np.trace(partial(rho,[i],n,d)@tau).real) for i in range(n)]
    lhs=ent(out); rhs=n*((1-lam)*ent(tau)+h(lam))-sum(h((1-lam)*t) for t in overlaps)
    return dict(n=n,d=d,lam=lam,spectrum=list(map(float,spec)),overlaps=overlaps,
                lhs=lhs,rhs=rhs,margin=lhs-rhs,
                psi_real=x.real.tolist(),psi_imag=x.imag.tolist())

def run():
    start=time.time(); structured=[]; optimized=[]; crosschecks=[]
    # Nonuniform Schmidt/GHZ weights, W and graph states; tau diagonal entails no
    # loss in the unrestricted optimization because the same unitary rotates inputs.
    for n,d in [(2,2),(2,3),(3,2),(3,3),(4,2)]:
        families=[]
        for a in [0.,1e-6,.01,.1,.5,.9,1.]:
            x=np.zeros(d**n,complex); x[0]=math.sqrt(a); x[-1]=math.sqrt(1-a)
            families.append(('ghz-'+str(a),x))
        w=np.zeros(d**n,complex)
        for i in range(n): w[d**i]=1/math.sqrt(n)
        families.append(('w',w))
        if d==2:
            x=np.ones(2**n,complex)/math.sqrt(2**n)
            for k in range(2**n):
                bits=[(k>>i)&1 for i in range(n)]
                x[k]*=(-1)**sum(bits[i]*bits[i+1] for i in range(n-1))
            families.append(('chain-graph',x))
        spectra=[np.ones(d)/d,np.eye(d)[0],np.eye(d)[-1]]
        for p in [.001,.01,.1,.5,.9,.99]:
            s=np.full(d,(1-p)/(d-1));s[0]=p;spectra.append(s)
        if d==3:spectra.append(np.array([.6,.3,.1]))
        for name,x in families:
            for spec in spectra:
                for lam in [.001,.01,.05,.2,.5,.8,.95,.99,.999]:
                    r=evaluate(x,spec,lam,n,d,check=True); r['family']=name
                    structured.append(r)
        print('structured',n,d,len(structured),flush=True)
    for n,d,starts in [(2,2,24),(2,3,12),(3,2,12)]:
        size=d**n
        def unpack(z):
            x=z[:size]+1j*z[size:2*size]
            logits=z[2*size:2*size+d];s=np.exp(logits-logits.max());s/=s.sum()
            return x,s,float(z[-1])
        def objective(z):
            x,s,lam=unpack(z)
            if np.linalg.norm(x)<1e-10:return 100.
            return evaluate(x,s,lam,n,d)['margin']
        for j in range(starts):
            z=RNG.normal(size=2*size+d+1);z[-1]=RNG.uniform(.02,.98)
            if j%3==0:z[2*size:2*size+d]*=5
            opt=minimize(objective,z,method='L-BFGS-B',bounds=[(-4,4)]*(2*size)+[(-12,12)]*d+[(.001,.999)],
                         options={'maxiter':180,'ftol':1e-12,'gtol':2e-7,'maxls':30})
            x,s,lam=unpack(opt.x);r=evaluate(x,s,lam,n,d,check=True)
            r.update(start=j,success=bool(opt.success),message=str(opt.message),nfev=int(opt.nfev),nit=int(opt.nit))
            optimized.append(r)
        print('optimized',n,d,'min',min(r['margin'] for r in optimized if r['n']==n and r['d']==d),flush=True)
    # Random mixed inputs validate independent branch and sequential constructions,
    # beyond the pure inputs in the conjecture.
    for n,d in [(2,2),(2,3),(3,2)]:
        for j in range(10):
            z=RNG.normal(size=(d**n,d**n))+1j*RNG.normal(size=(d**n,d**n))
            rho=z@z.conj().T;rho/=np.trace(rho);spec=RNG.dirichlet(np.ones(d));lam=RNG.uniform()
            err=float(np.max(np.abs(branch_output(rho,np.diag(spec),lam,n,d)-sequential(rho,np.diag(spec),lam,n,d))))
            assert err<2e-12;crosschecks.append(dict(n=n,d=d,error=err))
    result=dict(scope='Finite structured and local optimization falsifier search, no theorem or capacity claim',seed=2026102466,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),structured=structured,optimized=optimized,
                independent_channel_crosschecks=crosschecks,elapsed_seconds=time.time()-start)
    result['summary']=dict(structured_cases=len(structured),optimization_runs=len(optimized),
        structured_minimum=min(r['margin'] for r in structured),optimized_minimum=min(r['margin'] for r in optimized),
        numerical_negative_cases=sum(r['margin'] < -2e-10 for r in structured+optimized),
        optimizer_nonconvergences=sum(not r['success'] for r in optimized))
    (ROOT/'entropy-falsifier66-results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result['summary']),flush=True)

if __name__=='__main__':run()
