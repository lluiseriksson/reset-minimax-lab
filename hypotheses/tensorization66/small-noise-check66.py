from pathlib import Path
import json,hashlib,importlib.util,math
import numpy as np

ROOT=Path(__file__).parent
spec=importlib.util.spec_from_file_location('falsifier66',ROOT/'entropy-falsifier66.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
rng=np.random.default_rng(2026102467)
rows=[]
for n,d in [(2,2),(2,3),(3,2)]:
    for sample in range(4):
        v=rng.normal(size=d**n)+1j*rng.normal(size=d**n);v/=np.linalg.norm(v)
        tau_spec=rng.dirichlet(np.ones(d));tau=np.diag(tau_spec);psi=np.outer(v,v.conj())
        marg=[m.partial(psi,[i],n,d) for i in range(n)]
        delta=sum(float(np.trace(tau@(p-p@p)).real) for p in marg)
        a=n-sum(float(np.trace(tau@p@p).real) for p in marg)
        # Independent global expectation verifies the Schmidt-overlap identity.
        expectations=[]
        for i in range(n):
            rest=[j for j in range(n) if j!=i];order=[i]+rest;p=[order.index(j) for j in range(n)]
            block=np.kron(tau,m.partial(psi,rest,n,d)).reshape([d]*(2*n))
            out=block.transpose(p+[n+j for j in p]).reshape(d**n,d**n)
            ex=float(np.vdot(v,out@v).real);expectations.append(ex)
            assert abs(ex-np.trace(tau@marg[i]@marg[i]).real)<1e-12
        scales=[]
        for q in [1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]:
            r=m.evaluate(v,tau_spec,1-q,n,d,check=True)
            scales.append(dict(q=q,gap=r['margin'],ratio=r['margin']/(q*math.log(1/q)),
                               order_q_remainder=(r['margin']-delta*q*math.log(1/q))/q))
        rows.append(dict(n=n,d=d,sample=sample,delta=delta,entropy_coefficient=a,
                         replacement_expectations=expectations,scales=scales))
orthogonal=[]
for d,tau_spec,offset in [(3,[1,0,0],1),(4,[.3,.7,0,0],2)]:
    v=np.zeros(d*d,complex);v[offset*d+offset]=1/math.sqrt(2);v[(offset+1)*d+offset+1]=1/math.sqrt(2)
    for q in [.001,.01,.2,.5,.8,.99,.999]:
        r=m.evaluate(v,np.array(tau_spec),1-q,2,d,check=True)
        expected=2*q*(1-q)*math.log(2)
        assert abs(r['margin']-expected)<2e-12
        orthogonal.append(dict(d=d,q=q,gap=r['margin'],expected=expected,error=r['margin']-expected))
out=dict(scope='Floating-point checks of a candidate pointwise asymptotic lemma, not proof',
         source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=rows,orthogonal_support_controls=orthogonal)
(ROOT/'small-noise-check66-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(dict(asymptotic_states=len(rows),scales_per_state=6,orthogonal_support_controls=len(orthogonal),
                      maximum_orthogonal_control_error=max(abs(r['error']) for r in orthogonal))))
