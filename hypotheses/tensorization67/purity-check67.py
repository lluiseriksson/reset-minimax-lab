from pathlib import Path
import importlib.util,json,hashlib,math
import numpy as np
C=Path(__file__).parent
spec=importlib.util.spec_from_file_location('falsifier66',C/'entropy-falsifier66.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
data=json.loads((C/'entropy-falsifier66-results.json').read_text())
rows=[]
for group in ['structured','optimized']:
    for j,r in enumerate(data[group]):
        n,d,q=r['n'],r['d'],1-r['lam']
        if n*q>0.5:continue
        x=np.array(r['psi_real'])+1j*np.array(r['psi_imag']);rho=np.outer(x,x.conj());tau=np.diag(r['spectrum'])
        a=sum(1-float(np.trace(tau@m.partial(rho,[i],n,d)@m.partial(rho,[i],n,d)).real) for i in range(n))
        lower=m.h(a*q*(1-q)**(2*n-1));margin=r['lhs']-lower
        assert margin>=-2e-10,(group,j,margin)
        rows.append(dict(group=group,index=j,n=n,d=d,q=q,a=a,output_entropy=r['lhs'],purity_entropy_lower=lower,margin=margin))
out=dict(scope='Finite diagnostic of an unaudited candidate lower bound; not proof',
 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 candidate_sha256=hashlib.sha256((C/'purity-bound67-candidate.txt').read_bytes()).hexdigest(),
 checked_cases=len(rows),minimum_margin=min(r['margin'] for r in rows),rows=rows)
(C/'purity-check67-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:out[k] for k in ['scope','checked_cases','minimum_margin']}))
