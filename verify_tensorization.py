"""Replay an unproved entropy conjecture's falsifier search in a fresh directory."""
from pathlib import Path
import hashlib,json,platform,shutil,subprocess,sys,tempfile
import numpy,scipy

root=Path(__file__).resolve().parent
source=root/'hypotheses'/'tensorization66'
manifest=json.loads((source/'manifest.json').read_text(encoding='utf-8'))
for item in manifest['files']:
    raw=(source/item['path']).read_bytes()
    assert len(raw)==item['bytes'] and hashlib.sha256(raw).hexdigest()==item['sha256'],item['path']
(root/'runs').mkdir(exist_ok=True)
dest=Path(tempfile.mkdtemp(prefix='tensorization66-',dir=root/'runs'))
for name in ['entropy-falsifier66.py','small-noise-check66.py']:
    shutil.copyfile(source/name,dest/name)
for script in ['entropy-falsifier66.py','small-noise-check66.py']:
    with (dest/(script+'.log')).open('w',encoding='utf-8') as log:
        result=subprocess.run([sys.executable,'-B',str(dest/script)],cwd=dest,stdout=log,stderr=subprocess.STDOUT)
    assert result.returncode==0,(script,result.returncode,str(dest))
data=json.loads((dest/'entropy-falsifier66-results.json').read_text(encoding='utf-8'))
s=data['summary']
assert s['structured_cases']==3627 and s['optimization_runs']==48
# Convergence and finite minima are diagnostics; do not compare optimizer minima
# bit-for-bit or mistake their positive value for a global lower bound.
assert s['numerical_negative_cases']==0,('Counterexample candidate requires inspection',str(dest))
assert len(data['independent_channel_crosschecks'])==30
checks=json.loads((dest/'small-noise-check66-results.json').read_text(encoding='utf-8'))
assert len(checks['rows'])==12 and len(checks['orthogonal_support_controls'])==14
packet=dict(scope='Finite floating-point falsifier replay, not theorem verification',
            python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,
            run=str(dest),summary=s,independent_mixed_input_checks=30,asymptotic_states=12,orthogonal_controls=14)
(dest/'replay-summary.json').write_text(json.dumps(packet,indent=2)+'\n',encoding='utf-8')
print(json.dumps(packet,indent=2))
