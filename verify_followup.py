"""Replay finite diagnostics adjacent to separately audited partial entropy lemmas."""
from pathlib import Path
import hashlib,json,platform,shutil,subprocess,sys,tempfile
root=Path(__file__).resolve().parent;source=root/'hypotheses/tensorization67'
manifest=json.loads((source/'manifest.json').read_text(encoding='utf-8'))
for item in manifest['files']:
    raw=(source/item['path']).read_bytes()
    assert len(raw)==item['bytes'] and hashlib.sha256(raw).hexdigest()==item['sha256'],item['path']
(root/'runs').mkdir(exist_ok=True)
dest=Path(tempfile.mkdtemp(prefix='followup67-',dir=root/'runs'))
for script in ['entropy-falsifier66.py','purity-bound67-candidate.txt','purity-check67.py','ghz-falsifier68.py']:
    shutil.copyfile(source/script,dest/script)
# The bound diagnostic reuses the preserved, separately manifested 66 rows.
old_manifest=json.loads((root/'hypotheses/tensorization66/manifest.json').read_text(encoding='utf-8'))
old_item=next(i for i in old_manifest['files'] if i['path']=='entropy-falsifier66-results.json')
old_raw=(root/'hypotheses/tensorization66/entropy-falsifier66-results.json').read_bytes()
assert len(old_raw)==old_item['bytes'] and hashlib.sha256(old_raw).hexdigest()==old_item['sha256']
shutil.copyfile(root/'hypotheses/tensorization66/entropy-falsifier66-results.json',dest/'entropy-falsifier66-results.json')
for script in ['purity-check67.py','ghz-falsifier68.py']:
    with (dest/(script+'.log')).open('w',encoding='utf-8') as log:
        result=subprocess.run([sys.executable,'-B',str(dest/script)],cwd=dest,stdout=log,stderr=subprocess.STDOUT)
    assert result.returncode==0,(script,result.returncode,str(dest))
purity=json.loads((dest/'purity-check67-results.json').read_text(encoding='utf-8'))
ghz=json.loads((dest/'ghz-falsifier68-results.json').read_text(encoding='utf-8'))
assert purity['checked_cases']==1418 and purity['minimum_margin']>=-2e-10
assert ghz['summary']['cases']==2880 and ghz['summary']['optimization_runs']==36
assert ghz['summary']['negative_cases']==0
assert len(ghz['crosschecks'])==32
packet=dict(scope='Finite numerical follow-up, not proof/novelty/importance certification',
            platform=platform.platform(),native_files_verified=len(manifest['files']),
            purity_cases=purity['checked_cases'],GHZ_summary=ghz['summary'],run=str(dest))
(dest/'summary.json').write_text(json.dumps(packet,indent=2)+'\n',encoding='utf-8')
print(json.dumps(packet,indent=2))
