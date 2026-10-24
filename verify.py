"""Replay finite saddle diagnostics without changing preserved results."""
import hashlib,json,sys,platform
from pathlib import Path
import numpy,scipy
from finite_game50 import solve

root=Path(__file__).resolve().parent
manifest=json.loads((root/'manifest.json').read_text())
for item in manifest['files']:
    assert hashlib.sha256((root/item['path']).read_bytes()).hexdigest()==item['sha256'],item['path']
saved=json.loads((root/'finite-game-results50.json').read_text())
rows=[];guess=None;previous=None
for expected in saved['rows']:
    if expected['c']!=previous:guess=None
    row,guess=solve(expected['E'],expected['c'],guess,xmax=expected['xmax'])
    previous=expected['c'];rows.append(row)
    assert row['passed_numerical_checks'],row
    for key in ('finite_entropy_game_value_bits','input_scaled_mean','reset_scaled_mean','reset_mass_at_zero'):
        assert abs(row[key]-expected[key])<2e-7,(expected['E'],expected['c'],key,row[key],expected[key])
    print(f"E={row['E']} c={row['c']}: verified",flush=True)
out=root/'runs';out.mkdir(exist_ok=True)
packet={'status':'finite_floating_point_replay_only','python':sys.version,'platform':platform.platform(),'numpy':numpy.__version__,'scipy':scipy.__version__,'rows':rows,'all_pass':True}
(out/'verification.json').write_text(json.dumps(packet,indent=2)+'\n',encoding='utf-8')
print('All 28 finite diagnostics passed; no infinite-limit or quantum-capacity certification.')
