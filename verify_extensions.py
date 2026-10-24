"""Replay three frozen finite diagnostic groups; this is not a theorem checker."""
from pathlib import Path
import hashlib,json,math,os,platform,shutil,subprocess,sys,tempfile

ROOT=Path(__file__).resolve().parent
EXT=ROOT/'extensions'

def load(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

def finite(value):
    if isinstance(value,float):
        assert math.isfinite(value), 'Nonfinite numerical output'
    elif isinstance(value,dict):
        for item in value.values(): finite(item)
    elif isinstance(value,list):
        for item in value: finite(item)

def main():
    for entry in load(EXT/'manifest.json')['files']:
        path=EXT/entry['path']
        assert path.stat().st_size==entry['bytes']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256'],str(path)
    runs=ROOT/'runs';runs.mkdir(exist_ok=True)
    replay=Path(tempfile.mkdtemp(prefix='extensions-',dir=runs))
    summary={'scope':'Finite floating-point consistency only; no infinite-limit, novelty, formal or interval certification','platform':platform.platform(),'python':sys.version,'manifest_sha256':hashlib.sha256((EXT/'manifest.json').read_bytes()).hexdigest(),'groups':{}}
    configurations=[('holevo56','checks56.py','checks56-results.json'),('regular57','check57.py','results/checks57.json'),('quartic58','checks58.py','checks58-results.json')]
    for group,script,output in configurations:
        work=replay/group;shutil.copytree(EXT/group,work)
        env=dict(os.environ,PYTHONUTF8='1',PYTHONWARNINGS='error')
        with (work/'replay.log').open('w',encoding='utf-8') as log:
            subprocess.run([sys.executable,'-X','utf8',str(work/script)],cwd=work,env=env,check=True,stdout=log,stderr=subprocess.STDOUT)
        actual=load(work/output);finite(actual);s=actual['summary']
        if group=='holevo56':
            assert s['quantum_case_count']==48 and s['atomic_case_count']==60
            assert s['max_quantum_identity_error']<2e-10 and s['max_atomic_identity_error']<2e-10
            assert max(s['max_fano_violation'],s['max_jensen_violation'],s['max_flag_violation'],s['max_saddle_violation'])<2e-10
            invalid=[r for r in actual['atomic_finite_sums'] if r['ground_KKT_margin']<0]
            expected=load(EXT/group/'failures56.json')['finite_ground_condition_failure_count']
            assert len(invalid)==expected==14, 'Do not conceal failed ground certificates'
            summary['groups'][group]={'quantum_cases':48,'atomic_cases':60,'invalid_finite_ground_certificates_preserved':len(invalid),'max_identity_residual':max(s['max_quantum_identity_error'],s['max_atomic_identity_error'])}
        elif group=='regular57':
            assert s['cases']==24 and s['warnings']==0
            assert s['maximum_exact_identity_residual']<2e-9
            assert max(row['max_exact_identity_residual'] for row in actual['rows'])<2e-9
            summary['groups'][group]={'cases':24,'max_identity_residual':s['maximum_exact_identity_residual'],'initial_underflow_failure_preserved':True,'uniform_rate_claim':False}
        else:
            assert s['profiles_count']==54 and s['contrasts_count']==18 and s['all_positive_domain']
            assert float(s['max_moment_residual'])<1e-25 and float(s['max_C1_closed_formula_error'])<1e-40
            low=float(s['min_quartic_ratio_at_epsilon_0025']);high=float(s['max_quartic_ratio_at_epsilon_0025'])
            assert .97<low<=high<1.03
            summary['groups'][group]={'profiles':54,'contrasts':18,'quartic_ratio_range_at_epsilon_0025':[low,high],'max_moment_residual':s['max_moment_residual'],'interval_certification':False}
    summary['passed_stated_checks']=True
    data=json.dumps(summary,indent=2,allow_nan=False)+'\n'
    (replay/'verification.json').write_text(data,encoding='utf-8',newline='\n')
    (runs/'extensions-verification.json').write_text(data,encoding='utf-8',newline='\n')
    print(data)

if __name__=='__main__':main()
