"""Exact rational interval certificate for one finite concave-convex game. Stdlib only."""
from pathlib import Path
from fractions import Fraction as F
import argparse, json, hashlib, re
root=Path(__file__).parent
S=10**60;N=80
def ensure(condition, message):
    if not condition:
        raise ValueError(message)
def integer_text(value):
    ensure(isinstance(value,str) and re.fullmatch(r"[1-9][0-9]*",value) is not None, "positive integer string required")
    return int(value)
def rational_text(value):
    ensure(isinstance(value,str) and re.fullmatch(r"-?[0-9]+(?:/[1-9][0-9]*)?",value) is not None, "rational string required")
    return F(value)
def ceildiv(a,b):return -((-a)//b)
def unit_log(q):
    ensure(1<=q<=2, 'invalid certificate or interval')
    u=(q-1)/(q+1)
    lo=u.numerator*S//u.denominator;hi=ceildiv(u.numerator*S,u.denominator)
    u2lo=lo*lo//S;u2hi=ceildiv(hi*hi,S)
    pl,ph=lo,hi;sl=sh=0
    for j in range(N):
        sl+=pl//(2*j+1);sh+=ceildiv(ph,2*j+1)
        pl=pl*u2lo//S;ph=ceildiv(ph*u2hi,S)
    # ln(q)=2*atanh(u), u<=1/3; omitted tail <=9/(4(2N+1)3^(2N+1)).
    tail=ceildiv(9*S,4*(2*N+1)*3**(2*N+1))
    return F(2*sl,S),F(2*sh+tail,S)
L2=unit_log(F(2))
def scale(a,c):return (a[0]*c,a[1]*c) if c>=0 else (a[1]*c,a[0]*c)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def log(q):
    ensure(q>0, 'invalid certificate or interval')
    k=q.numerator.bit_length()-q.denominator.bit_length()
    r=q/(F(2)**k)
    while r<1:k-=1;r*=2
    while r>=2:k+=1;r/=2
    return add(unit_log(r),scale(L2,k))
def divide_pos(a,b):
    ensure(b[0]>0, 'invalid certificate or interval')
    z=[x/y for x in a for y in b]
    return min(z),max(z)
def validate(cert):
    ensure(isinstance(cert,dict), 'certificate object required')
    Q=integer_text(cert['denominator'])
    ensure(type(cert['E']) is int and type(cert['dimension']) is int, 'integer metadata required')
    ensure(type(cert['p']) is list and type(cert['t']) is list, 'probability arrays required')
    p=list(map(integer_text,cert['p']));t=list(map(integer_text,cert['t']))
    ensure(cert['E']==8 and cert['c']=='1/20' and cert['lambda']=='3/5', 'invalid certificate or interval')
    ensure(len(p)==len(t)==cert['dimension']==321 and min(p+t)>0, 'invalid certificate or interval')
    ensure(sum(p)==sum(t)==Q, 'invalid certificate or interval')
    ensure(sum(i*p[i] for i in range(321))==8*Q, 'invalid certificate or interval')
    ensure(5*sum(i*t[i] for i in range(321))==2*Q, 'invalid certificate or interval')
    return [F(v,Q) for v in p],[F(v,Q) for v in t]
def certify(cert, *, enabled=False):
    ensure(enabled is True, "explicit enabled=True required")
    p,t=validate(cert);lam=F(3,5);q=1-lam;c=F(1,20)
    B=rational_text(cert['input_slope']);b=rational_text(cert['reset_slope'])
    ensure(B>=0 and b<=0, 'dual slope signs required')
    m=[lam*a+q*z for a,z in zip(p,t)]
    lp=list(map(log,p));lt=list(map(log,t));lm=list(map(log,m))
    numerator=add(scale(log(lam),lam),scale(log(q),q))
    gp=[];gt=[]
    for i in range(321):
        numerator=add(numerator,scale(lp[i],-lam*p[i]))
        numerator=add(numerator,scale(lm[i],-m[i]))
        numerator=add(numerator,scale(lt[i],q*t[i]))
        gp.append(divide_pos(scale(add(add(lp[i],lm[i]),(F(2),F(2))),-lam),L2))
        gt.append(divide_pos(scale(add(lt[i],scale(lm[i],F(-1))),q),L2))
    phi=divide_pos(numerator,L2)
    A=max(gp[i][1]-B*F(i,8) for i in range(321))
    a=min(gt[i][0]-b*F(i,8) for i in range(321))
    dotp=(sum(p[i]*gp[i][0] for i in range(321)),sum(p[i]*gp[i][1] for i in range(321)))
    dott=(sum(t[i]*gt[i][0] for i in range(321)),sum(t[i]*gt[i][1] for i in range(321)))
    upper=phi[1]+A+B-dotp[0]
    lower=phi[0]+a+b*c-dott[1]
    ensure(lower<=upper, 'invalid certificate or interval')
    # Public compact interval, rounded outwards to 1e-10 bits.
    D=10**10;L=lower.numerator*D//lower.denominator;U=ceildiv(upper.numerator*D,upper.denominator)
    ensure(F(L,D)<=lower<=upper<=F(U,D), 'outward rounding failed')
    return {'lower':str(lower),'upper':str(upper),'compact_lower':str(F(L,D)),'compact_upper':str(F(U,D)),'compact_lower_decimal':f'{L//D}.{L%D:010d}','compact_upper_decimal':f'{U//D}.{U%D:010d}','width':str(upper-lower),'dual_input_gap_upper':str(A+B-dotp[0]),'dual_reset_gap_upper':str(dott[1]-a-b*c),'phi_interval':list(map(str,phi)),'all_642_dual_inequalities_certified':True,'log_fixed_scale':str(S),'atanh_terms':N}

CASE = Path(__file__).with_name("case_e8.json")
EXPECTED_SHA256 = "fde590462bae8b16d085ebe966156acc223b6fc40c8a0a0eaed45e9c28c9bff5"

def load_case():
    raw=CASE.read_bytes()
    ensure(hashlib.sha256(raw).hexdigest()==EXPECTED_SHA256, "stored certificate hash mismatch")
    return json.loads(raw)

def main():
    parser=argparse.ArgumentParser(description="Opt-in exact enclosure of one finite entropy game; no optimizer.")
    parser.add_argument("--enable", action="store_true", help="verify the fixed E8 rational certificate")
    args=parser.parse_args()
    if not args.enable:
        parser.error("--enable is required")
    result=certify(load_case(), enabled=True)
    ensure(F(result["lower"])>=F("4.9836437909") and F(result["upper"])<=F("4.9836437910"), "published enclosure not reproduced")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
