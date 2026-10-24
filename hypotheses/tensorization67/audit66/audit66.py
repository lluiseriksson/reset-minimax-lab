"""Targeted audit checks, authored without reading the candidate's test scripts.

These finite-precision checks support, but do not replace, the accompanying proof.
They do not certify the all-q conjecture. Every sampled tuple and spectrum is saved.
"""
from __future__ import annotations
import hashlib
import itertools
import json
import platform
import sys
import traceback
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parent
RNG = np.random.default_rng(660024)
TOL = 2e-10
FAILURES = []


def serial(z):
    z = np.asarray(z)
    return {"real": z.real.tolist(), "imag": z.imag.tolist()}


def pure(v):
    return np.outer(v, v.conj())


def unit(v):
    return v / np.linalg.norm(v)


def random_vector(m):
    return unit(RNG.normal(size=m) + 1j * RNG.normal(size=m))


def random_density(d, rank):
    z = RNG.normal(size=(d, rank)) + 1j * RNG.normal(size=(d, rank))
    z = z @ z.conj().T
    return z / np.trace(z)


def h(t):
    if t < -TOL or t > 1 + TOL:
        raise ValueError(f"binary entropy outside domain: {t}")
    t = min(1.0, max(0.0, float(t)))
    if t == 0 or t == 1:
        return 0.0
    return -t * np.log(t) - (1 - t) * np.log1p(-t)


def entropy(rho):
    ev = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    if ev.min() < -TOL:
        raise ValueError(f"negative eigenvalue {ev.min()}")
    pos = ev[ev > 0]
    return float(-np.dot(pos, np.log(pos)))


def marginal(rho, keep, n, d):
    """Partial trace by tensor-axis contraction; kept indices are ascending."""
    z = rho.reshape([d] * (2 * n))
    remaining = n
    for i in reversed(range(n)):
        if i not in keep:
            z = np.trace(z, axis1=i, axis2=i + remaining)
            remaining -= 1
    m = d ** len(keep)
    return z.reshape(m, m)


def branch(rho, keep, tau, n, d):
    """Construct a retained/replaced branch via Kronecker factors and permutation."""
    keep = tuple(keep)
    missing = tuple(i for i in range(n) if i not in keep)
    z = marginal(rho, keep, n, d)
    for _ in missing:
        z = np.kron(z, tau)
    order = keep + missing
    perm = [order.index(i) for i in range(n)]
    z = z.reshape([d] * (2 * n))
    return z.transpose(perm + [n + i for i in perm]).reshape(d**n, d**n)


def all_branches(rho, tau, n, d):
    return [(keep, branch(rho, keep, tau, n, d))
            for k in range(n + 1) for keep in itertools.combinations(range(n), k)]


def subset_channel(branches, q, n):
    return sum((1-q)**len(t) * q**(n-len(t)) * b for t, b in branches)


def sequential_channel(rho, tau, q, n, d):
    """Independent channel construction via one-site superoperator tensor contraction."""
    z = rho.reshape([d] * (2 * n))
    # C[a,b,c,e] sends input matrix index (c,e) to output index (a,b).
    identity = np.eye(d)
    c = ((1-q) * np.einsum('ac,be->abce', identity, identity)
         + q * np.einsum('ab,ce->abce', tau, identity))
    for i in range(n):
        contracted = np.tensordot(c, z, axes=([2, 3], [i, n+i]))
        remaining = [j for j in range(2*n) if j not in (i, n+i)]
        new_order = [i, n+i] + remaining
        z = contracted.transpose([new_order.index(j) for j in range(2*n)])
    return z.reshape(d**n, d**n)


def case_record(name, v, tau, n, d, zero_expected=False, e_sites=None):
    rho = pure(v)
    branches = all_branches(rho, tau, n, d)
    one = [marginal(rho, (i,), n, d) for i in range(n)]
    t = np.array([np.trace(tau @ x).real for x in one])
    overlaps = np.array([np.trace(tau @ x @ x).real for x in one])
    delta = float(np.sum(t-overlaps))
    a = float(n - sum(overlaps))
    A = sum(branch(rho, tuple(j for j in range(n) if j != i), tau, n, d)
            for i in range(n)) - n*rho
    P = np.eye(d**n) - rho
    mu = np.linalg.eigvalsh(P @ A @ P)
    mup = mu[mu > 1e-12]
    tsafe = t[t > 0]
    c_linear = float(delta - np.dot(mup, np.log(mup)) - n*entropy(tau)
                     - np.dot(tsafe, np.log(tsafe)))
    identities = {
        "a_trace_error": float(abs(np.trace(P @ A @ P).real-a)),
        "minimum_compression_eigenvalue": float(mu.min()),
        "schmidt_identity_max_error": float(max(
            abs(np.vdot(v, branch(rho, tuple(j for j in range(n) if j != i),
                                   tau, n, d) @ v).real-overlaps[i])
            for i in range(n))),
    }
    if zero_expected and abs(delta) > TOL:
        FAILURES.append({"name": name, "type": "delta_zero", "delta": delta})
    rows = []
    for q in (0., 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, .01, .1, .25, .5, .75, .9, .99, 1.):
        out = subset_channel(branches, q, n)
        alternative = sequential_channel(rho, tau, q, n, d)
        channel_error = float(np.max(np.abs(out-alternative)))
        g = float(entropy(out) - n*(q*entropy(tau)+h(q)) + sum(h(q*x) for x in t))
        row = {"q": q, "G": g, "channel_crosscheck_error": channel_error,
               "output_spectrum": np.linalg.eigvalsh(out).tolist()}
        if 0 < q < 1:
            row["G_div_qlog"] = g / (q*np.log(1/q))
            row["linear_remainder"] = (g-delta*q*np.log(1/q))/q
            row["linear_coefficient_prediction"] = c_linear
        if channel_error > TOL:
            FAILURES.append({"name": name, "q": q, "type": "channel_crosscheck",
                             "error": channel_error})
        if zero_expected and g < -TOL:
            FAILURES.append({"name": name, "q": q, "type": "claim_b_margin", "G": g})
        rows.append(row)
    if e_sites is not None and len(e_sites) == n:
        max_orthogonality = max(abs(np.trace(b @ c))
                                for j, (_, b) in enumerate(branches)
                                for _, c in branches[:j])
        identities["maximum_distinct_branch_trace_overlap"] = float(max_orthogonality)
        errors = []
        for row in rows:
            q = row["q"]
            et = sum((1-q)**len(k)*q**(n-len(k))*entropy(marginal(rho,k,n,d))
                     for k, _ in branches)
            errors.append(abs(row["G"]-et))
        identities["orthogonal_entropy_identity_max_error"] = float(max(errors))
        if max_orthogonality > TOL or max(errors) > TOL:
            FAILURES.append({"name": name, "type": "orthogonal_branch_identity",
                             "overlap": float(max_orthogonality), "error": float(max(errors))})
    return {"name": name, "n": n, "d": d, "v": serial(v), "tau": serial(tau),
            "delta_zero_by_construction": zero_expected, "delta": delta,
            "t": t.tolist(), "a": a, "PAP_spectrum": mu.tolist(),
            "identities": identities, "rows": rows}


def main():
    cases = []
    cases.append(case_record("dimension_one", np.ones(1), np.ones((1,1)), 1, 1, True))
    cases.append(case_record("unchanged_rank_one_input", np.array([1.,0.,0.,0.]), pure(np.array([1.,0.])), 2, 2, True))
    cases.append(case_record("orthogonal_product_with_quadratic_branch", np.array([1.,0.,0.,0.]), pure(np.array([0.,1.])), 2, 2, True))
    tau = pure(np.array([np.sqrt(.9), 1j*np.sqrt(.1)]))
    cases.append(case_record("pure_input_error_above_half", np.array([1.,0.]), tau, 1, 2, True))
    tau = np.array([[.6,.13+.09j],[.13-.09j,.4]])
    cases.append(case_record("product_noncommuting_reset", np.kron(random_vector(2), random_vector(2)), tau, 2, 2, True))
    cases.append(case_record("bell_pure_reset", unit(np.array([1.,0.,0.,1.])), pure(np.array([1.,0.])), 2, 2))
    ghz = np.zeros(8, complex); ghz[[0,7]] = [np.sqrt(.31), 1j*np.sqrt(.69)]
    cases.append(case_record("ghz_mixed_reset", ghz, tau, 3, 2))
    w = np.zeros(8, complex); w[[1,2,4]] = unit(np.array([1.,1j,.7]))
    cases.append(case_record("w_rank_one_reset", w, pure(random_vector(2)), 3, 2))
    tau4 = np.zeros((4,4), complex); tau4[2:,2:] = np.array([[.63,.16j],[-.16j,.37]])
    pair = np.zeros(16, complex); pair[[0,5]] = [np.sqrt(.3), 1j*np.sqrt(.7)]
    cases.append(case_record("orthogonal_entangled_pair", pair, tau4, 2, 4, True, [0,1]))
    phi = unit(np.array([1.,1j,2.,-.7j]))
    cases.append(case_record("orthogonal_pair_plus_pure_site", np.kron(pair,phi), tau4, 3, 4, True, [0,1]))
    ghz4 = np.zeros(64, complex); ghz4[[0,21]] = [np.sqrt(.4), np.sqrt(.6)]
    cases.append(case_record("orthogonal_entangled_ghz", ghz4, tau4, 3, 4, True, [0,1,2]))
    z = RNG.normal(size=(4,4)) + 1j*RNG.normal(size=(4,4))
    U, _ = np.linalg.qr(z)
    cases.append(case_record("rotated_orthogonal_pair", np.kron(U,U)@pair, U@tau4@U.conj().T, 2, 4, True, [0,1]))
    for eps in (1e-2,1e-6):
        v = unit(np.array([np.sqrt(1-eps),0.,0.,np.sqrt(eps)]))
        cases.append(case_record(f"near_product_schmidt_{eps}", v, pure(np.array([1.,0.])), 2, 2))
    for k in range(36):
        n,d = ((2,2),(2,3),(3,2))[k%3]
        rank = 1 + ((k//3) % d)
        cases.append(case_record(f"random_{k:02d}", random_vector(d**n), random_density(d,rank), n,d))
    summary = {
        "case_count": len(cases), "channel_evaluations": sum(len(x["rows"]) for x in cases),
        "delta_zero_cases": sum(x["delta_zero_by_construction"] for x in cases),
        "failure_count": len(FAILURES),
        "max_channel_error": max(r["channel_crosscheck_error"] for x in cases for r in x["rows"]),
        "max_schmidt_error": max(x["identities"]["schmidt_identity_max_error"] for x in cases),
        "max_a_trace_error": max(x["identities"]["a_trace_error"] for x in cases),
        "minimum_delta_zero_margin": min(r["G"] for x in cases if x["delta_zero_by_construction"] for r in x["rows"]),
        "max_linear_remainder_error_at_q_1e_minus_7": max(abs(x["rows"][1]["linear_remainder"]-x["rows"][1]["linear_coefficient_prediction"]) for x in cases),
        "tolerance": TOL,
    }
    result = {"scope": "Targeted finite-precision checks; no proof or all-q conjecture certificate",
              "python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
              "seed": 660024, "summary": summary, "failures": FAILURES, "cases": cases}
    (BASE/"numeric-evidence.json").write_text(json.dumps(result,indent=2,allow_nan=False),encoding="utf-8")
    print(json.dumps(summary,indent=2))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        with (BASE/"execution-failures.txt").open("a",encoding="utf-8") as f:
            f.write(traceback.format_exc()+"\n")
        raise
