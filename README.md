# Reset minimax laboratory

Original numerical diagnostics for a mean-constrained, concave–convex entropy game motivated by identity/replacement quantum channels. The default replay solves 28 finite games and checks normalization, energy constraints, stationary identities and independent linear-program best responses in both variables. It compares numerical values against the preserved execution.

```powershell
python -m pip install -r requirements.txt
python verify.py
```

Python 3.12; NumPy 2.2.6; SciPy 1.18.0. Linux and Windows are CI targets. Inspect actual workflow results before asserting cross-platform verification.

## Scope

The calculation maximizes in input probabilities and minimizes in reset probabilities, with separate mean constraints, lambda=0.6, seven reset/input budget ratios, input budgets 8,32,128,512 and level cutoff 40 times the input budget. Values are in bits. The distributions have full support on the finite grid. Log probabilities are retained to handle underflow in remote tails.

This repository does **not** certify an infinite-dimensional quantum capacity, an optimal limiting profile, novelty, a formal proof, an interval enclosure, or an AIRR importance grade. The analytical reduction and continuum theorem require separate proofs. A linear-program tangent check supplies finite numerical evidence; it is not an interval certificate.

## Preserved failure

`finite_game50-initial.py` and `finite-game-results50-initial.json` retain the first execution. At input budget 8 and ratio 0.05, the optimizer reported success despite an unacceptable energy residual. `finite_game50.py` uses a fixed bounded multistart policy and rejects such residuals. The corrected output retains attempted starts, including failures. Do not discard the initial files or infer validity from the solver flag alone.

`finite-game-results50.json` is the original corrected execution. `verify.py` writes a new result under ignored `runs/`, leaving all historical files unchanged. It exits unsuccessfully on constraint or best-response failures or disagreement with the saved values. `manifest.json` binds the immutable scientific files.

## Spectral and single-letter extensions

The native programs and original outputs under `extensions/` are separate finite checks for a combined manuscript in preparation. They include single-letter Holevo games, multiplicity-weighted spectral counts with slowly varying factors, and high-precision comparisons near equal budgets. They do not determine regularized classical capacity.

```powershell
python -m pip install -r requirements-extensions.txt
python verify_extensions.py
```

The replay copies each group into a new ignored `runs/` directory before executing, so historical inputs and outputs stay intact. `extensions/manifest.json` binds the native bytes. The latest summary is `runs/extensions-verification.json`; complete outputs and logs remain in the individual replay directory.

The Holevo group retains **14 invalid finite ground-state KKT certificates among60 grids**. The verification expects and reports those failures; it does not claim all60 are saddles. Their eventual validity is an analytical claim outside this program. The regular-variation group preserves its first underflow/NaN diagnostic and the corrected code. Only the corrected script runs by default. The quartic group uses65-digit arithmetic but no interval enclosure. No score or prior model report is packaged here.

The root-authored `input*.md` sketches are included because the native scripts record their hashes; those hashes do not certify their mathematical claims. All code, sketches and outputs were produced with substantial AI assistance and a shared filesystem. Labels such as independent in historical script headers do not establish human or statistical independence. The new workflow must actually pass on each platform before cross-platform reproduction is claimed.

## Optional oscillator larger checks

`continuum_check51.py` evaluates explicit continuum integrals and six larger finite grids (up to 327681 levels); the original output is `continuum-results51.json`. This can take several minutes per case and is excluded from default CI. It is a diagnostic of convergence, not its proof. Its formulas do not establish uniqueness of all optimizers or a physical phase transition.

The scientific code was developed with AI assistance. Human peer review and formal verification are not claimed. References and manuscript assessments belong to the knowledge tree; they are deliberately separate from this finite computation.

Related context: https://github.com/lluiseriksson/physmath-knowledge-tree .
