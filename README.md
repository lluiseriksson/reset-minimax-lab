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

## Optional larger checks

`continuum_check51.py` evaluates explicit continuum integrals and six larger finite grids (up to 327681 levels); the original output is `continuum-results51.json`. This can take several minutes per case and is excluded from default CI. It is a diagnostic of convergence, not its proof. Its formulas do not establish uniqueness of all optimizers or a physical phase transition.

The scientific code was developed with AI assistance. Human peer review and formal verification are not claimed. References and manuscript assessments belong to the knowledge tree; they are deliberately separate from this finite computation.

Related context: https://github.com/lluiseriksson/physmath-knowledge-tree .
