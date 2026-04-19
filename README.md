# Admissibility Constraints and Structural Saturation

### Interactive Mathematical Appendix to Paper 4 of the Admissibility Physics Framework

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18604845.svg)](https://doi.org/10.5281/zenodo.18604845) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Ethan-Brooke/APF-Paper-4-Admissibility-Constraints-Field-Content/blob/main/APF_Reviewer_Walkthrough.ipynb)

[Interactive Derivation DAG](https://ethan-brooke.github.io/APF-Paper-4-Admissibility-Constraints-Field-Content/) · [Theorem Map](#theorem-mapping-table) · [Reviewers' Guide](REVIEWERS_GUIDE.md) · [The full APF corpus](#the-full-apf-corpus) · [Citation](#citation)

> **AI agents:** start with [`START_HERE.md`](START_HERE.md) — operational checklist that loads the framework context in 5–10 minutes. The corpus inventory and full file map are in [`ai_context/repo_map.json`](ai_context/repo_map.json).

---

## Why this codebase exists

Field content and fermion spectrum of the Standard Model: 45 fermions in three generations, the 1-of-4680 admissible template, mass-ratio structure, Cabibbo angle, CKM matrix structure. Currently PDF-only; .tex source pending.

This repository is the executable proof.

The codebase is a faithful subset of the canonical APF codebase v6.9 (frozen 2026-04-18; 355 verify_all checks, 342 bank-registered theorems across 19 modules + `apf/standalone/`). Each theorem in the manuscript traces to a named `check_*` function in `apf/core.py`, which can be called independently and returns a structured result.

The codebase requires Python 3.8+ and NumPy / SciPy (some numerical lemmas use them; see `pyproject.toml`).

## How to verify

Three paths, in order of increasing friction:

**1. Colab notebook — zero install.** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Ethan-Brooke/APF-Paper-4-Admissibility-Constraints-Field-Content/blob/main/APF_Reviewer_Walkthrough.ipynb) Every key theorem is derived inline, with annotated cells you can inspect and modify. Run all cells — the full verification takes under a minute.

**2. Browser — zero install.** Open the [Interactive Derivation DAG](https://ethan-brooke.github.io/APF-Paper-4-Admissibility-Constraints-Field-Content/). Explore the dependency graph. Hover any node for its mathematical statement, key result, and shortest derivation chain to A1. Click **Run Checks** to watch all theorems verify in topological order.

**3. Local execution.**

```bash
git clone https://github.com/Ethan-Brooke/APF-Paper-4-Admissibility-Constraints-Field-Content.git
cd APF-Paper-4-Admissibility-Constraints-Field-Content
pip install -e .
python run_checks.py
```

Expected output:

```
      Paper 4 (Admissibility Constraints and Structural Saturation): 40 passed, 0 failed, 40 total — verified in <minutes>
```

**4. Individual inspection.**

```python
from apf.bank import get_check
r = get_check('name')()
print(r['key_result'])
```

For reviewers, a [dedicated guide](REVIEWERS_GUIDE.md) walks through the logical architecture, the structural assumptions, and the anticipated objections.

---

## Theorem mapping table

This table maps every result in the manuscript to its executable verification.

| Check | Type | Summary |
|-------|------|---------|
| `name` | Other |  |
| `check_Regime_R` | Other | Regime_R: PLEC Well-Posedness under R1..R4 [P]. |
| `check_Theorem_R` | Theorem | Theorem_R: Representation Requirements from Admissibility. |
| `check_L_gauge_template_uniqueness` | Lemma | L_gauge_template_uniqueness: SU(N_c)×SU(2)×U(1) is the Unique Gauge Template. |
| `check_L_count` | Lemma | L_count: Capacity Counting ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â 1 structural enforcement channel = 1 unit. |
| `check_L_nc` | Lemma | L_nc: Non-Closure from Admissibility Physics + Locality. |
| `check_T_M` | Theorem | T_M: Interface Monogamy. |
| `check_L_area_scaling` | Lemma |  |
| `check_Regime_exit_Type_IV` | Other | Regime_exit_Type_IV: Loss of Smooth or Local Structure [P]. |
| `check_T4F` | Theorem | T4F: Flavor-Capacity Saturation. |
| `check_T7` | Theorem | T7: Generation Bound N_gen = 3 [P]. |
| `check_T_BH_information` | Theorem | T_BH_information: Black Hole Information Preservation [P]. |
| `check_L_BH_page_curve_capacity` | Lemma | L_BH_page_curve_capacity: Page Curve from Capacity Counting [P]. |
| `check_T_gauge` | Theorem | T_gauge: SU(3)*SU(2)*U(1) from Capacity Budget. |
| `check_L_anomaly_free` | Lemma | L_anomaly_free: Gauge Anomaly Cancellation Cross-Check [P]. |
| `check_L_Witten_parity` | Lemma |  |
| `check_T_field` | Theorem | T_field: SM Fermion Content -- Exhaustive Derivation. |
| `check_L_saturation_partition` | Lemma | L_saturation_partition: Type-Count Partition is Saturation-Independent [P]. |
| `check_T_capacity_ladder` | Theorem | T_capacity_ladder: Capacity Charges from Budget [P]. |
| `check_L_FN_ladder_uniqueness` | Lemma | L_FN_ladder_uniqueness: q_B = (7,4,0) is Unique Cost-Minimal Partition [P]. |
| `check_T25a` | Theorem | T25a: Overlap Bounds from Interface Monogamy. |
| `check_T25b` | Theorem | T25b: Overlap Bound from Saturation. |
| `check_T27d` | Theorem | T27d: gamma_2/gamma_1 = d + 1/d from Representation Principles. |
| `check_T23` | Theorem | T23: Fixed-Point Formula for sin^2theta_W. |
| `check_T_sin2theta` | Theorem | T_sin2theta: Weinberg Angle -- structurally derived from fixed point. |
| `check_L_capacity_per_dimension` | Lemma | L_capacity_per_dimension: Neutrino d_1 = x^(q_B1/d_Y) [P]. |
| `check_L_Yukawa_bilinear` | Lemma | L_Yukawa_bilinear: Yukawa Coupling Is Bilinear in Generation Amplitudes [P]. |
| `check_T_mass_ratios` | Theorem | T_mass_ratios: Six Charged Fermion Mass Ratios from Zero Parameters [P]. |
| `check_L_mass_from_capacity` | Lemma | L_mass_from_capacity: Complete Mass Matrix Derivation — Zero FN Imports [P]. |
| `check_L_c_FN_gap` | Lemma | L_c_FN_gap: NNLO Coefficient c = x^Δq from FN Charge Gap [P]. |
| `check_T_CKM` | Theorem | T_CKM: Zero-Parameter CKM Matrix Prediction [P]. |
| `check_L_gen_path` | Lemma | L_gen_path: Generation Graph Is a Path [P]. |
| `check_L_CP_channel` | Lemma | L_CP_channel: Channel Asymmetry Enables CP Violation [P \| L_H_curv, T_q_Higgs]. |
| `check_T_PMNS` | Theorem | T_PMNS: Zero-Parameter PMNS Neutrino Mixing Matrix [P]. |
| `check_T_PMNS_CP` | Theorem | T_PMNS_CP: Leptonic CP Violation Vanishes Exactly [P]. |
| `check_T_nu_ordering` | Theorem | T_nu_ordering: Normal Neutrino Mass Ordering [P]. |
| `check_L_dm2_hierarchy` | Lemma | L_dm2_hierarchy: Neutrino Mass-Squared Splitting Ratio [P]. |
| `check_L_mbb_prediction` | Lemma | L_mbb_prediction: Neutrinoless Double Beta Decay Effective Mass [P]. |
| `check_L_DUNE_response` | Lemma | L_DUNE_response: APF δ_PMNS Prediction vs DUNE/Hyper-K Sensitivity [P]. |
| `check_L_cosmological_constant` | Lemma |  |

All check functions reside in `apf/core.py`. Every function listed above can be called independently and returns a structured result including its logical dependencies and the mathematical content it verifies.

---

## The derivation chain

```
  Level 0: name · Regime_R · Theorem_R · L_gauge_template_uniqueness · L_count · L_nc · T_M · L_area_scaling · Regime_exit_Type_IV · T4F · T7 · T_BH_information · L_BH_page_curve_capacity · T_gauge · L_anomaly_free · L_Witten_parity · T_field · L_saturation_partition · T_capacity_ladder · L_FN_ladder_uniqueness · T25a · T25b · T27d · T23 · T_sin2theta · L_capacity_per_dimension · L_Yukawa_bilinear · T_mass_ratios · L_mass_from_capacity · L_c_FN_gap · T_CKM · L_gen_path · L_CP_channel · T_PMNS · T_PMNS_CP · T_nu_ordering · L_dm2_hierarchy · L_mbb_prediction · L_DUNE_response · L_cosmological_constant
```

The [interactive DAG](https://ethan-brooke.github.io/APF-Paper-4-Admissibility-Constraints-Field-Content/) shows the full graph with hover details and animated verification.

---

## Repository structure

```
├── README.md                              ← you are here
├── START_HERE.md                          ← AI operational checklist; read-first for AI agents
├── REVIEWERS_GUIDE.md                     ← physics-first walkthrough for peer reviewers
├── interactive_dag.html                   ← interactive D3.js derivation DAG (also served at docs/ via GitHub Pages)
├── repo_map.json                          ← machine-readable map of this repo (root copy of ai_context/repo_map.json)
├── theorems.json                          ← theorem catalog (root copy of ai_context/theorems.json)
├── derivation_graph.json                  ← theorem DAG as JSON (root copy of ai_context/derivation_graph.json)
├── ai_context/                            ← AI onboarding pack (corpus map, theorems, glossary, etc.)
│   ├── AGENTS.md                          ← authoritative entry point for AI agents
│   ├── FRAMEWORK_OVERVIEW.md              ← APF in 5 minutes
│   ├── GLOSSARY.md                        ← axioms, PLEC primitives, epistemic tags
│   ├── AUDIT_DISCIPLINE.md                ← engagement posture for critique/proposal
│   ├── OPEN_PROBLEMS.md                   ← catalog of open problems + verdicts
│   ├── repo_map.json                      ← machine-readable map of this repo
│   ├── theorems.json                      ← machine-readable theorem catalog
│   ├── derivation_graph.json              ← theorem DAG as JSON
│   └── wiki/                              ← bundled APF wiki (concepts, papers, codebase)
├── apf/
│   ├── core.py                            ← 40 theorem check functions
│   ├── apf_utils.py                       ← exact arithmetic + helpers
│   └── bank.py                            ← registry and runner
├── docs/
│   └── index.html                         ← interactive derivation DAG (GitHub Pages)
├── APF_Reviewer_Walkthrough.ipynb         ← Colab notebook
├── run_checks.py                          ← convenience entry point
├── pyproject.toml                         ← package metadata
├── zenodo.json                            ← archival metadata
├── Paper_4_Admissibility_Constraints_Field_Content_v2.0.tex                ← the paper
├── Paper_4_Admissibility_Constraints_Field_Content_Supplement_v1.0.tex                ← Technical Supplement

└── LICENSE                                ← MIT
```

---

## What this paper derives and what it does not

**Derived:** (see Theorem mapping table above)

**Not derived here:** Specific results outside this paper's scope live in companion papers — see the corpus table below for the full 9-paper series.

---

## Citation

```bibtex
@software{apf-paper4,
  title   = {Admissibility Constraints and Structural Saturation},
  author  = {Brooke, Ethan},
  year    = {2026},
  doi     = {10.5281/zenodo.18604845},
  url     = {https://github.com/Ethan-Brooke/APF-Paper-4-Admissibility-Constraints-Field-Content}
}
```

For the full citation lineage (concept-DOI vs version-DOI, related identifiers, bibtex for all corpus papers), see [`ai_context/CITING.md`](ai_context/CITING.md).

---

## The full APF corpus

This repository is **one paper-companion** in a 9-paper series. Each paper has its own companion repo following this same layout. The full corpus, with canonical references:

| # | Title | Zenodo DOI | GitHub repo | Status |
|---|---|---|---|---|
| 0 | What Physics Permits | [10.5281/zenodo.18605692](https://doi.org/10.5281/zenodo.18605692) | [`APF-Paper-0-What-Physics-Permits`](https://github.com/Ethan-Brooke/APF-Paper-0-What-Physics-Permits) | public |
| 1 | The Enforceability of Distinction | [10.5281/zenodo.18604678](https://doi.org/10.5281/zenodo.18604678) | [`APF-Paper-1-The-Enforceability-of-Distinction`](https://github.com/Ethan-Brooke/APF-Paper-1-The-Enforceability-of-Distinction) | public |
| 2 | The Structure of Admissible Physics | [10.5281/zenodo.18604839](https://doi.org/10.5281/zenodo.18604839) | [`APF-Paper-2-The-Structure-of-Admissible-Physics`](https://github.com/Ethan-Brooke/APF-Paper-2-The-Structure-of-Admissible-Physics) | public |
| 3 | Ledgers | [10.5281/zenodo.18604844](https://doi.org/10.5281/zenodo.18604844) | [`APF-Paper-3-Ledgers-Entropy-Time-Cost`](https://github.com/Ethan-Brooke/APF-Paper-3-Ledgers-Entropy-Time-Cost) | public |
| 4 | Admissibility Constraints and Structural Saturation **(this repo)** | [10.5281/zenodo.18604845](https://doi.org/10.5281/zenodo.18604845) | [`APF-Paper-4-Admissibility-Constraints-Field-Content`](https://github.com/Ethan-Brooke/APF-Paper-4-Admissibility-Constraints-Field-Content) | public |
| 5 | Quantum Structure from Finite Enforceability | [10.5281/zenodo.18604861](https://doi.org/10.5281/zenodo.18604861) | [`APF-Paper-5-Quantum-Structure-Hilbert-Born`](https://github.com/Ethan-Brooke/APF-Paper-5-Quantum-Structure-Hilbert-Born) | public |
| 6 | Dynamics and Geometry as Optimal Admissible Reallocation | [10.5281/zenodo.18604874](https://doi.org/10.5281/zenodo.18604874) | [`APF-Paper-6-Dynamics-Geometry-Spacetime-Gravity`](https://github.com/Ethan-Brooke/APF-Paper-6-Dynamics-Geometry-Spacetime-Gravity) | public |
| 7 | Action, Internalization, and the Lagrangian | [10.5281/zenodo.18604875](https://doi.org/10.5281/zenodo.18604875) | [`APF-Paper-7-Action-Internalization-Lagrangian`](https://github.com/Ethan-Brooke/APF-Paper-7-Action-Internalization-Lagrangian) | public |
| 13 | The Minimal Admissibility Core | [10.5281/zenodo.18614663](https://doi.org/10.5281/zenodo.18614663) | [`APF-Paper-13-The-Minimal-Admissibility-Core`](https://github.com/Ethan-Brooke/APF-Paper-13-The-Minimal-Admissibility-Core) | public |
| — | Canonical codebase (v6.9) | [10.5281/zenodo.18604548](https://doi.org/10.5281/zenodo.18604548) | [`APF-Codebase`](https://github.com/Ethan-Brooke/APF-Codebase) | pending |

The canonical computational engine — the full bank of 342 theorems across 19 modules — is the **APF Codebase** ([Zenodo](https://doi.org/10.5281/zenodo.18604548)). Every per-paper repo is a faithful subset of that engine.

---

## License

MIT. See [LICENSE](LICENSE).

---

*Generated by the APF `create-repo` skill on 2026-04-19. Codebase snapshot: v6.9 (frozen 2026-04-18; 355 verify_all checks, 342 bank-registered theorems, 48 quantitative predictions).*
