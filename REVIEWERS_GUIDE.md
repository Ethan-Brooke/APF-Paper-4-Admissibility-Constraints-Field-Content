# Reviewers' Guide — APF Paper 4: Admissibility Constraints and Structural Saturation

This guide is the physics-first walkthrough for peer reviewers of Admissibility Constraints and Structural Saturation. It addresses the structural assumptions, anticipated objections, and falsification surfaces of the paper's argument. The intended audience is a peer reviewer with a background in physics or mathematics.

The guide complements three other verification routes:
- The [executable codebase](apf/) (run `python run_checks.py` for the full 40-theorem verification);
- The [Colab notebook](APF_Reviewer_Walkthrough.ipynb) (zero-install, theorem-by-theorem cells with prose + code);
- The [interactive derivation DAG](https://ethan-brooke.github.io/Admissibility-Constraints-Field-Content/) (hover-for-details, click-for-dependencies, animated verification).

---

## The derivation in 40 steps

The full chain from the foundational commitments to this paper's results:

**Step 1 — name** (`name`)

  


**Step 2 — Regime_R** (`check_Regime_R`)

  Regime_R: PLEC Well-Posedness under R1..R4 [P].


**Step 3 — Theorem_R** (`check_Theorem_R`)

  Theorem_R: Representation Requirements from Admissibility.


**Step 4 — L_gauge_template_uniqueness** (`check_L_gauge_template_uniqueness`)

  L_gauge_template_uniqueness: SU(N_c)×SU(2)×U(1) is the Unique Gauge Template.


**Step 5 — L_count** (`check_L_count`)

  L_count: Capacity Counting ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â 1 structural enforcement channel = 1 unit.


**Step 6 — L_nc** (`check_L_nc`)

  L_nc: Non-Closure from Admissibility Physics + Locality.


**Step 7 — T_M** (`check_T_M`)

  T_M: Interface Monogamy.


**Step 8 — L_area_scaling** (`check_L_area_scaling`)

  


**Step 9 — Regime_exit_Type_IV** (`check_Regime_exit_Type_IV`)

  Regime_exit_Type_IV: Loss of Smooth or Local Structure [P].


**Step 10 — T4F** (`check_T4F`)

  T4F: Flavor-Capacity Saturation.


**Step 11 — T7** (`check_T7`)

  T7: Generation Bound N_gen = 3 [P].


**Step 12 — T_BH_information** (`check_T_BH_information`)

  T_BH_information: Black Hole Information Preservation [P].


**Step 13 — L_BH_page_curve_capacity** (`check_L_BH_page_curve_capacity`)

  L_BH_page_curve_capacity: Page Curve from Capacity Counting [P].


**Step 14 — T_gauge** (`check_T_gauge`)

  T_gauge: SU(3)*SU(2)*U(1) from Capacity Budget.


**Step 15 — L_anomaly_free** (`check_L_anomaly_free`)

  L_anomaly_free: Gauge Anomaly Cancellation Cross-Check [P].


**Step 16 — L_Witten_parity** (`check_L_Witten_parity`)

  


**Step 17 — T_field** (`check_T_field`)

  T_field: SM Fermion Content -- Exhaustive Derivation.


**Step 18 — L_saturation_partition** (`check_L_saturation_partition`)

  L_saturation_partition: Type-Count Partition is Saturation-Independent [P].


**Step 19 — T_capacity_ladder** (`check_T_capacity_ladder`)

  T_capacity_ladder: Capacity Charges from Budget [P].


**Step 20 — L_FN_ladder_uniqueness** (`check_L_FN_ladder_uniqueness`)

  L_FN_ladder_uniqueness: q_B = (7,4,0) is Unique Cost-Minimal Partition [P].


**Step 21 — T25a** (`check_T25a`)

  T25a: Overlap Bounds from Interface Monogamy.


**Step 22 — T25b** (`check_T25b`)

  T25b: Overlap Bound from Saturation.


**Step 23 — T27d** (`check_T27d`)

  T27d: gamma_2/gamma_1 = d + 1/d from Representation Principles.


**Step 24 — T23** (`check_T23`)

  T23: Fixed-Point Formula for sin^2theta_W.


**Step 25 — T_sin2theta** (`check_T_sin2theta`)

  T_sin2theta: Weinberg Angle -- structurally derived from fixed point.


**Step 26 — L_capacity_per_dimension** (`check_L_capacity_per_dimension`)

  L_capacity_per_dimension: Neutrino d_1 = x^(q_B1/d_Y) [P].


**Step 27 — L_Yukawa_bilinear** (`check_L_Yukawa_bilinear`)

  L_Yukawa_bilinear: Yukawa Coupling Is Bilinear in Generation Amplitudes [P].


**Step 28 — T_mass_ratios** (`check_T_mass_ratios`)

  T_mass_ratios: Six Charged Fermion Mass Ratios from Zero Parameters [P].


**Step 29 — L_mass_from_capacity** (`check_L_mass_from_capacity`)

  L_mass_from_capacity: Complete Mass Matrix Derivation — Zero FN Imports [P].


**Step 30 — L_c_FN_gap** (`check_L_c_FN_gap`)

  L_c_FN_gap: NNLO Coefficient c = x^Δq from FN Charge Gap [P].


**Step 31 — T_CKM** (`check_T_CKM`)

  T_CKM: Zero-Parameter CKM Matrix Prediction [P].


**Step 32 — L_gen_path** (`check_L_gen_path`)

  L_gen_path: Generation Graph Is a Path [P].


**Step 33 — L_CP_channel** (`check_L_CP_channel`)

  L_CP_channel: Channel Asymmetry Enables CP Violation [P | L_H_curv, T_q_Higgs].


**Step 34 — T_PMNS** (`check_T_PMNS`)

  T_PMNS: Zero-Parameter PMNS Neutrino Mixing Matrix [P].


**Step 35 — T_PMNS_CP** (`check_T_PMNS_CP`)

  T_PMNS_CP: Leptonic CP Violation Vanishes Exactly [P].


**Step 36 — T_nu_ordering** (`check_T_nu_ordering`)

  T_nu_ordering: Normal Neutrino Mass Ordering [P].


**Step 37 — L_dm2_hierarchy** (`check_L_dm2_hierarchy`)

  L_dm2_hierarchy: Neutrino Mass-Squared Splitting Ratio [P].


**Step 38 — L_mbb_prediction** (`check_L_mbb_prediction`)

  L_mbb_prediction: Neutrinoless Double Beta Decay Effective Mass [P].


**Step 39 — L_DUNE_response** (`check_L_DUNE_response`)

  L_DUNE_response: APF δ_PMNS Prediction vs DUNE/Hyper-K Sensitivity [P].


**Step 40 — L_cosmological_constant** (`check_L_cosmological_constant`)

  



Each step is a single theorem. The dependencies form a directed acyclic graph (DAG) with no cycles; the [interactive DAG](https://ethan-brooke.github.io/Admissibility-Constraints-Field-Content/) shows multiple derivation paths to the same result, which is structural redundancy: many results are over-determined and would survive even if individual derivation steps were weakened.

---

## Pre-empting the operational-witness objection

A common objection to executable mathematical appendices is that the test-case witnesses (specific small examples used to verify the theorem) "smuggle in" the very structure being derived. Each test-case witness in this paper is the smallest concrete instance of a structure already proved in the abstract argument. Witnesses are constructed *after* the abstract claim is established, never before. The reviewer can re-derive each witness from the abstract derivation chain, not the other way around. If a witness appears to smuggle structure in, the corresponding `check_*` function will fail when the witness is replaced with an explicit countermodel.

---

## The structural assumptions

Beyond Axiom A1, this paper relies on the four PLEC components and a small set of structural regularity conditions:

- **A1 (Finite Enforcement Capacity)** — at any interface, the total enforcement cost across maintained distinctions is bounded.
- **MD (Minimum Distinction)** — every physical distinction has strictly positive enforcement cost (positive floor μ\* > 0).
- **A2 (Argmin selection)** — the realized configuration is the argmin of cost over the admissible set.
- **BW (Boundary Weight / Non-degeneracy)** — the cost spectrum is rich enough that the argmin is generically informative.

The four are pairwise logically independent (proved in Paper 1's Technical Supplement §1 with explicit countermodels). Together they define the variational structure named the **Principle of Least Enforcement Cost (PLEC)**: *reality is the minimum-cost expression of distinction compatible with finite capacity.*

_No paper-specific assumptions beyond the PLEC framework._

---

## The identifications

Several mappings between enforcement concepts and standard mathematical structures are used in the paper. They are **motivated, not deduced**, and reviewers should pay particular attention to their justification:

1. **Distinction ↔ binary partition of admissible configurations.** Spencer-Brown's distinction primitive (drawing the line between *this* and *not-this*) is identified with a binary measurement / partition at an enforcement interface. Motivated, not deduced.

2. **Cost ↔ resource expenditure.** The framework's $\varepsilon(d)$ functional is identified with the operational resource (energy / mass / barrier height) that maintains a distinction in any specific laboratory instance. This is a multi-resource correspondence; specific identifications appear in domain-specific papers (Paper 7 for action / Lagrangian, Paper 6 for spacetime metric, etc.).

3. **Capacity bound ↔ physical scarcity.** $C(\Gamma) < \infty$ at every interface is identified with the universally observed pattern that any specific physical system supports finitely many distinguishable states (Bekenstein--Hawking, Landauer, Ashby). Motivated by these traditions; not derived from them.

4. **Argmin ↔ realization.** The realized configuration at an interface is identified with the $\arg\min$ of total enforcement cost over the admissible set. This is a *locator*, not a process: nothing is doing the optimization; the argmin is what realization *means* under the framework.

5. **Hilbert space ↔ minimal representation of the noncommutative enforcement algebra.** Per the Paper 1 GNS argument, complex Hilbert space appears as the *unique* representational target for the algebraic structure that finite enforcement forces. Not assumed; selected via Frobenius.

These identifications are conventional rather than discovered: they are the choice that makes the framework's vocabulary and standard mathematics align. None is internally derivable.

---

## What is *not* proved in this paper

To prevent scope inflation, the following are flagged explicitly as *outside* the scope of Admissibility Constraints Field Content:

- **Specific gauge-group derivations** belong to Paper 2; this paper assumes them as imports if cited.
- **Specific particle content / generation count** belongs to Papers 2 and 4; not derived here.
- **Quantitative cosmological observables** beyond what is explicitly cited belong to Papers 3 (entropy / horizon) and 6 (geometry / DESI).
- **Quantum-gravity backreaction** is out of scope for any single paper in the current APF series; it is a future direction.
- **Numerical mass values** (absolute scales for $m_t$, $m_b$, $m_\nu$) are open problems noted in Paper 4 and the Engine paper; not within this paper's scope.
- **Spacetime dimension** (D = 4) is structurally derived in Paper 6 from Lovelock uniqueness; this paper assumes it where used.

If a reviewer concludes that the paper claims any of the above without supplying a proof for it, the reviewer is correct that the paper does not deliver that claim — those claims belong to other papers in the series and are explicitly flagged as such.

---

## How to falsify: attack surfaces

Each falsifier below targets a specific theorem or structural assumption. The corresponding code change in `apf/core.py` is also identified; reviewers can modify the codebase and re-run `python run_checks.py` to test each surface directly.

| # | Surface | Code-level test |
|---|---------|------------------|
| 1 | **Numerical disagreement.** A predicted observable disagrees with experiment beyond the framework's stated tolerance. | Modify the corresponding `check_*` to use the published experimental value; observe failure. |
| 2 | **Structural redundancy collapse.** Removing one PLEC component (A1, MD, A2, BW) leaves the derivation chain intact. | Comment out the test for the removed component in `apf/core.py`; observe other downstream checks failing. |
| 3 | **Reconstruction-program parity.** A standard quantum reconstruction (Hardy / CDP / Masanes--Müller) reaches the same conclusion with a strictly weaker assumption set. | Extract the GPT axioms used by the comparator; supply them as inputs to the relevant `apf/core.py` functions. |
| 4 | **Composition / locality break.** A multi-interface test with one interface deliberately violating $L_{\rm loc}$ does not produce the expected falsification mode. | Modify `check_L_loc` countermodel; rerun and observe expected vs. actual behavior. |
| 5 | **Cost-functional uniqueness fails.** An alternative cost functional satisfies all framework axioms equally well. | Replace the cost functional in `apf_utils.py`; observe whether downstream checks still pass. |
| 6 | **Scope-creep test.** A claim attributed to this paper is shown to actually require an unstated assumption. | Trace the claim's `\coderef`s through the bank dependency chain; identify any check that exits the paper's named scope. |

This is a structured threat model. If any of the surfaces fails empirically, the paper falsifies on that specific point.

---

## Reading the code

The codebase has three files in `apf/`:

- **`apf/core.py`** — the 40 theorem check functions for this paper. Each function constructs a mathematical witness, verifies the theorem's claim, and returns a structured result with name, dependencies, status, and key result.
- **`apf/apf_utils.py`** — exact arithmetic utilities (mostly `Fraction`-based; numpy/scipy where required by specific numerical lemmas).
- **`apf/bank.py`** — registry of all check functions in this repo, plus the `run_all()` runner.

Execution model: `run_checks.py` calls `bank.run_all()`, which iterates over every registered check and aggregates pass/fail/error counts. Individual checks can be invoked via `apf.bank.get_check('T_Born')()`.

---

## Scalar-to-matrix boundary

A characteristic feature of the APF program is that the early derivations use only finite sets and exact rational arithmetic (no matrices, no complex numbers, no linear algebra). Matrices first appear at the GNS construction (T2 in Paper 1), as the *minimal representation* of structure that the earlier scalar-only theorems prove must exist. This paper inherits the scalar-to-matrix transition from Paper 1's $T_2$. No new derivations of matrix structure occur in this paper; matrices appear where Paper 1 already established them.

This stratification is a deliberate methodological commitment: matrices are derived as representations of an already-proved abstract structure, not assumed as the substrate of the framework.

---

## The complex-field justification

The complex numbers as the unique admissible amplitude field is not a postulate but a derived selection (Paper 1 T2c, with proved exclusions $P_{\rm tom}$ and $P_{\rm cls}$ ruling out $\mathbb{R}$ and $\mathbb{H}$ respectively). This paper does not re-derive the complex-field selection; it inherits from Paper 1 ($T_{2c}$, with proved exclusions $P_{\rm tom}$ and $P_{\rm cls}$). For the original derivation, see Paper 1's [REVIEWERS_GUIDE.md](https://github.com/Ethan-Brooke/The-Enforceability-of-Distinction/blob/main/REVIEWERS_GUIDE.md).

---

## Citation and Zenodo

This repository is the executable mathematical appendix to APF Paper 4. The canonical archival deposit is at [https://doi.org/10.5281/zenodo.18604845](https://doi.org/10.5281/zenodo.18604845) (DOI: 10.5281/zenodo.18604845).

---

*Generated by the APF `create-repo` skill on 2026-04-19. Codebase snapshot: v6.9 (frozen 2026-04-18; 355 verify_all checks, 342 bank-registered theorems, 48 quantitative predictions).*
