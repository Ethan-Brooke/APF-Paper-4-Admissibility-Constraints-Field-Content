# Claims Ledger — Paper 4

| # | Claim | Status | Proof location | Code check | Failure mode |
|---|---|---|---|---|---|
| 1 | $T_M$ biconditional | nontrivial | Supp §T_M | `check_T_M_biconditional` | counterexample to biconditional direction |
| 2 | 1-of-4680 filter yields 45-fermion content | nontrivial (load-bearing) | Supp §T_{4F} | `check_T_field`, `check_L_4680_filter` | alternative candidate passes filter |
| 3 | $\sin^2 \theta_W = 3/13$ | nontrivial | Supp §$T_{25a/b} + T_{27d} + T_{23}$ | `check_T_sin2theta` | gate $S_0$ ambiguity |
| 4 | $m_c$ at 2.6% (structural limit) | structural | Supp §Schur | `check_L_Schur_chain_mc` | alternative Schur route |
| 5 | $m_s/m_b$ mass ratio | structural + match | Supp §mass | `check_T_mass_ratios` | PDG disagreement |
| 6 | $N_{\rm gen} = 3$ via capacity ladder | nontrivial | Supp §N_gen | `check_T_N_gen_3` | 4-generation admissible alternative |
| 7 | CKM + Cabibbo structure | nontrivial | Supp §CKM | `check_T_CKM` | measurement disagreement |
| 8 | PMNS structure | nontrivial | Supp §PMNS | `check_T_PMNS` | DUNE measurement disagreement |
| 9 | $C_{\rm ext}$ captures dark sector structurally | structural | Main §IV | `check_L_C_ext` | dark sector count wrong |
| 10 | $\Lambda = $ capacity residual | structural | Main §IV | `check_L_Lambda_residual` | non-residual $\Lambda$ alternative |
| 11 | Falsifier F1-F5 | empirical | Main §IV | — | any one confirmed null |

## Attack surface priority

Claims 2, 3, 4, 6. Claim 2 is the spine of Paper 4 — 1-of-4680 is what turns the gauge template into the observed SM.

---

*35 bank-registered checks verify this paper (full v2.0 Supplement rebuild).*
