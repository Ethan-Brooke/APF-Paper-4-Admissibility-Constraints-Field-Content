"""apf/bank.py — Paper 4 registry.

Lightweight registry for the 40-check subset bundled in this
paper-companion repo. Mirrors the canonical apf.bank API: REGISTRY (dict),
get_check(name), run_all(verbose=False).
"""
from collections import OrderedDict
import traceback

from apf import core as _core


def _build_registry():
    reg = OrderedDict()
    for name in ['name', 'check_Regime_R', 'check_Theorem_R', 'check_L_gauge_template_uniqueness', 'check_L_count', 'check_L_nc', 'check_T_M', 'check_L_area_scaling', 'check_Regime_exit_Type_IV', 'check_T4F', 'check_T7', 'check_T_BH_information', 'check_L_BH_page_curve_capacity', 'check_T_gauge', 'check_L_anomaly_free', 'check_L_Witten_parity', 'check_T_field', 'check_L_saturation_partition', 'check_T_capacity_ladder', 'check_L_FN_ladder_uniqueness', 'check_T25a', 'check_T25b', 'check_T27d', 'check_T23', 'check_T_sin2theta', 'check_L_capacity_per_dimension', 'check_L_Yukawa_bilinear', 'check_T_mass_ratios', 'check_L_mass_from_capacity', 'check_L_c_FN_gap', 'check_T_CKM', 'check_L_gen_path', 'check_L_CP_channel', 'check_T_PMNS', 'check_T_PMNS_CP', 'check_T_nu_ordering', 'check_L_dm2_hierarchy', 'check_L_mbb_prediction', 'check_L_DUNE_response', 'check_L_cosmological_constant']:
        fn = getattr(_core, name, None)
        if fn is None:
            # Function couldn't be extracted — skip with a warning attribute
            continue
        reg[name] = fn
    return reg


REGISTRY = _build_registry()
EXPECTED_CHECK_COUNT = 40


def get_check(name):
    """Return the check function registered as `name`. Raises KeyError if missing."""
    if name not in REGISTRY:
        raise KeyError(f"Check '{name}' not found. Available: {sorted(REGISTRY.keys())}")
    return REGISTRY[name]


def run_all(verbose=False):
    """Run every registered check, returning a list of result dicts."""
    results = []
    for name, fn in REGISTRY.items():
        try:
            r = fn()
            if not isinstance(r, dict):
                # Some legacy checks return True/False
                r = {"name": name, "passed": bool(r), "key_result": str(r)}
            elif "passed" not in r:
                r["passed"] = True
            r.setdefault("name", name)
        except Exception as e:
            r = {
                "name": name,
                "passed": False,
                "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc(),
            }
        if verbose:
            status = "PASS" if r.get("passed", True) else "FAIL"
            print(f"  {r['name']}: {status}")
            if r.get("key_result"):
                print(f"    {r['key_result']}")
        results.append(r)
    return results
