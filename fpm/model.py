"""Versioned run-spec contract for the Quarantine Games simulator.

Mission-specific design code must emit this schema; it must not create a new
simulator. Topology is a graph artifact plus metadata. Run Spec V4 adds exact
certificate identity and boundary semantics at the raw-kernel trust boundary.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping
import hashlib, json, math, string
from decimal import Decimal, InvalidOperation

SCHEMA_VERSION = "FPM_RUN_SPEC_V6"
ANALYSIS_ROLES = {"primary", "trace_diagnostic", "benchmark", "calibration", "smoke_test"}

FINGERPRINT_VERSION = "FPM_TSV_ID_V3"
_FINGERPRINT_EXCLUDE = {
    "run_fingerprint", "run_index", "run_id", "graph_path",
    "graph_certificate_path", "fingerprint_version",
}
_INT_ID_FIELDS = {"population","n0","n_occ","n_gc","graph_seed","q_init_seed","i_init_seed","dynamic_seed","signal_seed","information_radius","event_limit"}
_FLOAT_ID_FIELDS = {"beta_n","beta_a","beta_q","beta_mix_a","q0","i0","tau","i_c","fermi_k","omega","delta","information_alpha","behavior_bin_width","tmax"}

def semantic_identity_payload(row: Mapping[str,str]) -> dict:
    """Canonical typed identity projection for future run specs.

    Human-readable labels, filesystem paths and run ordering are excluded. Graph bytes
    enter through graph_sha256; stochastic replicate identity enters through the four
    domain-separated seeds and graph_seed.
    """
    payload={}
    for key in REQUIRED_COLUMNS:
        if key in _FINGERPRINT_EXCLUDE: continue
        value=row[key]
        if key in _INT_ID_FIELDS: payload[key]=int(value)
        elif key in _FLOAT_ID_FIELDS: payload[key]=float(value)
        else: payload[key]=value
    return payload

def semantic_run_fingerprint(row: Mapping[str,str]) -> str:
    """Return the versioned semantic SHA-256 identity for a future run-spec row."""
    return hashlib.sha256(canonical_json(semantic_identity_payload(row))).hexdigest()
REQUIRED_COLUMNS = (
    "schema_version","analysis_role","run_fingerprint","run_index","run_id","stage",
    "topology","graph_path","graph_sha256","graph_certificate_path",
    "graph_certificate_sha256","boundary_condition","population","n0","n_occ","n_gc",
    "graph_seed","q_init_seed","i_init_seed","dynamic_seed","contact_mode","beta_n","beta_a","beta_q","beta_mix_a","beta_mixing_rule",
    "profile","arm","protocol","q0","i0","tau","i_c","fermi_k",
    "omega","delta","seed_policy","behavior_mode","q_placement","information_mode",
    "information_radius","information_alpha","signal_seed","instrumentation_mode",
    "trace_mode","behavior_bin_width","event_limit","tmax",
)


def canonical_json(payload: object) -> bytes:
    """Canonical JSON bytes used as the identity surface for run fingerprints."""
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False).encode()


def stable_seed(text: str) -> int:
    """Map an arbitrary identity string to a nonzero 64-bit master seed."""
    return 1 + int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "big") % (2**64 - 2)


def derive_seed(master: int, domain: str) -> int:
    """Derive a deterministic domain-separated 64-bit seed."""
    return 1 + int.from_bytes(hashlib.sha256(f"{master}|{domain}".encode()).digest()[:8], "big") % (2**64 - 2)


def normalize_tau_aliases(row: Mapping[str, str]) -> dict[str, str]:
    """Normalize the behavioral-rate alias at ingestion only.

    New persisted V5 artifacts use ``tau`` only. Legacy ``chi`` is accepted as
    an input alias. If both are supplied they must denote exactly the same finite
    decimal value; approximate equality is deliberately not used.
    """
    out=dict(row)
    has_tau="tau" in out and out["tau"]!=""
    has_chi="chi" in out and out["chi"]!=""
    if not has_tau and not has_chi:
        return out
    if has_tau and has_chi:
        try:
            tau=Decimal(out["tau"]); chi=Decimal(out["chi"])
        except InvalidOperation as exc:
            raise ValueError("invalid tau/chi alias value") from exc
        if not tau.is_finite() or not chi.is_finite() or tau != chi:
            raise ValueError("conflicting tau and legacy chi values")
    elif has_chi:
        out["tau"]=out["chi"]
    out.pop("chi",None)
    return out


def normalize_fermi_k_aliases(row: Mapping[str, str]) -> dict[str, str]:
    """Normalize historical Fermi-noise ``epsilon`` to current ``fermi_k``.

    New persisted V6 artifacts use ``fermi_k`` only. If both names are supplied,
    they must denote exactly the same finite decimal value.
    """
    out=dict(row)
    has_new="fermi_k" in out and out["fermi_k"]!=""
    has_old="epsilon" in out and out["epsilon"]!=""
    if not has_new and not has_old:
        return out
    if has_new and has_old:
        try:
            new=Decimal(out["fermi_k"]); old=Decimal(out["epsilon"])
        except InvalidOperation as exc:
            raise ValueError("invalid fermi_k/epsilon alias value") from exc
        if not new.is_finite() or not old.is_finite() or new != old:
            raise ValueError("conflicting fermi_k and legacy epsilon values")
    elif has_old:
        out["fermi_k"]=out["epsilon"]
    out.pop("epsilon",None)
    return out


def run_fingerprint(full_spec: Mapping) -> str:
    """Return SHA-256(canonical run spec)."""
    return hashlib.sha256(canonical_json(dict(full_spec))).hexdigest()


@dataclass(frozen=True)
class RuntimeContract:
    """Stable runtime switches. Experimental backends are explicit and opt-in."""
    pool_backend: str = "fenwick_exact"
    information_mode: str = "global"
    information_radius: int = 1
    instrumentation_mode: str = "on"
    trace_mode: str = "off"

    def validate(self) -> None:
        if self.pool_backend not in {"fenwick_exact", "uniform_regular"}:
            raise ValueError(f"unknown pool backend: {self.pool_backend}")
        if self.information_mode not in {"global", "local", "local_shuffled", "mixed"}:
            raise ValueError(f"unknown information mode: {self.information_mode}")
        if not 1 <= int(self.information_radius) <= 16:
            raise ValueError("information radius must be in [1,16]")
        if self.instrumentation_mode not in {"on", "off"}:
            raise ValueError("instrumentation_mode must be on/off")
        if self.trace_mode not in {"on", "off"}:
            raise ValueError("trace_mode must be on/off")


def _finite_float(row: Mapping[str,str], key: str) -> float:
    value=float(row[key])
    if not math.isfinite(value): raise ValueError(f"non-finite {key}")
    return value


def validate_design_row(row: Mapping[str, str], *, stable_surface: bool = True, allow_legacy_identity: bool = False) -> None:
    """Validate one corrected run-spec row without silently coercing semantics.

    `stable_surface=True` enforces the M20B-certified branch: normalized contact,
    uniform seeding, global information and fixed-count independent Q/I streams.
    Experimental/generalized branches require separate audit authority.
    """
    missing=[c for c in REQUIRED_COLUMNS if c not in row or row[c]==""]
    if missing: raise ValueError(f"missing required run-spec fields: {missing}")
    if row["schema_version"] != SCHEMA_VERSION: raise ValueError("unsupported schema_version")
    if row["analysis_role"] not in ANALYSIS_ROLES: raise ValueError("unknown analysis_role")
    hexchars=set(string.hexdigits)
    if (len(row["run_fingerprint"]) != 64 or len(row["graph_sha256"]) != 64
        or len(row["graph_certificate_sha256"]) != 64
        or any(ch not in hexchars for ch in row["run_fingerprint"])
        or any(ch not in hexchars for ch in row["graph_sha256"])
        or any(ch not in hexchars for ch in row["graph_certificate_sha256"])):
        raise ValueError("fingerprint/hash must be 64-character hexadecimal SHA-256")
    fp_version=row.get("fingerprint_version","")
    if not fp_version and not allow_legacy_identity:
        raise ValueError("legacy opaque run identity requires allow_legacy_identity=True")
    if fp_version:
        if fp_version != FINGERPRINT_VERSION:
            raise ValueError("unsupported fingerprint_version")
        identity_strings={k: str(row[k]) for k in REQUIRED_COLUMNS
                          if k not in _FINGERPRINT_EXCLUDE | _INT_ID_FIELDS | _FLOAT_ID_FIELDS}
        if any(not value.isascii() for value in identity_strings.values()):
            raise ValueError("FPM_TSV_ID_V3 semantic identity strings must be ASCII")
        expected=semantic_run_fingerprint(row)
        if row["run_fingerprint"].lower()!=expected:
            raise ValueError("semantic run_fingerprint mismatch")
    n=int(row["population"]); n0=int(row["n0"]); nocc=int(row["n_occ"]); ngc=int(row["n_gc"])
    if n < 2 or not (n0 >= nocc >= ngc == n): raise ValueError("invalid population accounting")
    lam=_finite_float(row,"beta_n"); q0=_finite_float(row,"q0"); i0=_finite_float(row,"i0")
    if lam <= 0 or not 0 <= q0 <= 1 or not 0 <= i0 <= 1: raise ValueError("invalid epidemiological/init parameter")
    bq=_finite_float(row,"beta_q"); ba=_finite_float(row,"beta_a"); a=_finite_float(row,"beta_mix_a")
    if row["beta_mixing_rule"] != "INTERPOLATED": raise ValueError("corrected V3 stable surface requires INTERPOLATED")
    if not (bq>0 and lam>ba>bq and 0<a<1): raise ValueError("beta hierarchy/mixing contract violated")
    expected_ba=bq+a*(lam-bq)
    if abs(ba-expected_ba)>1e-12*max(1,abs(expected_ba)): raise ValueError("beta_a resolver mismatch")
    if _finite_float(row,"tau") < 0: raise ValueError("negative tau")
    for key in ("i_c","fermi_k","omega","delta","behavior_bin_width","tmax"):
        if _finite_float(row,key) <= 0: raise ValueError(f"non-positive {key}")
    seeds=[int(row[k]) for k in ("q_init_seed","i_init_seed","dynamic_seed","signal_seed")]
    if len(set(seeds)) != 4:
        raise ValueError("Q_INIT, I_INIT, DYNAMICS and SIGNAL seeds must differ")
    if int(row["event_limit"]) <= 0: raise ValueError("event_limit must be positive")
    if row["trace_mode"] == "on" and row["analysis_role"] != "trace_diagnostic":
        raise ValueError("trace output must be segregated as trace_diagnostic")
    spatial={"honeycomb","truncated_square_488","4.8.8","square","kagome","triangular"}
    zero_sha="0"*64
    if row["topology"] in spatial:
        if row["graph_certificate_path"] != row["graph_path"]+".periodic.json":
            raise ValueError("spatial certificate path must be graph_path + .periodic.json")
        if row["graph_certificate_sha256"] == zero_sha:
            raise ValueError("spatial certificate identity cannot be the null sentinel")
        if row["boundary_condition"] != "periodic_torus":
            raise ValueError("spatial topology requires periodic_torus")
    else:
        expected_boundary="well_mixed" if row["topology"]=="well_mixed" else "nonspatial_network"
        if (row["graph_certificate_path"] != "-"
            or row["graph_certificate_sha256"] != zero_sha
            or row["boundary_condition"] != expected_boundary):
            raise ValueError("nonspatial graph certificate sentinel/boundary mismatch")
    if stable_surface:
        if row["contact_mode"] != "normalized": raise ValueError("stable V4 surface currently certifies normalized only")
        if row["seed_policy"] != "uniform": raise ValueError("stable V4 surface currently certifies uniform seeds only")
        if row["q_placement"] != "none": raise ValueError("stable V4 surface currently certifies fixed-fraction Q only")
        if row["information_mode"] not in {"global","local","local_shuffled","mixed"}: raise ValueError("unsupported information mode")
