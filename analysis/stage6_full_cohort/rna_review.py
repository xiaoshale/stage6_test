"""Make a capture-level RNA QC/seed-identity review pack; never filter source data.

Input JSON: {"records": [{"path": "capture.h5ad", "resource": "GSE...",
"library_id": "...", "donor_id": "...", "assay": "scRNA|snRNA",
"fraction": null, "source_label_columns": [], "capture_identity_status": "..."}]}.
Paths may be relative to the manifest. Each checkpoint must contain one capture.
Outputs preserve source obs, feature axes and every profile. All identities remain
provisional. No integration, state discovery, disease test or final eligibility.
"""

import argparse
from datetime import datetime, timezone
import gc
import hashlib
from importlib import metadata, util
import json
from pathlib import Path
import re
import time
import traceback
import warnings

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc


SEED = 0
MAD_MULTIPLIERS = (3, 4, 5)
PRIORS = (0.02, 0.05, 0.10)
PREFIX = "review_"
# Fixed normal-identity review seeds from the locked dictionary, not learned
# disease signatures. Anchors guard against ACTA2/PDPN/keratin-only assignments.
PROGRAMS = {
    "epithelial_unspecified": ("epithelial", "EPCAM KRT8 KRT18 KRT19 KRT17", "EPCAM KRT8 KRT18 KRT19 KRT17"),
    "basal_like": ("epithelial", "KRT5 TP63 KRT15", "KRT5 TP63"),
    "club_like": ("epithelial", "SCGB1A1 SCGB3A1", "SCGB1A1 SCGB3A1"),
    "mucous_secretory": ("epithelial", "MUC5AC MUC5B SPDEF AGR2", "MUC5AC MUC5B SPDEF"),
    "multiciliated": ("epithelial", "FOXJ1 TPPP3 PIFO RSPH1", "FOXJ1 PIFO RSPH1"),
    "AT1_comparator": ("epithelial", "AGER CAV1 CLIC5 RTKN2", "AGER CLIC5 RTKN2"),
    "AT2_comparator": ("epithelial", "SFTPC SFTPA1 SFTPA2 ABCA3 LAMP3", "SFTPC ABCA3"),
    "fibroblast_lineage": ("stromal", "DCN LUM COL1A1 COL1A2 PDGFRA", "DCN LUM PDGFRA"),
    "pericyte_like": ("stromal", "RGS5 CSPG4 MCAM PDGFRB", "RGS5 CSPG4"),
    "smooth_muscle_like": ("stromal", "MYH11 CNN1 TAGLN ACTA2 MYL9", "MYH11 CNN1"),
    "myeloid_background": ("immune", "PTPRC TYROBP LYZ FCER1G LST1", "PTPRC TYROBP FCER1G LST1"),
    "lymphoid_background": ("immune", "PTPRC CD3D CD3E CD79A MS4A1 NKG7 GNLY", "PTPRC CD3D CD3E CD79A MS4A1 NKG7 GNLY"),
    "vascular_endothelium": ("endothelial", "PECAM1 VWF CLDN5 KDR", "PECAM1 VWF CLDN5 KDR"),
    "lymphatic_guard": ("endothelial", "PROX1 FLT4 CCL21 PDPN", "PROX1 FLT4"),
    "mesothelial_guard": ("other", "MSLN WT1 KRT8 KRT18", "MSLN WT1"),
}
MARKERS = sorted({g for _, genes, _ in PROGRAMS.values() for g in genes.split()})
BOUNDARY = (
    "Human-review preparation only. Every source profile retained. Zero-count profiles are hard-fail flagged; "
    "other QC flags and all seed-marker identities require review. No expert PASS, final eligibility, "
    "disease comparison, state/niche discovery or target ranking is produced."
)


def log(event, **fields):
    print(json.dumps({"utc": datetime.now(timezone.utc).isoformat(), "event": event, **fields}), flush=True)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def axis_hash(values):
    digest = hashlib.sha256()
    for value in values:
        raw = str(value).encode("utf-8")
        digest.update(len(raw).to_bytes(8, "little"))
        digest.update(raw)
    return digest.hexdigest()


def matrix_hash(matrix):
    digest = hashlib.sha256(str(matrix.shape).encode())
    for array in (matrix.data, matrix.indices, matrix.indptr):
        digest.update(str(array.dtype).encode())
        raw = memoryview(np.ascontiguousarray(array)).cast("B")
        for start in range(0, len(raw), 4 * 1024**2):
            digest.update(raw[start:start + 4 * 1024**2])
    return digest.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def finite_number(value):
    return float(value) if np.isfinite(value) else None


def qc_metrics(counts, mt):
    total = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    detected = np.diff(counts.indptr)
    mitochondrial = np.asarray(counts[:, mt].sum(axis=1)).ravel().astype(np.float64)
    percent = np.divide(100 * mitochondrial, total, out=np.full(len(total), np.nan), where=total > 0)
    return pd.DataFrame({"total_counts": total, "detected_source_features": detected,
                         "mitochondrial_counts": mitochondrial, "pct_mt": percent,
                         "zero_counts": total == 0})


def mad_flags(metrics):
    positive = metrics["total_counts"].to_numpy() > 0
    logs = {}
    for name, source in (("log2_counts", "total_counts"), ("log2_features", "detected_source_features")):
        values = metrics[source].to_numpy(dtype=float)
        logs[name] = np.log2(values, out=np.full(len(values), np.nan), where=values > 0)
    logs["pct_mt"] = metrics["pct_mt"].to_numpy()
    flags, thresholds, problems = pd.DataFrame(index=metrics.index), [], []
    for metric, values in logs.items():
        usable = values[positive & np.isfinite(values)]
        median = np.median(usable) if len(usable) else np.nan
        mad = 1.4826 * np.median(np.abs(usable - median)) if len(usable) else np.nan
        stable = len(usable) >= 50 and np.isfinite(mad) and mad > 0
        if not stable:
            problems.append(f"{metric}: fewer than 50 valid profiles or zero/undefined MAD; manual threshold review")
        for n in MAD_MULTIPLIERS:
            lower, upper = median - n * mad, median + n * mad
            if metric != "pct_mt":
                flags[f"low_{metric}_{n}mad"] = positive & np.isfinite(values) & (values < lower) if np.isfinite(mad) and mad > 0 else False
            flags[f"high_{metric}_{n}mad"] = positive & np.isfinite(values) & (values > upper) if np.isfinite(mad) and mad > 0 else False
            thresholds.append({"metric": metric, "n_mads": n, "valid_profiles": len(usable),
                               "median": finite_number(median), "scaled_mad": finite_number(mad),
                               "lower": finite_number(lower) if metric != "pct_mt" and mad > 0 else None,
                               "upper": finite_number(upper) if mad > 0 else None,
                               "stratum_requires_manual_review": not stable})
    for n in MAD_MULTIPLIERS:
        columns = [c for c in flags if c.endswith(f"_{n}mad")]
        flags[f"any_{n}mad"] = flags[columns].any(axis=1)
        flags[f"joint_low_complexity_{n}mad"] = flags[f"low_log2_counts_{n}mad"] & flags[f"low_log2_features_{n}mad"]
    return flags, thresholds, problems


def prior_score(scores, prior, baseline=0.05):
    """Change Scrublet's rho with fixed neighbors: posterior odds scale by prior odds.

    The core score is q*rho/r / [1-rho-q*(1-rho-rho/r)]. Its odds are
    (q/[r*(1-q)]) * rho/(1-rho), so the likelihood factor cancels.
    This is conditional-prior sensitivity, not another fitted model.
    """
    factor = (prior / (1 - prior)) / (baseline / (1 - baseline))
    scores = np.asarray(scores, dtype=float)
    return factor * scores / (1 - scores + factor * scores)


def scrublet_review(counts, names, output):
    result = pd.DataFrame(index=names)
    for prior in PRIORS:
        tag = f"{prior:.2f}"
        result[f"doublet_score_prior_{tag}"] = np.nan
        result[f"doublet_flag_prior_{tag}"] = pd.array([pd.NA] * len(names), dtype="boolean")
    result["doublet_diagnostic_status"] = "not_scored"
    valid = np.flatnonzero(np.asarray((counts > 0).sum(axis=1)).ravel() >= 3)
    report = {"software_prior": 0.05, "study_doublet_rate": "unknown", "sensitivity_priors": list(PRIORS),
              "sim_doublet_ratio": 2, "seed": SEED, "mean_center": False, "normalize_variance": True,
              "embedding": "sparse TruncatedSVD; avoids default dense centered PCA",
              "use_approx_neighbors": False, "knn_dist_metric": "euclidean",
              "neighbor_method": "Exact sklearn nearest neighbors; explicit official Scrublet option avoids the automatic approximate-neighbor compatibility path",
              "threshold_method": "official Scrublet/skimage threshold_minimum; provisional human review",
              "diagnostic_input_profiles": len(valid), "excluded_from_algorithm_only": len(names) - len(valid),
              "sensitivity": "One fitted neighbor graph; Bayesian prior-odds score transformation at fixed baseline numeric threshold, not three refits.",
              "warnings": [], "status": "not_scored", "threshold": None}
    simulated = np.array([], dtype=float)
    if len(valid) < 50:
        report["status"] = "too_few_profiles_for_stable_doublet_review"
    else:
        work = ad.AnnData(X=counts[valid].copy(), obs=pd.DataFrame(index=names[valid]))
        n_features = int(np.asarray((work.X > 0).sum(axis=0)).ravel().astype(bool).sum())
        pcs = min(30, len(valid) - 2, n_features - 2)
        report["n_prin_comps"] = pcs
        auto_threshold = util.find_spec("skimage") is not None
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                sc.pp.scrublet(work, sim_doublet_ratio=2, expected_doublet_rate=0.05,
                               stdev_doublet_rate=0.02, n_prin_comps=pcs, mean_center=False,
                               normalize_variance=True, random_state=SEED,
                               use_approx_neighbors=False, knn_dist_metric="euclidean",
                               threshold=None if auto_threshold else 1.0, verbose=False)
            report["warnings"] = [str(w.message) for w in caught]
            base = work.obs["doublet_score"].to_numpy(dtype=float)
            simulated = np.asarray(work.uns["scrublet"]["doublet_scores_sim"], dtype=float)
            threshold = work.uns["scrublet"].get("threshold") if auto_threshold else None
            if threshold is not None and np.isfinite(threshold):
                report["threshold"] = float(threshold)
            report["status"] = "scored_threshold_requires_review" if report["threshold"] is not None else "scored_threshold_unresolved"
            result.loc[work.obs_names, "doublet_diagnostic_status"] = report["status"]
            for prior in PRIORS:
                tag = f"{prior:.2f}"
                transformed = prior_score(base, prior)
                result.loc[work.obs_names, f"doublet_score_prior_{tag}"] = transformed
                if report["threshold"] is not None:
                    result.loc[work.obs_names, f"doublet_flag_prior_{tag}"] = transformed > report["threshold"]
            if not auto_threshold:
                report["warnings"].append("skimage missing: threshold=1.0 used only to obtain scores; those temporary calls are discarded, final threshold/calls missing")
            pd.DataFrame({f"score_prior_{p:.2f}": prior_score(simulated, p) for p in PRIORS}).to_csv(
                output / "scrublet_simulated_scores.tsv.gz", sep="\t", index=False)
        except Exception as exc:
            report["status"] = "algorithm_failed_requires_review"
            report["error"] = repr(exc)
            report["traceback"] = traceback.format_exc()
        del work
        gc.collect()
    report["flag_counts"] = {f"{p:.2f}": int(result[f"doublet_flag_prior_{p:.2f}"].fillna(False).sum()) for p in PRIORS}
    report["unscored_profiles"] = int(result["doublet_score_prior_0.05"].isna().sum())
    result["doublet_score_available"] = result["doublet_score_prior_0.05"].notna()
    figure, axes = plt.subplots(1, 3, figsize=(12, 3.2))
    for axis, prior in zip(axes, PRIORS):
        observed = result[f"doublet_score_prior_{prior:.2f}"].dropna().to_numpy()
        if len(observed):
            axis.hist(observed, bins=75, density=True, alpha=.6, label="Observed")
            axis.hist(prior_score(simulated, prior), bins=75, density=True, histtype="step", label="Simulated")
            if report["threshold"] is not None:
                axis.axvline(report["threshold"], color="red", linestyle="--", label="Fixed baseline threshold")
            axis.legend(fontsize=6)
        else:
            axis.text(.05, .5, report["status"], wrap=True, fontsize=8)
        axis.set_title(f"Diagnostic prior {prior:.2f}")
        axis.set_xlabel("Scrublet score")
    figure.tight_layout()
    figure.savefig(output / "scrublet_scores.png", dpi=140)
    plt.close(figure)
    write_json(output / "scrublet.json", report)
    return result, report


def marker_evidence(counts, symbols, totals):
    rows = {marker: np.flatnonzero(symbols == marker) for marker in MARKERS}
    selected = np.array(sorted({int(i) for indices in rows.values() for i in indices}), dtype=int)
    block = counts[:, selected]
    # Only one dense vector per fixed marker; never densify cells x all genes.
    values = np.zeros((counts.shape[0], len(MARKERS)), dtype=np.float32)
    for j, marker in enumerate(MARKERS):
        if len(rows[marker]):
            values[:, j] = np.asarray(block[:, np.searchsorted(selected, rows[marker])].sum(axis=1)).ravel()
    detected = values > 0
    scale = np.divide(1e4, totals, out=np.zeros(len(totals)), where=totals > 0)
    logged = np.log1p(values * scale[:, None]).astype(np.float32)
    available = np.array([bool(len(rows[m])) for m in MARKERS])
    return logged, detected, available, rows


def provisional_identity(logged, detected, available):
    names = list(PROGRAMS)
    score = np.zeros((len(logged), len(names)), dtype=np.float32)
    support = np.zeros_like(score, dtype=np.int16)
    anchors = np.zeros_like(score, dtype=bool)
    coverage = []
    for j, (name, (broad, genes, anchor)) in enumerate(PROGRAMS.items()):
        ids = [MARKERS.index(g) for g in genes.split() if available[MARKERS.index(g)]]
        aids = [MARKERS.index(g) for g in anchor.split() if available[MARKERS.index(g)]]
        if ids:
            score[:, j] = logged[:, ids].mean(axis=1)
            support[:, j] = detected[:, ids].sum(axis=1)
        if aids:
            anchors[:, j] = detected[:, aids].any(axis=1)
        coverage.append({"program": name, "broad": broad, "seed_genes": genes,
                         "available_markers": len(ids), "seed_markers": len(genes.split()),
                         "available_anchors": len(aids)})
    guarded = score * anchors
    broad_names = ["epithelial", "stromal", "immune", "endothelial", "other"]
    broad_scores = np.column_stack([guarded[:, [j for j, v in enumerate(PROGRAMS.values()) if v[0] == b]].max(axis=1) for b in broad_names])
    order = np.argsort(broad_scores, axis=1)
    top = order[:, -1]
    second = order[:, -2]
    indices = np.arange(len(logged))
    top_value = broad_scores[indices, top]
    runner = broad_scores[indices, second]
    ratio = np.divide(runner, top_value, out=np.zeros(len(logged)), where=top_value > 0)
    conflict = (top_value > 0) & (runner > 0) & (ratio >= 0.8)
    best = guarded.argmax(axis=1)
    best_support = support[indices, best]
    broad = np.array([broad_names[i] for i in top], dtype=object)
    broad[(top_value == 0) | conflict] = "unknown"
    candidate = np.array([names[i] for i in best], dtype=object)
    candidate[top_value == 0] = "unknown"
    core = candidate.copy()
    core[(best_support < 2) | conflict | (top_value == 0)] = "unresolved"
    evidence = np.full(len(logged), "provisional_multi_seed", dtype=object)
    evidence[best_support < 2] = "weak_single_seed_review"
    evidence[conflict] = "competing_broad_programs_review"
    evidence[top_value == 0] = "insufficient_seed_evidence"
    result = pd.DataFrame({"provisional_broad_identity": broad, "provisional_core_identity": core,
                           "marker_candidate": candidate, "marker_support_n": best_support,
                           "broad_competitor_ratio": ratio, "identity_evidence": evidence,
                           "identity_status": "provisional_requires_human_review"})
    # The 0.8 ratio is only a transparent review-queue heuristic, never confidence calibration.
    return result, coverage


def source_comparison(obs, broad):
    columns = [c for c in ("source_broad_annotation", "source_coarse_cell_type", "source_population", "source_cellType", "source_cell_type") if c in obs]
    mapping = {"epithelial": "epithelial", "mesenchymal": "stromal", "stroma": "stromal", "stromal": "stromal",
               "endothelial": "endothelial", "immune": "immune"}
    if not columns:
        return "", pd.Series("source_broad_annotation_unavailable", index=obs.index)
    column = columns[0]
    source = obs[column].astype(str).str.lower().map(mapping)
    compare = np.full(len(obs), "not_comparable_or_unknown", dtype=object)
    valid = source.notna().to_numpy() & (broad.to_numpy() != "unknown")
    compare[valid] = np.where(source.to_numpy()[valid] == broad.to_numpy()[valid], "agrees", "discordant_review")
    return column, pd.Series(compare, index=obs.index)


def marker_groups(labels, logged, detected, available, scope):
    rows = []
    for group in pd.unique(labels.astype(str)):
        selection = labels.astype(str).to_numpy() == group
        for j, marker in enumerate(MARKERS):
            rows.append({"scope": scope, "group": group, "profiles": int(selection.sum()), "marker": marker,
                         "measurement": "exported_source_feature" if available[j] else "not_in_source_export_coverage_unknown",
                         "fraction_nonzero": float(detected[selection, j].mean()) if available[j] else None,
                         "mean_log1p_cp10k": float(logged[selection, j].mean()) if available[j] else None})
    return rows


def marker_plot(table, output, scope="provisional_core_identity", filename="identity_marker_review.png"):
    frame = pd.DataFrame(table)
    frame = frame[frame["scope"] == scope]
    mean = frame.pivot(index="group", columns="marker", values="mean_log1p_cp10k").reindex(columns=MARKERS)
    fraction = frame.pivot(index="group", columns="marker", values="fraction_nonzero").reindex(index=mean.index, columns=MARKERS)
    figure, axis = plt.subplots(figsize=(16, max(3, len(mean) * .38 + 1.6)))
    x, y = np.meshgrid(np.arange(len(MARKERS)), np.arange(len(mean)))
    points = axis.scatter(x.ravel(), y.ravel(), s=110 * fraction.fillna(0).to_numpy().ravel(),
                          c=mean.fillna(0).to_numpy().ravel(), cmap="viridis", vmin=0)
    axis.set_xticks(np.arange(len(MARKERS)), MARKERS, rotation=90, fontsize=6)
    axis.set_yticks(np.arange(len(mean)), mean.index, fontsize=7)
    axis.invert_yaxis()
    axis.set_title(f"{scope}: size=detection fraction; color=mean log1p(CP10K). Feature availability is in TSV.", fontsize=9)
    figure.colorbar(points, ax=axis, shrink=.6)
    figure.tight_layout()
    figure.savefig(output / filename, dpi=150)
    plt.close(figure)


def qc_plot(metrics, flags, output):
    figure, axes = plt.subplots(2, 3, figsize=(11, 6))
    for axis, column in zip(axes[0], ("total_counts", "detected_source_features", "pct_mt")):
        values = metrics[column].to_numpy()
        axis.hist(values[np.isfinite(values)], bins=70)
        axis.set_title(column, fontsize=9)
    rng = np.random.default_rng(SEED)
    selection = np.arange(len(metrics))
    if len(selection) > 12000:
        selection = np.sort(rng.choice(selection, 12000, replace=False))
    colors = np.where(flags["any_3mad"].to_numpy()[selection], "#c44e52", "#4c72b0")
    for axis, a, b in ((axes[1, 0], "total_counts", "detected_source_features"),
                       (axes[1, 1], "total_counts", "pct_mt"),
                       (axes[1, 2], "detected_source_features", "pct_mt")):
        axis.scatter(metrics[a].to_numpy()[selection], metrics[b].to_numpy()[selection], s=2, c=colors, alpha=.45, rasterized=True)
        axis.set_xscale("symlog", linthresh=1)
        if b != "pct_mt":
            axis.set_yscale("symlog", linthresh=1)
        axis.set_xlabel(a, fontsize=8)
        axis.set_ylabel(b, fontsize=8)
    figure.suptitle("Source profiles retained; red=any 3-MAD review flag, not deletion or expert failure", fontsize=10)
    figure.tight_layout()
    figure.savefig(output / "qc_review.png", dpi=150)
    plt.close(figure)


def resolve_symbols(data, record):
    field = record.get("symbol_field")
    if field is None:
        field = next((x for x in ("source_symbol", "source_gene_symbol", "gene_symbol", "source_gene_name", "gene_name") if x in data.var), None)
    if field is not None:
        if field not in data.var:
            raise ValueError(f"Symbol field absent: {field}")
        symbols = data.var[field].astype(str).to_numpy()
    else:
        symbols = data.var_names.astype(str).to_numpy()
        if np.mean(np.char.startswith(symbols.astype(str), "ENSG")) > .8:
            raise ValueError("Ensembl-only axis has no source symbol field; do not invent symbol mappings")
    return symbols, field or "var_names_source_symbols"


def process_record(record, output, script_hash):
    started = time.monotonic()
    path = Path(record["path"])
    fingerprint = sha256(path)
    summary_path = output / "summary.json"
    if summary_path.exists():
        old = json.loads(summary_path.read_text())
        if old.get("input_sha256") == fingerprint and old.get("script_sha256") == script_hash and old.get("record") == record:
            if all((output / name).is_file() and sha256(output / name) == digest for name, digest in old.get("artifacts", {}).items()):
                log("reuse_completed_capture", library=record["library_id"])
                return old
        raise ValueError(f"Existing output differs; select a new output directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    log("capture_start", resource=record["resource"], library=record["library_id"], path=str(path))
    data = ad.read_h5ad(path)
    if not sparse.issparse(data.X) or not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("Require sparse original X and unique storage axes")
    if any(str(c).startswith(PREFIX) for c in data.obs):
        raise ValueError("Source obs already has review_ fields; refusing to overwrite source metadata")
    for field, expected in (("donor_id", record["donor_id"]), ("library_id", record["library_id"])):
        if field in data.obs and set(data.obs[field].dropna().astype(str)) != {str(expected)}:
            raise ValueError(f"Record does not represent exactly one {field}: {expected}")
    for field in record.get("source_label_columns", []):
        if field not in data.obs:
            raise ValueError(f"Declared original label absent: {field}")
    counts = data.X.tocsr(copy=True)
    if counts.data.size and (not np.isfinite(counts.data).all() or counts.data.min() < 0 or not np.equal(counts.data, np.floor(counts.data)).all()):
        raise ValueError("Original UMI counts must be finite, nonnegative integers")
    counts.sum_duplicates()
    counts.eliminate_zeros()
    counts.sort_indices()
    before = matrix_hash(counts)
    symbols, symbol_field = resolve_symbols(data, record)
    if "mt" in data.var and pd.api.types.is_bool_dtype(data.var["mt"].dtype):
        mt = data.var["mt"].to_numpy(dtype=bool)
        mt_definition = "preserved source-checkpoint var.mt; see source provenance"
    else:
        mt = np.char.startswith(np.char.upper(symbols.astype(str)), "MT-")
        mt_definition = "source symbol MT- prefix fallback, not official chromosome remapping"
    metrics = qc_metrics(counts, mt)
    metrics.index = data.obs_names.copy()
    if not mt.any():
        metrics["pct_mt"] = np.nan
    flags, thresholds, problems = mad_flags(metrics)
    flags.index = data.obs_names.copy()
    if not mt.any():
        problems.append("No source mitochondrial rows recognized; MT percentage unmeasurable, not zero")
    doublets, doublet_report = scrublet_review(counts, data.obs_names, output)
    logged, detected, available, marker_rows = marker_evidence(counts, symbols, metrics["total_counts"].to_numpy())
    identity, program_coverage = provisional_identity(logged, detected, available)
    identity.index = data.obs_names.copy()
    source_column, comparison = source_comparison(data.obs, identity["provisional_broad_identity"])
    identity["source_broad_comparison"] = comparison
    core_groups = marker_groups(identity["provisional_core_identity"], logged, detected, available, "provisional_core_identity")
    group_rows = core_groups + marker_groups(identity["provisional_broad_identity"], logged, detected, available, "provisional_broad_identity")
    if source_column:
        group_rows += marker_groups(data.obs[source_column], logged, detected, available, "original_source_broad_annotation")
    for column in ("source_fine_annotation", "source_subcluster", "source_celltype"):
        if column in data.obs:
            group_rows += marker_groups(data.obs[column], logged, detected, available, f"original_{column}")
    pd.DataFrame(group_rows).to_csv(output / "group_marker_evidence.tsv.gz", sep="\t", index=False)
    pd.DataFrame(program_coverage).to_csv(output / "seed_program_coverage.tsv", sep="\t", index=False)
    pd.DataFrame([{"marker": m, "source_feature_rows": len(marker_rows[m]),
                   "available_in_source_axis": bool(available[j]), "duplicate_symbol_detection": "union; count intensity sums only for review, X unchanged"}
                  for j, m in enumerate(MARKERS)]).to_csv(output / "marker_axis_coverage.tsv", sep="\t", index=False)
    marker_plot(core_groups, output)
    if source_column:
        marker_plot(group_rows, output, "original_source_broad_annotation", "source_identity_marker_review.png")
    qc_plot(metrics, flags, output)
    obs = data.obs.copy()
    additions = pd.concat([metrics, flags, doublets, identity], axis=1)
    for key, value in (("resource", record["resource"]), ("donor_uid", f"{record['resource']}::{record['donor_id']}"),
                       ("library_id", record["library_id"]), ("assay", record["assay"]), ("fraction", record.get("fraction"))):
        additions[key] = "" if value is None else value
    additions["automated_disposition"] = np.where(metrics["zero_counts"], "HARD_FAIL_ZERO_COUNTS_RETAINED", "REVIEW_PENDING")
    additions["expert_qc_status"] = "not_adjudicated"
    additions["expert_identity_status"] = "not_adjudicated"
    additions["endpoint_eligibility"] = "not_adjudicated"
    for column in additions:
        obs[PREFIX + column] = additions[column]
    obs.to_csv(output / "obs_review.tsv.gz", sep="\t", index_label="source_obs_name")
    data.var.to_csv(output / "source_feature_axis.tsv.gz", sep="\t", index_label="source_var_name")
    restored = pd.read_csv(output / "obs_review.tsv.gz", sep="\t", usecols=[0], dtype=str, keep_default_na=False)
    if restored.iloc[:, 0].tolist() != data.obs_names.astype(str).tolist():
        raise ValueError("Exported obs axis does not exactly reproduce source order")
    if before != matrix_hash(counts) or fingerprint != sha256(path):
        raise ValueError("Source file or review count matrix changed during processing")
    source_fields = [str(c) for c in data.obs if "source" in str(c).lower() or "annot" in str(c).lower()]
    counts_by_identity = identity.groupby(["provisional_broad_identity", "provisional_core_identity"], dropna=False).size()
    coverage = []
    no_flags = ~flags["any_3mad"] & ~metrics["zero_counts"] & ~doublets["doublet_flag_prior_0.05"].fillna(False)
    scored = doublets["doublet_score_available"]
    for (broad, core), n in counts_by_identity.items():
        selected = (identity["provisional_broad_identity"] == broad) & (identity["provisional_core_identity"] == core)
        coverage.append({"resource": record["resource"], "donor_id": record["donor_id"],
                         "donor_uid": f"{record['resource']}::{record['donor_id']}", "library_id": record["library_id"],
                         "assay": record["assay"], "fraction": record.get("fraction"), "broad": broad, "core": core,
                         "all_retained_profiles": int(n), "no_automated_3mad_or_doublet_flag": int((selected & no_flags).sum()),
                         "doublet_unscored_profiles": int((selected & ~scored).sum()),
                         "no_3mad_flags_and_scored_doublet_negative": int((selected & no_flags & scored).sum()),
                         "interpretation": "Provisional review counts, not final identity/eligibility; unscored doublets remain disclosed"})
    pd.DataFrame(coverage).to_csv(output / "identity_coverage.tsv", sep="\t", index=False)
    summary = {"record": record, "status": "review_pack_generated_not_adjudicated", "boundary": BOUNDARY,
               "input_sha256": fingerprint, "script_sha256": script_hash,
               "source_shape": list(data.shape), "all_profiles_retained": len(obs) == data.n_obs,
               "source_counts_file_unchanged": True, "canonical_review_counts_sha256": before,
               "obs_axis_sha256": axis_hash(data.obs_names), "var_axis_sha256": axis_hash(data.var_names),
               "source_obs_fields_preserved": list(data.obs.columns), "source_annotation_fields": source_fields,
               "symbol_field": symbol_field, "mitochondrial_definition": mt_definition, "mt_feature_rows": int(mt.sum()),
               "total_counts": int(metrics["total_counts"].sum()), "hard_fail_zero_profiles": int(metrics["zero_counts"].sum()),
               "mad_review_flags": {str(n): int(flags[f"any_{n}mad"].sum()) for n in MAD_MULTIPLIERS},
               "threshold_review_problems": problems, "qc_stratum": "single input capture/assay/fraction; physical capture provenance may remain unresolved",
               "scrublet": doublet_report, "capture_identity_status": record.get("capture_identity_status", "unknown"),
               "doublet_capture_limit": "Diagnostics run per submitted capture/library proxy. Unresolved physical captures and sorted-out singlet parents limit interpretation; scores do not establish doublet rates or final exclusion.",
               "provisional_broad_counts": {str(k): int(v) for k, v in identity["provisional_broad_identity"].value_counts().items()},
               "source_comparison_field": source_column, "source_comparison_counts": {str(k): int(v) for k, v in comparison.value_counts().items()},
               "marker_scoring": "mean log1p(CP10K) over available fixed seeds with anchor guard; maximum program per broad identity; competitor/top >=0.8 queues conflict; weak single-seed profiles retain candidate but unresolved core",
               "marker_threshold_boundary": "0.8 competition and support counts are transparent provisional review heuristics, not validated identity cutoffs or expert assignments",
               "ambient_rna": "Original-count branch retained; filtered/uncertain background prohibits forced correction. Off-lineage seed coexpression is review evidence, not a contamination estimate.",
               "elapsed_seconds": time.monotonic() - started,
               "versions": {name: metadata.version(name) for name in ("scanpy", "anndata", "numpy", "scipy", "pandas", "matplotlib")}}
    pd.DataFrame(thresholds).to_csv(output / "mad_thresholds.tsv", sep="\t", index=False)
    summary["artifacts"] = {p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file() and p.name != "summary.json"}
    write_json(summary_path, summary)
    log("capture_complete", library=record["library_id"], profiles=data.n_obs,
        doublet_status=doublet_report["status"], seconds=round(summary["elapsed_seconds"], 2))
    del data, counts, obs, logged, detected
    gc.collect()
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.inputs.read_text())
    records = manifest["records"] if isinstance(manifest, dict) else manifest
    if not records:
        raise ValueError("No capture records")
    seen, paths = set(), set()
    for record in records:
        for field in ("path", "resource", "library_id", "donor_id", "assay"):
            if not record.get(field):
                raise ValueError(f"Missing record field {field}")
        if record["assay"] not in ("scRNA", "snRNA"):
            raise ValueError("Assay must be scRNA or snRNA")
        path = Path(record["path"])
        record["path"] = str((args.inputs.parent / path).resolve(strict=True)) if not path.is_absolute() else str(path.resolve(strict=True))
        uid = f"{record['resource']}::{record['library_id']}"
        if uid in seen or record["path"] in paths:
            raise ValueError("Duplicate capture identity or reused checkpoint; supply one physical/proxy capture per file")
        seen.add(uid)
        paths.add(record["path"])
    args.output.mkdir(parents=True, exist_ok=True)
    script_hash = sha256(Path(__file__))
    write_json(args.output / "seed_dictionary.json", {"programs": PROGRAMS, "markers": MARKERS,
               "boundary": "Fixed review priors from QC_ANNOTATION_STANDARD; not trained on reserved cohorts; no disease-state signatures"})
    summaries, failures, coverage = [], [], []
    for record in records:
        uid = f"{record['resource']}::{record['library_id']}"
        folder = re.sub(r"[^A-Za-z0-9_.-]", "_", uid) + "_" + hashlib.sha256(uid.encode()).hexdigest()[:8]
        output = args.output / folder
        try:
            summary = process_record(record, output, script_hash)
            summaries.append(summary)
            coverage.append(pd.read_csv(output / "identity_coverage.tsv", sep="\t", keep_default_na=False))
        except Exception as exc:
            failure = {"record": record, "status": "quarantined_processing_error", "error": repr(exc), "traceback": traceback.format_exc()}
            failures.append(failure)
            output.mkdir(parents=True, exist_ok=True)
            write_json(output / "ERROR.json", failure)
            log("capture_error", library=record["library_id"], error=repr(exc))
        write_json(args.output / "progress.json", {"completed_captures": len(summaries), "failed_captures": failures, "expected_captures": len(records)})
    if coverage:
        all_coverage = pd.concat(coverage, ignore_index=True)
        all_coverage.to_csv(args.output / "library_identity_coverage.tsv", sep="\t", index=False)
        keys = ["resource", "donor_uid", "donor_id", "assay", "fraction", "broad", "core"]
        donor = all_coverage.groupby(keys, dropna=False, observed=True)[["all_retained_profiles", "no_automated_3mad_or_doublet_flag", "doublet_unscored_profiles", "no_3mad_flags_and_scored_doublet_negative"]].sum().reset_index()
        donor["interpretation"] = "Same source donor/assay/fraction review sums; specimen/region compatibility and expert QC still pending, no pseudobulk produced"
        for n in (10, 20, 50):
            donor[f"review_profiles_at_least_{n}"] = donor["all_retained_profiles"] >= n
        donor.to_csv(args.output / "donor_identity_coverage.tsv", sep="\t", index=False)
    library_rows = []
    for summary in summaries:
        record = summary["record"]
        library_rows.append({**{k: record.get(k) for k in ("resource", "donor_id", "library_id", "assay", "fraction")},
                             "profiles_retained": summary["source_shape"][0], "source_feature_rows": summary["source_shape"][1],
                             "zero_hard_fail_retained": summary["hard_fail_zero_profiles"],
                             **{f"review_{n}mad": summary["mad_review_flags"][str(n)] for n in MAD_MULTIPLIERS},
                             "scrublet_status": summary["scrublet"]["status"],
                             "source_counts_unchanged": summary["source_counts_file_unchanged"]})
    pd.DataFrame(library_rows).to_csv(args.output / "library_review_summary.tsv", sep="\t", index=False)
    result = {"status": "review_pack_generated" if not failures else "incomplete_quarantined_captures",
              "boundary": BOUNDARY, "inputs_sha256": sha256(args.inputs), "script_sha256": script_hash,
              "expected_captures": len(records), "completed_captures": len(summaries), "failures": failures,
              "retained_profiles": sum(s["source_shape"][0] for s in summaries),
              "expert_qc_adjudicated": False, "expert_identity_adjudicated": False,
              "checkpoint_A_passed": False, "formal_discovery_started": False,
              "capture_summaries": summaries}
    write_json(args.output / "summary.json", result)
    lines = ["# Full-cohort RNA QC and identity review pack", "", BOUNDARY, "",
             f"Completed captures: {len(summaries)}/{len(records)}. Retained source profiles: {result['retained_profiles']:,}.", "",
             "Review library_review_summary.tsv and each capture's qc_review.png, mad_thresholds.tsv and scrublet_scores.png.",
             "Use identity_marker_review.png with group_marker_evidence.tsv.gz, original labels and obs_review.tsv.gz; identity hints are not expert labels.",
             "Where source identities exist, compare source_identity_marker_review.png against the provisional groups. Missing feature rows are coverage-unknown exports, not measured biological zero.",
             "No profile has an expert PASS. No 200-gene/25%-MT filter, doublet quota or ambient correction was applied.",
             "Scrublet 0.02/0.10 sensitivity changes prior odds at the fixed 0.05 numeric threshold with the same neighbors; inspect uncertainty, not just flags.",
             "Provisional group/donor counts do not establish an estimable disease contrast; final identity, region compatibility, QC and independence remain review gates.",
             "Source H5ADs and their complete axes remain unchanged; exports retain every original obs column and profile.", ""]
    (args.output / "REVIEW_README.md").write_text("\n".join(lines))
    log("review_pack_complete", completed=len(summaries), failed=len(failures), profiles=result["retained_profiles"])
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
