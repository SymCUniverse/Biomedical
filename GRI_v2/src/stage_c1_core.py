from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd


def sha256_file(path: Path, chunk: int = 16 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


def stable_seed(namespace: str, cancer: str, track: str, stratum: str,
                resample: int, null_type: str, replicate: int = 0) -> int:
    payload = f"{namespace}|{cancer}|{track}|{stratum}|{int(resample)}|{null_type}|{int(replicate)}".encode('utf-8')
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], 'big', signed=False) % (2**32)


def exact_sign_test_positive(values: Sequence[float]) -> tuple[int, int, int, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    pos = int(np.sum(vals > 0))
    neg = int(np.sum(vals < 0))
    ties = int(np.sum(vals == 0))
    n = pos + neg
    if n == 0:
        return pos, neg, ties, float('nan')
    p = sum(math.comb(n, k) for k in range(pos, n + 1)) / (2 ** n)
    return pos, neg, ties, float(p)


def bh_fdr(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    good = np.flatnonzero(np.isfinite(p))
    if len(good) == 0:
        return out
    vals = p[good]
    order = np.argsort(vals, kind='mergesort')
    ranked = vals[order]
    m = len(ranked)
    q = ranked * m / np.arange(1, m + 1, dtype=float)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    tmp = np.empty_like(q)
    tmp[order] = q
    out[good] = tmp
    return out


def parse_gmt(path: Path) -> dict[str, list[str]]:
    modules: dict[str, list[str]] = {}
    with Path(path).open('r', encoding='utf-8', errors='strict') as fh:
        for line in fh:
            parts = line.rstrip('\r\n').split('\t')
            if len(parts) >= 3:
                name = parts[0].strip()
                genes = [g.strip() for g in parts[2:] if g.strip()]
                if name and genes:
                    modules[name] = list(dict.fromkeys(genes))
    if len(modules) != 50 or any(not x.startswith('HALLMARK_') for x in modules):
        raise ValueError(f'expected exactly 50 HALLMARK modules; found {len(modules)}')
    return modules


def _rank_columns(x: np.ndarray) -> np.ndarray:
    a = np.asarray(x, dtype=float)
    if a.ndim == 1:
        a = a[:, None]
    out = np.full_like(a, np.nan, dtype=float)
    for j in range(a.shape[1]):
        s = pd.Series(a[:, j], dtype=float)
        good = s.notna()
        if int(good.sum()) >= 3:
            r = s[good].rank(method='average').to_numpy(dtype=float, copy=True)
            r -= r.mean()
            den = np.sqrt(np.sum(r * r))
            if den > 0:
                out[np.flatnonzero(good.to_numpy()), j] = r / den
    return out


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    a = np.asarray(x, float)
    b = np.asarray(y, float)
    good = np.isfinite(a) & np.isfinite(b)
    if int(good.sum()) < 3:
        return float('nan')
    ra = pd.Series(a[good]).rank(method='average').to_numpy(dtype=float)
    rb = pd.Series(b[good]).rank(method='average').to_numpy(dtype=float)
    if np.std(ra, ddof=1) <= 0 or np.std(rb, ddof=1) <= 0:
        return float('nan')
    return float(np.corrcoef(ra, rb)[0, 1])


def abs_spearman_columns(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if a.shape != b.shape:
        raise ValueError('paired Hallmark matrices must have identical shape')
    vals = []
    signed = []
    for j in range(a.shape[1]):
        r = spearman(a[:, j], b[:, j])
        signed.append(r)
        vals.append(abs(r) if np.isfinite(r) else np.nan)
    return np.asarray(vals, float), np.asarray(signed, float)


def same_hallmark_stat(a: np.ndarray, b: np.ndarray) -> tuple[float, np.ndarray]:
    vals, signed = abs_spearman_columns(a, b)
    good = np.isfinite(vals)
    if not np.any(good):
        return float('nan'), signed
    return float(np.median(vals[good])), signed


def modal_analysis(beta: np.ndarray, top_modes: int = 5, contribution_dtype=np.float32) -> dict:
    b = np.asarray(beta, dtype=float)
    if b.ndim != 2 or b.shape[0] != 30:
        raise ValueError(f'modal beta matrix must be 30 x p, got {b.shape}')
    if not np.isfinite(b).all():
        raise ValueError('modal beta matrix contains non-finite values')
    p = b.shape[1]
    if p < 2:
        raise ValueError('modal analysis requires at least two probes')
    x = b - b.mean(axis=0, keepdims=True)
    g = (x @ x.T) / float(p)
    vals, vecs = np.linalg.eigh(g)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    scale = float(np.sum(np.abs(vals)))
    tol = 1e-10 * scale
    if np.any(vals < -tol):
        raise ValueError(f'modal eigenvalue negative beyond tolerance: min={vals.min()} tol={tol}')
    vals = np.where(vals < 0, 0.0, vals)
    vals29 = vals[:29]
    total = float(vals29.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError('zero total modal variance')
    q = vals29 / total
    nz = q > 0
    h = float(-np.sum(q[nz] * np.log(q[nz])))
    hnorm = h / math.log(29.0)
    reff = math.exp(h)
    rpr = float(1.0 / np.sum(q * q))
    s_spec = float(1.0 - hnorm)
    s_pr = float(1.0 - rpr / 29.0)
    tm = min(int(top_modes), 5)
    top_u = vecs[:, :tm].copy()
    for k in range(tm):
        j = int(np.argmax(np.abs(top_u[:, k])))
        if top_u[j, k] < 0:
            top_u[:, k] *= -1.0
    contrib = np.empty((tm, p), dtype=contribution_dtype)
    for k in range(tm):
        lam = float(vals[k])
        if lam <= 0:
            contrib[k] = 0
            continue
        c = ((x.T @ top_u[:, k]) ** 2) / (float(p) * lam)
        s = float(np.sum(c))
        if not np.isfinite(s) or abs(s - 1.0) > 1e-6:
            raise ValueError(f'mode contribution sum drift: {s}')
        contrib[k] = c.astype(contribution_dtype, copy=False)
    return {
        'eigenvalues_29': vals29,
        'q_29': q,
        'q1': float(q[0]),
        'h_norm': float(hnorm),
        'effective_rank': float(reff),
        'participation_ratio': float(rpr),
        's_spec': s_spec,
        's_pr': s_pr,
        'top_u': top_u,
        'contributions': contrib,
        'centered': x,
    }


def marginal_permute_columns(beta: np.ndarray, seed: int) -> np.ndarray:
    b = np.asarray(beta, dtype=float)
    rng = np.random.default_rng(int(seed))
    keys = rng.random(b.shape)
    order = np.argsort(keys, axis=0, kind='mergesort')
    cols = np.arange(b.shape[1])[None, :]
    return b[order, cols]


def linear_cka_from_centered(x: np.ndarray, y: np.ndarray) -> float:
    a = np.asarray(x, float)
    b = np.asarray(y, float)
    if a.shape[0] != b.shape[0]:
        raise ValueError('CKA sample mismatch')
    ka = a @ a.T
    kb = b @ b.T
    n = a.shape[0]
    h = np.eye(n) - np.ones((n, n), dtype=float) / float(n)
    kac = h @ ka @ h
    kbc = h @ kb @ h
    den = float(np.linalg.norm(kac, 'fro') * np.linalg.norm(kbc, 'fro'))
    if den <= 0 or not np.isfinite(den):
        return float('nan')
    return float(np.sum(kac * kbc) / den)


def top_principal_angles(x: np.ndarray, y: np.ndarray, k: int = 5) -> np.ndarray:
    a = np.asarray(x, float)
    b = np.asarray(y, float)
    ua, _, _ = np.linalg.svd(a, full_matrices=False)
    ub, _, _ = np.linalg.svd(b, full_matrices=False)
    kk = min(int(k), ua.shape[1], ub.shape[1])
    s = np.linalg.svd(ua[:, :kk].T @ ub[:, :kk], compute_uv=False)
    s = np.clip(s, 0.0, 1.0)
    return np.arccos(s)


def methylation_pc1(gene_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    m = np.asarray(gene_matrix, dtype=float)
    if m.ndim != 2 or m.shape[0] != 30 or not np.isfinite(m).all():
        raise ValueError('methylation Hallmark matrix invalid')
    xc = m - m.mean(axis=0, keepdims=True)
    if np.allclose(xc, 0):
        raise ValueError('methylation Hallmark has zero variation')
    u, s, vt = np.linalg.svd(xc, full_matrices=False)
    eig = u[:, 0] * s[0]
    load = vt[0].copy()
    mean_state = np.mean(m, axis=1)
    if np.std(eig, ddof=1) > 0 and np.std(mean_state, ddof=1) > 0:
        pr = float(np.corrcoef(eig, mean_state)[0, 1])
        if pr < 0:
            eig = -eig; load = -load
        elif pr == 0:
            j = int(np.argmax(np.abs(load)))
            if load[j] < 0:
                eig = -eig; load = -load
    else:
        j = int(np.argmax(np.abs(load)))
        if load[j] < 0:
            eig = -eig; load = -load
    return eig, load


def rna_hallmark_pc1(raw: np.ndarray, minimum_genes: int = 15) -> tuple[np.ndarray, int, float]:
    x = np.asarray(raw, dtype=float)
    if x.ndim != 2 or x.shape[0] != 30:
        raise ValueError('RNA Hallmark matrix must be 30 x genes')
    finite = np.isfinite(x)
    counts = finite.sum(axis=0)
    required = max(20, int(math.ceil(0.95 * x.shape[0])))
    sd = np.nanstd(x, axis=0, ddof=1)
    valid = (counts >= required) & np.isfinite(sd) & (sd > 0)
    if int(valid.sum()) < int(minimum_genes):
        raise ValueError(f'RNA Hallmark has fewer than {minimum_genes} eligible genes')
    xv = x[:, valid]
    mean = np.nanmean(xv, axis=0)
    sdv = np.nanstd(xv, axis=0, ddof=1)
    z = (xv - mean) / sdv
    imputed = ~np.isfinite(z)
    z = np.where(np.isfinite(z), z, 0.0)
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    eig = u[:, 0] * s[0]
    mean_state = np.mean(z, axis=1)
    if np.std(eig, ddof=1) > 0 and np.std(mean_state, ddof=1) > 0:
        if float(np.corrcoef(eig, mean_state)[0, 1]) < 0:
            eig = -eig
    return eig, int(valid.sum()), float(imputed.sum() / imputed.size if imputed.size else 0.0)


def residual_project(x: np.ndarray, purity: np.ndarray, leukocyte: np.ndarray) -> np.ndarray:
    a = np.asarray(x, dtype=float)
    p = np.asarray(purity, dtype=float)
    l = np.asarray(leukocyte, dtype=float)
    if len(p) != a.shape[0] or len(l) != a.shape[0] or not (np.isfinite(p).all() and np.isfinite(l).all()):
        raise ValueError('composition covariates incomplete')
    z = np.column_stack([np.ones(a.shape[0]), p, l])
    rz = np.eye(a.shape[0]) - z @ np.linalg.pinv(z.T @ z) @ z.T
    return rz @ a


def topology_similarity(meth: np.ndarray, rna: np.ndarray) -> float:
    a = np.asarray(meth, float)
    b = np.asarray(rna, float)
    if a.shape != b.shape or a.shape[1] < 3:
        return float('nan')
    ra = np.empty((a.shape[1], a.shape[1]), float)
    rb = np.empty_like(ra)
    for i in range(a.shape[1]):
        for j in range(a.shape[1]):
            ra[i, j] = spearman(a[:, i], a[:, j])
            rb[i, j] = spearman(b[:, i], b[:, j])
    iu = np.triu_indices(a.shape[1], 1)
    return spearman(np.abs(ra[iu]), np.abs(rb[iu]))


def build_covariates(patient_split: pd.DataFrame, purity_path: Path, leuk_path: Path) -> pd.DataFrame:
    stage = patient_split[['cancer_type', 'participant_root', 'stage_a_sample_id']].copy()
    stage['sample_root'] = stage['stage_a_sample_id'].astype(str).str.extract(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z])', expand=False)
    a = pd.read_csv(purity_path, sep='\t', dtype=str)
    a['participant_root'] = a['sample'].str.extract(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})', expand=False)
    a['sample_root'] = a['sample'].str.extract(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z])', expand=False)
    a['sample_type'] = a['sample'].str.extract(r'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-([0-9]{2})', expand=False)
    a['value'] = pd.to_numeric(a['purity'], errors='coerce')
    a = a[(a['call status'].str.lower() == 'called') & a['value'].notna() & (a['sample_type'] == '01')].copy()
    if a['participant_root'].duplicated().any():
        raise ValueError('ABSOLUTE eligible primary table not unique by participant')
    l = pd.read_csv(leuk_path, sep='\t', header=None, names=['source_cancer', 'sample', 'value'])
    l['participant_root'] = l['sample'].str.extract(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})', expand=False)
    l['sample_root'] = l['sample'].str.extract(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z])', expand=False)
    l['sample_type'] = l['sample'].str.extract(r'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-([0-9]{2})', expand=False)
    l['value'] = pd.to_numeric(l['value'], errors='coerce')
    l = l[l['value'].notna() & (l['sample_type'] == '01')].copy()
    l = l.groupby(['participant_root', 'sample_root'], as_index=False).agg(value=('value', 'median'))
    def attach(src: pd.DataFrame, label: str):
        roots = src.groupby('sample_root').size()
        if (roots > 1).any():
            raise ValueError(f'duplicate exact roots in {label}')
        root_map = src.set_index('sample_root')['value'].to_dict()
        pc = src.groupby('participant_root').size()
        unique = set(pc[pc == 1].index)
        pat_map = src[src['participant_root'].isin(unique)].set_index('participant_root')['value'].to_dict()
        vals = []
        for r in stage.itertuples(index=False):
            if r.sample_root in root_map:
                vals.append(float(root_map[r.sample_root]))
            elif r.participant_root in pat_map:
                vals.append(float(pat_map[r.participant_root]))
            else:
                vals.append(float('nan'))
        return vals
    stage['absolute_purity'] = attach(a, 'purity')
    stage['leukocyte_fraction'] = attach(l, 'leukocyte')
    stage['composition_complete_case'] = np.isfinite(stage['absolute_purity']) & np.isfinite(stage['leukocyte_fraction'])
    return stage
