#!/usr/bin/env python3
"""
🔄 Wine CycloPep — De novo cyclic peptide detection from Bruker timsTOF PASEF
═══════════════════════════════════════════════════════════════════════════════
Reads a Bruker .d folder directly, identifies cyclic peptide candidates
by their characteristic MS1/MS2 signatures, and generates energy-minimized
3D conformers for each candidate using RDKit MMFF force field.

Cyclic peptide detection strategy
──────────────────────────────────
  1. MS1 mass filter  — precursor neutral mass matches sum of n residues
                        (NO +H₂O, the defining feature of head-to-tail cyclics)
  2. Residue loss ions — [M+H]⁺ − residue_mass peaks present in MS2
                        (each ring-opening produces a diagnostic bn−1 ion)
  3. bn ion coverage  — all rotations of the candidate sequence are checked;
                        cyclic peptides show bn ions from multiple start positions
  4. Absence of y1   — no free C-terminus means no y1 immonium series
  5. Score            — weighted sum of all evidence types → confidence tier

3D structure output
───────────────────
  • RDKit ETKDGv3 distance geometry embedding
  • MMFF94 force-field minimization (50 conformers, lowest energy selected)
  • PDB file per candidate, SDF with all conformers optional
  • Compatible with PyMOL, UCSF ChimeraX, VMD

Usage
─────
  pip install alphatims rdkit

  # Full de novo search
  python wine_cyclopep.py --d path/to/sample.d

  # Only short peptides 2–6 residues (fastest)
  python wine_cyclopep.py --d path/to/sample.d --min_aa 2 --max_aa 6

  # With known sequence hints (restricts search space)
  python wine_cyclopep.py --d path/to/sample.d --hints VAAG,IAA,GGF

  # Quick test
  python wine_cyclopep.py --d path/to/sample.d --max_spectra 500

Author : Pol Giménez-Gil — ISVV, Université de Bordeaux
ORCID  : 0000-0002-7720-3733
Series : Step 1 WinePeptidome → Step 2 WineStructure → Step 3 WineSeq
         → Step 3b WineCycloPep (here)
"""

from __future__ import annotations

import argparse
import csv
import itertools
import logging
import sys
import time
from pathlib import Path
from collections import defaultdict

import numpy as np

# ── optional imports ──────────────────────────────────────────────────────────
try:
    from alphatims.bruker import TimsTOF
    HAS_ALPHATIMS = True
except ImportError:
    HAS_ALPHATIMS = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("wine_cyclopep")

# ══════════════════════════════════════════════════════════════════════════════
#  CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

PROTON = 1.007276
H2O    = 18.010565

RESIDUE_MASS: dict[str, float] = {
    "A": 71.03711,  "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
    "M": 131.04049, "F": 147.06841, "P":  97.05276, "S":  87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V":  99.06841,
}

# RDKit sidechain SMILES fragments (attached at Cα)
AA_SIDECHAIN: dict[str, str | None] = {
    "G": "",           "A": "C",          "V": "CC(C)",
    "L": "CC(C)C",     "I": "C(CC)C",     "F": "Cc1ccccc1",
    "W": "Cc1c[nH]c2ccccc12",             "M": "CCSC",
    "S": "CO",         "T": "C(O)C",      "C": "CS",
    "Y": "Cc1ccc(O)cc1",                  "H": "Cc1c[nH]cn1",
    "D": "CC(=O)O",    "E": "CCC(=O)O",   "N": "CC(=O)N",
    "Q": "CCC(=O)N",   "K": "CCCCN",      "R": "CCCNC(=N)N",
    "P": None,          # proline: special ring — skip in de novo SMILES builder
}

# Canonical sorted key for deduplication of rotational isomers
def canonical_cyclic(seq: str) -> str:
    """Smallest rotation + its reverse — gives one canonical key per cyclic sequence."""
    rotations = [seq[i:] + seq[:i] for i in range(len(seq))]
    rev = seq[::-1]
    rev_rotations = [rev[i:] + rev[:i] for i in range(len(rev))]
    return min(rotations + rev_rotations)

# ══════════════════════════════════════════════════════════════════════════════
#  MASS FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def cyclic_neutral_mass(sequence: str) -> float:
    """Neutral mass of a head-to-tail cyclic peptide (no H₂O)."""
    return sum(RESIDUE_MASS.get(aa, 0.0) for aa in sequence)

def linear_neutral_mass(sequence: str) -> float:
    return sum(RESIDUE_MASS.get(aa, 0.0) for aa in sequence) + H2O

def precursor_neutral(mz: float, charge: int) -> float:
    return mz * charge - charge * PROTON

# ══════════════════════════════════════════════════════════════════════════════
#  DE NOVO CYCLIC CANDIDATE SEARCH
# ══════════════════════════════════════════════════════════════════════════════

def _all_compositions(
    target_mass: float,
    residues: list[tuple[str, float]],
    n_min: int,
    n_max: int,
    tol: float,
    current: list[str],
    current_mass: float,
    results: list[list[str]],
) -> None:
    """Recursive mass-directed composition finder (branch-and-bound)."""
    n = len(current)
    if n >= n_min and abs(current_mass - target_mass) <= tol:
        results.append(current[:])

    if n >= n_max or current_mass > target_mass + tol:
        return

    last = current[-1] if current else None
    for aa, mass in residues:
        # Prune: only extend in lexicographic order of AA to avoid redundant combos
        if last and aa < last:
            continue
        new_mass = current_mass + mass
        if new_mass > target_mass + tol:
            continue
        current.append(aa)
        _all_compositions(target_mass, residues, n_min, n_max, tol,
                          current, new_mass, results)
        current.pop()


def compositions_for_mass(
    target_mass: float,
    n_min: int = 2,
    n_max: int = 9,
    tol: float = 0.02,
    allowed_aa: str | None = None,
) -> list[str]:
    """
    Find all amino acid compositions (sorted) whose summed residue mass
    matches target_mass within tol Da.
    Returns list of composition strings (e.g. 'AAI', 'AV', 'GGF').
    """
    if allowed_aa:
        residues = [(aa, RESIDUE_MASS[aa]) for aa in allowed_aa if aa in RESIDUE_MASS]
    else:
        residues = sorted(RESIDUE_MASS.items(), key=lambda x: x[0])

    results: list[list[str]] = []
    _all_compositions(target_mass, residues, n_min, n_max, tol,
                      [], 0.0, results)
    return ["".join(r) for r in results]


def all_permutations(composition: str) -> list[str]:
    """All unique sequences for a given amino acid composition."""
    seen = set()
    out = []
    for perm in itertools.permutations(composition):
        s = "".join(perm)
        key = canonical_cyclic(s)
        if key not in seen:
            seen.add(key)
            out.append(s)
    return out

# ══════════════════════════════════════════════════════════════════════════════
#  CYCLIC MS2 FRAGMENT FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def bn_ions_all_rotations(sequence: str, charge: int = 1) -> list[float]:
    """
    Generate all bn ions for a cyclic peptide.
    For a cyclic peptide of n residues, each of the n ring-opening positions
    produces a series of (n-1) bn ions. Total: n*(n-1) ions (with overlap).
    Returns flat list of unique m/z values.
    """
    n = len(sequence)
    mz_set: set[float] = set()
    for start in range(n):
        rot = sequence[start:] + sequence[:start]
        b_mass = 0.0
        for aa in rot[:-1]:  # n-1 bn ions per rotation
            b_mass += RESIDUE_MASS.get(aa, 0.0)
            mz = round((b_mass + PROTON * charge) / charge, 5)
            mz_set.add(mz)
    return sorted(mz_set)


def residue_loss_ions(sequence: str, cyclic_mh: float) -> dict[str, float]:
    """
    Diagnostic loss ions: [M+H]⁺ − residue_mass for each unique residue.
    These are the first ions formed by ring-opening (loss of one residue
    from the cyclic [M+H]⁺ to give a linear bn-1 fragment).
    """
    losses = {}
    for aa in set(sequence):
        loss_mz = round(cyclic_mh - RESIDUE_MASS[aa], 5)
        losses[f"loss_{aa}"] = loss_mz
    return losses


def score_spectrum_vs_cyclic(
    mz_obs: np.ndarray,
    int_obs: np.ndarray,
    sequence: str,
    tol: float = 0.02,
) -> dict:
    """
    Score a MS2 spectrum against a cyclic peptide candidate.

    Scoring dimensions
    ──────────────────
    1. bn_coverage       — fraction of theoretical bn ions matched
    2. loss_coverage     — fraction of residue-loss ions matched
    3. intensity_score   — matched intensity / total intensity
    4. y1_absence        — penalty if a y1 ion is present (linear indicator)
    5. a1_immonium       — immonium ions present (cyclic-neutral)

    Returns dict with per-dimension scores and composite score (0–1).
    """
    if len(mz_obs) == 0:
        return {"composite": 0.0}

    total_int = float(np.sum(int_obs))
    cyclic_mh = cyclic_neutral_mass(sequence) + PROTON

    # ── bn ions ──────────────────────────────────────────────────────────────
    bn_theo = bn_ions_all_rotations(sequence, charge=1)
    bn_matched = sum(
        1 for mz in bn_theo if np.any(np.abs(mz_obs - mz) <= tol)
    )
    bn_coverage = bn_matched / len(bn_theo) if bn_theo else 0.0

    # ── residue loss ions ────────────────────────────────────────────────────
    loss_ions = residue_loss_ions(sequence, cyclic_mh)
    loss_matched = sum(
        1 for mz in loss_ions.values() if np.any(np.abs(mz_obs - mz) <= tol)
    )
    loss_coverage = loss_matched / len(loss_ions) if loss_ions else 0.0

    # ── intensity score (matched peaks' contribution) ────────────────────────
    all_theo = bn_theo + list(loss_ions.values())
    matched_int = 0.0
    for mz in all_theo:
        hits = np.where(np.abs(mz_obs - mz) <= tol)[0]
        if len(hits) > 0:
            matched_int += float(np.max(int_obs[hits]))
    intensity_score = matched_int / total_int if total_int > 0 else 0.0

    # ── y1 absence check (y1 = C-terminal residue + H2O + H) ─────────────────
    # If y1 ions are strongly present → more likely linear
    y1_penalty = 0.0
    for aa in set(sequence):
        y1_mz = RESIDUE_MASS[aa] + H2O + PROTON
        if np.any(np.abs(mz_obs - y1_mz) <= tol):
            y1_penalty += 0.1
    y1_absence = max(0.0, 1.0 - y1_penalty)

    # ── immonium ions ─────────────────────────────────────────────────────────
    # Immonium = residue_mass - 28 (CO loss) + H
    immonium_theo = {
        aa: round(RESIDUE_MASS[aa] - 27.9949 + PROTON, 5)
        for aa in set(sequence)
    }
    imm_matched = sum(
        1 for mz in immonium_theo.values() if np.any(np.abs(mz_obs - mz) <= tol)
    )
    immonium_score = imm_matched / len(immonium_theo) if immonium_theo else 0.0

    # ── composite score ───────────────────────────────────────────────────────
    composite = (
        0.35 * bn_coverage +
        0.25 * loss_coverage +
        0.20 * intensity_score +
        0.10 * y1_absence +
        0.10 * immonium_score
    )

    return {
        "composite":       round(composite, 4),
        "bn_coverage":     round(bn_coverage, 3),
        "loss_coverage":   round(loss_coverage, 3),
        "intensity_score": round(intensity_score, 3),
        "y1_absence":      round(y1_absence, 3),
        "immonium_score":  round(immonium_score, 3),
        "bn_matched":      bn_matched,
        "bn_total":        len(bn_theo),
        "loss_matched":    loss_matched,
        "loss_total":      len(loss_ions),
    }


def confidence_tier(score: float) -> str:
    if score >= 0.60:
        return "HIGH"
    if score >= 0.35:
        return "MEDIUM"
    if score >= 0.15:
        return "LOW"
    return "VERY_LOW"

# ══════════════════════════════════════════════════════════════════════════════
#  3D CONFORMER GENERATION (RDKit)
# ══════════════════════════════════════════════════════════════════════════════

def build_cyclic_smiles(sequence: str) -> str | None:
    """
    Build head-to-tail cyclic peptide SMILES.
    Ring: O=C1-[Cα(sc)]-N-C(=O)-[Cα(sc)]-N-...-C(=O)-N-1
    Returns None for sequences containing proline (P) — ring-within-ring
    requires separate handling.
    """
    if "P" in sequence:
        return None  # proline excluded (future extension)

    parts = []
    for i, aa in enumerate(sequence):
        sc = AA_SIDECHAIN.get(aa)
        if sc is None:
            return None
        if i == 0:
            parts.append(f"[C@@H]({sc})" if sc else "C")
        else:
            parts.append(f"NC(=O)[C@@H]({sc})" if sc else "NC(=O)C")

    return "O=C1" + "".join(parts) + "N1"



# ── xTB optional import ───────────────────────────────────────────────────────
try:
    from tblite.interface import Calculator as _XTBCalculator
    from scipy.optimize import minimize as _scipy_minimize
    HAS_XTB = True
except ImportError:
    HAS_XTB = False

ANGSTROM_TO_BOHR = 1.8897259886
BOHR_TO_ANGSTROM = 1.0 / ANGSTROM_TO_BOHR
HARTREE_TO_KCAL  = 627.509474
KT_298           = 0.593   # kcal/mol at 298 K


# ─── MMFF energy minimization ─────────────────────────────────────────────────

def _mmff_minimize(mol_h: "Chem.Mol", cid: int) -> float | None:
    """Minimize one conformer in-place. Returns energy or None."""
    ff_props = AllChem.MMFFGetMoleculeProperties(mol_h)
    if ff_props:
        ff = AllChem.MMFFGetMoleculeForceField(mol_h, ff_props, confId=cid)
    else:
        ff = AllChem.UFFGetMoleculeForceField(mol_h, confId=cid)
    if ff is None:
        return None
    ff.Minimize()
    return ff.CalcEnergy()


# ─── Statistical conformer selection (energy + RMSD clustering) ───────────────

def _cluster_conformers(
    mol_h: "Chem.Mol",
    cids: list[int],
    kt_window: float = 2.0,
    rmsd_thresh: float = 0.5,
) -> list[tuple[float, int, float]]:
    """
    Return the statistically significant, structurally diverse subset of
    conformers. No arbitrary rank cutoff.

    Algorithm
    ─────────
    1. MMFF-minimize all conformers → energy list
    2. Energy window filter: keep only conformers within kt_window × kT(298K)
       of the global minimum. This is the Boltzmann-relevant population at RT;
       conformers outside this window contribute < e^(-kt_window) ≈ 1–14% of
       the partition function and are thermodynamically negligible.
    3. RMSD diversity filter: within the energy-selected pool, discard any
       conformer whose heavy-atom RMSD to an already-kept conformer is below
       rmsd_thresh Å. This removes structurally redundant frames while
       preserving genuinely different ring puckers, cis/trans amide isomers,
       and sidechain rotamers.
    4. Boltzmann weight: w_i = exp(−ΔE_i / kT) / Z, where ΔE_i is relative
       to the ensemble minimum.

    Returns list of (energy_kcal, confId, boltzmann_weight), sorted by energy.
    """
    energies: list[tuple[float, int]] = []
    for cid in cids:
        e = _mmff_minimize(mol_h, cid)
        if e is not None:
            energies.append((e, cid))

    if not energies:
        return []

    energies.sort(key=lambda x: x[0])
    e_min = energies[0][0]
    e_cut  = e_min + kt_window * KT_298

    # Step 2: energy window
    pool = [(e, cid) for e, cid in energies if e <= e_cut]

    # Step 3: RMSD diversity filter
    kept: list[tuple[float, int]] = [pool[0]]
    for e, cid in pool[1:]:
        diverse = True
        for _, k_cid in kept:
            try:
                rmsd = AllChem.GetBestRMS(mol_h, mol_h, k_cid, cid)
            except Exception:
                rmsd = 0.0
            if rmsd < rmsd_thresh:
                diverse = False
                break
        if diverse:
            kept.append((e, cid))

    # Step 4: Boltzmann weights
    e_arr   = [e for e, _ in kept]
    e_rel   = [e - e_min for e in e_arr]
    boltz   = [2.718281828 ** (-de / KT_298) for de in e_rel]
    z       = sum(boltz)
    weights = [b / z for b in boltz]

    return [(e, cid, round(w, 4)) for (e, cid), w in zip(kept, weights)]


# ─── GFN2-xTB geometry optimization ──────────────────────────────────────────

def _xtb_optimize(
    mol_h: "Chem.Mol",
    cid:   int,
    method: str = "GFN2-xTB",
    max_steps: int = 300,
    grad_tol:  float = 1e-3,
) -> tuple["Chem.Mol", float, bool, int]:
    """
    Refine one conformer with GFN2-xTB using L-BFGS (scipy) + tblite.

    The MMFF conformer is used as the starting geometry. L-BFGS uses analytic
    gradients from tblite, converging in ~20–50 single-point evaluations per
    conformer — ~0.2–2 s per structure depending on ring size.

    Returns (optimized_mol_copy, energy_kcal, converged, n_grad_calls).
    """
    from copy import deepcopy

    mol_opt = deepcopy(mol_h)
    conf    = mol_opt.GetConformer(cid)
    numbers = np.array([int(a.GetAtomicNum()) for a in mol_opt.GetAtoms()])

    x0 = np.array([
        [conf.GetAtomPosition(i).x,
         conf.GetAtomPosition(i).y,
         conf.GetAtomPosition(i).z]
        for i in range(mol_opt.GetNumAtoms())
    ]).flatten() * ANGSTROM_TO_BOHR

    call_count = [0]

    def energy_grad(x: np.ndarray) -> tuple[float, np.ndarray]:
        call_count[0] += 1
        pos  = x.reshape(-1, 3)
        calc = _XTBCalculator(method, numbers, pos)
        calc.set("verbosity", 0)
        res  = calc.singlepoint()
        return float(res.get("energy")), res.get("gradient").flatten()

    result = _scipy_minimize(
        energy_grad, x0,
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": max_steps, "gtol": grad_tol, "ftol": 1e-12},
    )

    pos_opt = result.x.reshape(-1, 3) * BOHR_TO_ANGSTROM
    # Write optimised coords back into a new conformer (confId 0)
    from copy import deepcopy as _dc
    mol_out = _dc(mol_h)
    # Remove all conformers, add a single new one
    mol_out.RemoveAllConformers()
    from rdkit.Chem import rdchem
    new_conf = rdchem.Conformer(mol_out.GetNumAtoms())
    from rdkit.Geometry import Point3D
    for i, (x, y, z) in enumerate(pos_opt):
        new_conf.SetAtomPosition(i, Point3D(x, y, z))
    mol_out.AddConformer(new_conf, assignId=True)

    energy_kcal = result.fun * HARTREE_TO_KCAL
    return mol_out, energy_kcal, result.success, call_count[0]


# ─── Main conformer pipeline ──────────────────────────────────────────────────

def generate_conformers(
    sequence:    str,
    n_embed:     int   = 200,
    kt_window:   float = 2.0,
    rmsd_thresh: float = 0.5,
    use_xtb:     bool  = True,
    random_seed: int   = 42,
) -> dict:
    """
    Full statistical conformer pipeline for a cyclic peptide sequence.

    Stage 1 — Embedding
        ETKDGv3 generates n_embed starting geometries (distance geometry,
        experimentally derived torsion preferences, macrocycle torsions for
        rings ≥ 8 residues).

    Stage 2 — MMFF94 minimization of all embeddings
        Force-field minimization of every conformer provides a reliable
        energy surface for the statistical selection step.

    Stage 3 — Statistical selection (energy + RMSD)
        Energy window: keep conformers within kt_window × kT(298K) of the
        global minimum. Default 2.0 kT retains ~87% of the equilibrium
        population at room temperature; conformers outside this window are
        thermodynamically negligible at physiological conditions.
        RMSD filter: remove structurally redundant frames (< rmsd_thresh Å
        heavy-atom RMSD to any already-selected conformer).

    Stage 4 — GFN2-xTB refinement (optional, requires tblite)
        Each selected MMFF conformer is re-optimized with GFN2-xTB using
        L-BFGS (analytic gradients). xTB captures hydrogen bonding, proper
        amide geometry, and intramolecular electrostatics at ~1/1000 the
        cost of DFT while being significantly more accurate than MMFF94
        for polar cyclic peptides.

    Parameters
    ──────────
    sequence    : one-letter AA sequence
    n_embed     : conformers to generate (more = better sampling for large rings)
    kt_window   : energy window in multiples of kT(298K) = 0.593 kcal/mol
    rmsd_thresh : minimum RMSD (Å) between kept conformers
    use_xtb     : whether to refine with GFN2-xTB if available
    random_seed : ETKDGv3 seed for reproducibility

    Returns
    ───────
    dict with keys:
        smiles           : SMILES string
        conformers       : list of dicts (one per kept conformer)
            confId       : RDKit conformer ID (pre-xTB)
            mmff_energy  : MMFF energy (kcal/mol)
            delta_e_mmff : ΔE from MMFF global min (kcal/mol)
            boltzmann_w  : Boltzmann weight in MMFF ensemble
            xtb_energy   : GFN2-xTB energy (kcal/mol) or None
            xtb_converged: bool
            pdb_block    : PDB string of final geometry
        n_embedded       : total conformers embedded
        n_selected       : conformers passing statistical filter
        global_min_kcal  : lowest MMFF energy
        xtb_available    : bool
        status           : 'ok' | 'no_smiles' | 'embed_failed' | 'proline'
    """
    result_base = {
        "smiles":          None,
        "conformers":      [],
        "n_embedded":      0,
        "n_selected":      0,
        "global_min_kcal": None,
        "xtb_available":   HAS_XTB,
        "status":          "ok",
    }

    if not HAS_RDKIT:
        return {**result_base, "status": "no_rdkit"}

    smiles = build_cyclic_smiles(sequence)
    if smiles is None:
        return {**result_base, "status": "proline" if "P" in sequence else "no_smiles"}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {**result_base, "smiles": smiles, "status": "invalid_smiles"}

    result_base["smiles"] = smiles
    mol_h = Chem.AddHs(mol)

    # ── Stage 1: Embedding ────────────────────────────────────────────────────
    params = AllChem.ETKDGv3()
    params.randomSeed      = random_seed
    params.numThreads      = 0
    params.useSmallRingTorsions   = True
    params.useMacrocycleTorsions  = len(sequence) >= 8

    cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=n_embed, params=params))
    if not cids:
        params2          = AllChem.ETDG()
        params2.randomSeed = random_seed
        cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=20, params=params2))
    if not cids:
        return {**result_base, "smiles": smiles, "status": "embed_failed"}

    result_base["n_embedded"] = len(cids)

    # ── Stage 2+3: MMFF minimization + statistical selection ─────────────────
    selected = _cluster_conformers(
        mol_h, cids,
        kt_window=kt_window,
        rmsd_thresh=rmsd_thresh,
    )

    if not selected:
        return {**result_base, "smiles": smiles, "status": "cluster_failed"}

    result_base["n_selected"]      = len(selected)
    result_base["global_min_kcal"] = round(selected[0][0], 4)

    # ── Stage 4: xTB refinement + PDB export per conformer ───────────────────
    conformer_list = []
    do_xtb = use_xtb and HAS_XTB

    for rank, (mmff_e, cid, boltz_w) in enumerate(selected, start=1):
        entry: dict = {
            "rank":         rank,
            "confId":       cid,
            "mmff_energy":  round(mmff_e, 4),
            "delta_e_mmff": round(mmff_e - selected[0][0], 4),
            "boltzmann_w":  boltz_w,
            "xtb_energy":   None,
            "xtb_delta_e":  None,
            "xtb_converged": None,
            "xtb_calls":    None,
            "pdb_block":    None,
        }

        if do_xtb:
            try:
                mol_xtb, xtb_e, conv, n_calls = _xtb_optimize(mol_h, cid)
                entry["xtb_energy"]    = round(xtb_e, 4)
                entry["xtb_converged"] = conv
                entry["xtb_calls"]     = n_calls
                pdb = Chem.MolToPDBBlock(mol_xtb)
            except Exception as exc:
                log.warning(f"xTB failed for {sequence} conf {rank}: {exc}")
                pdb = Chem.MolToPDBBlock(mol_h, confId=cid)
        else:
            pdb = Chem.MolToPDBBlock(mol_h, confId=cid)

        entry["pdb_block"] = pdb
        conformer_list.append(entry)

    # Re-rank by xTB energy if available
    if do_xtb and any(c["xtb_energy"] is not None for c in conformer_list):
        conformer_list.sort(key=lambda c: c["xtb_energy"] if c["xtb_energy"] is not None else 1e9)
        for i, c in enumerate(conformer_list, start=1):
            c["rank"] = i
            if c["xtb_energy"] is not None:
                c["xtb_delta_e"] = round(c["xtb_energy"] - conformer_list[0]["xtb_energy"], 4)

    result_base["conformers"] = conformer_list
    return result_base

# ══════════════════════════════════════════════════════════════════════════════
#  METAL CHELATION PREDICTION
# ══════════════════════════════════════════════════════════════════════════════

# Donor atom weights per residue and metal (0 = no donation, 3 = strong)
# Based on Irving-Williams series and coordination chemistry of short peptides:
# Cu²⁺ > Zn²⁺ > Fe²⁺ for N/S donors; Fe²⁺ > Cu²⁺ for carboxylate O-donors
CHELATION_WEIGHTS: dict[str, dict[str, float]] = {
    #         Fe2+   Cu2+   Zn2+
    "H": {"Fe": 2.0, "Cu": 3.0, "Zn": 3.0},   # imidazole Nε2 — strong
    "C": {"Fe": 1.0, "Cu": 3.0, "Zn": 2.5},   # thiol S — very strong Cu/Zn
    "D": {"Fe": 2.0, "Cu": 1.5, "Zn": 1.5},   # carboxylate O — good Fe
    "E": {"Fe": 2.0, "Cu": 1.5, "Zn": 1.5},   # carboxylate O
    "Y": {"Fe": 1.0, "Cu": 2.0, "Zn": 1.0},   # phenolate O
    "N": {"Fe": 0.5, "Cu": 1.0, "Zn": 0.5},   # amide O/N — weak
    "Q": {"Fe": 0.5, "Cu": 1.0, "Zn": 0.5},   # amide O/N — weak
    "S": {"Fe": 0.3, "Cu": 0.5, "Zn": 0.3},   # hydroxyl O — very weak
    "T": {"Fe": 0.3, "Cu": 0.5, "Zn": 0.3},   # hydroxyl O — very weak
    # Backbone C=O — present in all residues, weak but additive
    "_backbone": {"Fe": 0.2, "Cu": 0.2, "Zn": 0.2},
}

DONOR_LABEL: dict[str, str] = {
    "H": "His-Nε2", "C": "Cys-S", "D": "Asp-COO⁻", "E": "Glu-COO⁻",
    "Y": "Tyr-OH",  "N": "Asn-CO", "Q": "Gln-CO",
    "S": "Ser-OH",  "T": "Thr-OH",
}

def chelation_tier(score: float) -> str:
    if score >= 0.60:
        return "HIGH"
    if score >= 0.35:
        return "MEDIUM"
    if score >= 0.10:
        return "LOW"
    return "NONE"


def predict_chelation(
    sequence: str,
    metals: list[str] | None = None,
) -> dict:
    """
    Rule-based metal chelation predictor for cyclic peptides.

    Scoring logic
    ─────────────
    1. Sum donor atom weights for each metal across all residues.
    2. Apply a geometric bonus when ≥ 2 donors are within 3 residues of
       each other in the cyclic ring (favours 5-/6-membered chelate rings).
    3. Normalise to [0, 1] using empirical maximum for a perfect chelator
       (e.g. cyclo-HCHC would score ~1.0 for Cu).

    Parameters
    ----------
    sequence : str  — one-letter AA sequence
    metals   : list of metal symbols to predict for (default: Fe, Cu, Zn)

    Returns
    -------
    dict with per-metal scores, donor list, chelate ring sizes, tier, fenton flag.
    """
    if metals is None:
        metals = ["Fe", "Cu", "Zn"]

    n = len(sequence)
    donors = []
    for i, aa in enumerate(sequence):
        if aa in CHELATION_WEIGHTS:
            donors.append((i, aa, CHELATION_WEIGHTS[aa]))
    # Backbone C=O: every residue contributes weakly
    backbone_w = CHELATION_WEIGHTS["_backbone"]

    scores: dict[str, float] = {}
    for metal in metals:
        raw = 0.0
        # Sidechain donors
        for _, aa, w in donors:
            raw += w.get(metal, 0.0)
        # Backbone contribution (all n residues)
        raw += backbone_w.get(metal, 0.0) * n

        # Geometric bonus: pairs of strong donors close in sequence
        # In a cyclic ring, residues i and j are "close" if min(|i-j|, n-|i-j|) ≤ 3
        bonus = 0.0
        strong_donors = [(i, aa) for i, aa, w in donors if w.get(metal, 0.0) >= 1.5]
        for a_idx in range(len(strong_donors)):
            for b_idx in range(a_idx + 1, len(strong_donors)):
                i_pos, _ = strong_donors[a_idx]
                j_pos, _ = strong_donors[b_idx]
                ring_dist = min(abs(i_pos - j_pos), n - abs(i_pos - j_pos))
                if ring_dist <= 3:
                    # 5-membered chelate (dist=2) > 6-membered (dist=3)
                    if ring_dist == 2:
                        bonus += 0.20
                    elif ring_dist == 3:
                        bonus += 0.12
                    else:
                        bonus += 0.05

        raw += bonus

        # Normalise: empirical max ≈ 3 strong donors + 2 backbone + 0.4 bonus = ~8
        MAX_SCORE = 8.0
        scores[metal] = round(min(raw / MAX_SCORE, 1.0), 4)

    # Donor atom labels for report
    donor_labels = [DONOR_LABEL.get(aa, aa) for _, aa, _ in donors]

    # Chelate ring sizes (distances between strong donor pairs in ring)
    chelate_rings: list[int] = []
    strong_all = [(i, aa) for i, aa, w in donors if any(w.get(m, 0) >= 1.5 for m in metals)]
    for a_idx in range(len(strong_all)):
        for b_idx in range(a_idx + 1, len(strong_all)):
            i_pos, _ = strong_all[a_idx]
            j_pos, _ = strong_all[b_idx]
            rd = min(abs(i_pos - j_pos), n - abs(i_pos - j_pos))
            # chelate ring size = ring_dist + 2 (the metal closes the ring)
            chelate_rings.append(rd + 2)

    chelate_rings = sorted(set(chelate_rings))

    # Best overall score
    best_score = max(scores.values()) if scores else 0.0

    # Fenton risk reduction: needs both Fe AND Cu chelation
    fenton_reduction = (
        scores.get("Fe", 0.0) >= 0.35 and scores.get("Cu", 0.0) >= 0.35
    )

    return {
        "sequence":            sequence,
        "donor_residues":      ", ".join(donor_labels) if donor_labels else "backbone only",
        "n_donors":            len(donors),
        **{f"chelation_score_{m}": scores.get(m, 0.0) for m in metals},
        "chelate_ring_sizes":  str(chelate_rings) if chelate_rings else "none",
        "best_chelation_score": round(best_score, 4),
        "chelation_tier":      chelation_tier(best_score),
        "fenton_risk_reduction": fenton_reduction,
    }

# ══════════════════════════════════════════════════════════════════════════════
#  STATISTICAL VALIDATION: DECOY SCORING / FDR / F-STATISTIC
# ══════════════════════════════════════════════════════════════════════════════

SCORE_DIMS = ["bn_coverage", "loss_coverage", "intensity_score",
              "y1_absence", "immonium_score"]
SCORE_WEIGHTS = np.array([0.35, 0.25, 0.20, 0.10, 0.10])
MIN_DECOYS_FOR_PVAL = 20  # below this, p-value is unreliable


def generate_decoys(sequence: str, n_decoys: int = 100) -> list[str]:
    """
    Generate decoy sequences by shuffling residue order, preserving composition
    (and therefore mass — identical precursor mass → matched null distribution).

    For short sequences with few unique permutations (e.g. IAA → only 3 rotations),
    isobaric substitution is used: I↔L share nominal mass 113.084 Da, so replacing
    one I with L gives a valid null that differs only in structural isomerism.

    Returns a list of unique decoy sequences, possibly fewer than n_decoys for
    short or compositionally degenerate sequences.
    """
    import itertools as _it
    import random as _rand

    n = len(sequence)
    # All cyclic equivalents of the original (not valid decoys)
    excluded: set[str] = set()
    for i in range(n):
        excluded.add(sequence[i:] + sequence[:i])
        rev = sequence[::-1]
        excluded.add(rev[i:] + rev[:i])

    # All unique permutations of the composition
    all_perms: set[str] = set(
        "".join(p) for p in _it.permutations(sequence)
    )
    valid = all_perms - excluded

    # Isobaric fallback for degenerate cases
    if len(valid) < MIN_DECOYS_FOR_PVAL:
        pseudo: set[str] = set()
        for i, aa in enumerate(sequence):
            sub = "L" if aa == "I" else ("I" if aa == "L" else None)
            if sub:
                s = sequence[:i] + sub + sequence[i+1:]
                if s not in excluded:
                    pseudo.add(s)
                    # also add permutations of the substituted sequence
                    for p in _it.permutations(s):
                        ps = "".join(p)
                        if ps not in excluded:
                            pseudo.add(ps)
        valid |= pseudo

    result = list(valid)
    if len(result) > n_decoys:
        result = _rand.sample(result, n_decoys)
    return result


def _score_vector(candidate: dict) -> np.ndarray:
    """Extract the 5-dimensional score vector from a candidate dict."""
    return np.array([candidate.get(d, 0.0) for d in SCORE_DIMS], dtype=float)


def compute_empirical_pvalue(
    target_score: float,
    decoy_scores: np.ndarray,
) -> float:
    """
    One-tailed empirical p-value: P(decoy ≥ target_score | H₀).
    Laplace (add-1) pseudocount prevents p = 0 for unseen extremes.
    """
    if len(decoy_scores) == 0:
        return np.nan
    return float((np.sum(decoy_scores >= target_score) + 1) / (len(decoy_scores) + 1))


def bh_fdr(
    p_values: np.ndarray,
    alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Benjamini-Hochberg (1995) FDR correction.

    Assumes independence or positive dependence (PRDS) between tests —
    satisfied here because all candidates are scored against the same
    global decoy pool.

    Returns
    ───────
    (adj_p, rejected) — BH-adjusted p-values and boolean rejection mask.
    """
    n = len(p_values)
    if n == 0:
        return np.array([]), np.array([], dtype=bool)

    order  = np.argsort(p_values)
    ranked = np.arange(1, n + 1, dtype=float)

    # Step-up: find the largest k such that p_(k) ≤ k/m × α
    bh_thresh  = ranked / n * alpha
    below      = p_values[order] <= bh_thresh
    k_max      = int(np.max(np.where(below)[0])) if np.any(below) else -1

    rejected          = np.zeros(n, dtype=bool)
    if k_max >= 0:
        rejected[order[:k_max + 1]] = True

    # BH-adjusted p-values (step-up minimum)
    adj = np.minimum(1.0, p_values[order] * n / ranked)
    for i in range(n - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])

    adj_out          = np.zeros(n)
    adj_out[order]   = adj
    return adj_out, rejected


def storey_qvalue(
    p_values: np.ndarray,
    lambdas: np.ndarray | None = None,
) -> tuple[np.ndarray, float]:
    """
    Storey & Tibshirani (2003) q-value with π₀ estimation.

    π₀ = proportion of truly null hypotheses, estimated from the tail of the
    p-value distribution. Produces less conservative FDR estimates than BH
    when many nulls are present (typical in metabolomics/peptidomics screens).

    Returns (q_values, pi0_estimate).
    """
    n = len(p_values)
    if n == 0:
        return np.array([]), 1.0

    if lambdas is None:
        lambdas = np.arange(0.05, 0.90, 0.05)

    # π₀ estimation: use the plateau of π̂₀(λ) in the tail
    pi0_hats = np.array([
        np.sum(p_values > lam) / (n * (1.0 - lam))
        for lam in lambdas
    ])
    # Conservative: take the minimum of the last 5 estimates (right tail)
    pi0 = float(np.clip(np.min(pi0_hats[-5:]) if len(pi0_hats) >= 5
                        else pi0_hats[-1], 0.0, 1.0))

    order  = np.argsort(p_values)
    ranked = np.arange(1, n + 1, dtype=float)

    q = p_values[order] * n * pi0 / ranked
    # Enforce monotonicity (step-up)
    for i in range(n - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])

    q_out          = np.zeros(n)
    q_out[order]   = np.minimum(1.0, q)
    return q_out, pi0


def hotellings_t2(
    X_target: np.ndarray,
    X_decoy:  np.ndarray,
) -> dict:
    """
    Hotelling's T² test: multivariate analogue of the t-test.

    Tests H₀: μ_target = μ_decoy in the 5-dimensional score space
    (bn_coverage, loss_coverage, intensity_score, y1_absence, immonium_score).

    Converted to an exact F-statistic under the assumption that both groups
    follow multivariate normal distributions (approximately satisfied by the
    Central Limit Theorem for n > 30).

    F = T² × (n₁+n₂−p−1) / ((n₁+n₂−2) × p)  ~  F(p, n₁+n₂−p−1)

    Returns dict with T2, F, df1, df2, p_value, per-dimension MW-U tests.
    """
    from scipy.stats import f as _f_dist, mannwhitneyu as _mwu

    n1, p = X_target.shape
    n2    = X_decoy.shape[0]

    if n1 < 5 or n2 < 5:
        return {
            "T2": np.nan, "F_stat": np.nan,
            "df1": p, "df2": np.nan, "p_value_multivariate": np.nan,
            "note": "insufficient_samples",
        }

    mu1 = X_target.mean(axis=0)
    mu2 = X_decoy.mean(axis=0)
    S1  = np.cov(X_target.T)
    S2  = np.cov(X_decoy.T)
    Sp  = ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)

    try:
        Sp_inv = np.linalg.pinv(Sp)
    except np.linalg.LinAlgError:
        Sp_inv = np.eye(p)

    diff = mu1 - mu2
    T2   = float((n1 * n2 / (n1 + n2)) * float(diff @ Sp_inv @ diff))
    df1  = p
    df2  = n1 + n2 - p - 1
    F    = T2 * (n1 + n2 - p - 1) / ((n1 + n2 - 2) * p)
    p_mv = float(1.0 - _f_dist.cdf(F, df1, df2)) if df2 > 0 else np.nan

    # Per-dimension Mann-Whitney U (non-parametric, one-tailed)
    dim_results = {}
    mw_pvals = []
    for i, dim in enumerate(SCORE_DIMS):
        try:
            u_stat, p_mw = _mwu(
                X_target[:, i], X_decoy[:, i], alternative="greater"
            )
            dim_results[f"mw_U_{dim}"]   = round(float(u_stat), 1)
            dim_results[f"mw_p_{dim}"]   = round(float(p_mw), 6)
            mw_pvals.append(p_mw)
        except Exception:
            dim_results[f"mw_U_{dim}"] = np.nan
            dim_results[f"mw_p_{dim}"] = np.nan
            mw_pvals.append(1.0)

    # BH correction across the 5 dimension tests
    mw_arr       = np.array(mw_pvals)
    adj_mw, _    = bh_fdr(mw_arr, alpha=0.05)
    for i, dim in enumerate(SCORE_DIMS):
        dim_results[f"mw_padj_{dim}"] = round(float(adj_mw[i]), 6)

    return {
        "n_targets":              n1,
        "n_decoys":               n2,
        "T2":                     round(T2, 3),
        "F_stat":                 round(F, 3),
        "df1":                    df1,
        "df2":                    df2,
        "p_value_multivariate":   round(p_mv, 8) if not np.isnan(p_mv) else np.nan,
        **dim_results,
    }


def run_statistics(
    candidates:   list[dict],
    n_decoys_per: int   = 100,
    fdr_alpha:    float = 0.05,
) -> tuple[list[dict], dict]:
    """
    Full statistical validation pipeline.

    Strategy
    ────────
    Builds a GLOBAL decoy score pool by scoring each spectrum against
    shuffled (same-composition) versions of its best candidate. This mirrors
    target-decoy competition in shotgun proteomics (Elias & Gygi 2007) and
    avoids the problem of insufficient permutations for short sequences.

    Steps
    ─────
    1. For each candidate, generate up to n_decoys_per shuffled sequences
       and score them against the same MS2 spectrum.
    2. Pool ALL decoy composite scores into one global null distribution.
    3. Compute per-candidate empirical p-value against the global pool.
    4. Apply BH-FDR and Storey q-value correction across all candidates.
    5. Run Hotelling T² comparing the 5D score vectors of all targets
       vs all decoys (global multivariate test).

    Parameters
    ──────────
    candidates   : list of candidate dicts (from the main scoring loop),
                   each must contain the 5 score dimensions and the raw
                   mz_array / int_array for decoy rescoring (passed via
                   '_mz' and '_int' keys — stripped before output).
    n_decoys_per : number of decoy shuffles per candidate
    fdr_alpha    : FDR threshold for rejection

    Returns
    ───────
    (enriched_candidates, global_stats_dict)
    """
    if not candidates:
        return candidates, {}

    log.info(
        f"Running statistical validation "
        f"(n_decoys_per={n_decoys_per}, FDR α={fdr_alpha}) …"
    )

    # ── Step 1: score decoys per candidate ───────────────────────────────────
    all_decoy_composites: list[float] = []
    target_vectors:       list[np.ndarray] = []
    decoy_vectors:        list[np.ndarray] = []

    for cand in candidates:
        seq     = cand["sequence"]
        mz_arr  = cand.pop("_mz",  np.array([]))
        int_arr = cand.pop("_int", np.array([]))

        decoy_seqs = generate_decoys(seq, n_decoys=n_decoys_per)
        cand["n_decoys_generated"] = len(decoy_seqs)

        target_vectors.append(_score_vector(cand))

        for d_seq in decoy_seqs:
            d_sc = score_spectrum_vs_cyclic(mz_arr, int_arr, d_seq)
            all_decoy_composites.append(d_sc["composite"])
            decoy_vectors.append(_score_vector(d_sc))

    decoy_pool = np.array(all_decoy_composites)
    log.info(
        f"  Decoy pool: {len(decoy_pool)} scores "
        f"(mean={decoy_pool.mean():.3f}, "
        f"sd={decoy_pool.std():.3f})"
    )

    # ── Step 2: empirical p-values ────────────────────────────────────────────
    p_vals = np.array([
        compute_empirical_pvalue(c["composite"], decoy_pool)
        for c in candidates
    ])

    # ── Step 3: BH-FDR ───────────────────────────────────────────────────────
    adj_bh, rejected = bh_fdr(p_vals, alpha=fdr_alpha)

    # ── Step 4: Storey q-value ────────────────────────────────────────────────
    q_storey, pi0 = storey_qvalue(p_vals)

    # ── Attach per-candidate results ──────────────────────────────────────────
    for i, cand in enumerate(candidates):
        cand["p_value_empirical"]  = round(float(p_vals[i]), 6)
        cand["p_adj_BH"]           = round(float(adj_bh[i]), 6)
        cand["q_storey"]           = round(float(q_storey[i]), 6)
        cand["rejected_BH"]        = bool(rejected[i])
        cand["n_decoys_in_pool"]   = int(len(decoy_pool))
        # Flag if decoy pool was too small for this candidate
        cand["stat_note"] = (
            "low_decoy_count" if cand["n_decoys_generated"] < MIN_DECOYS_FOR_PVAL
            else "ok"
        )

    # ── Step 5: Hotelling T² ─────────────────────────────────────────────────
    X_target = np.array(target_vectors)
    X_decoy  = np.vstack(decoy_vectors) if decoy_vectors else np.empty((0, 5))

    global_stats = hotellings_t2(X_target, X_decoy)
    global_stats["pi0_storey"]         = round(pi0, 4)
    global_stats["fdr_alpha"]          = fdr_alpha
    global_stats["n_rejected_BH"]      = int(rejected.sum())
    global_stats["n_q_below_alpha"]    = int((q_storey <= fdr_alpha).sum())
    global_stats["decoy_pool_size"]     = int(len(decoy_pool))
    global_stats["decoy_mean_score"]    = round(float(decoy_pool.mean()), 4)
    global_stats["decoy_sd_score"]      = round(float(decoy_pool.std()),  4)

    log.info(
        f"  BH-FDR:    {global_stats['n_rejected_BH']}/{len(candidates)} "
        f"rejected at α={fdr_alpha}"
    )
    log.info(
        f"  Storey q:  {global_stats['n_q_below_alpha']}/{len(candidates)} "
        f"below α={fdr_alpha}  (π₀={pi0:.3f})"
    )
    log.info(
        f"  Hotelling: T²={global_stats['T2']}, "
        f"F({global_stats['df1']},{global_stats['df2']})="
        f"{global_stats['F_stat']}, p={global_stats['p_value_multivariate']}"
    )

    return candidates, global_stats

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def run(args: argparse.Namespace) -> None:

    if not HAS_ALPHATIMS:
        log.error("alphatims not installed. Run: pip install alphatims")
        sys.exit(1)
    if not HAS_RDKIT:
        log.warning("RDKit not found — 3D conformers will be skipped. "
                    "Install with: conda install -c conda-forge rdkit")

    d_path = Path(args.d)
    if not d_path.exists():
        log.error(f".d folder not found: {d_path}")
        sys.exit(1)

    out_dir = Path(args.out) if args.out else (
        d_path.parent / (d_path.stem + "_cyclopep")
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_dir = out_dir / "pdb_structures"
    pdb_dir.mkdir(exist_ok=True)

    log.info("═══ Wine CycloPep — de novo cyclic peptide detection ═══")
    log.info(f"Input   : {d_path}")
    log.info(f"Output  : {out_dir}")
    log.info(f"AA range: {args.min_aa}–{args.max_aa} residues")

    # ── Load .d ───────────────────────────────────────────────────────────────
    t0 = time.time()
    log.info("Loading .d folder …")
    data = TimsTOF(str(d_path))
    log.info(
        f"Loaded {time.time()-t0:.1f}s | "
        f"Frames: {data.frame_max_index} | "
        f"DDA precursors: {data.precursor_max_index - 1}"
    )

    # ── Extract MS2 ───────────────────────────────────────────────────────────
    log.info("Extracting MS2 spectra …")
    spectrum_indptr, spectrum_tof_idx, spectrum_intensities = data.index_precursors(
        centroiding_window=args.centroid,
        keep_n_most_abundant_peaks=args.top_peaks,
    )

    precursors  = data.precursors
    frames_df   = data.frames
    frame_rt    = dict(zip(frames_df.index, frames_df["Time"].values))
    mono_mzs    = precursors["MonoisotopicMz"].values.astype(float)
    avg_mzs     = precursors["AverageMz"].values.astype(float)
    charges     = precursors["Charge"].values.astype(int)
    parent_ids  = precursors["Parent"].values.astype(int)
    scan_nums   = precursors["ScanNumber"].values.astype(int)
    rt_seconds  = np.array([frame_rt.get(pid, 0.0) for pid in parent_ids])

    n_spectra = data.precursor_max_index - 1
    n_process = min(args.max_spectra, n_spectra) if args.max_spectra else n_spectra
    log.info(f"Processing {n_process}/{n_spectra} spectra")

    # ── Optional allowed AA set ───────────────────────────────────────────────
    allowed_aa = args.allowed_aa.upper().replace(",", "").replace(" ", "") if args.allowed_aa else None
    if allowed_aa:
        log.info(f"Restricting to AA: {allowed_aa}")

    # ── Optional hints (known sequences to prioritize) ────────────────────────
    hint_sequences: list[str] = []
    if args.hints:
        hint_sequences = [h.strip().upper() for h in args.hints.split(",")]
        log.info(f"Hint sequences: {hint_sequences}")

    # ══════════════════════════════════════════════════════════════════════════
    #  SPECTRUM LOOP
    # ══════════════════════════════════════════════════════════════════════════
    all_candidates: list[dict] = []
    seen_canonical: dict[str, list] = defaultdict(list)  # canonical_seq → spectra

    log.info("Searching for cyclic peptide candidates …")
    n_searched = 0

    for idx in range(1, n_process + 1):
        start = int(spectrum_indptr[idx])
        end   = int(spectrum_indptr[idx + 1])
        if end == start:
            continue

        tof_slice = spectrum_tof_idx[start:end]
        int_slice = spectrum_intensities[start:end]
        mz_slice  = data.mz_values[tof_slice]

        if len(mz_slice) < 3:
            continue

        i        = idx - 1
        mono_mz  = float(mono_mzs[i]) if mono_mzs[i] > 0 else float(avg_mzs[i])
        charge   = int(charges[i]) if charges[i] > 0 else 1
        rt_s     = float(rt_seconds[i])
        scan     = int(scan_nums[i])
        mob_val  = float(data.mobility_values[scan]) if scan < len(data.mobility_values) else 0.0

        neutral = precursor_neutral(mono_mz, charge)

        # ── Mass filter: is this neutral mass compatible with a cyclic peptide? ─
        # Cyclic peptide mass range for n_min to n_max residues
        mass_min = RESIDUE_MASS["G"] * args.min_aa        # lightest possible
        mass_max = min(args.mass_max,
                       RESIDUE_MASS["W"] * args.max_aa)   # user ceiling or AA ceiling
        if not (mass_min <= neutral <= mass_max):
            continue

        n_searched += 1

        # ── Find all compositions matching this mass ───────────────────────────
        compositions = compositions_for_mass(
            target_mass=neutral,
            n_min=args.min_aa,
            n_max=args.max_aa,
            tol=args.tol,
            allowed_aa=allowed_aa,
        )

        if not compositions:
            continue

        # ── For each composition, try all unique cyclic sequences ──────────────
        best_score   = 0.0
        best_hit     = None

        for comp in compositions:
            seqs = all_permutations(comp)

            # Prioritize hint sequences
            priority = [s for s in seqs if s in hint_sequences]
            rest     = [s for s in seqs if s not in hint_sequences]
            ordered  = priority + rest

            # Cap: very long compositions → too many permutations
            if len(ordered) > args.max_perms:
                ordered = ordered[: args.max_perms]

            for seq in ordered:
                sc = score_spectrum_vs_cyclic(mz_slice, int_slice, seq, tol=args.tol)
                if sc["composite"] > best_score:
                    best_score = sc["composite"]
                    best_hit   = {
                        "spectrum_idx":    idx,
                        "rt_s":            round(rt_s, 2),
                        "rt_min":          round(rt_s / 60, 3),
                        "precursor_mz":    round(mono_mz, 5),
                        "charge":          charge,
                        "neutral_mass":    round(neutral, 5),
                        "ion_mobility":    round(mob_val, 4),
                        "sequence":        seq,
                        "composition":     "".join(sorted(seq)),
                        "n_residues":      len(seq),
                        "theo_cyclic_mass":round(cyclic_neutral_mass(seq), 5),
                        "mass_error_da":   round(neutral - cyclic_neutral_mass(seq), 5),
                        "mass_error_ppm":  round(
                            (neutral - cyclic_neutral_mass(seq)) / cyclic_neutral_mass(seq) * 1e6, 2
                        ),
                        **sc,
                        "confidence":      confidence_tier(sc["composite"]),
                        "n_fragments":     len(mz_slice),
                    }

        if best_hit and best_score >= args.min_score:
            can_key = canonical_cyclic(best_hit["sequence"])
            best_hit["canonical_key"] = can_key
            # Store raw spectrum arrays for decoy rescoring (stripped before TSV write)
            best_hit["_mz"]  = mz_slice.copy()
            best_hit["_int"] = int_slice.copy()
            all_candidates.append(best_hit)
            seen_canonical[can_key].append(idx)

    log.info(
        f"Spectra searched: {n_searched} | "
        f"Candidates found: {len(all_candidates)} | "
        f"Unique cyclic sequences: {len(seen_canonical)}"
    )

    if not all_candidates:
        log.warning(
            "No cyclic peptide candidates found. "
            "Try lowering --min_score or widening --tol."
        )
        return

    # ── Sort by composite score ───────────────────────────────────────────────
    all_candidates.sort(key=lambda x: x["composite"], reverse=True)

    # ── Statistical validation (FDR / F-stat) ────────────────────────────────
    all_candidates, global_stats = run_statistics(
        all_candidates,
        n_decoys_per = args.n_decoys,
        fdr_alpha    = args.fdr_alpha,
    )

    # ── Write candidate TSV (strip internal _mz/_int arrays) ─────────────────
    tsv_path = out_dir / "cyclic_candidates.tsv"
    fields = [k for k in all_candidates[0].keys() if not k.startswith("_")]
    with open(tsv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(all_candidates)
    log.info(f"Candidates TSV → {tsv_path.name} ({len(all_candidates)} entries)")

    # ── Write global statistics report ───────────────────────────────────────
    stat_path = out_dir / "statistical_report.tsv"
    with open(stat_path, "w", newline="") as f:
        w2 = csv.writer(f, delimiter="\t")
        w2.writerow(["metric", "value"])
        for k, v in global_stats.items():
            w2.writerow([k, v])
    log.info(f"Statistics report → {stat_path.name}")

    # ── Unique sequences summary ──────────────────────────────────────────────
    # Best hit per unique canonical sequence
    best_per_seq: dict[str, dict] = {}
    for c in all_candidates:
        k = c["canonical_key"]
        if k not in best_per_seq or c["composite"] > best_per_seq[k]["composite"]:
            best_per_seq[k] = c

    uniq_path = out_dir / "unique_cyclic_sequences.tsv"
    uniq_rows = []
    for k, c in sorted(best_per_seq.items(), key=lambda x: x[1]["composite"], reverse=True):
        uniq_rows.append({
            "canonical_key":    k,
            "best_sequence":    c["sequence"],
            "n_residues":       c["n_residues"],
            "neutral_mass":     c["neutral_mass"],
            "best_score":       c["composite"],
            "confidence":       c["confidence"],
            "p_value_empirical":c.get("p_value_empirical", ""),
            "p_adj_BH":         c.get("p_adj_BH", ""),
            "q_storey":         c.get("q_storey", ""),
            "rejected_BH":      c.get("rejected_BH", ""),
            "stat_note":        c.get("stat_note", ""),
            "n_spectra":        len(seen_canonical[k]),
            "mass_error_ppm":   c["mass_error_ppm"],
            "bn_coverage":      c["bn_coverage"],
            "loss_coverage":    c["loss_coverage"],
        })

    with open(uniq_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(uniq_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(uniq_rows)
    log.info(f"Unique sequences  → {uniq_path.name} ({len(uniq_rows)} unique)")

    # ── Metal chelation prediction ────────────────────────────────────────────
    metals = [m.strip() for m in args.metals.split(",")] if args.metals else ["Fe", "Cu", "Zn"]
    log.info(f"Predicting metal chelation for: {metals} …")
    chelation_rows = []
    for row in uniq_rows:
        ch = predict_chelation(row["best_sequence"], metals=metals)
        chelation_rows.append(ch)

    chel_path = out_dir / "chelation_report.tsv"
    if chelation_rows:
        with open(chel_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(chelation_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(chelation_rows)
        n_high_chel = sum(1 for r in chelation_rows if r["chelation_tier"] == "HIGH")
        n_fenton    = sum(1 for r in chelation_rows if r["fenton_risk_reduction"])
        log.info(
            f"Chelation report  → {chel_path.name} | "
            f"HIGH: {n_high_chel} | Fenton-active: {n_fenton}"
        )

    # ── b/y ion detail for top candidates ────────────────────────────────────
    ion_path = out_dir / "bn_ions_detail.tsv"
    ion_rows = []
    for c in all_candidates[:50]:  # top 50 spectra
        seq = c["sequence"]
        # Reload spectrum for this index
        idx   = c["spectrum_idx"]
        start = int(spectrum_indptr[idx])
        end   = int(spectrum_indptr[idx + 1])
        tof_slice = spectrum_tof_idx[start:end]
        int_slice = spectrum_intensities[start:end]
        mz_slice  = data.mz_values[tof_slice]

        bn_theo = bn_ions_all_rotations(seq, charge=1)
        for bn_mz in bn_theo:
            hits = np.where(np.abs(mz_slice - bn_mz) <= args.tol)[0]
            matched = len(hits) > 0
            obs_mz  = float(mz_slice[hits[0]]) if matched else None
            ion_rows.append({
                "spectrum_idx": idx,
                "sequence":     seq,
                "ion_type":     "bn_cyclic",
                "theo_mz":      bn_mz,
                "obs_mz":       round(obs_mz, 5) if obs_mz else "",
                "matched":      matched,
                "delta_da":     round(obs_mz - bn_mz, 4) if obs_mz else "",
            })

        loss_ions = residue_loss_ions(seq, cyclic_neutral_mass(seq) + PROTON)
        for label, loss_mz in loss_ions.items():
            hits = np.where(np.abs(mz_slice - loss_mz) <= args.tol)[0]
            matched = len(hits) > 0
            obs_mz  = float(mz_slice[hits[0]]) if matched else None
            ion_rows.append({
                "spectrum_idx": idx,
                "sequence":     seq,
                "ion_type":     label,
                "theo_mz":      loss_mz,
                "obs_mz":       round(obs_mz, 5) if obs_mz else "",
                "matched":      matched,
                "delta_da":     round(obs_mz - loss_mz, 4) if obs_mz else "",
            })

    if ion_rows:
        with open(ion_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(ion_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(ion_rows)
        log.info(f"bn ion detail     → {ion_path.name}")

    # ── FASTA of candidates ───────────────────────────────────────────────────
    fasta_path = out_dir / "cyclic_candidates.fasta"
    with open(fasta_path, "w") as f:
        for row in uniq_rows:
            seq = row["best_sequence"]
            f.write(
                f">cyclo_{seq} score={row['best_score']} "
                f"mass={row['neutral_mass']} n_spectra={row['n_spectra']} "
                f"confidence={row['confidence']}\n"
            )
            f.write(seq + "\n")
    log.info(f"Candidates FASTA  → {fasta_path.name}")

    # ── 3D CONFORMER GENERATION (statistical + xTB) ───────────────────────────
    if HAS_RDKIT:
        do_xtb = HAS_XTB and not args.no_xtb
        log.info(
            f"Generating 3D conformers "
            f"({'MMFF → GFN2-xTB' if do_xtb else 'MMFF only'}) … "
            f"[kT-window={args.kt_window}, RMSD={args.rmsd_thresh} Å]"
        )
        n_ok   = 0
        n_skip = 0
        conformer_log  = []   # one row per sequence (summary)
        conformer_rows = []   # one row per conformer (detail)

        for row in uniq_rows:
            seq  = row["best_sequence"]
            name = f"cyclo_{seq}"
            conf_subdir = pdb_dir / name
            conf_subdir.mkdir(exist_ok=True)

            result = generate_conformers(
                seq,
                n_embed     = args.n_embed,
                kt_window   = args.kt_window,
                rmsd_thresh = args.rmsd_thresh,
                use_xtb     = do_xtb,
                random_seed = 42,
            )

            log_entry = {
                "sequence":       seq,
                "n_residues":     row["n_residues"],
                "neutral_mass":   row["neutral_mass"],
                "smiles":         result["smiles"] or "N/A",
                "status":         result["status"],
                "n_embedded":     result["n_embedded"],
                "n_selected":     result["n_selected"],
                "global_min_mmff": result["global_min_kcal"] or "",
                "xtb_used":       do_xtb and result["status"] == "ok",
            }

            if result["status"] == "ok" and result["conformers"]:
                for conf_entry in result["conformers"]:
                    rank    = conf_entry["rank"]
                    pdb_out = conf_subdir / f"{name}_conf{rank:02d}.pdb"
                    pdb_out.write_text(conf_entry["pdb_block"])

                    method_tag = "xTB" if conf_entry["xtb_energy"] is not None else "MMFF"
                    final_e    = conf_entry["xtb_energy"] or conf_entry["mmff_energy"]
                    delta_e    = conf_entry["xtb_delta_e"] if conf_entry["xtb_energy"] else conf_entry["delta_e_mmff"]

                    conformer_rows.append({
                        "sequence":     seq,
                        "rank":         rank,
                        "pdb_file":     str(pdb_out.relative_to(out_dir)),
                        "method":       method_tag,
                        "energy_kcal":  final_e,
                        "delta_e_kcal": delta_e,
                        "boltzmann_w":  conf_entry["boltzmann_w"],
                        "xtb_converged":conf_entry["xtb_converged"],
                        "mmff_energy":  conf_entry["mmff_energy"],
                    })

                log_entry["n_pdb_written"] = len(result["conformers"])
                log_entry["best_energy"]   = (
                    result["conformers"][0]["xtb_energy"]
                    or result["conformers"][0]["mmff_energy"]
                )
                n_ok += 1
            else:
                log_entry["n_pdb_written"] = 0
                log_entry["best_energy"]   = ""
                n_skip += 1

            conformer_log.append(log_entry)

        # Write conformer summary TSV
        conf_tsv = out_dir / "conformer_log.tsv"
        if conformer_log:
            with open(conf_tsv, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(conformer_log[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(conformer_log)

        # Write per-conformer detail TSV
        conf_detail_tsv = out_dir / "conformer_detail.tsv"
        if conformer_rows:
            with open(conf_detail_tsv, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(conformer_rows[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(conformer_rows)

        total_pdbs = sum(r.get("n_pdb_written", 0) for r in conformer_log)
        log.info(
            f"3D conformers: {n_ok} sequences OK, {n_skip} skipped | "
            f"{total_pdbs} PDB files → {pdb_dir.name}/"
        )
    else:
        log.warning("RDKit not available — skipping 3D conformer generation")

    # ── Pipeline summary ──────────────────────────────────────────────────────
    summary_path = out_dir / "pipeline_summary.txt"
    with open(summary_path, "w") as f:
        f.write("Wine CycloPep — Pipeline Summary\n")
        f.write("=" * 55 + "\n\n")
        f.write(f"Input .d folder        : {d_path.name}\n")
        f.write(f"Sample name            : {data.sample_name}\n")
        f.write(f"Total DDA precursors   : {n_spectra}\n")
        f.write(f"Spectra processed      : {n_process}\n")
        f.write(f"Spectra in mass range  : {n_searched}\n")
        f.write(f"Total candidates       : {len(all_candidates)}\n")
        f.write(f"Unique cyclic sequences: {len(uniq_rows)}\n\n")
        f.write("Confidence distribution:\n")
        for tier in ["HIGH", "MEDIUM", "LOW", "VERY_LOW"]:
            n = sum(1 for r in uniq_rows if r["confidence"] == tier)
            f.write(f"  {tier:<10}: {n}\n")
        f.write("\nTop 10 candidates:\n")
        f.write(f"  {'Sequence':<15} {'Score':>7} {'Mass':>10} {'nSp':>5} {'Conf'}\n")
        f.write("  " + "-" * 50 + "\n")
        for row in uniq_rows[:10]:
            f.write(
                f"  {row['best_sequence']:<15} "
                f"{row['best_score']:>7.4f} "
                f"{row['neutral_mass']:>10.4f} "
                f"{row['n_spectra']:>5} "
                f"{row['confidence']}\n"
            )
        f.write(f"\nElapsed: {time.time()-t0:.1f}s\n")
        f.write(f"\nNext step: open PDB files in PyMOL / ChimeraX\n")
        f.write(f"  pymol {pdb_dir}/*.pdb\n")
        f.write(f"\nFor MS confirmation: validate top hits with CycloBranch\n")
        f.write(f"  https://ms.biomed.cas.cz/cyclobranch/\n")

    log.info(f"Pipeline summary  → {summary_path.name}")
    log.info(f"Done in {time.time()-t0:.1f}s — outputs in: {out_dir}/")


# ══════════════════════════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════════════════════════

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Wine CycloPep — de novo cyclic peptide detection from .d",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--d", required=True,
                   help="Bruker timsTOF .d folder")
    p.add_argument("--out", default=None,
                   help="Output directory")
    p.add_argument("--min_aa", type=int, default=2,
                   help="Minimum number of residues")
    p.add_argument("--max_aa", type=int, default=9,
                   help="Maximum number of residues")
    p.add_argument("--mass_max", type=float, default=1300.0,
                   help="Hard neutral mass ceiling in Da (default 1300)")
    p.add_argument("--tol", type=float, default=0.02,
                   help="Fragment mass tolerance (Da)")
    p.add_argument("--min_score", type=float, default=0.10,
                   help="Minimum composite score to report a candidate")
    p.add_argument("--centroid", type=int, default=3,
                   help="Centroiding window for MS2 extraction")
    p.add_argument("--top_peaks", type=int, default=50,
                   help="Top N MS2 peaks to retain per spectrum")
    p.add_argument("--max_spectra", type=int, default=None,
                   help="Limit to first N spectra (for testing)")
    p.add_argument("--max_perms", type=int, default=200,
                   help="Max sequence permutations to score per composition")
    p.add_argument("--allowed_aa", type=str, default=None,
                   help="Restrict AA alphabet, e.g. GAILVFWM for non-polar search")
    p.add_argument("--hints", type=str, default=None,
                   help="Comma-separated sequence hints to prioritize, e.g. VAAG,IAA")
    p.add_argument("--n_embed", type=int, default=200,
                   help="Conformers to embed per sequence (more = better sampling)")
    p.add_argument("--kt_window", type=float, default=2.0,
                   help="Energy window in multiples of kT(298K)=0.593 kcal/mol "
                        "for conformer selection (default 2.0 = ~87%% of Boltzmann population)")
    p.add_argument("--rmsd_thresh", type=float, default=0.5,
                   help="Minimum RMSD (Å) between selected conformers (diversity filter)")
    p.add_argument("--no_xtb", action="store_true",
                   help="Skip GFN2-xTB refinement even if tblite is installed")
    p.add_argument("--metals", type=str, default="Fe,Cu,Zn",
                   help="Comma-separated metals for chelation prediction (default: Fe,Cu,Zn)")
    p.add_argument("--n_decoys", type=int, default=100,
                   help="Decoy shuffles per candidate for empirical p-value estimation")
    p.add_argument("--fdr_alpha", type=float, default=0.05,
                   help="FDR significance threshold for BH correction")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args)
