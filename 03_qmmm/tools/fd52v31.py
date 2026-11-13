#!/usr/bin/env python3
"""fd52v3.py  QM/MM partial-Hessian (PHVA) by external central finite differences (CP2K 2024.1).

Configuration is mandatory: environment variable FD_CFG=<json> (see fd_qm52v3.json).
Subcommands (run in 07_cp2k_qmmm):
  validate  TEMPLATE                        template audit only (read-only, exit 0 = OK)
  check     OUT [--forces F]                parse an existing CP2K output exactly like 'work' (read-only;
                                            a multi-block forces file is read at its LAST block via mmap)
  prep      SET --template T --restart R --wfn W --mode test|full [--natoms N --must LIST] [--order-from SET0]
            [--h 0.01] [--center-grid] [--dry-run]
  work      SET [--hrt SEC]                 inside a job; exit 0 all done, 2 failure, 3 stopped/incomplete
  status    SET [--release-stale]           read-only unless --release-stale (renames, never deletes)
  test      SET                             verdict for --mode test; exit 0 only on PASS
  hess      SET                             only when ALL tasks of SET are complete and valid; exit 0 on success
  thermo    R.json TS.json [--T 310.15]     refuses mismatched or invalid inputs
  kick      SET --hess fd_SET/hess.json --name NAME [--amp 0.10]
  dimer     NAME --template GEO.inp --start S.restart --wfn W.wfn --a A.restart --b B.restart
  mf        T,RC,dE ...                     mean force along RC (K = 0.10 a.u.)
Nothing here deletes files. Existing outputs are never overwritten.
"""
import sys, os, re, json, math, time, shutil, hashlib, gzip, argparse, subprocess, glob, mmap

VERSION = "fd52v3.1"
NATOM = 95456
BOHR = 0.529177210903
HA2KCAL = 627.509474
FREQ = 5140.48714
K1 = 2.8591440e-3
RGAS = 1.98720425864e-3
C2 = 1.438776877
MASS = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "P": 30.973762,
        "S": 32.065, "MN": 54.938045, "MG": 24.305, "NA": 22.98977, "CL": 35.453}
AX = "xyz"
EXIT_FAIL, EXIT_STOP = 2, 3


def die(msg, code=1):
    sys.stderr.write("ERROR: " + msg + "\n"); sys.exit(code)


def rng(t):
    out = []
    for x in re.split(r"[,\s]+", str(t).strip()):
        if not x or x.upper() == "LIST": continue
        if ".." in x: a_, b_ = x.split(".."); out += list(range(int(a_), int(b_) + 1))
        else: out.append(int(x))
    return out


# ---------------- configuration ----------------
CFGP = os.environ.get("FD_CFG")
if not CFGP:
    die("FD_CFG is not set (fd52v3 has no built-in atom numbers)")
CFG = json.load(open(CFGP))
QM = rng(CFG["qm"])
if len(QM) != len(set(QM)) or len(QM) != CFG["nqm"]:
    die("QM list: %d entries, %d unique, config nqm %d" % (len(QM), len(set(QM)), CFG["nqm"]))
P, OW, O3P = CFG["P"], CFG["OW"], CFG["O3P"]
LINK_ATOM = CFG["LINK_ATOM"]
FIXED_FILE, PSF_FILE = CFG["fixed"], CFG["psf"]
SPIN_RANGE = tuple(CFG["spin_range"])
GATE = tuple(CFG.get("base_gate", [1.0e-3, 1.0e-3]))          # mobile, QM max force component
REPRO = tuple(CFG.get("repro_tol", [1.0e-7, 1.0e-6]))          # |dE| Ha, max|dF| Ha/bohr
NOISE_TOL = CFG.get("sym_noise_tol", 5.0e-5)          # even-order (anharmonic+noise) indicator, informative
ASYM_TOL = CFG.get("asym_tol", 0.02)                   # max|H-H^T| / max|H| : numerical-noise criterion
CFG_SHA = hashlib.sha256(open(CFGP, "rb").read()).hexdigest()
for q in (P, OW, O3P, LINK_ATOM):
    if q not in QM: die("atom %d of the config is not in the QM list" % q)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def fixed_set():
    s = set(rng(open(FIXED_FILE).read()))
    bad = [q for q in QM if q in s]
    if bad: die("QM atoms in fixed list %s: %s" % (FIXED_FILE, bad))
    return s


# ---------------- queue (never guesses) ----------------
def queue_job_ids():
    """returns set of job ids (possibly empty). Raises on any qstat failure or unknown format."""
    p = subprocess.run(["qstat", "-u", os.environ.get("USER", "")], stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True)
    if p.returncode != 0:
        raise RuntimeError("qstat exit %d: %s" % (p.returncode, p.stderr.strip()))
    out = p.stdout
    if not out.strip():
        return set()                       # Grid Engine prints nothing when the user has no jobs
    if "job-ID" not in out:
        raise RuntimeError("unexpected qstat output:\n" + out[:500])
    return set(l.split()[0] for l in out.splitlines()[2:] if l.split())


# ---------------- input parsing with section paths ----------------
def sections(text):
    """yield (path_tuple, line_index, stripped_line) for every line; path = open sections at that line"""
    stack = []
    for i, raw in enumerate(text.splitlines()):
        l = raw.split("!")[0].split("#")[0].strip() if not raw.strip().startswith("@") else raw.strip()
        if not l:
            yield tuple(stack), i, ""; continue
        if re.match(r"(?i)^&END\b", l):
            yield tuple(stack), i, l
            if not stack: raise ValueError("unbalanced &END at line %d" % (i + 1))
            stack.pop(); continue
        if l.startswith("&"):
            stack.append(l[1:].split()[0].upper())
            yield tuple(stack), i, l; continue
        yield tuple(stack), i, l
    if stack: raise ValueError("unclosed sections: %s" % stack)


def expand_includes(text, base=".", depth=0):
    if depth > 5: raise ValueError("include depth > 5")
    out = []
    for raw in text.splitlines():
        m = re.match(r"^\s*@INCLUDE\s+['\"]?([^'\"\s]+)['\"]?", raw, re.I)
        if m:
            f = os.path.join(base, m.group(1))
            if not os.path.isfile(f): raise ValueError("@INCLUDE file missing: " + f)
            out.append(expand_includes(open(f).read(), base, depth + 1))
        else:
            out.append(raw)
    return "\n".join(out)


def kw(text, path, key):
    """values of KEY lines whose section path ends with PATH (tuple)"""
    vals = []
    n = len(path)
    for pth, i, l in sections(text):
        if pth[-n:] == path and l and not l.startswith("&") and re.match(r"(?i)^%s\b" % key, l):
            vals.append((i, l.split(None, 1)[1] if len(l.split(None, 1)) > 1 else ""))
    return vals


def truthy(v):
    return v.strip().upper() in ("", "T", "TRUE", ".TRUE.", "ON", "YES")


def audit_template(path, center_grid=False):
    """returns rendered text with placeholders; dies with a list of problems"""
    raw = open(path).read()
    probs = []
    try:
        full = expand_includes(raw)
        list(sections(raw)); list(sections(full))
    except ValueError as e:
        die("template %s: %s" % (path, e))
    if re.search(r"(?mi)^\s*&(COLLECTIVE|RESTRAINT)\b", full): probs.append("&COLLECTIVE/&RESTRAINT present")
    rt = kw(full, ("GLOBAL",), "RUN_TYPE")
    if len(rt) != 1 or rt[0][1].split()[0].upper() != "ENERGY_FORCE": probs.append("GLOBAL/RUN_TYPE must be exactly one ENERGY_FORCE: %s" % rt)
    if len(kw(full, ("GLOBAL",), "PROJECT")) + len(kw(full, ("GLOBAL",), "PROJECT_NAME")) != 1: probs.append("exactly one GLOBAL/PROJECT required")
    inner = kw(full, ("FORCE_EVAL", "DFT", "SCF"), "EPS_SCF")
    outer = kw(full, ("FORCE_EVAL", "DFT", "SCF", "OUTER_SCF"), "EPS_SCF")
    for lab, v in (("SCF/EPS_SCF", inner), ("SCF/OUTER_SCF/EPS_SCF", outer)):
        if len(v) != 1: probs.append("%s must be given explicitly exactly once (found %d)" % (lab, len(v)))
        elif float(v[0][1].split()[0].replace("D", "E").replace("d", "e")) > 1.0001e-8: probs.append("%s = %s (> 1.0E-8)" % (lab, v[0][1]))
    g = kw(full, ("FORCE_EVAL", "DFT", "SCF"), "SCF_GUESS")
    if len(g) != 1 or g[0][1].split()[0].upper() != "RESTART": probs.append("SCF/SCF_GUESS RESTART required: %s" % g)
    if len(kw(full, ("FORCE_EVAL", "DFT"), "WFN_RESTART_FILE_NAME")) != 1: probs.append("exactly one DFT/WFN_RESTART_FILE_NAME required")
    er = dict((k, kw(full, ("EXT_RESTART",), k)) for k in ("RESTART_FILE_NAME", "RESTART_DEFAULT", "RESTART_POS", "RESTART_CONSTRAINT"))
    if len(er["RESTART_FILE_NAME"]) != 1: probs.append("exactly one EXT_RESTART/RESTART_FILE_NAME required")
    pos_ok = (er["RESTART_POS"] and truthy(er["RESTART_POS"][-1][1])) or (not er["RESTART_POS"] and (not er["RESTART_DEFAULT"] or truthy(er["RESTART_DEFAULT"][-1][1])))
    if not pos_ok: probs.append("EXT_RESTART does not restart positions")
    if not er["RESTART_CONSTRAINT"] or truthy(er["RESTART_CONSTRAINT"][-1][1]): probs.append("EXT_RESTART/RESTART_CONSTRAINT .FALSE. required explicitly")
    # QM atoms
    qmset = []
    for pth, i, l in sections(full):
        if len(pth) >= 2 and pth[-1] == "QM_KIND" and "QMMM" in pth and re.match(r"(?i)^MM_INDEX\b", l):
            qmset += rng(l.split(None, 1)[1])
    if sorted(qmset) != sorted(QM):
        probs.append("QMMM/QM_KIND MM_INDEX set (%d atoms) differs from config QM (%d): only-in-input %s only-in-config %s" % (
            len(qmset), len(QM), sorted(set(qmset) - set(QM))[:10], sorted(set(QM) - set(qmset))[:10]))
    # PSF
    cf = kw(full, ("FORCE_EVAL", "SUBSYS", "TOPOLOGY"), "CONN_FILE_NAME")
    if len(cf) != 1 or os.path.basename(cf[0][1].split()[0].strip("'\"")) != os.path.basename(PSF_FILE): probs.append("TOPOLOGY/CONN_FILE_NAME must be %s: %s" % (PSF_FILE, cf))
    # fixed atoms include must be the configured file
    fx_inc = [m.group(1) for m in re.finditer(r"@INCLUDE\s+['\"]?([^'\"\s]+)", raw, re.I)]
    for pth, i, l in sections(raw):
        if "FIXED_ATOMS" in pth and l.upper().startswith("@INCLUDE"):
            f = re.split(r"\s+", l, maxsplit=1)[1].strip("'\"")
            if os.path.basename(f) != os.path.basename(FIXED_FILE): probs.append("FIXED_ATOMS includes %s, config says %s" % (f, FIXED_FILE))
    if not any("FIXED_ATOMS" in pth for pth, _, _ in sections(raw)): probs.append("no MOTION/CONSTRAINT/FIXED_ATOMS section")
    # forces print section: FORCE_EVAL/PRINT/FORCES exactly once
    fsec = [(i, l) for pth, i, l in sections(raw) if pth == ("FORCE_EVAL", "PRINT", "FORCES") and l.startswith("&") and not l.upper().startswith("&END")]
    if len(fsec) != 1: probs.append("FORCE_EVAL/PRINT/FORCES section must appear exactly once (found %d)" % len(fsec))
    elif len(fsec[0][1].split()) > 1 and not truthy(fsec[0][1].split()[1]): probs.append("&FORCES is switched off")
    cg = kw(full, ("FORCE_EVAL", "QMMM"), "CENTER_GRID")
    if cg and truthy(cg[-1][1]) != center_grid: probs.append("template CENTER_GRID %s disagrees with --center-grid %s" % (cg[-1][1], center_grid))
    if probs:
        die("template %s rejected:\n  - " % path + "\n  - ".join(probs))
    # ---- render (on the un-expanded text) ----
    lines = raw.splitlines()
    fixes = {}
    for pth, i, l in sections(raw):
        if pth == ("GLOBAL",) and re.match(r"(?i)^PROJECT(_NAME)?\b", l): fixes[i] = "  PROJECT @PROJ@"
        if pth == ("EXT_RESTART",) and re.match(r"(?i)^RESTART_FILE_NAME\b", l): fixes[i] = "  RESTART_FILE_NAME @RST@"
        if pth[-2:] == ("FORCE_EVAL", "DFT") and re.match(r"(?i)^WFN_RESTART_FILE_NAME\b", l): fixes[i] = "    WFN_RESTART_FILE_NAME @WFN@"
        if pth == ("FORCE_EVAL", "PRINT", "FORCES") and re.match(r"(?i)^(FILENAME|FORMAT|NDIGITS|UNIT)\b", l): fixes[i] = None
    for k in ("PROJECT", "RESTART_FILE_NAME", "WFN_RESTART_FILE_NAME"):
        pass
    out = []
    for i, l in enumerate(lines):
        if i in fixes:
            if fixes[i] is not None: out.append(fixes[i])
            continue
        out.append(l)
        if i == fsec[0][0]:
            out.append("      FILENAME =@FRC@")
            out.append("      NDIGITS 12")
    s = "\n".join(out) + "\n"
    if center_grid and not cg:
        s, c = re.subn(r"(?mi)^([ \t]*&QMMM[ \t]*)$", r"\1\n    CENTER_GRID .TRUE.", s, count=1)
        if c != 1: die("&QMMM line not found for CENTER_GRID")
    for ph in ("@PROJ@", "@RST@", "@WFN@", "@FRC@"):
        if s.count(ph) != 1: die("internal: placeholder %s appears %d times" % (ph, s.count(ph)))
    return s


# ---------------- restart &COORD ----------------
def coord_block(path):
    with open(path) as f: L = f.readlines()
    st = [i for i, l in enumerate(L) if re.match(r"^\s*&COORD\b", l, re.I)]
    if len(st) != 1: die("%s: &COORD sections = %d" % (path, len(st)))
    i0 = st[0] + 1; unit = 1.0
    while True:
        t = L[i0].split()
        if not t: i0 += 1; continue
        k = t[0].upper()
        if k == "UNIT":
            u = t[1].lower()
            if u in ("angstrom", "ang"): unit = 1.0
            elif u in ("bohr", "au", "a.u."): unit = BOHR
            else: die("unknown COORD UNIT " + u)
            i0 += 1; continue
        if k == "SCALED":
            if len(t) == 1 or truthy(t[1]): die("SCALED coordinates not supported")
            i0 += 1; continue
        break
    n = 0
    while not re.match(r"^\s*&END\b", L[i0 + n], re.I): n += 1
    if n != NATOM: die("%s: %d atom lines in &COORD (expected %d)" % (path, n, NATOM))
    return L, i0, unit


def atom_xyz(line):
    t = line.split(); return [float(t[1]), float(t[2]), float(t[3])]


def displaced_line(line, axis, dA):
    m = re.match(r"^(\s*\S+)(\s+)(\S+)(\s+)(\S+)(\s+)(\S+)(.*)$", line.rstrip("\n"))
    if not m: die("cannot parse COORD line: " + line)
    g = list(m.groups()); xyz = [float(g[2]), float(g[4]), float(g[6])]; xyz[axis] += dA
    g[2], g[4], g[6] = ["%.15E" % v for v in xyz]
    return "".join(g) + "\n"


def all_xyz(path):
    L, i0, u = coord_block(path)
    return [[c * u for c in atom_xyz(L[i0 + k])] for k in range(NATOM)]


# ---------------- output parsing (strict) ----------------
HDR = "ATOMIC FORCES in [a.u.]"
ROW = re.compile(r"^\s*(\d+)\s+(\d+)\s+([A-Za-z]{1,3})\s+(\S+)\s+(\S+)\s+(\S+)\s*$")


def parse_force_block(block):
    """block starts at the header. Strict: all NATOM ids 1..NATOM once, in order, finite, SUM line present."""
    lines = block.splitlines()
    if not lines or HDR not in lines[0]: return None, "header '%s' missing" % HDR
    F, why, nxt, ended = {}, [], 1, False
    for l in lines[1:]:
        if not l.strip() or l.lstrip().startswith("#"): continue
        if "SUM OF ATOMIC FORCES" in l: ended = True; break
        m = ROW.match(l)
        if not m: why.append("unparsable row: %r" % l[:80]); break
        a = int(m.group(1))
        if a != nxt: why.append("atom id %d where %d expected (duplicate/missing/out of order)" % (a, nxt)); break
        try:
            v = [float(m.group(j)) for j in (4, 5, 6)]
        except ValueError:
            why.append("non-numeric force at atom %d" % a); break
        if not all(math.isfinite(x) for x in v): why.append("non-finite force at atom %d" % a); break
        F[a] = v; nxt += 1
    if not why and not ended: why.append("table not terminated by SUM OF ATOMIC FORCES")
    if not why and len(F) != NATOM: why.append("%d rows (expected %d)" % (len(F), NATOM))
    return (F if not why else None), "; ".join(why)


def last_force_block_file(path):
    """mmap the file, take the LAST header; returns (block_text, n_blocks)"""
    with open(path, "rb") as f:
        if os.path.getsize(path) == 0: return None, 0
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            pos = mm.rfind(HDR.encode())
            if pos < 0: return None, 0
            n, p = 0, mm.find(HDR.encode())
            while p >= 0 and n < 100000: n += 1; p = mm.find(HDR.encode(), p + 1)
            end = mm.find(b"SUM OF ATOMIC FORCES", pos)
            stop = mm.find(b"\n", end) if end >= 0 else len(mm)
            return mm[pos:(stop if stop >= 0 else len(mm))].decode(errors="replace"), n
        finally:
            mm.close()


E_RE = re.compile(r"ENERGY\| Total FORCE_EVAL \( QMMM \) energy \[a\.u\.\]:\s*(-?[\d.]+(?:[EeDd][-+]?\d+)?)")


def parse_out(text):
    e = E_RE.findall(text)
    s = re.findall(r"SCF run converged in\s+(\d+)\s+steps", text)
    sp = re.findall(r"Integrated absolute spin density\s*:\s*(-?[\d.]+)", text)
    return dict(ended="PROGRAM ENDED AT" in text, notconv=text.count("SCF run NOT converged"),
                nscf=int(s[-1]) if s else -1, nE=len(e), E=float(e[-1]) if e else None,
                spin=float(sp[-1]) if sp else None, walltime_hit="exceeded requested execution time" in text)


def evaluate_output(out_path, frc_path, fixed):
    """strict judgement used by work/check: returns (record dict, list_of_problems)"""
    why = []
    text = open(out_path, errors="replace").read() if os.path.exists(out_path) else ""
    r = parse_out(text)
    if not r["ended"]: why.append("no PROGRAM ENDED")
    if r["walltime_hit"]: why.append("WALLTIME exceeded")
    if r["notconv"]: why.append("SCF NOT converged x%d" % r["notconv"])
    if r["nscf"] < 0: why.append("no 'SCF run converged' line")
    if r["nE"] != 1: why.append("%d energy lines (expected exactly 1 for ENERGY_FORCE)" % r["nE"])
    if r["spin"] is None: why.append("no spin density line")
    elif not (SPIN_RANGE[0] <= r["spin"] <= SPIN_RANGE[1]): why.append("spin %.4f outside %s" % (r["spin"], SPIN_RANGE))
    F = None
    if not os.path.exists(frc_path): why.append("forces file %s not written" % frc_path)
    else:
        blk, nb = last_force_block_file(frc_path)
        if nb != 1: why.append("forces file has %d blocks (expected 1)" % nb)
        if blk is None: why.append("no force block")
        else:
            F, w = parse_force_block(blk)
            if F is None: why.append("forces: " + w)
    if F is not None:
        nz = max(abs(c) for q in fixed for c in F[q]) if fixed else 0.0
        r["fmax_qm"] = max(abs(c) for q in QM for c in F[q])
        r["fmax_mob"] = max(abs(c) for q in F if q not in fixed for c in F[q])
        r["fmax_fixed"] = nz
        r["F_qm"] = [F[q][k] for q in QM for k in range(3)]
    return r, why


# ---------------- prep ----------------
def setdir(name): return "fd_" + name


def tag_of(name): return re.sub(r"[^A-Za-z0-9]", "", name)


def cmd_validate(a):
    audit_template(a.template, a.center_grid); fixed_set()
    print("template OK:", a.template, " config", CFGP, CFG_SHA[:12]); return 0


def cmd_prep(a):
    sd = setdir(a.set)
    if os.path.exists(sd): die(sd + " already exists (never overwritten)")
    for other in glob.glob("fd_*/manifest.json"):
        if tag_of(json.load(open(other))["set"]) == tag_of(a.set): die("project tag collision with " + other)
    for f in (a.template, a.restart, a.wfn, FIXED_FILE, PSF_FILE):
        if not os.path.isfile(f): die("missing " + f)
    body = audit_template(a.template, a.center_grid)
    fixed_set()
    L, i0, unit = coord_block(a.restart)
    xyz = {q: [c * unit for c in atom_xyz(L[i0 + q - 1])] for q in QM}
    h = a.h
    tasks = [(0, 0, -1, 0, 0.0, "base")]
    def add(atom, ax, hh, tag):
        for sg in (+1, -1): tasks.append((len(tasks), atom, ax, sg, hh, tag))
    order = None
    if a.mode == "test":
        if a.natoms or a.must or a.order_from: die("--natoms/--must/--order-from are for --mode full")
        v = [xyz[O3P][k] - xyz[P][k] for k in range(3)]; axP = max(range(3), key=lambda k: abs(v[k]))
        add(P, axP, h, "P_bondaxis"); add(O3P, axP, h, "O3p_bondaxis"); add(P, axP, h / 2, "P_half")
        cen = [sum(xyz[q][k] for q in QM) / len(QM) for k in range(3)]
        axL = max(range(3), key=lambda k: abs(xyz[LINK_ATOM][k] - cen[k]))
        add(LINK_ATOM, axL, h, "link_atom")
        ext = [(fn(QM, key=lambda q: xyz[q][k]), k) for k in range(3) for fn in (min, max)]
        ext = [e for e in ext if e[0] not in (P, O3P, LINK_ATOM)] or ext
        add(ext[0][0], ext[0][1], h, "extremal_center")
        tasks.append((len(tasks), 0, -1, 0, 0.0, "repeat"))      # independent repeat of the base point
    else:
        if a.order_from:
            if a.natoms or a.must: die("--order-from fixes the atom set; do not combine with --natoms/--must")
            m0 = json.load(open(os.path.join(setdir(a.order_from), "manifest.json")))
            if m0["cfg_sha"] != CFG_SHA: die("config differs from %s" % a.order_from)
            if abs(m0["h_bohr"] - h) > 1e-15 or m0["center_grid"] != a.center_grid: die("h/center_grid differ from %s" % a.order_from)
            if m0["rendered_sha"] != hashlib.sha256(body.encode()).hexdigest(): die("rendered template differs from %s" % a.order_from)
            order = m0["order"]
        else:
            if not a.natoms: die("--mode full needs --natoms N (and normally --must) or --order-from")
            mid = [(xyz[P][k] + xyz[OW][k]) / 2 for k in range(3)]
            byd = sorted(QM, key=lambda q: sum((xyz[q][k] - mid[k]) ** 2 for k in range(3)))
            must = rng(a.must or "")
            miss = [q for q in must if q not in QM]
            if miss: die("--must atoms not in QM: %s" % miss)
            if len(set(must)) > a.natoms: die("--must has more atoms than --natoms")
            sel = list(dict.fromkeys(must))
            for q in byd:
                if len(sel) >= a.natoms: break
                if q not in sel: sel.append(q)
            order = [q for q in byd if q in sel]
        for q in order:
            for k in range(3): add(q, k, h, "full")
    print("restart:", a.restart, " unit->A", unit, " atom lines", NATOM)
    for q in (P, O3P, OW): print("  atom %6d: %s" % (q, L[i0 + q - 1].rstrip()))
    if order: print("atom set (%d, displacement order): %s" % (len(order), order))
    print("tasks:", len(tasks), " h(bohr)", h, " mode", a.mode, " CENTER_GRID", a.center_grid)
    if a.dry_run: print("dry-run: nothing written"); return 0
    for d in ("", "inp", "out", "res", "claims"): os.makedirs(os.path.join(sd, d))
    shutil.copy2(a.restart, os.path.join(sd, "base.restart")); shutil.copy2(a.wfn, os.path.join(sd, "base.wfn"))
    os.chmod(os.path.join(sd, "base.restart"), 0o444); os.chmod(os.path.join(sd, "base.wfn"), 0o444)
    open(os.path.join(sd, "template.inp"), "w").write(body)
    with open(os.path.join(sd, "tasks.tsv"), "w") as f:
        for t in tasks: f.write("%d\t%d\t%d\t%d\t%.10f\t%s\n" % t)
    man = dict(version=VERSION, set=a.set, created=time.strftime("%F %T"), cfg=CFGP, cfg_sha=CFG_SHA,
               template=a.template, template_sha=sha(a.template), rendered_sha=hashlib.sha256(body.encode()).hexdigest(),
               restart=a.restart, restart_sha=sha(a.restart), wfn=a.wfn, wfn_sha=sha(a.wfn), h_bohr=h, mode=a.mode,
               center_grid=a.center_grid, ntasks=len(tasks), order=order, qm=QM, qm_xyz_A={str(q): xyz[q] for q in QM},
               gate=(tuple(CFG.get("path_gate", [0.05, 0.05])) if a.path_point else GATE), path_point=a.path_point, spin_range=SPIN_RANGE)
    json.dump(man, open(os.path.join(sd, "manifest.json"), "w"), indent=1)
    print("written:", sd); return 0


# ---------------- work ----------------
def load_tasks(sd):
    T = []
    for l in open(os.path.join(sd, "tasks.tsv")):
        i, at, ax, sg, hh, tag = l.rstrip("\n").split("\t"); T.append((int(i), int(at), int(ax), int(sg), float(hh), tag))
    return T


def write_atomic(path, text):
    tmp = path + ".tmp%d" % os.getpid()
    with open(tmp, "w") as f: f.write(text)
    os.rename(tmp, path)


def cmd_work(a):
    sd = setdir(a.set)
    man = json.load(open(os.path.join(sd, "manifest.json")))
    if man["cfg_sha"] != CFG_SHA: die("FD_CFG differs from the one used at prep")
    body = open(os.path.join(sd, "template.inp")).read()
    if hashlib.sha256(body.encode()).hexdigest() != man["rendered_sha"]: die("template.inp changed after prep")
    if sha(os.path.join(sd, "base.restart")) != man["restart_sha"] or sha(os.path.join(sd, "base.wfn")) != man["wfn_sha"]:
        die("base.restart/base.wfn changed after prep")
    L, i0, unit = coord_block(os.path.join(sd, "base.restart"))
    fx = fixed_set()
    gate = tuple(man["gate"])
    cmd = os.environ.get("CP2K_CMD", "mpirun -np 4 cp2k.psmp")
    job = os.environ.get("JOB_ID", "nojob") + "." + os.environ.get("SGE_TASK_ID", "0")
    t_start = time.time(); worst = 0.0; nfail = 0; stopped = False
    tag0 = tag_of(a.set)
    PASS, FAIL, STOP = (os.path.join(sd, x) for x in ("BASE_PASS", "BASE_FAIL", "STOP"))
    for (i, at, ax, sg, hh, tag) in load_tasks(sd):
        rid = "d%04d" % i
        if os.path.exists(STOP) or os.path.exists(FAIL):
            print("STOP/BASE_FAIL present -> stop"); stopped = True; break
        if os.path.exists(os.path.join(sd, "res", rid + ".txt")):
            if i == 0 and not os.path.exists(PASS) and not os.path.exists(FAIL):   # crash between result and verdict
                d0 = read_res(sd)[0]; fm, fq = float(d0["fmax_mob"]), float(d0["fmax_qm"])
                msg = "base (re-judged): mobile %.3e (gate %.1e) QM %.3e (gate %.1e)" % (fm, gate[0], fq, gate[1])
                write_atomic(FAIL if (fm > gate[0] or fq > gate[1]) else PASS, msg + "\n")
                if os.path.exists(FAIL): print("BASE_FAIL written"); break
            continue
        if i > 0:
            while not os.path.exists(PASS):
                if os.path.exists(FAIL) or os.path.exists(STOP): break
                if time.time() - t_start > a.hrt - 900: print("time guard while waiting for BASE_PASS"); stopped = True; break
                time.sleep(20)
            if stopped: break
            if not os.path.exists(PASS): print("base not passed -> stop"); stopped = True; break
        if time.time() - t_start + max(1.5 * worst, 300) > a.hrt - 300:
            print("time guard: stop before", rid); stopped = True; break
        try: os.mkdir(os.path.join(sd, "claims", rid))
        except FileExistsError: continue
        write_atomic(os.path.join(sd, "claims", rid, "info"), "%s %s %s\n" % (job, os.uname()[1], time.strftime("%F %T")))
        proj = "fd%s_%s" % (tag0, rid)
        out = os.path.join(sd, "out", proj + ".out"); frc = os.path.join(sd, "out", proj + ".forces")
        if os.path.exists(out) or os.path.exists(frc) or glob.glob(proj + "-*"):
            write_atomic(os.path.join(sd, "res", rid + ".fail"), "stale files exist for %s\n" % proj); nfail += 1; continue
        rst = os.path.join(sd, "inp", proj + ".restart")
        with open(rst, "w") as f:
            for k, l in enumerate(L):
                f.write(displaced_line(l, ax, sg * hh * BOHR / unit) if (at and k == i0 + at - 1) else l)
        wfn = os.path.join(sd, "inp", proj + "_in.wfn")
        shutil.copy2(os.path.join(sd, "base.wfn"), wfn); os.chmod(wfn, 0o644)
        inp = os.path.join(sd, "inp", proj + ".inp")
        open(inp, "w").write(body.replace("@PROJ@", proj).replace("@RST@", rst).replace("@WFN@", wfn).replace("@FRC@", frc))
        t0 = time.time()
        with open(os.path.join(sd, "out", proj + ".launcher.log"), "w") as lg:
            rc = subprocess.call(cmd.split() + ["-i", inp, "-o", out], stdout=lg, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
        wall = time.time() - t0; worst = max(worst, wall)
        for g in glob.glob(proj + "-*"): shutil.move(g, os.path.join(sd, "out", os.path.basename(g)))
        r, why = evaluate_output(out, frc, fx)
        if rc != 0: why.insert(0, "exit %d" % rc)
        if why:
            write_atomic(os.path.join(sd, "res", rid + ".fail"), "; ".join(why) + "\n"); nfail += 1
            print(rid, "FAIL", why, flush=True)
            if i == 0: write_atomic(FAIL, "base evaluation failed: %s\n" % "; ".join(why)); break
            continue
        head = "id=%d atom=%d axis=%d sign=%d h=%.10f tag=%s E=%.12f nscf=%d spin=%.6f wall=%.1f job=%s fmax_qm=%.6e fmax_mob=%.6e fmax_fixed=%.3e" % (
            i, at, ax, sg, hh, tag, r["E"], r["nscf"], r["spin"], wall, job, r["fmax_qm"], r["fmax_mob"], r["fmax_fixed"])
        write_atomic(os.path.join(sd, "res", rid + ".txt"), head + "\n" + " ".join("%.12e" % v for v in r["F_qm"]) + "\n")
        for big in (rst, frc):
            with open(big, "rb") as fi, gzip.open(big + ".gz", "wb") as fo: shutil.copyfileobj(fi, fo)
            os.remove(big)
        print(rid, "ok E=%.10f nscf=%d wall=%.0fs" % (r["E"], r["nscf"], wall), flush=True)
        if i == 0:
            msg = "base: max force component mobile %.3e (gate %.1e)  QM %.3e (gate %.1e)" % (r["fmax_mob"], gate[0], r["fmax_qm"], gate[1])
            print(msg, flush=True)
            if r["fmax_mob"] > gate[0] or r["fmax_qm"] > gate[1]:
                write_atomic(FAIL, "gate failed: " + msg + "\n"); print("BASE_FAIL written"); break
            write_atomic(PASS, msg + "\n")
    if nfail or os.path.exists(FAIL): return EXIT_FAIL
    if stopped: return EXIT_STOP
    T = load_tasks(sd)
    if all(os.path.exists(os.path.join(sd, "res", "d%04d.txt" % t[0])) for t in T): return 0
    return EXIT_STOP            # other workers still running or tasks left


# ---------------- status ----------------
def read_res(sd):
    R = {}
    for p in glob.glob(os.path.join(sd, "res", "d*.txt")):
        l1, l2 = open(p).read().splitlines()[:2]
        d = dict(kv.split("=", 1) for kv in l1.split()); d["F"] = [float(x) for x in l2.split()]
        if len(d["F"]) != 3 * len(QM): die("%s: %d force values (expected %d)" % (p, len(d["F"]), 3 * len(QM)))
        R[int(d["id"])] = d
    return R


def cmd_status(a):
    sd = setdir(a.set); T = load_tasks(sd); R = read_res(sd)
    fails = sorted(glob.glob(os.path.join(sd, "res", "*.fail")))
    claims = [os.path.basename(c) for c in glob.glob(os.path.join(sd, "claims", "d*")) if os.path.isdir(c) and ".released" not in c]
    try:
        running = queue_job_ids(); qok = True
    except Exception as e:
        print("qstat unavailable (%s): stale detection skipped" % e); running, qok = set(), False
    stale = []
    for c in claims:
        rid = int(c[1:])
        if rid in R or os.path.exists(os.path.join(sd, "res", c + ".fail")): continue
        job = open(os.path.join(sd, "claims", c, "info")).read().split()[0].split(".")[0]
        if qok and job not in running: stale.append(c)
    for flag in ("BASE_PASS", "BASE_FAIL", "STOP"):
        if os.path.exists(os.path.join(sd, flag)): print("  %s: %s" % (flag, open(os.path.join(sd, flag)).read().strip()))
    print("%s: tasks %d  done %d  failed %d  claimed-unfinished %d  stale %d" % (
        sd, len(T), len(R), len(fails), len([c for c in claims if int(c[1:]) not in R]), len(stale)))
    if R:
        w = [float(d["wall"]) for d in R.values()]
        print("  wall per eval: mean %.0f s  max %.0f s ; SCF steps %d-%d ; spin %.3f-%.3f" % (
            sum(w) / len(w), max(w), min(int(d["nscf"]) for d in R.values()), max(int(d["nscf"]) for d in R.values()),
            min(float(d["spin"]) for d in R.values()), max(float(d["spin"]) for d in R.values())))
    for f in fails: print("  FAIL", os.path.basename(f), open(f).read().strip())
    for c in stale: print("  STALE", c)
    if a.release_stale:
        if not qok: die("cannot release without a working qstat")
        for c in stale + [os.path.basename(f)[:-5] for f in fails]:
            ts = time.strftime("%H%M%S")
            src = os.path.join(sd, "claims", c)
            if os.path.isdir(src): os.rename(src, src + ".released." + ts)
            ff = os.path.join(sd, "res", c + ".fail")
            if os.path.exists(ff): os.rename(ff, ff + "." + ts)
            for g in glob.glob(os.path.join(sd, "out", "fd*_%s*" % c)) + glob.glob(os.path.join(sd, "inp", "fd*_%s*" % c)):
                os.rename(g, g + ".attempt" + ts)
            print("  released", c)
    return 0


# ---------------- Hessian helpers ----------------
def columns(R, T):
    F0 = R[0]["F"]; cols, sym, pairs = {}, {}, {}
    for (i, at, ax, sg, hh, tag) in T:
        if i and at and i in R: pairs.setdefault((at, ax, round(hh, 10)), {})[sg] = R[i]["F"]
    for key, d in pairs.items():
        if 1 in d and -1 in d:
            hh = key[2]
            cols[key] = [-(d[1][r] - d[-1][r]) / (2 * hh) for r in range(len(F0))]
            sym[key] = max(abs((d[1][r] + d[-1][r]) / 2 - F0[r]) for r in range(len(F0)))
    return cols, sym


def ridx(atom, ax): return QM.index(atom) * 3 + ax


def completeness(sd, T, R):
    probs = []
    miss = [t[0] for t in T if t[0] not in R]
    if miss: probs.append("missing results %s%s" % (miss[:10], " ..." if len(miss) > 10 else ""))
    fails = glob.glob(os.path.join(sd, "res", "*.fail"))
    if fails: probs.append("%d .fail files" % len(fails))
    for f in ("STOP", "BASE_FAIL"):
        if os.path.exists(os.path.join(sd, f)): probs.append(f + " present")
    if not os.path.exists(os.path.join(sd, "BASE_PASS")): probs.append("BASE_PASS absent")
    bad = [k for k, d in R.items() if not (SPIN_RANGE[0] <= float(d["spin"]) <= SPIN_RANGE[1])]
    if bad: probs.append("spin outside range in %s" % bad)
    return probs


def cmd_test(a):
    sd = setdir(a.set); T = load_tasks(sd); R = read_res(sd); man = json.load(open(os.path.join(sd, "manifest.json")))
    if man["mode"] != "test": die("not a test set")
    probs = completeness(sd, T, R)
    if 0 not in R:
        print("VERDICT: FAIL (base not finished)"); return EXIT_FAIL
    h = man["h_bohr"]; cols, sym = columns(R, T)
    print("E0 = %.12f  nscf %s  base max force: QM %s mobile %s" % (float(R[0]["E"]), R[0]["nscf"], R[0]["fmax_qm"], R[0]["fmax_mob"]))
    walls = [float(d["wall"]) for k, d in R.items() if k]
    if walls: print("wall per displaced evaluation: mean %.0f s, max %.0f s" % (sum(walls) / len(walls), max(walls)))
    need = {"P_bondaxis", "O3p_bondaxis", "P_half", "link_atom", "extremal_center", "repeat"}
    have = set(t[5] for t in T if t[0] in R)
    if not need <= have: probs.append("test points missing: %s" % sorted(need - have))
    flagged = []
    print("%-24s %12s %12s" % ("column (atom,axis,h)", "max|sym|", "H_diag"))
    for key in sorted(cols):
        print("%-24s %12.3e %12.5f" % ("%d,%s,%.4f" % (key[0], AX[key[1]], key[2]), sym[key], cols[key][ridx(key[0], key[1])]))
        if sym[key] > NOISE_TOL: flagged.append(key)
    full = [k for k in cols if k[0] == P and abs(k[2] - h) < 1e-12]
    half = [k for k in cols if k[0] == P and abs(k[2] - h / 2) < 1e-12]
    o = [k for k in cols if k[0] == O3P and abs(k[2] - h) < 1e-12]
    if flagged and not (full and half): probs.append("cannot classify large even-order deviation without the P h/2 column")
    if full and half:
        cf, ch = cols[full[0]], cols[half[0]]
        rel = max(abs(x - y) for x, y in zip(cf, ch)) / max(abs(x) for x in cf)
        print("P column h vs h/2: max rel diff %.3e ; sym(h/2)/sym(h) = %.2f (0.25 curvature, ~1 noise)" % (rel, sym[half[0]] / max(sym[full[0]], 1e-30)))
        if rel > 0.02: probs.append("h vs h/2 rel diff %.3e > 0.02" % rel)
        ratio = sym[half[0]] / max(sym[full[0]], 1e-30)
        if flagged and ratio > 0.6:
            probs.append("even-order deviation does not scale as h^2 (ratio %.2f): numerical noise in %s" % (ratio, flagged))
        elif flagged:
            print("INFO: |(F+ + F-)/2 - F0| above %.0e at %s; h^2 scaling (ratio %.2f) shows anharmonicity, not noise (central differences unaffected)" % (NOISE_TOL, flagged, ratio))
    if full and o:
        ax = full[0][1]; hpo = cols[full[0]][ridx(O3P, ax)]; hop = cols[o[0]][ridx(P, ax)]
        asym = abs(hpo - hop) / (0.5 * (abs(hpo) + abs(hop)))
        print("H[P,O3'] = %.6f  H[O3',P] = %.6f  rel asym %.3e" % (hpo, hop, asym))
        if asym > 0.03: probs.append("Hessian asymmetry %.3e > 0.03" % asym)
    rep = [t[0] for t in T if t[5] == "repeat" and t[0] in R]
    if rep:
        d = R[rep[0]]
        dE = abs(float(d["E"]) - float(R[0]["E"])); dF = max(abs(x - y) for x, y in zip(d["F"], R[0]["F"]))
        print("repeat of base: |dE| %.3e Ha  max|dF| %.3e Ha/bohr (tol %.0e / %.0e)" % (dE, dF, REPRO[0], REPRO[1]))
        if dE > REPRO[0] or dF > REPRO[1]: probs.append("base not reproducible")
    for p in probs: print("  PROBLEM:", p)
    if probs: print("VERDICT: FAIL"); return EXIT_FAIL
    print("VERDICT: PASS"); return 0


def jacobi(A):
    n = len(A); A = [r[:] for r in A]; V = [[float(i == j) for j in range(n)] for i in range(n)]
    for sweep in range(100):
        if sum(A[i][j] ** 2 for i in range(n) for j in range(i + 1, n)) < 1e-22: break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(A[p][q]) < 1e-14: continue
                th = (A[q][q] - A[p][p]) / (2 * A[p][q])
                t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1)); c = 1 / math.sqrt(t * t + 1); s = t * c
                for k in range(n):
                    akp, akq = A[k][p], A[k][q]; A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = A[p][k], A[q][k]; A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k][p], V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq
    w = [A[i][i] for i in range(n)]; idx = sorted(range(n), key=lambda i: w[i])
    return [w[i] for i in idx], [[V[k][i] for k in range(n)] for i in idx]


def eig(H):
    try:
        import numpy as np
        w, v = np.linalg.eigh(np.array(H)); return [float(x) for x in w], [[float(t) for t in v[:, i]] for i in range(len(w))]
    except ImportError:
        return jacobi(H)


def cmd_hess(a):
    sd = setdir(a.set); T = load_tasks(sd); R = read_res(sd); man = json.load(open(os.path.join(sd, "manifest.json")))
    if man["mode"] != "full": die("hess needs a --mode full set")
    probs = completeness(sd, T, R)
    if len(T) != 1 + 6 * len(man["order"]): probs.append("task list size %d != 1+6*%d" % (len(T), len(man["order"])))
    if probs: die("hess refused for %s:\n  - " % sd + "\n  - ".join(probs), EXIT_FAIL)
    h = man["h_bohr"]; cols, sym = columns(R, T); A = man["order"]
    for q in A:
        for k in range(3):
            if (q, k, round(h, 10)) not in cols: die("column %d/%s missing" % (q, AX[k]), EXIT_FAIL)
    noise = max(sym[(q, k, round(h, 10))] for q in A for k in range(3))   # informative (anharmonic + noise)
    psf = open(PSF_FILE).read().splitlines()
    i = next(k for k, l in enumerate(psf) if "!NATOM" in l)
    el = {}
    for l in psf[i + 1:i + 1 + NATOM]:
        t = l.split()
        if int(t[0]) in A: el[int(t[0])] = (t[4], float(t[7]))
    m = []
    for q in A:
        name, pm = el[q]; e = "MN" if name.upper().startswith("MN") else name[0].upper(); sm = MASS[e]
        if abs(pm - sm) > 0.05: print("WARN atom %d %s PSF mass %.4f vs standard %.4f (standard used)" % (q, name, pm, sm))
        m += [sm] * 3
    n = 3 * len(A); idx = [ridx(q, k) for q in A for k in range(3)]
    H = [[cols[(A[j // 3], j % 3, round(h, 10))][idx[i_]] for j in range(n)] for i_ in range(n)]
    asym = max(abs(H[x][y] - H[y][x]) for x in range(n) for y in range(x + 1, n)); hmax = max(abs(H[x][y]) for x in range(n) for y in range(n))
    if asym / hmax > ASYM_TOL: die("Hessian asymmetry %.2e (relative %.3f) > %.3f" % (asym, asym / hmax, ASYM_TOL), EXIT_FAIL)
    Hm = [[0.5 * (H[x][y] + H[y][x]) / math.sqrt(m[x] * m[y]) for y in range(n)] for x in range(n)]
    proj = None
    if a.project_path:
        # reaction-path projection at a non-stationary path point: remove the mass-weighted energy-gradient direction
        g = [-R[0]["F"][idx[i_]] / math.sqrt(m[i_]) for i_ in range(n)]
        gn_ = math.sqrt(sum(t * t for t in g))
        if gn_ < 1e-5: die("energy gradient at the base point is ~0 (%.2e): projection undefined (stationary point? use normal hess)" % gn_, EXIT_FAIL)
        u = [t / gn_ for t in g]
        Pm = [[(1.0 if x == y else 0.0) - u[x] * u[y] for y in range(n)] for x in range(n)]
        PH = [[sum(Pm[x][k] * Hm[k][y] for k in range(n)) for y in range(n)] for x in range(n)]
        Hm = [[sum(PH[x][k] * Pm[k][y] for k in range(n)) for y in range(n)] for x in range(n)]
        proj = dict(grad_norm_mw=gn_, base_fmax_qm=float(R[0]["fmax_qm"]))
    w, V = eig(Hm)
    if proj:
        jz = max(range(n), key=lambda j: abs(sum(V[j][i_] * u[i_] for i_ in range(n))))
        proj["dropped_eigenvalue"] = w[jz]; proj["dropped_overlap"] = abs(sum(V[jz][i_] * u[i_] for i_ in range(n)))
        w = [x for j, x in enumerate(w) if j != jz]; V = [v for j, v in enumerate(V) if j != jz]; n_modes = n - 1
    nu = [math.copysign(FREQ * math.sqrt(abs(x)), x) for x in w]
    xyz = man["qm_xyz_A"]
    def unit(v): s = math.sqrt(sum(x * x for x in v)); return [x / s for x in v]
    uPW = unit([xyz[str(P)][k] - xyz[str(OW)][k] for k in range(3)]); uPO = unit([xyz[str(P)][k] - xyz[str(O3P)][k] for k in range(3)])
    g = {P: [uPW[k] - uPO[k] for k in range(3)], OW: [-uPW[k] for k in range(3)], O3P: [uPO[k] for k in range(3)]}
    gv = [g.get(q, [0, 0, 0])[k] for q in A for k in range(3)]; gn = math.sqrt(sum(t * t for t in gv))
    modes = []
    for j in range(min(6, len(w))):
        x = [V[j][i_] / math.sqrt(m[i_]) for i_ in range(n)]; s = math.sqrt(sum(t * t for t in x)); x = [t / s for t in x]
        cos = sum(x[i_] * gv[i_] for i_ in range(n)) / gn if gn else 0.0
        wt = sorted(((sum(x[3 * a_ + k] ** 2 for k in range(3)), A[a_]) for a_ in range(len(A))), reverse=True)[:5]
        modes.append(dict(nu=nu[j], cos_RC=cos, top=[(q, round(v, 3)) for v, q in wt], vec=x))
    res = dict(version=VERSION, set=a.set, cfg_sha=CFG_SHA, E0=float(R[0]["E"]), atoms=A, h=h, center_grid=man["center_grid"],
               rendered_sha=man["rendered_sha"], restart_sha=man["restart_sha"], asym_max=asym, H_max=hmax, sym_noise_max=noise,
               base_fmax_qm=float(R[0]["fmax_qm"]), base_fmax_mob=float(R[0]["fmax_mob"]),
               nu=nu, n_negative=sum(1 for x in nu if x < 0), lowest_modes=modes, projected=bool(proj), projection=proj,
               asym_rel=asym / hmax)
    outp = os.path.join(sd, "hess_path.json" if proj else "hess.json")
    if os.path.exists(outp): die(outp + " exists (not overwritten)")
    json.dump(res, open(outp, "w"), indent=1)
    print("%s: atoms %d  max|H-H^T| %.2e (relative %.4f; max|H| %.3f)  even-order max %.2e  negative eigenvalues %d%s" % (
        a.set, len(A), asym, asym / hmax, hmax, noise, res["n_negative"], "  [path-projected: 1 mode removed, overlap %.3f, eigenvalue %.2e]" % (proj["dropped_overlap"], proj["dropped_eigenvalue"]) if proj else ""))
    for md in modes: print("  nu %9.1f cm-1  cos(RC) %+.2f  top %s" % (md["nu"], md["cos_RC"], md["top"]))
    print("written", outp); return 0


def gvib(nus, T, cut):
    kT = RGAS * T; G = Z = Hv = S = 0.0
    for v in nus:
        v = max(v, cut) if cut else v; x = C2 * v / T
        Z += K1 * v / 2; G += K1 * v / 2 + kT * math.log(1 - math.exp(-x))
        Hv += K1 * v / 2 + K1 * v / (math.exp(x) - 1); S += RGAS * (x / (math.exp(x) - 1) - math.log(1 - math.exp(-x)))
    return G, Z, Hv, S


def cmd_thermo(a):
    r = json.load(open(a.r)); t = json.load(open(a.ts))
    for k in ("version", "cfg_sha", "atoms", "h", "center_grid", "rendered_sha"):
        if r.get(k) != t.get(k): die("R and TS differ in %s" % k, EXIT_FAIL)
    if r["n_negative"]: die("reactant has %d negative eigenvalue(s); not removed. Inspect first." % r["n_negative"], EXIT_FAIL)
    if r.get("projected"): die("reactant Hessian must not be path-projected", EXIT_FAIL)
    if a.path:
        if not t.get("projected"): die("--path needs a path-projected Hessian of the path maximum (hess --project-path)", EXIT_FAIL)
        if t["n_negative"]: die("path-projected Hessian still has %d negative eigenvalue(s); not removed" % t["n_negative"], EXIT_FAIL)
        nr = r["nu"]; nt = list(t["nu"]); dE = (t["E0"] - r["E0"]) * HA2KCAL
        print("*** ESTIMATE: path maximum (restrained-path point, gradient direction projected), TS NOT VERIFIED ***")
        print("T = %.2f K  atoms %d  modes R %d / path-max %d (1 removed by projection)" % (a.T, len(r["atoms"]), len(nr), len(nt)))
        print("dE_path_max = %.3f kcal/mol" % dE)
        print("%-12s %8s %8s %8s %8s %12s" % ("low-mode", "dZPE", "dHvib", "-TdS", "dGvib", "dG_path_est"))
        for cut, lab in ((0, "harmonic"), (50, "raise<50"), (100, "raise<100")):
            Gr, Zr, Hr, Sr = gvib(nr, a.T, cut); Gt, Zt, Ht, St = gvib(nt, a.T, cut)
            print("%-12s %8.3f %8.3f %8.3f %8.3f %12.3f" % (lab, Zt - Zr, Ht - Hr, -a.T * (St - Sr), Gt - Gr, dE + Gt - Gr))
        return 0
    if t.get("projected"): die("path-projected Hessian given without --path", EXIT_FAIL)
    if t["n_negative"] != 1: die("TS has %d negative eigenvalues (need exactly 1)" % t["n_negative"], EXIT_FAIL)
    if abs(t["lowest_modes"][0]["cos_RC"]) <= 0.5: print("WARNING: |cos(imag mode, grad RC)| = %.2f <= 0.5" % abs(t["lowest_modes"][0]["cos_RC"]))
    nr = r["nu"]; nt = [v for v in t["nu"] if v > 0]; dE = (t["E0"] - r["E0"]) * HA2KCAL
    print("T = %.2f K  atoms %d  modes R %d / TS %d (+1 imaginary %.1f i cm-1, cos_RC %+.2f)" % (a.T, len(r["atoms"]), len(nr), len(nt), -t["nu"][0], t["lowest_modes"][0]["cos_RC"]))
    print("base residual forces (max comp): R QM %.2e mob %.2e | TS QM %.2e mob %.2e" % (r["base_fmax_qm"], r["base_fmax_mob"], t["base_fmax_qm"], t["base_fmax_mob"]))
    print("dE_double-dagger = %.3f kcal/mol (EPS_SCF 1E-8 base points)" % dE)
    print("%-12s %8s %8s %8s %8s %10s" % ("low-mode", "dZPE", "dHvib", "-TdS", "dGvib", "dG_PHVA"))
    for cut, lab in ((0, "harmonic"), (50, "raise<50"), (100, "raise<100")):
        Gr, Zr, Hr, Sr = gvib(nr, a.T, cut); Gt, Zt, Ht, St = gvib(nt, a.T, cut)
        print("%-12s %8.3f %8.3f %8.3f %8.3f %10.3f" % (lab, Zt - Zr, Ht - Hr, -a.T * (St - Sr), Gt - Gr, dE + Gt - Gr))
    print("modes < 100 cm-1: R %d, TS %d" % (sum(v < 100 for v in nr), sum(v < 100 for v in nt))); return 0


# ---------------- kick / check / dimer / mf ----------------
def cmd_kick(a):
    sd = setdir(a.set); man = json.load(open(os.path.join(sd, "manifest.json"))); hj = json.load(open(a.hess))
    if hj.get("set") != a.set or hj.get("restart_sha") != man["restart_sha"] or hj.get("cfg_sha") != CFG_SHA:
        die("hessian %s does not belong to set %s (set/restart/config mismatch)" % (a.hess, a.set))
    if sha(os.path.join(sd, "base.restart")) != man["restart_sha"]: die("base.restart changed")
    md = hj["lowest_modes"][0]
    if md["nu"] >= 0: die("lowest mode is not imaginary (%.1f)" % md["nu"])
    A = hj["atoms"]; v = md["vec"]
    big = max(math.sqrt(sum(v[3 * i + k] ** 2 for k in range(3))) for i in range(len(A)))
    L, i0, unit = coord_block(os.path.join(sd, "base.restart"))
    outs = ["%s_%s.restart" % (a.name, t) for t in ("p", "m")]
    for o in outs:
        if os.path.exists(o): die(o + " exists")
    for sign, o in zip((1, -1), outs):
        LL = list(L)
        for i, q in enumerate(A):
            for k in range(3): LL[i0 + q - 1] = displaced_line(LL[i0 + q - 1], k, sign * a.amp * v[3 * i + k] / big / unit)
        open(o, "w").writelines(LL)
        print("written", o, "(max atom displacement %.3f A along %.1f i cm-1 mode)" % (a.amp, -md["nu"]))
    return 0


def cmd_check(a):
    fx = fixed_set()
    text = open(a.out, errors="replace").read(); r = parse_out(text)
    print("out: ended %s  SCF-not-converged %d  nscf(last) %d  energy lines %d  E %s  |spin| %s (range %s)  walltime-hit %s" % (
        r["ended"], r["notconv"], r["nscf"], r["nE"], r["E"], r["spin"], SPIN_RANGE, r["walltime_hit"]))
    src = a.forces or a.out
    blk, nb = last_force_block_file(src)
    print("forces source %s: %d block(s); using the LAST one" % (src, nb))
    if blk is None: die("no force block found", EXIT_FAIL)
    F, why = parse_force_block(blk)
    if F is None: die("force block invalid: " + why, EXIT_FAIL)
    print("force table complete: %d rows, ids 1..%d in order, finite, SUM line present" % (len(F), NATOM))
    print("max component QM %.3e  mobile %.3e  fixed %.3e" % (
        max(abs(c) for q in QM for c in F[q]), max(abs(c) for q in F if q not in fx for c in F[q]), max(abs(c) for q in fx for c in F[q])))
    print("P(%d) %s  Ow(%d) %s" % (P, F[P], OW, F[OW]))
    if nb > 1: print("NOTE: multi-block file; energy/forces correspondence is only guaranteed for one-block per-task files")
    return 0


def remove_block(s, name):
    out, depth, skip = [], 0, False
    for l in s.splitlines(True):
        if re.match(r"^\s*&%s\b" % name, l, re.I) and not skip: skip = True; depth = 0
        if skip:
            if re.match(r"^\s*&END\b", l, re.I): depth -= 1
            elif re.match(r"^\s*&", l): depth += 1
            if depth == 0: skip = False
            continue
        out.append(l)
    return "".join(out)


def cmd_dimer(a):
    for f in (a.name + ".inp", a.name + "_dimer_vector.inc", a.name + "_in.wfn"):
        if os.path.exists(f): die(f + " exists")
    fx = fixed_set(); XA, XB = all_xyz(a.a), all_xyz(a.b)
    d = [[XB[i][k] - XA[i][k] for k in range(3)] for i in range(NATOM)]
    fl = sorted(fx); t = [sum(d[i - 1][k] for i in fl) / len(fl) for k in range(3)]
    dev = max(abs(d[i - 1][k] - t[k]) for i in fl for k in range(3))
    print("A->B translation (fixed-atom mean) %s  max deviation %.2e A" % (["%.6f" % x for x in t], dev))
    if dev > 1e-4: die("fixed atoms are not rigidly translated")
    v = [0.0] * (3 * NATOM)
    for q in QM:
        for k in range(3): v[3 * (q - 1) + k] = d[q - 1][k] - t[k]
    nrm = math.sqrt(sum(x * x for x in v))
    if nrm < 1e-6: die("zero dimer vector")
    with open(a.name + "_dimer_vector.inc", "w") as f:
        for i in range(0, len(v), 6): f.write(" ".join("%.10E" % (x / nrm) for x in v[i:i + 6]) + "\n")
    s = remove_block(open(a.template).read(), "COLLECTIVE")
    if re.search(r"(?mi)^\s*&(COLLECTIVE|RESTRAINT)\b", s): die("restraint still present")
    s, c = re.subn(r"(?m)^(\s*)PROJECT(_NAME)?\s+\S+.*$", r"\g<1>PROJECT " + a.name, s)
    if c != 1: die("PROJECT")
    s, c = re.subn(r"(?mi)^(\s*)RESTART_FILE_NAME\s+\S+.*$", r"\g<1>RESTART_FILE_NAME " + a.start, s)
    if c != 1: die("RESTART_FILE_NAME")
    s, c = re.subn(r"(?mi)^(\s*)WFN_RESTART_FILE_NAME\s+\S+.*$", r"\g<1>WFN_RESTART_FILE_NAME " + a.name + "_in.wfn", s)
    if c != 1: die("WFN_RESTART_FILE_NAME")
    if not re.search(r"(?mi)^\s*RESTART_CONSTRAINT\s+(F|\.FALSE\.|FALSE)", s): die("RESTART_CONSTRAINT FALSE not in template")
    s = re.sub(r"(?mi)^(\s*WALLTIME\s+)\S+", r"\g<1>%d" % a.walltime, s)
    for k_ in ("RESTART_DIMER", "RESTART_COUNTERS"):
        if re.search(r"(?mi)^\s*%s\b" % k_, s): s = re.sub(r"(?mi)^(\s*%s\s+)\S+" % k_, r"\g<1>F", s)
        else:
            s, c = re.subn(r"(?mi)^([ \t]*&EXT_RESTART[ \t]*)$", r"\1\n    %s F" % k_, s, count=1)
            if c != 1: die("&EXT_RESTART not found")
    s = re.sub(r"(?mi)^(\s*RUN_TYPE\s+)\S+", r"\g<1>GEO_OPT", s)
    blk = ("  &GEO_OPT\n    TYPE TRANSITION_STATE\n    OPTIMIZER CG\n    MAX_ITER %d\n    MAX_FORCE 4.5E-4\n    RMS_FORCE 3.0E-4\n"
           "    MAX_DR 3.0E-3\n    RMS_DR 1.5E-3\n    &CG\n      MAX_STEEP_STEPS 0\n      &LINE_SEARCH\n        TYPE 2PNT\n      &END LINE_SEARCH\n    &END CG\n"
           "    &TRANSITION_STATE\n      METHOD DIMER\n      &DIMER\n        DR [angstrom] 0.01\n        ANGLE_TOLERANCE [deg] 5.0\n        INTERPOLATE_GRADIENT T\n"
           "        &ROT_OPT\n          OPTIMIZER CG\n          MAX_ITER 10\n          &CG\n            MAX_STEEP_STEPS 0\n            &LINE_SEARCH\n              TYPE 2PNT\n"
           "            &END LINE_SEARCH\n          &END CG\n        &END ROT_OPT\n        &DIMER_VECTOR\n@INCLUDE '%s'\n        &END DIMER_VECTOR\n"
           "      &END DIMER\n    &END TRANSITION_STATE\n  &END GEO_OPT\n") % (a.maxiter, a.name + "_dimer_vector.inc")
    m = re.search(r"(?mis)^[ \t]*&GEO_OPT\b.*?^[ \t]*&END\s+GEO_OPT\b[^\n]*\n", s)
    if not m: die("&GEO_OPT ... &END GEO_OPT not found")
    s = s[:m.start()] + blk + s[m.end():]
    shutil.copy2(a.wfn, a.name + "_in.wfn"); open(a.name + ".inp", "w").write(s)
    print("written", a.name + ".inp", a.name + "_dimer_vector.inc"); return 0


def cmd_mf(a):
    pts = [tuple(float(x) for x in p.split(",")) for p in a.pts]; f = HA2KCAL / BOHR
    s = [2 * a.K * (T - R) / BOHR * f for (T, R, E) in pts]
    print("%7s %9s %9s %12s" % ("T", "RC", "dE", "dE/dRC"))
    for (T, R, E), sl in zip(pts, s): print("%7.3f %9.5f %9.3f %12.2f" % (T, R, E, sl))
    for j in range(1, len(pts)):
        (T0, R0, E0), (T1, R1, E1) = pts[j - 1], pts[j]
        integ = 0.5 * (s[j - 1] + s[j]) * (R1 - R0); curv = (s[j] - s[j - 1]) / (R1 - R0) * BOHR ** 2 / HA2KCAL
        print("  %.3f->%.3f  int %.2f  dE %.2f  diff %+.2f   d(slope)/dRC %+.4f Ha/bohr^2 (stability needs > -%.2f)" % (R0, R1, integ, E1 - E0, integ - (E1 - E0), curv, 2 * a.K))
    for j in range(1, len(pts)):
        if s[j - 1] < 0 <= s[j]:
            (T0, R0, _), (T1, R1, _) = pts[j - 1], pts[j]; g0, g1 = T0 - R0, T1 - R1
            print("BRACKETED between RC %.3f and %.3f -> spring-zero target T* = %.4f A" % (R0, R1, T0 - g0 * (T1 - T0) / (g1 - g0))); return 0
    if len(pts) >= 2:
        (T0, R0, _), (T1, R1, _) = pts[-2], pts[-1]; d = (s[-1] - s[-2]) / (R1 - R0)
        print(("not bracketed; |slope| decreasing -> slope zero predicted at RC ~ %.3f A" % (R1 - s[-1] / d)) if d < 0 else
              "not bracketed; |slope| still INCREASING toward smaller RC")
    return 0


def main():
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="c")
    p = sp.add_parser("validate"); p.add_argument("template"); p.add_argument("--center-grid", action="store_true")
    p = sp.add_parser("check"); p.add_argument("out"); p.add_argument("--forces")
    p = sp.add_parser("prep"); p.add_argument("set"); p.add_argument("--template", required=True); p.add_argument("--restart", required=True)
    p.add_argument("--wfn", required=True); p.add_argument("--mode", choices=("test", "full"), required=True)
    p.add_argument("--h", type=float, default=0.01); p.add_argument("--natoms", type=int); p.add_argument("--must")
    p.add_argument("--order-from"); p.add_argument("--center-grid", action="store_true"); p.add_argument("--dry-run", action="store_true")
    p.add_argument("--path-point", action="store_true", help="base is a non-stationary path point: force gate relaxed to path_gate of the config")
    p = sp.add_parser("work"); p.add_argument("set"); p.add_argument("--hrt", type=float, default=float(os.environ.get("FD_HRT", "14400")))
    p = sp.add_parser("status"); p.add_argument("set"); p.add_argument("--release-stale", action="store_true")
    p = sp.add_parser("test"); p.add_argument("set")
    p = sp.add_parser("hess"); p.add_argument("set"); p.add_argument("--project-path", action="store_true",
        help="path-maximum estimate: project out the energy-gradient direction at the base point")
    p = sp.add_parser("thermo"); p.add_argument("r"); p.add_argument("ts"); p.add_argument("--T", type=float, default=310.15)
    p.add_argument("--path", action="store_true", help="second file is a path-projected path maximum (estimate, not a TS)")
    p = sp.add_parser("kick"); p.add_argument("set"); p.add_argument("--hess", required=True); p.add_argument("--name", required=True); p.add_argument("--amp", type=float, default=0.10)
    p = sp.add_parser("dimer"); p.add_argument("name"); p.add_argument("--template", required=True); p.add_argument("--start", required=True)
    p.add_argument("--wfn", required=True); p.add_argument("--a", required=True); p.add_argument("--b", required=True)
    p.add_argument("--maxiter", type=int, default=60); p.add_argument("--walltime", type=int, default=19800)
    p = sp.add_parser("mf"); p.add_argument("pts", nargs="+"); p.add_argument("--K", type=float, default=0.10)
    a = ap.parse_args()
    if not a.c: ap.print_help(); return 1
    return {"validate": cmd_validate, "check": cmd_check, "prep": cmd_prep, "work": cmd_work, "status": cmd_status,
            "test": cmd_test, "hess": cmd_hess, "thermo": cmd_thermo, "kick": cmd_kick, "dimer": cmd_dimer, "mf": cmd_mf}[a.c](a)


if __name__ == "__main__":
    sys.exit(main() or 0)
