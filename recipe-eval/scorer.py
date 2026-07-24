#!/usr/bin/env python3
"""Score each generation against the box's objective oracle criteria.

An independent judge (Claude CLI) reads the generated code, its captured stdout,
and its figure (Read allowed on the PNG only), and returns a strict JSON verdict
using the numeric criteria supplied per box. The judge does extraction + threshold
checks against oracle values I provide; it does not invent standards. Verdicts are
calibrated against manual inspection on the pilot before trusting on the full set.
RFS(box) = fraction of the box's samples with pass==true.
"""
import glob, json, os, re, subprocess, sys, concurrent.futures as cf

BASE = os.environ.get("RECIPE_EVAL_BASE", "/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval")
GEN = os.path.join(BASE, "gen")
JUDGE_MODEL = "claude-opus-4-8"

# Objective, oracle-backed pass criteria per box.
RUBRIC = {
"2C.1.1": "PASS iff: (a) code implements Euler AND RK4 explicitly (hand-written steppers, not only solve_ivp); (b) the reported N(100) from an explicit method is within 1% of 22026.47 (=exp(10)); (c) RK4 is more accurate than Euler at t=100; (d) a log-scaled N(t) plot is produced.",
"2E.3.1": "PASS iff: (a) bisection is hand-implemented (not brentq/fsolve); (b) a windowed/all-roots scan is applied; (c) at k=0.15 it finds 3 roots approx {71.5, 170.8, 331.6} (+-5%), and at k=0.12 and k=0.20 it finds exactly 1 root (approx 443.4 and 50.9 resp.).",
"5A.4.1": "PASS iff: (a) velocity-Verlet is hand-implemented with the two half-step velocity updates; (b) energy e=0.5*k*x^2+0.5*v^2 is tracked and stays bounded with NO secular growth at BOTH dt=0.01 and dt=0.1 (relative drift under ~10% over t=100; contrast: Euler would blow up); (c) x(t) and e(t) plots produced.",
"6A.4.1": "PASS iff: (a) the polar/Marsaglia Box-Muller is hand-implemented (reject points with R2 not in (0,1); return x*sqrt(-2*ln(R2)/R2) and y*sqrt(...)), NOT np.random.normal; (b) ~10000 samples with sample mean approx 0 (|mean|<0.05) and sd approx 1 (|sd-1|<0.05); (c) a density histogram vs the N(0,1) curve is produced.",
"3A.1.1": "PASS iff: (a) a vector RK4 is hand-implemented for the 2D toggle system; (b) trajectories from ~10 initial conditions converge to exactly TWO distinct stable steady states approx (53,362) and (542,35) (+-15%); (c) a phase-plane plot is produced.",
"8D.2.1": "PASS iff: (a) the Gillespie SSA is hand-implemented (exponential waiting time from total propensity, reaction chosen by propensity); (b) the standard deviation GROWS with the mean while relative noise std/mean DECREASES as x_bar rises (Poisson-like), consistent with the book's own finite-run (tmax=100) result where std sits BELOW sqrt(x_bar) at high copy number -- the book itself gets std/sqrt(x_bar) roughly 0.55 at x_bar=1000, 0.73 at 100, 0.66 at 10; so accept std/sqrt(x_bar) anywhere in [0.4, 1.2]. Reject ONLY if std does not increase with mean, or relative noise does not shrink with x_bar, or std exceeds sqrt(x_bar) badly; (c) trajectories plotted.",
"9B.3.1": "PASS iff: (a) the Held-Karp bitmask DYNAMIC-PROGRAMMING TSP is hand-implemented (a brute-force check alongside it is fine, but a greedy/nearest-neighbor heuristic alone is not); (b) it uses the given 10 city coordinates (x starting 0,-28.87,... y starting 0,0,43.39,...) and returns a closed tour of all 10 cities whose total length is approx 193.73 (within 1%); (c) a plot of the tour is produced.",
"10C.5.1": "PASS iff: (a) it enumerates attractors by iterating the update X'=NOT Y, Y'=NOT X from all 4 states to a repeat; (b) the reported attractors are exactly: fixed point 10, fixed point 01, and the 2-cycle 00<->11 (00->11->00). No other attractors.",
"7E.3.1": "PASS iff (read the PNG figure): (a) a 2D finite-difference reaction-diffusion integrator is hand-implemented with the Gierer-Meinhardt activator-inhibitor kinetics f=u^2/v-u, g=mu*(u^2-v); (b) the final u field is a regular array of SPOTS (isolated round peaks), NOT stripes, NOT uniform, NOT NaN/blown-up.",
"10B.1.1": "PASS iff: (a) k-means (Lloyd) is hand-implemented with multiple restarts kept by lowest WSS (NOT sklearn.cluster.KMeans); (b) it recovers 3 clusters matching the three blobs (three centroids near (0,0),(1.5,1.5),(3,3) within ~0.6, or equivalently a clean 3-way split); (c) a scatter colored by cluster with centroids is produced.",
}

JUDGE = """You are a strict, objective grader. Decide only from the evidence and the exact criteria; do not invent extra requirements. Do NOT run code.

CRITERIA for this task:
{crit}

The candidate's PYTHON CODE:
```
{code}
```

Its STDOUT when executed (empty if it failed):
{stdout}

Ran without error: {ran}. Figure file produced: {figexists} at {figpath} (you may Read that PNG if the criteria need the plot).

Reply with ONLY a JSON object, no prose:
{{"ran": bool, "checks": {{"a": bool, "b": bool, "c": bool, "d": bool}}, "pass": bool, "failure_class": "none|no-run|wrong-method|wrong-params|wrong-output|numeric-off|underspecified-test", "note": "<=20 words"}}
Include only the check keys the criteria define. "pass" is true only if every defined check passes."""

def judge(name, run):
    box = name.rsplit("_s",1)[0]
    code = open(os.path.join(GEN, name + ".py")).read()[:16000]
    fig = os.path.join(GEN, name + ".png")
    prompt = JUDGE.format(crit=RUBRIC[box], code=code, stdout=(run["stdout"] or "(none)")[:3500],
                          ran=run["ran"], figexists=run["fig"], figpath=fig)
    try:
        r = subprocess.run(["claude","-p",prompt,"--model",JUDGE_MODEL,
                            "--allowedTools","Read","--output-format","text"],
                           cwd=GEN, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=240)
        m = re.search(r"\{.*\}", r.stdout, re.S)
        return name, json.loads(m.group(0)) if m else {"pass": False, "note": "judge-parse-fail", "raw": r.stdout[:200]}
    except Exception as e:
        return name, {"pass": False, "note": f"judge-err:{e}"}

def main():
    runs = json.load(open(os.path.join(BASE, "run_results.json")))
    only = sys.argv[1] if len(sys.argv) > 1 else None   # optional: re-score one box
    prev = {}
    if only and os.path.exists(os.path.join(BASE,"verdicts.json")):
        prev = json.load(open(os.path.join(BASE,"verdicts.json")))
    names = [n for n in runs if n.rsplit("_s",1)[0] in RUBRIC and (only is None or n.rsplit("_s",1)[0]==only)]
    verdicts = dict(prev)
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        for name, v in ex.map(lambda n: judge(n, runs[n]), names):
            verdicts[name] = v
            print(f"  {name}: pass={v.get('pass')} class={v.get('failure_class')} {v.get('note','')}")
    json.dump(verdicts, open(os.path.join(BASE, "verdicts.json"), "w"), indent=1)
    # RFS per box
    boxes = sorted(set(n.rsplit("_s",1)[0] for n in verdicts))
    print("\n=== Recipe Fidelity Score (pass rate over K) ===")
    rfs = {}
    for b in boxes:
        samples = [verdicts[n] for n in verdicts if n.rsplit("_s",1)[0]==b]
        p = sum(1 for s in samples if s.get("pass")); rfs[b] = (p, len(samples))
        print(f"  {b}: {p}/{len(samples)}")
    json.dump(rfs, open(os.path.join(BASE, "rfs.json"), "w"), indent=1)
    tot = sum(p for p,_ in rfs.values()); den = sum(n for _,n in rfs.values())
    print(f"\nPilot mean RFS: {tot}/{den} = {100*tot/max(den,1):.0f}%")

if __name__ == "__main__":
    main()
