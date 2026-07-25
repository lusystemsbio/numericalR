#!/usr/bin/env python3
"""Score each generation against the box's oracle criteria.

Criterion per box = a hand-written OVERRIDE if present (for soft/statistical/
visual boxes where the verification prose diverges from the book's actual
output, or where a tighter numeric check is wanted), otherwise a DEFAULT built
from the box's own Method + Verification + Show fields (the verification is the
acceptance test). An independent judge (Claude CLI, Read allowed only on the
figure) returns a strict JSON verdict. RFS(box) = passes / K.

Judge model defaults to Haiku (cheap; the task is extract-and-threshold);
validate against Opus on the pilot before trusting on the full set:
  JUDGE_MODEL=claude-opus-4-8 python3 scorer.py <box>   # re-score one box
"""
import glob, json, os, re, subprocess, sys, concurrent.futures as cf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import prompts

BASE = os.environ.get("RECIPE_EVAL_BASE", os.path.dirname(os.path.abspath(__file__)))
GEN = os.path.join(BASE, "gen")
JUDGE_MODEL = os.environ.get("JUDGE_MODEL", "claude-haiku-4-5-20251001")
JUDGE_WORKERS = int(os.environ.get("JUDGE_WORKERS", "3"))   # gentler concurrency

def valid_verdict(v):
    """A usable prior verdict (not a usage-limit / call failure)."""
    return (isinstance(v, dict) and "pass" in v
            and "session limit" not in str(v.get("raw", ""))
            and v.get("note") not in ("judge-err", "judge-parse-fail"))
FIELDS = prompts.all_boxes()   # current box fields (match the prompts used)

# Hand-written, book-grounded criteria for boxes whose verification prose is
# looser/different from the book's real output, or that are purely visual.
OVERRIDES = {
"2C.1.1": "PASS iff: (a) Euler AND RK4 are hand-implemented (not only solve_ivp); (b) the reported N(100) from an explicit method is within 1% of 22026.47 (=exp(10)); (c) RK4 beats Euler at t=100; (d) a log-scaled N(t) plot is produced.",
"2E.3.1": "PASS iff: (a) bisection is hand-implemented (not brentq/fsolve); (b) a windowed/all-roots scan is applied; (c) at k=0.15 it finds 3 roots ~{71.5,170.8,331.6} (+-5%), and at k=0.12 and 0.20 exactly 1 root (~443.4 and ~50.9).",
"5A.4.1": "PASS iff: (a) velocity-Verlet hand-implemented with two half-step velocity updates; (b) energy e=0.5*k*x^2+0.5*v^2 stays bounded with NO secular growth at BOTH dt=0.01 and 0.1 (relative drift <~10% over t=100); (c) x(t) and e(t) plots produced.",
"6A.4.1": "PASS iff: (a) polar/Marsaglia Box-Muller hand-implemented (reject R2 not in (0,1); return x*sqrt(-2*ln(R2)/R2) etc.), NOT np.random.normal; (b) ~10000 samples with mean ~0 (|mean|<0.05) and sd ~1 (|sd-1|<0.05); (c) density histogram vs N(0,1).",
"3A.1.1": "PASS iff: (a) vector RK4 hand-implemented for the 2D toggle; (b) trajectories from ~10 ICs converge to exactly TWO distinct stable states ~(53,362) and (542,35) (+-15%); (c) phase-plane plot.",
"8D.2.1": "PASS iff: (a) Gillespie SSA hand-implemented (exp waiting time from total propensity, reaction by propensity); (b) std GROWS with mean while relative noise std/mean DECREASES as x_bar rises, consistent with the book's own finite-run result where std sits BELOW sqrt(x_bar) at high copy number (book gets std/sqrt(x_bar) ~0.55 at x_bar=1000, ~0.73 at 100, ~0.66 at 10; accept std/sqrt(x_bar) in [0.4,1.2]). Reject only if std does not grow with mean, or relative noise does not shrink, or std far exceeds sqrt(x_bar); (c) trajectories plotted.",
"9B.3.1": "PASS iff: (a) Held-Karp bitmask DP hand-implemented (a brute-force check alongside is fine; a greedy heuristic alone is not); (b) uses the given 10 coordinates (x starting 0,-28.87,...) and returns a closed 10-city tour of length ~193.73 (within 1%); (c) a tour plot.",
"10C.5.1": "PASS iff: (a) attractors enumerated by iterating X'=NOT Y, Y'=NOT X from all 4 states; (b) reported attractors are exactly fixed points 10 and 01 and the 2-cycle 00<->11. No others.",
"7E.3.1": "PASS iff (read the PNG): (a) a 2D FD reaction-diffusion integrator with Gierer-Meinhardt activator-inhibitor kinetics f=u^2/v-u, g=mu*(u^2-v); (b) the final u field is a regular array of isolated round SPOTS, NOT stripes/labyrinth, NOT uniform, NOT NaN.",
"10B.1.1": "PASS iff: (a) k-means (Lloyd) hand-implemented with restarts kept by lowest WSS (NOT sklearn KMeans); (b) recovers 3 clusters matching the blobs (centroids near (0,0),(1.5,1.5),(3,3) within ~0.6, or a clean 3-way split); (c) scatter colored by cluster with centroids.",
# --- additional book-grounded overrides for Part 7 pattern morphology ---
"7E.2.1": "PASS iff (read the PNG): (a) a 2D FD reaction-diffusion integrator with Gierer-Meinhardt substrate-depletion kinetics f=u^2*v-u, g=mu*(1-u^2*v); (b) the final u field is a labyrinth of winding STRIPES (connected ridges), NOT isolated spots, NOT uniform, NOT NaN.",
"7C.2.1": "PASS iff (read the PNG): (a) a multi-component 1D FD reaction-diffusion integrator with the substrate-depletion kinetics; (b) the pattern-forming case (d=0.1) shows a STATIONARY spatially-periodic pattern (regular peaks in u across x); reject if the u profile stays flat/uniform or is NaN. (Reproducing all three regimes is a bonus, not required.)",
"7C.3.1": "PASS iff (read the PNG): (a) a multi-component 1D FD reaction-diffusion integrator with activator-inhibitor kinetics f=u^2/v-u, g=mu*(u^2-v); (b) the final state shows a stationary spatially-periodic pattern (regular peaks) with u and v peaks roughly IN PHASE; reject if flat/uniform or NaN.",
"7D.2.1": "PASS iff (read the PNG): (a) a 2D FD integrator for a cAMP field coupled to discrete excitable cells (inactive->excited->refractory); (b) the field shows an organized traveling/spiral wave structure (curved wavefronts), NOT random noise, NOT uniform, NOT NaN.",
}

JUDGE = """You are a strict, objective grader. Decide only from the evidence and the exact criteria; invent no extra requirements. Do NOT run code.

CRITERIA:
{crit}

CANDIDATE PYTHON CODE:
```
{code}
```

ITS STDOUT (empty if it failed):
{stdout}

Ran without error: {ran}. Figure produced: {figexists} at {figpath} (Read that PNG only if the criteria mention the plot/figure).

Reply with ONLY a JSON object, no prose:
{{"ran": bool, "pass": bool, "failure_class": "none|no-run|wrong-method|wrong-params|wrong-output|numeric-off|underspecified-test", "note": "<=20 words"}}
"pass" is true only if EVERY lettered check in the criteria holds."""

def criterion(box):
    if box in OVERRIDES:
        return OVERRIDES[box]
    f = FIELDS.get(box, {})
    parts = ["PASS iff the candidate: (a) uses the stated method -- " + f.get("Method","(any correct method)"),
             "(b) produces a result matching this expected outcome (numeric claims within ~5%, "
             "counts/sets exactly): " + f.get("Verification","(correct result)")]
    if f.get("Show"):
        parts.append("(c) produces the requested output: " + f.get("Show"))
    return " ; ".join(parts) + ". Treat the expected outcome as the acceptance test; ignore cross-references to other sections."

def judge(name, run):
    box = name.rsplit("_s",1)[0]
    code = open(os.path.join(GEN, name + ".py")).read()[:16000]
    fig = os.path.join(GEN, name + ".png")
    prompt = JUDGE.format(crit=criterion(box), code=code, stdout=(run["stdout"] or "(none)")[:3500],
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
    only = sys.argv[1] if len(sys.argv) > 1 else None
    vpath = os.path.join(BASE, "verdicts.json")
    prev = json.load(open(vpath)) if os.path.exists(vpath) else {}
    verdicts = dict(prev)
    # judge a box's samples if requested; else only samples lacking a valid verdict (resume)
    if only:
        names = [n for n in runs if n.rsplit("_s",1)[0] == only]
    else:
        names = [n for n in runs if n not in prev or not valid_verdict(prev[n])]
    print(f"judging {len(names)} generations with {JUDGE_MODEL} (x{JUDGE_WORKERS}); "
          f"{len(prev)-len([n for n in prev if not valid_verdict(prev[n])])} already valid", flush=True)
    with cf.ThreadPoolExecutor(max_workers=JUDGE_WORKERS) as ex:
        for name, v in ex.map(lambda n: judge(n, runs[n]), names):
            verdicts[name] = v
            print(f"  {name}: pass={v.get('pass')} {v.get('failure_class','')} {v.get('note','')}", flush=True)
    json.dump(verdicts, open(vpath, "w"), indent=1)
    boxes = sorted(set(n.rsplit("_s",1)[0] for n in verdicts),
                   key=lambda b:(int(re.match(r'(\d+)',b)[1]), b))
    rfs = {}
    for b in boxes:
        s = [verdicts[n] for n in verdicts if n.rsplit("_s",1)[0]==b]
        rfs[b] = (sum(1 for x in s if x.get("pass")), len(s))
    json.dump(rfs, open(os.path.join(BASE, "rfs.json"), "w"), indent=1)
    tot = sum(p for p,_ in rfs.values()); den = sum(n for _,n in rfs.values())
    print(f"\nmean RFS: {tot}/{den} = {100*tot/max(den,1):.0f}%")

if __name__ == "__main__":
    main()
