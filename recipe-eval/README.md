# Recipe-box optimization baseline

`original_boxes.{json,md}` snapshot all 116 recipe boxes in their
**pre-optimization** state (the working tree before the cold-generation
evaluation began; 7E.3.1 and 9B.3.1 are their pre-pilot-fix originals).

- `.json` — structured: per box, its file, the six fields, the raw callout, and
  the generated context-free Python prompt.
- `.md` — human-readable: each box's full text plus its generated prompt.

Kept so recipe edits made during optimization can be diffed against the originals.
The live originals also remain in git history.

## Evaluation harness state (paused)

`generate.py`, `run.py`, `scorer.py`, `prompts.py` are the cold-generation
evaluation harness (see the `recipe-eval` memory). `gen/*.py` holds the
context-free generations produced so far. The full generation was paused
partway (see `gen_status.json`); the scripts are cached, so nothing regenerates.

### Resume generation
The harness reads its working dir from `$RECIPE_EVAL_BASE` (defaults to a
session temp dir). To resume against this committed copy:

```sh
export RECIPE_EVAL_BASE="$PWD/recipe-eval"
python3 recipe-eval/generate.py recipe-eval/all_prompts.json   # skips cached, does the rest
python3 recipe-eval/run.py                                     # execute all
python3 recipe-eval/scorer.py                                  # judge -> rfs.json
```

Generation and judging call the `claude` CLI (context-free, tools off), which
consumes usage, so run when you have budget.
