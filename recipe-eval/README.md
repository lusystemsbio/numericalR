# Recipe-box optimization baseline

`original_boxes.{json,md}` snapshot all 116 recipe boxes in their
**pre-optimization** state (the working tree before the cold-generation
evaluation began; 7E.3.1 and 9B.3.1 are their pre-pilot-fix originals).

- `.json` — structured: per box, its file, the six fields, the raw callout, and
  the generated context-free Python prompt.
- `.md` — human-readable: each box's full text plus its generated prompt.

Kept so recipe edits made during optimization can be diffed against the originals.
The live originals also remain in git history.
