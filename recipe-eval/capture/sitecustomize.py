"""Auto-imported (via PYTHONPATH) before any script runs. Monkeypatches
matplotlib so every plotting call records the actual data it was handed; at exit
the captured series are dumped to $CAPTURE_OUT as JSON. Used to capture BOTH the
book's reference plots (ground truth) and each generation's plots, so scoring
compares plotted-data to plotted-data instead of parsing screen output.

Never raises into the host script: all hooks are best-effort.
"""
import os, json, atexit

_CAP = []   # list of {"kind":..., "data":...}

def _arr(x):
    try:
        import numpy as np
        a = np.asarray(x, dtype=float)
        if a.size == 0 or a.size > 200000:
            return None
        return a.ravel().tolist() if a.ndim == 1 else a.tolist()
    except Exception:
        return None

def _install():
    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.axes import Axes
    except Exception:
        return

    def wrap(name, recorder):
        orig = getattr(Axes, name, None)
        if orig is None: return
        def patched(self, *a, **k):
            try: recorder(a, k)
            except Exception: pass
            return orig(self, *a, **k)
        setattr(Axes, name, patched)

    def xy(a, k):
        # plot/step/loglog/semilog: (y) or (x,y) [+ fmt], possibly repeated
        nums = [x for x in a if not isinstance(x, str)]
        if len(nums) >= 2 and _arr(nums[0]) and _arr(nums[1]):
            _CAP.append({"kind":"xy", "x":_arr(nums[0]), "y":_arr(nums[1])})
        elif len(nums) == 1 and _arr(nums[0]) is not None:
            _CAP.append({"kind":"y", "y":_arr(nums[0])})
    def scat(a, k):
        if len(a) >= 2:
            _CAP.append({"kind":"scatter", "x":_arr(a[0]), "y":_arr(a[1])})
    def hist(a, k):
        if a: _CAP.append({"kind":"hist", "data":_arr(a[0])})
    def img(a, k):
        if a: _CAP.append({"kind":"image", "z":_arr(a[0])})
    def cont(a, k):
        z = a[-1] if a else None
        _CAP.append({"kind":"contour", "z":_arr(z)})
    def bar(a, k):
        if len(a) >= 2: _CAP.append({"kind":"bar", "x":_arr(a[0]), "y":_arr(a[1])})

    for nm in ("plot","step","loglog","semilogx","semilogy"): wrap(nm, xy)
    wrap("scatter", scat); wrap("hist", hist)
    wrap("imshow", img); wrap("pcolormesh", img); wrap("pcolor", img)
    wrap("contour", cont); wrap("contourf", cont)
    wrap("bar", bar); wrap("barh", bar)

    @atexit.register
    def _dump():
        out = os.environ.get("CAPTURE_OUT")
        if out:
            try: json.dump(_CAP, open(out, "w"))
            except Exception: pass

_install()
