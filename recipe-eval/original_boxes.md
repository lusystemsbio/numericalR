# Original recipe boxes (pre-optimization baseline)

116 boxes. Working-tree state before the recipe-eval optimization; 7E.3.1 and 9B.3.1 are their pre-pilot-fix originals.

## Recipe 1A.3.1  (01-intro/01a-getting-started.qmd)

::: {.callout-note title="Recipe 1A.3.1"}
**Objective.** Label each number in a list as "even" or "odd".

**Method.** Loop over the numbers, test each with the modulo (remainder) operator,
and print its label.

**Test.** The list 4, 7, 10, 13.

**Show.** The printed even/odd label for each number in the list.

**Verification.** 4 and 10 print as even; 7 and 13 print as odd.
:::

**Generated prompt (Python):** In Python, label each number in a list as "even" or "odd". Use loop over the numbers, test each with the modulo (remainder) operator, and print its label. Test it on The list 4, 7, 10, 13. Produce the printed even/odd label for each number in the list. As a separate check, confirm that 4 and 10 print as even; 7 and 13 print as odd. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 1A.5.1  (01-intro/01a-getting-started.qmd)

::: {.callout-note title="Recipe 1A.5.1"}
**Objective.** Simulate Conway's Game of Life and follow a glider across the grid.

**Model.** A grid of cells, each alive (1) or dead (0), updated synchronously by three rules applied to each cell's eight neighbors: a live cell with 2 or 3 live neighbors survives; a dead cell with exactly 3 live neighbors becomes alive; every other cell is dead next generation. The grid edges wrap around (a periodic boundary), so it is a torus.

**Method.** Count each cell's live neighbors with wraparound, then apply the three rules to all cells at once.

**Test.** A 10x10 grid seeded with a glider at cells (1,2), (2,3), (3,1), (3,2), (3,3); run for several generations.

**Show.** The grid at several successive generations, with the glider visible.

**Verification.** The glider keeps its five-cell shape and drifts one step diagonally every four generations; on reaching an edge it reappears on the opposite side.
:::

**Generated prompt (Python):** In Python, simulate Conway's Game of Life and follow a glider across the grid. The model is a grid of cells, each alive (1) or dead (0), updated synchronously by three rules applied to each cell's eight neighbors: a live cell with 2 or 3 live neighbors survives; a dead cell with exactly 3 live neighbors becomes alive; every other cell is dead next generation. The grid edges wrap around (a periodic boundary), so it is a torus. Use count each cell's live neighbors with wraparound, then apply the three rules to all cells at once. Test it on A 10x10 grid seeded with a glider at cells (1,2), (2,3), (3,1), (3,2), (3,3); run for several generations. Produce the grid at several successive generations, with the glider visible. As a separate check, confirm that the glider keeps its five-cell shape and drifts one step diagonally every four generations; on reaching an edge it reappears on the opposite side. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 1C.3.1  (01-intro/01c-numerical-methods.qmd)

::: {.callout-note title="Recipe 1C.3.1"}
**Objective.** Compute the square root of a positive number by iteration, without a built-in square-root function.

**Model.** The square root of `a > 0` is the positive solution of `x^2 = a`.

**Method.** The Babylonian iteration `x_next = (x + a/x) / 2`, from any positive guess, repeated until `|x_next^2 - a| < epsilon`.

**Test.** `a = 2`; starting guesses 1, 100, 200; tolerance `epsilon = 1e-6`.

**Show.** The successive iterates and the final estimate for each starting guess.

**Verification.** Every guess converges to about 1.41421; the result squared is within the tolerance of 2, and a far starting guess just needs a few more iterations.
:::

**Generated prompt (Python):** In Python, compute the square root of a positive number by iteration, without a built-in square-root function. The model is the square root of `a > 0` is the positive solution of `x^2 = a`. Use the Babylonian iteration `x_next = (x + a/x) / 2`, from any positive guess, repeated until `|x_next^2 - a| < epsilon`. Test it on `a = 2`; starting guesses 1, 100, 200; tolerance `epsilon = 1e-6`. Produce the successive iterates and the final estimate for each starting guess. As a separate check, confirm that every guess converges to about 1.41421; the result squared is within the tolerance of 2, and a far starting guess just needs a few more iterations. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2A.1.1  (02-odes/02a-modeling-gene-circuits.qmd)

::: {.callout-note title="Recipe 2A.1.1"}
**Objective.** Plot the time course of a constitutively expressed gene from several initial conditions and watch them settle to a common steady state.

**Model.** A gene transcribed at a constant rate g with linear degradation at rate k: `dX/dt = g - k*X`; its exact solution is `X(t) = g/k + (X0 - g/k)*exp(-k*t)`.

**Method.** Direct evaluation of the exact solution over time (no numerical integration).

**Test.** `g = 50` nM/min, `k = 0.1` per min; initial values `X0 = 300, 400, 500, 600, 700` nM; t = 0 to 80 min.

**Show.** A plot of X(t) versus time, one curve per initial condition, with the steady state `g/k` marked.

**Verification.** Every curve begins at its own X0 and approaches the steady state `g/k = 500` nM; curves above 500 decay to it and those below rise.
:::

**Generated prompt (Python):** In Python, plot the time course of a constitutively expressed gene from several initial conditions and watch them settle to a common steady state. The model is a gene transcribed at a constant rate g with linear degradation at rate k: `dX/dt = g - k*X`; its exact solution is `X(t) = g/k + (X0 - g/k)*exp(-k*t)`. Use direct evaluation of the exact solution over time (no numerical integration). Test it on `g = 50` nM/min, `k = 0.1` per min; initial values `X0 = 300, 400, 500, 600, 700` nM; t = 0 to 80 min. Produce a plot of X(t) versus time, one curve per initial condition, with the steady state `g/k` marked. As a separate check, confirm that every curve begins at its own X0 and approaches the steady state `g/k = 500` nM; curves above 500 decay to it and those below rise. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2B.1.1  (02-odes/02b-numerical-integration.qmd)

::: {.callout-note title="Recipe 2B.1.1"}
**Objective.** Integrate the gene-expression ODE from an initial value and compare the numerical solution to the exact one at two step sizes.

**Model.** A constitutively expressed gene: transcription at a constant rate g with linear degradation at rate k, `dX/dt = g - k*X` (exact solution `X(t) = g/k + (X0 - g/k)*exp(-k*t)`).

**Method.** The Euler method, `X(t+dt) = X(t) + (g - k*X(t))*dt`, stepped from t = 0.

**Test.** `g = 50`, `k = 0.1`, `X0 = 300`; t = 0 to 80; dt = 1 and dt = 0.1.

**Show.** A plot of the Euler solution at dt = 1 and dt = 0.1 overlaid on the exact solution.

**Verification.** At dt = 0.1 the numerical curve tracks the exact solution closely (max error well under 1 nM); at dt = 1 there is a visible gap near t = 10; both reach `g/k = 500`.
:::

**Generated prompt (Python):** In Python, integrate the gene-expression ODE from an initial value and compare the numerical solution to the exact one at two step sizes. The model is a constitutively expressed gene: transcription at a constant rate g with linear degradation at rate k, `dX/dt = g - k*X` (exact solution `X(t) = g/k + (X0 - g/k)*exp(-k*t)`). Use the Euler method, `X(t+dt) = X(t) + (g - k*X(t))*dt`, stepped from t = 0. Test it on `g = 50`, `k = 0.1`, `X0 = 300`; t = 0 to 80; dt = 1 and dt = 0.1. Produce a plot of the Euler solution at dt = 1 and dt = 0.1 overlaid on the exact solution. As a separate check, confirm that at dt = 0.1 the numerical curve tracks the exact solution closely (max error well under 1 nM); at dt = 1 there is a visible gap near t = 10; both reach `g/k = 500`. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2B.4.1  (02-odes/02b-numerical-integration.qmd)

::: {.callout-note title="Recipe 2B.4.1"}
**Objective.** Implement a second-order Heun integrator for a one-variable ODE.

**Model.** A constitutively expressed gene: transcription at a constant rate g with linear degradation at rate k, `dX/dt = g - k*X` (exact solution `X(t) = g/k + (X0 - g/k)*exp(-k*t)`).

**Method.** Heun's method: take an Euler step to a predicted endpoint, evaluate the slope there, and advance by the average of the starting and predicted slopes.

**Test.** `g = 50`, `k = 0.1`, `X0 = 300`.

**Show.** A plot of the Heun solution overlaid on the exact solution.

**Verification.** At the same step size it is closer to the exact solution than Euler; halving the step cuts its error by roughly a factor of 4.
:::

**Generated prompt (Python):** In Python, implement a second-order Heun integrator for a one-variable ODE. The model is a constitutively expressed gene: transcription at a constant rate g with linear degradation at rate k, `dX/dt = g - k*X` (exact solution `X(t) = g/k + (X0 - g/k)*exp(-k*t)`). Use heun's method: take an Euler step to a predicted endpoint, evaluate the slope there, and advance by the average of the starting and predicted slopes. Test it on `g = 50`, `k = 0.1`, `X0 = 300`. Produce a plot of the Heun solution overlaid on the exact solution. As a separate check, confirm that at the same step size it is closer to the exact solution than Euler; halving the step cuts its error by roughly a factor of 4. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2B.5.1  (02-odes/02b-numerical-integration.qmd)

::: {.callout-note title="Recipe 2B.5.1"}
**Objective.** Implement the second-order Runge-Kutta (midpoint) integrator.

**Model.** A constitutively expressed gene: transcription at a constant rate g with linear degradation at rate k, `dX/dt = g - k*X` (exact solution `X(t) = g/k + (X0 - g/k)*exp(-k*t)`).

**Method.** RK2: take a trial half-step to the interval midpoint and use the slope there for the full step.

**Test.** `g = 50`, `k = 0.1`, `X0 = 300`.

**Show.** A plot of the RK2 solution overlaid on the exact solution.

**Verification.** Accuracy comparable to Heun and better than Euler at the same step size, checked against the exact solution.
:::

**Generated prompt (Python):** In Python, implement the second-order Runge-Kutta (midpoint) integrator. The model is a constitutively expressed gene: transcription at a constant rate g with linear degradation at rate k, `dX/dt = g - k*X` (exact solution `X(t) = g/k + (X0 - g/k)*exp(-k*t)`). Use rK2: take a trial half-step to the interval midpoint and use the slope there for the full step. Test it on `g = 50`, `k = 0.1`, `X0 = 300`. Produce a plot of the RK2 solution overlaid on the exact solution. As a separate check, confirm that accuracy comparable to Heun and better than Euler at the same step size, checked against the exact solution. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2B.6.1  (02-odes/02b-numerical-integration.qmd)

::: {.callout-note title="Recipe 2B.6.1"}
**Objective.** Implement the fourth-order Runge-Kutta integrator.

**Model.** A constitutively expressed gene: transcription at a constant rate g with linear degradation at rate k, `dX/dt = g - k*X` (exact solution `X(t) = g/k + (X0 - g/k)*exp(-k*t)`).

**Method.** RK4: combine four slope evaluations across the step (start, twice at the midpoint, end) with weights 1, 2, 2, 1.

**Test.** `g = 50`, `k = 0.1`, `X0 = 300`.

**Show.** A plot of the RK4 solution overlaid on the exact solution.

**Verification.** Far more accurate than RK2 or Euler at the same step size; halving the step cuts its error by roughly a factor of 16.
:::

**Generated prompt (Python):** In Python, implement the fourth-order Runge-Kutta integrator. The model is a constitutively expressed gene: transcription at a constant rate g with linear degradation at rate k, `dX/dt = g - k*X` (exact solution `X(t) = g/k + (X0 - g/k)*exp(-k*t)`). Use rK4: combine four slope evaluations across the step (start, twice at the midpoint, end) with weights 1, 2, 2, 1. Test it on `g = 50`, `k = 0.1`, `X0 = 300`. Produce a plot of the RK4 solution overlaid on the exact solution. As a separate check, confirm that far more accurate than RK2 or Euler at the same step size; halving the step cuts its error by roughly a factor of 16. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2B.8.1  (02-odes/02b-numerical-integration.qmd)

::: {.callout-note title="Recipe 2B.8.1"}
**Objective.** Integrate a stiff version of the model where the forward Euler method blows up, using an implicit step.

**Model.** A constitutively expressed gene with a large degradation rate (stiff), `dX/dt = g - k*X`; the exact solution still decays to `g/k`.

**Method.** Backward (implicit) Euler, `X(t+dt) = X(t) + (g - k*X(t+dt))*dt`; since the right-hand side is linear this rearranges to `X_next = (X + dt*g) / (1 + dt*k)`.

**Test.** Stiff parameters `g = 50`, `k = 10`, `X0 = 3`; t = 0 to 4; dt = 0.2.

**Show.** A plot comparing forward Euler and backward Euler at dt = 0.2 against the exact solution.

**Verification.** Forward Euler at dt = 0.2 oscillates and diverges, while backward Euler stays stable and tracks the exact solution toward `g/k = 5`.
:::

**Generated prompt (Python):** In Python, integrate a stiff version of the model where the forward Euler method blows up, using an implicit step. The model is a constitutively expressed gene with a large degradation rate (stiff), `dX/dt = g - k*X`; the exact solution still decays to `g/k`. Use backward (implicit) Euler, `X(t+dt) = X(t) + (g - k*X(t+dt))*dt`; since the right-hand side is linear this rearranges to `X_next = (X + dt*g) / (1 + dt*k)`. Test it on Stiff parameters `g = 50`, `k = 10`, `X0 = 3`; t = 0 to 4; dt = 0.2. Produce a plot comparing forward Euler and backward Euler at dt = 0.2 against the exact solution. As a separate check, confirm that forward Euler at dt = 0.2 oscillates and diverges, while backward Euler stays stable and tracks the exact solution toward `g/k = 5`. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2C.1.1  (02-odes/02c-bacterial-growth.qmd)

::: {.callout-note title="Recipe 2C.1.1"}
**Objective.** Model exponential bacterial growth and reproduce the exact solution numerically.

**Model.** A population growing at a constant per-capita rate r: `dN/dt = r*N`, with exact solution `N(t) = N0*exp(r*t)`.

**Method.** Explicit integration (Euler and RK4) and a general ODE solver.

**Test.** r = 0.1, N0 = 1; t = 0 to 100; step dt = 0.1.

**Show.** A log-scaled plot of N(t) for Euler, RK4, and the exact solution on one set of axes.

**Verification.** On a log-scaled y axis the exact solution is a straight line; at large t, RK4 tracks it far more closely than Euler; the final value matches `N(100) = exp(10)`, about 22026.
:::

**Generated prompt (Python):** In Python, model exponential bacterial growth and reproduce the exact solution numerically. The model is a population growing at a constant per-capita rate r: `dN/dt = r*N`, with exact solution `N(t) = N0*exp(r*t)`. Use explicit integration (Euler and RK4) and a general ODE solver. Test it on r = 0.1, N0 = 1; t = 0 to 100; step dt = 0.1. Produce a log-scaled plot of N(t) for Euler, RK4, and the exact solution on one set of axes. As a separate check, confirm that on a log-scaled y axis the exact solution is a straight line; at large t, RK4 tracks it far more closely than Euler; the final value matches `N(100) = exp(10)`, about 22026. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2C.2.1  (02-odes/02c-bacterial-growth.qmd)

::: {.callout-note title="Recipe 2C.2.1"}
**Objective.** Plot the exact logistic-growth curve as a target to check simulations against.

**Model.** Logistic growth at rate r with carrying capacity B: `dN/dt = r*N*(1 - N/B)`, exact solution `N(t) = N0*B / (N0 + (B - N0)*exp(-r*t))`.

**Method.** Direct evaluation of the exact solution over time.

**Test.** r = 0.1, B = 100, N0 = 1; t = 0 to 100.

**Show.** A plot of the logistic curve N(t) rising to the carrying capacity B.

**Verification.** The curve rises sigmoidally from N0 and levels off at the carrying capacity `B = 100`.
:::

**Generated prompt (Python):** In Python, plot the exact logistic-growth curve as a target to check simulations against. The model is logistic growth at rate r with carrying capacity B: `dN/dt = r*N*(1 - N/B)`, exact solution `N(t) = N0*B / (N0 + (B - N0)*exp(-r*t))`. Use direct evaluation of the exact solution over time. Test it on r = 0.1, B = 100, N0 = 1; t = 0 to 100. Produce a plot of the logistic curve N(t) rising to the carrying capacity B. As a separate check, confirm that the curve rises sigmoidally from N0 and levels off at the carrying capacity `B = 100`. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2D.1.1  (02-odes/02d-effective-potential.qmd)

::: {.callout-note title="Recipe 2D.1.1"}
**Objective.** Decide whether a gene circuit's steady state is stable or unstable.

**Model.** A self-inhibiting gene: basal transcription plus a repressive Hill function, minus linear degradation, `f(X) = g0 + g1/(1 + (X/Xth)^n) - k*X`.

**Method.** Linear stability: a steady state (`f(X) = 0`) is stable if the slope `df/dX` there is negative, estimated with a central finite difference.

**Test.** `g0 = 10, g1 = 60, Xth = 200, n = 4, k = 0.1`; the single steady state is near 250 nM.

**Show.** The steady-state value and the sign of `df/dX` there, with the stability verdict.

**Verification.** `df/dX` at the steady state (about 250 nM) is negative, so the state is stable.
:::

**Generated prompt (Python):** In Python, decide whether a gene circuit's steady state is stable or unstable. The model is a self-inhibiting gene: basal transcription plus a repressive Hill function, minus linear degradation, `f(X) = g0 + g1/(1 + (X/Xth)^n) - k*X`. Use linear stability: a steady state (`f(X) = 0`) is stable if the slope `df/dX` there is negative, estimated with a central finite difference. Test it on `g0 = 10, g1 = 60, Xth = 200, n = 4, k = 0.1`; the single steady state is near 250 nM. Produce the steady-state value and the sign of `df/dX` there, with the stability verdict. As a separate check, confirm that `df/dX` at the steady state (about 250 nM) is negative, so the state is stable. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2D.3.1  (02-odes/02d-effective-potential.qmd)

::: {.callout-note title="Recipe 2D.3.1"}
**Objective.** Compute and plot the effective potential of a one-variable gene circuit to see its stable and unstable states as valleys and peaks.

**Model.** A self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`.

**Method.** The effective potential `U(X) = -integral of f(x) dx` with `U(0) = 0`, accumulated by the trapezoidal rule `U(X+dx) = U(X) - (f(X) + f(X+dx))/2 * dx`.

**Test.** `g0 = 10, g1 = 45, Xth = 200, n = 4`, at `k = 0.15` (and again `k = 0.2` and `k = 0.1`).

**Show.** A plot of the effective potential U(X) at k = 0.15, 0.2, and 0.1.

**Verification.** At k = 0.15 the potential has two basins (near 100 and 300 nM, the stable states) split by a barrier near 200 nM (the unstable state); at k = 0.2 and k = 0.1 it has a single well.
:::

**Generated prompt (Python):** In Python, compute and plot the effective potential of a one-variable gene circuit to see its stable and unstable states as valleys and peaks. The model is a self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`. Use the effective potential `U(X) = -integral of f(x) dx` with `U(0) = 0`, accumulated by the trapezoidal rule `U(X+dx) = U(X) - (f(X) + f(X+dx))/2 * dx`. Test it on `g0 = 10, g1 = 45, Xth = 200, n = 4`, at `k = 0.15` (and again `k = 0.2` and `k = 0.1`). Produce a plot of the effective potential U(X) at k = 0.15, 0.2, and 0.1. As a separate check, confirm that at k = 0.15 the potential has two basins (near 100 and 300 nM, the stable states) split by a barrier near 200 nM (the unstable state); at k = 0.2 and k = 0.1 it has a single well. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2E.2.1  (02-odes/02e-bifurcation.qmd)

::: {.callout-note title="Recipe 2E.2.1"}
**Objective.** Trace the bifurcation curve: at each value of the parameter k, find every steady state and mark it stable or unstable.

**Model.** A self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter.

**Method.** Sample k on a grid; at each k find all roots of `f(X, k) = 0` with a library root-finder, classify each by the sign of `df/dX`, and order the scattered points with a nearest-neighbor walk.

**Test.** `g0 = 10, g1 = 45, Xth = 200, n = 4`; k swept across a range (bistable near k = 0.15).

**Show.** A plot of steady-state X versus k, points colored by stability (the bifurcation curve).

**Verification.** The result is an S-shaped curve: over a middle k range there are three steady states (two stable, one unstable), collapsing to one outside it.
:::

**Generated prompt (Python):** In Python, trace the bifurcation curve: at each value of the parameter k, find every steady state and mark it stable or unstable. The model is a self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter. Use sample k on a grid; at each k find all roots of `f(X, k) = 0` with a library root-finder, classify each by the sign of `df/dX`, and order the scattered points with a nearest-neighbor walk. Test it on `g0 = 10, g1 = 45, Xth = 200, n = 4`; k swept across a range (bistable near k = 0.15). Produce a plot of steady-state X versus k, points colored by stability (the bifurcation curve). As a separate check, confirm that the result is an S-shaped curve: over a middle k range there are three steady states (two stable, one unstable), collapsing to one outside it. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2E.3.1  (02-odes/02e-bifurcation.qmd)

::: {.callout-note title="Recipe 2E.3.1"}
**Objective.** Find steady states with a homemade bracketing root-finder, and recover all roots by scanning windows.

**Model.** A self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter.

**Method.** Bisection: on an interval where f changes sign, repeatedly halve it, keeping the half that still brackets the root; apply it to each small window with a sign change to find all roots.

**Test.** `g0 = 10, g1 = 45, Xth = 200, n = 4` on [0, 600], at k = 0.12, 0.15, and 0.2.

**Show.** The roots found on the whole interval versus by the windowed scan, at each k.

**Verification.** On the whole interval, bisection returns the single root at monostable k (0.12, 0.2) but only one of the three at bistable k = 0.15; the windowed scan recovers all three.
:::

**Generated prompt (Python):** In Python, find steady states with a homemade bracketing root-finder, and recover all roots by scanning windows. The model is a self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter. Use bisection: on an interval where f changes sign, repeatedly halve it, keeping the half that still brackets the root; apply it to each small window with a sign change to find all roots. Test it on `g0 = 10, g1 = 45, Xth = 200, n = 4` on [0, 600], at k = 0.12, 0.15, and 0.2. Produce the roots found on the whole interval versus by the windowed scan, at each k. As a separate check, confirm that on the whole interval, bisection returns the single root at monostable k (0.12, 0.2) but only one of the three at bistable k = 0.15; the windowed scan recovers all three. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2E.4.1  (02-odes/02e-bifurcation.qmd)

::: {.callout-note title="Recipe 2E.4.1"}
**Objective.** Find roots with the false-position bracketing method.

**Model.** A self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter.

**Method.** False position: like bisection, but replace the midpoint with the point where the line through the interval endpoints crosses zero, `x_new = (xmin*f2 - xmax*f1) / (f2 - f1)`; drive it over windows to find all roots.

**Test.** `g0 = 10, g1 = 45, Xth = 200, n = 4` on [0, 600] at k = 0.15.

**Show.** The roots found on the whole interval versus by the windowed scan at k = 0.15.

**Verification.** On the whole interval at k = 0.15 it returns only one of the three roots (possibly a different one than bisection); the windowed scan finds all three, often in fewer iterations.
:::

**Generated prompt (Python):** In Python, find roots with the false-position bracketing method. The model is a self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter. Use false position: like bisection, but replace the midpoint with the point where the line through the interval endpoints crosses zero, `x_new = (xmin*f2 - xmax*f1) / (f2 - f1)`; drive it over windows to find all roots. Test it on `g0 = 10, g1 = 45, Xth = 200, n = 4` on [0, 600] at k = 0.15. Produce the roots found on the whole interval versus by the windowed scan at k = 0.15. As a separate check, confirm that on the whole interval at k = 0.15 it returns only one of the three roots (possibly a different one than bisection); the windowed scan finds all three, often in fewer iterations. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2E.5.1  (02-odes/02e-bifurcation.qmd)

::: {.callout-note title="Recipe 2E.5.1"}
**Objective.** Trace the full S-shaped bifurcation curve, including the unstable branch, by following steady states with ODE integration rather than root-finding.

**Model.** A self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter.

**Method.** From a starting state, integrate `dX/dt = f(X, k)` to a stable steady state, nudge k, and continue from the previous state; reverse the k-sweep at a large jump, and integrate `dX/dt = -f(X, k)` to capture the unstable branch as an attractor.

**Test.** `g0 = 10, g1 = 45, Xth = 200, n = 4`, sweeping k across the bistable range.

**Show.** A plot of steady-state X versus k, including the unstable branch.

**Verification.** The traced points reproduce the S-curve of 2E.2, now including the unstable middle branch.
:::

**Generated prompt (Python):** In Python, trace the full S-shaped bifurcation curve, including the unstable branch, by following steady states with ODE integration rather than root-finding. The model is a self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter. Use from a starting state, integrate `dX/dt = f(X, k)` to a stable steady state, nudge k, and continue from the previous state; reverse the k-sweep at a large jump, and integrate `dX/dt = -f(X, k)` to capture the unstable branch as an attractor. Test it on `g0 = 10, g1 = 45, Xth = 200, n = 4`, sweeping k across the bistable range. Produce a plot of steady-state X versus k, including the unstable branch. As a separate check, confirm that the traced points reproduce the S-curve of 2E.2, now including the unstable middle branch. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2F.2.1  (02-odes/02f-bifurcation-curves.qmd)

::: {.callout-note title="Recipe 2F.2.1"}
**Objective.** Trace the bifurcation curve as a connected, ordered set of steady states.

**Model.** A self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter.

**Method.** Treat `z = f(k, X)` as a surface over the (k, X) plane and take its zero-level contour on a grid.

**Test.** `g0 = 10, g1 = 45, Xth = 200, n = 4`, evaluated over a grid of (k, X).

**Show.** A plot of the zero contour of f(k, X) in the (k, X) plane (the bifurcation curve).

**Verification.** The zero contour is the same S-shaped curve as in Part 2E, returned as connected line segments rather than a point cloud; one contour captures the whole connected curve.
:::

**Generated prompt (Python):** In Python, trace the bifurcation curve as a connected, ordered set of steady states. The model is a self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter. Use treat `z = f(k, X)` as a surface over the (k, X) plane and take its zero-level contour on a grid. Test it on `g0 = 10, g1 = 45, Xth = 200, n = 4`, evaluated over a grid of (k, X). Produce a plot of the zero contour of f(k, X) in the (k, X) plane (the bifurcation curve). As a separate check, confirm that the zero contour is the same S-shaped curve as in Part 2E, returned as connected line segments rather than a point cloud; one contour captures the whole connected curve. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2F.5.1  (02-odes/02f-bifurcation-curves.qmd)

::: {.callout-note title="Recipe 2F.5.1"}
**Objective.** Follow the bifurcation curve `X(k)` with a predictor-corrector method.

**Model.** A self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter.

**Method.** Numerical continuation: from a known point on `f(X, k) = 0`, predict along the tangent `dX/dk = -(df/dk) / (df/dX)`, then correct back onto the curve with a root solve at the new k.

**Test.** `g0 = 10, g1 = 45, Xth = 200, n = 4`, continued across a range of k.

**Show.** A plot of the continued branch X(k), showing where it stalls at the fold.

**Verification.** It follows a branch accurately away from folds but stalls at a fold, where `dX/dk` diverges because k is multivalued there.
:::

**Generated prompt (Python):** In Python, follow the bifurcation curve `X(k)` with a predictor-corrector method. The model is a self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter. Use numerical continuation: from a known point on `f(X, k) = 0`, predict along the tangent `dX/dk = -(df/dk) / (df/dX)`, then correct back onto the curve with a root solve at the new k. Test it on `g0 = 10, g1 = 45, Xth = 200, n = 4`, continued across a range of k. Produce a plot of the continued branch X(k), showing where it stalls at the fold. As a separate check, confirm that it follows a branch accurately away from folds but stalls at a fold, where `dX/dk` diverges because k is multivalued there. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 2F.6.1  (02-odes/02f-bifurcation-curves.qmd)

::: {.callout-note title="Recipe 2F.6.1"}
**Objective.** Trace the entire S-curve, through the folds, without stalling.

**Model.** A self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter.

**Method.** Arc-length continuation: parameterize by arc length s instead of k, splitting each step into `dk` and `dX` with `ds^2 = dk^2 + dX^2` and `dX = h*dk` (h = dX/dk), so the step stays finite where `dX/dk` diverges.

**Test.** `g0 = 10, g1 = 45, Xth = 200, n = 4`, across the full bistable range.

**Show.** A plot of the full S-shaped curve traced through both folds.

**Verification.** The curve is traced smoothly through both folds, recovering the full S-shape (both stable branches and the unstable middle) where plain k-continuation stalled.
:::

**Generated prompt (Python):** In Python, trace the entire S-curve, through the folds, without stalling. The model is a self-activating gene: basal transcription plus an excitatory Hill function, minus linear degradation, `f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X`, with k the control parameter. Use arc-length continuation: parameterize by arc length s instead of k, splitting each step into `dk` and `dX` with `ds^2 = dk^2 + dX^2` and `dX = h*dk` (h = dX/dk), so the step stays finite where `dX/dk` diverges. Test it on `g0 = 10, g1 = 45, Xth = 200, n = 4`, across the full bistable range. Produce a plot of the full S-shaped curve traced through both folds. As a separate check, confirm that the curve is traced smoothly through both folds, recovering the full S-shape (both stable branches and the unstable middle) where plain k-continuation stalled. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3A.1.1  (03-phase-plane/03a-nullclines.qmd)

::: {.callout-note title="Recipe 3A.1.1"}
**Objective.** Simulate a genetic toggle switch from many initial conditions and see which steady states the trajectories reach.

**Model.** Genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`.

**Method.** A vector form of RK4, with the state and each stage carried as 2-vectors.

**Test.** `gX0=5, gX1=50, Yth=100, nY=4, kX=0.1`; `gY0=4, gY1=40, Xth=150, nX=4, kY=0.12`; ten random initial (X, Y) in [0, 600].

**Show.** A phase-plane plot of the ten trajectories converging to the two stable steady states.

**Verification.** Every trajectory settles on one of two distinct stable steady states (different X and Y levels), showing the switch is bistable.
:::

**Generated prompt (Python):** In Python, simulate a genetic toggle switch from many initial conditions and see which steady states the trajectories reach. The model is genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`. Use a vector form of RK4, with the state and each stage carried as 2-vectors. Test it on `gX0=5, gX1=50, Yth=100, nY=4, kX=0.1`; `gY0=4, gY1=40, Xth=150, nX=4, kY=0.12`; ten random initial (X, Y) in [0, 600]. Produce a phase-plane plot of the ten trajectories converging to the two stable steady states. As a separate check, confirm that every trajectory settles on one of two distinct stable steady states (different X and Y levels), showing the switch is bistable. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3A.4.1  (03-phase-plane/03a-nullclines.qmd)

::: {.callout-note title="Recipe 3A.4.1"}
**Objective.** Draw the X- and Y-nullclines of the toggle switch, whose intersections are the steady states.

**Model.** Genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`.

**Method.** Solve each nullcline in closed form by separation of variables: rearrange `fX(X,Y) = 0` for X and `fY(X,Y) = 0` for Y, sweeping the other variable.

**Test.** The toggle switch with `gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`.

**Show.** A phase-plane plot of the X- and Y-nullclines with their three crossings.

**Verification.** The two nullclines cross three times, matching the two stable steady states and one unstable state seen by simulation.
:::

**Generated prompt (Python):** In Python, draw the X- and Y-nullclines of the toggle switch, whose intersections are the steady states. The model is genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`. Use solve each nullcline in closed form by separation of variables: rearrange `fX(X,Y) = 0` for X and `fY(X,Y) = 0` for Y, sweeping the other variable. Test it on The toggle switch with `gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`. Produce a phase-plane plot of the X- and Y-nullclines with their three crossings. As a separate check, confirm that the two nullclines cross three times, matching the two stable steady states and one unstable state seen by simulation. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3A.5.1  (03-phase-plane/03a-nullclines.qmd)

::: {.callout-note title="Recipe 3A.5.1"}
**Objective.** Draw the nullclines with a general method that does not need separation of variables.

**Model.** Genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`.

**Method.** Treat each nullcline as the zero-level contour of the surface `Z = fX(X,Y)` (and `Z = fY(X,Y)`), evaluated on a grid and traced by a contour routine.

**Test.** The toggle switch (`gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`) on a grid over the (X, Y) plane.

**Show.** A phase-plane plot of the nullclines drawn as zero contours, with their crossings.

**Verification.** The zero contours reproduce the same nullclines as separation of variables, and their crossings match the steady states.
:::

**Generated prompt (Python):** In Python, draw the nullclines with a general method that does not need separation of variables. The model is genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`. Use treat each nullcline as the zero-level contour of the surface `Z = fX(X,Y)` (and `Z = fY(X,Y)`), evaluated on a grid and traced by a contour routine. Test it on The toggle switch (`gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`) on a grid over the (X, Y) plane. Produce a phase-plane plot of the nullclines drawn as zero contours, with their crossings. As a separate check, confirm that the zero contours reproduce the same nullclines as separation of variables, and their crossings match the steady states. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3A.6.1  (03-phase-plane/03a-nullclines.qmd)

::: {.callout-note title="Recipe 3A.6.1"}
**Objective.** Trace a nullcline by following it as a curve, when neither separation nor a grid contour is convenient.

**Model.** Genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`.

**Method.** Arc-length numerical continuation (predictor-corrector): treat Y as the control parameter, follow `X(Y)` along `fX(X,Y) = 0`, stepping by arc length and correcting back with Newton's method using the Jacobian.

**Test.** The toggle switch's X-nullcline `fX(X,Y) = 0`, with `gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`.

**Show.** A phase-plane plot of the X-nullcline traced by continuation.

**Verification.** The traced curve matches the nullcline found by the other two methods.
:::

**Generated prompt (Python):** In Python, trace a nullcline by following it as a curve, when neither separation nor a grid contour is convenient. The model is genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`. Use arc-length numerical continuation (predictor-corrector): treat Y as the control parameter, follow `X(Y)` along `fX(X,Y) = 0`, stepping by arc length and correcting back with Newton's method using the Jacobian. Test it on The toggle switch's X-nullcline `fX(X,Y) = 0`, with `gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`. Produce a phase-plane plot of the X-nullcline traced by continuation. As a separate check, confirm that the traced curve matches the nullcline found by the other two methods. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3B.2.1  (03-phase-plane/03b-steady-states-2d.qmd)

::: {.callout-note title="Recipe 3B.2.1"}
**Objective.** Find every steady state of the toggle switch as the intersections of its nullclines.

**Model.** Genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`.

**Method.** Build the two nullclines, then find where they cross with a segment-intersection test; a faster version first narrows to segments that change sign before testing pairs.

**Test.** The toggle-switch nullclines with `gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`.

**Show.** The three steady-state coordinates, marked where the nullclines cross.

**Verification.** It returns three steady states (two stable, one unstable); the fast version agrees with the exhaustive all-pairs search while testing far fewer pairs.
:::

**Generated prompt (Python):** In Python, find every steady state of the toggle switch as the intersections of its nullclines. The model is genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`. Use build the two nullclines, then find where they cross with a segment-intersection test; a faster version first narrows to segments that change sign before testing pairs. Test it on The toggle-switch nullclines with `gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`. Produce the three steady-state coordinates, marked where the nullclines cross. As a separate check, confirm that it returns three steady states (two stable, one unstable); the fast version agrees with the exhaustive all-pairs search while testing far fewer pairs. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3B.3.1  (03-phase-plane/03b-steady-states-2d.qmd)

::: {.callout-note title="Recipe 3B.3.1"}
**Objective.** Classify each steady state of a two-variable system as stable or unstable.

**Model.** Genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`.

**Method.** Linearize about the steady state and take the eigenvalues of the 2x2 Jacobian (computed numerically); the state is stable when both eigenvalues have negative real part.

**Test.** The three steady states of the toggle switch (`gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`).

**Show.** The Jacobian eigenvalues and the stability label for each steady state.

**Verification.** The two outer states are stable and the middle one is an unstable saddle, matching the bistable simulations.
:::

**Generated prompt (Python):** In Python, classify each steady state of a two-variable system as stable or unstable. The model is genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`. Use linearize about the steady state and take the eigenvalues of the 2x2 Jacobian (computed numerically); the state is stable when both eigenvalues have negative real part. Test it on The three steady states of the toggle switch (`gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`). Produce the Jacobian eigenvalues and the stability label for each steady state. As a separate check, confirm that the two outer states are stable and the middle one is an unstable saddle, matching the bistable simulations. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3C.2.1  (03-phase-plane/03c-chemostat.qmd)

::: {.callout-note title="Recipe 3C.2.1"}
**Objective.** Simulate a chemostat (population N and substrate C) and see it reach either washout or coexistence.

**Model.** Population N grows on nutrient by Michaelis-Menten uptake and washes out at the dilution rate, while substrate C is consumed and fed in at a scaled rate: `dN/dt = a1*(C/(C+1))*N - N`, `dC/dt = -(C/(C+1))*N - C + a2`.

**Method.** The generic multi-variable RK4 from Part 3A, integrating the chemostat ODEs with the feed and dilution parameters.

**Test.** `a1 = 2`, `a2 = 5` (giving a washout state (0, 5) and a coexistence state (8, 1)); from N(0) = 0 and from a tiny seed N(0) = 0.01.

**Show.** A phase-plane plot of the two trajectories, one reaching washout (0, 5) and one coexistence (8, 1).

**Verification.** Starting at N = 0 relaxes to the washout state (0, 5); a tiny seed instead ends at the coexistence state (8, 1), first drifting toward washout then peeling away.
:::

**Generated prompt (Python):** In Python, simulate a chemostat (population N and substrate C) and see it reach either washout or coexistence. The model is population N grows on nutrient by Michaelis-Menten uptake and washes out at the dilution rate, while substrate C is consumed and fed in at a scaled rate: `dN/dt = a1*(C/(C+1))*N - N`, `dC/dt = -(C/(C+1))*N - C + a2`. Use the generic multi-variable RK4 from Part 3A, integrating the chemostat ODEs with the feed and dilution parameters. Test it on `a1 = 2`, `a2 = 5` (giving a washout state (0, 5) and a coexistence state (8, 1)); from N(0) = 0 and from a tiny seed N(0) = 0.01. Produce a phase-plane plot of the two trajectories, one reaching washout (0, 5) and one coexistence (8, 1). As a separate check, confirm that starting at N = 0 relaxes to the washout state (0, 5); a tiny seed instead ends at the coexistence state (8, 1), first drifting toward washout then peeling away. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3D.2.1  (03-phase-plane/03d-predator-prey.qmd)

::: {.callout-note title="Recipe 3D.2.1"}
**Objective.** Simulate the Lotka-Volterra predator-prey model and see its oscillations and closed orbits.

**Model.** The Lotka-Volterra predator-prey model: prey N grow and are eaten by predators P, which grow on the prey and die off, `dN/dt = N*(a - b*P)`, `dP/dt = P*(c*N - d)`.

**Method.** Integrate the two-variable system with the generic RK4, from several initial conditions.

**Test.** `a = 1, b = 0.03, c = 0.02, d = 1`; starting counts including (N, P) = (30, 10), (40, 20), (30, 25), (20, 40); integrated 50 time units at dt = 0.01.

**Show.** A time-series plot of N and P, and a phase-plane plot of the nested closed orbits.

**Verification.** The two populations oscillate out of phase in time, and in the phase plane the trajectories form a family of nested closed loops.
:::

**Generated prompt (Python):** In Python, simulate the Lotka-Volterra predator-prey model and see its oscillations and closed orbits. The model is the Lotka-Volterra predator-prey model: prey N grow and are eaten by predators P, which grow on the prey and die off, `dN/dt = N*(a - b*P)`, `dP/dt = P*(c*N - d)`. Use integrate the two-variable system with the generic RK4, from several initial conditions. Test it on `a = 1, b = 0.03, c = 0.02, d = 1`; starting counts including (N, P) = (30, 10), (40, 20), (30, 25), (20, 40); integrated 50 time units at dt = 0.01. Produce a time-series plot of N and P, and a phase-plane plot of the nested closed orbits. As a separate check, confirm that the two populations oscillate out of phase in time, and in the phase plane the trajectories form a family of nested closed loops. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3D.6.1  (03-phase-plane/03d-predator-prey.qmd)

::: {.callout-note title="Recipe 3D.6.1"}
**Objective.** Fit the Lotka-Volterra parameters to observed data by minimizing the error.

**Model.** The Lotka-Volterra predator-prey model: prey N grow and are eaten by predators P, which grow on the prey and die off, `dN/dt = N*(a - b*P)`, `dP/dt = P*(c*N - d)`.

**Method.** Simulate the model for trial parameters, measure the sum of squared differences from the data, and adjust the parameters with a local optimizer (Newton or BFGS) from a reasonable starting guess.

**Test.** The lynx-hare yearly data, from a hand-picked starting guess for a, b, c, d.

**Show.** The fitted parameters and the fitted trajectory over the lynx-hare data, in time and in the phase plane.

**Verification.** The fitted trajectory follows both the time courses and the phase-plane loop of the data.
:::

**Generated prompt (Python):** In Python, fit the Lotka-Volterra parameters to observed data by minimizing the error. The model is the Lotka-Volterra predator-prey model: prey N grow and are eaten by predators P, which grow on the prey and die off, `dN/dt = N*(a - b*P)`, `dP/dt = P*(c*N - d)`. Use simulate the model for trial parameters, measure the sum of squared differences from the data, and adjust the parameters with a local optimizer (Newton or BFGS) from a reasonable starting guess. Test it on The lynx-hare yearly data, from a hand-picked starting guess for a, b, c, d. Produce the fitted parameters and the fitted trajectory over the lynx-hare data, in time and in the phase plane. As a separate check, confirm that the fitted trajectory follows both the time courses and the phase-plane loop of the data. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3D.7.1  (03-phase-plane/03d-predator-prey.qmd)

::: {.callout-note title="Recipe 3D.7.1"}
**Objective.** Fit the Lotka-Volterra parameters by linear regression, using that the model is linear in its parameters.

**Model.** The Lotka-Volterra predator-prey model: prey N grow and are eaten by predators P, which grow on the prey and die off, `dN/dt = N*(a - b*P)`, `dP/dt = P*(c*N - d)`.

**Method.** Replace the time derivatives with centered finite differences at each data point, turning the rate equations into a linear system, and solve it by ordinary least squares (no intercept).

**Test.** The lynx-hare yearly data.

**Show.** The regression-fitted parameters and the fitted trajectory over the data.

**Verification.** It returns parameters directly and deterministically; on the sparse yearly data the fit is rougher than error minimization but far cheaper.
:::

**Generated prompt (Python):** In Python, fit the Lotka-Volterra parameters by linear regression, using that the model is linear in its parameters. The model is the Lotka-Volterra predator-prey model: prey N grow and are eaten by predators P, which grow on the prey and die off, `dN/dt = N*(a - b*P)`, `dP/dt = P*(c*N - d)`. Use replace the time derivatives with centered finite differences at each data point, turning the rate equations into a linear system, and solve it by ordinary least squares (no intercept). Test it on The lynx-hare yearly data. Produce the regression-fitted parameters and the fitted trajectory over the data. As a separate check, confirm that it returns parameters directly and deterministically; on the sparse yearly data the fit is rougher than error minimization but far cheaper. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3E.1.1  (03-phase-plane/03e-bifurcation-2d.qmd)

::: {.callout-note title="Recipe 3E.1.1"}
**Objective.** Trace a bifurcation diagram: the steady states of the toggle switch versus a control parameter, colored by stability.

**Model.** Genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`.

**Method.** At each parameter value, build the nullclines (separation of variables), find their intersections (steady states), and classify each by its Jacobian eigenvalues; sweep the parameter and plot steady-state X against it.

**Test.** The toggle switch (`gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`), sweeping the X production rate gX1 from 0 to 100.

**Show.** A plot of steady-state X versus the control parameter, points colored by stability.

**Verification.** Stable and unstable branches appear and merge at bifurcation points, so the number and stability of steady states change with the parameter.
:::

**Generated prompt (Python):** In Python, trace a bifurcation diagram: the steady states of the toggle switch versus a control parameter, colored by stability. The model is genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`. Use at each parameter value, build the nullclines (separation of variables), find their intersections (steady states), and classify each by its Jacobian eigenvalues; sweep the parameter and plot steady-state X against it. Test it on The toggle switch (`gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`), sweeping the X production rate gX1 from 0 to 100. Produce a plot of steady-state X versus the control parameter, points colored by stability. As a separate check, confirm that stable and unstable branches appear and merge at bifurcation points, so the number and stability of steady states change with the parameter. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3F.2.1  (03-phase-plane/03f-separatrix.qmd)

::: {.callout-note title="Recipe 3F.2.1"}
**Objective.** Trace the separatrix dividing the toggle switch's two basins of attraction.

**Model.** Genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`.

**Method.** Integrate the time-reversed system `dX/dt = -f(X)` from points just beside the saddle; the reversed flow carries them outward along the separatrix.

**Test.** The toggle switch (`gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`), seeding ten points near the saddle (about (190, 127)).

**Show.** A phase-plane plot of the separatrix together with sample forward trajectories.

**Verification.** The reversed trajectories trace a curve that no forward trajectory from random starts ever crosses, confirming it is the boundary between the two basins.
:::

**Generated prompt (Python):** In Python, trace the separatrix dividing the toggle switch's two basins of attraction. The model is genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`. Use integrate the time-reversed system `dX/dt = -f(X)` from points just beside the saddle; the reversed flow carries them outward along the separatrix. Test it on The toggle switch (`gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`), seeding ten points near the saddle (about (190, 127)). Produce a phase-plane plot of the separatrix together with sample forward trajectories. As a separate check, confirm that the reversed trajectories trace a curve that no forward trajectory from random starts ever crosses, confirming it is the boundary between the two basins. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3G.2.1  (03-phase-plane/03g-effective-potential-2d.qmd)

::: {.callout-note title="Recipe 3G.2.1"}
**Objective.** Compute an effective potential for the two-variable toggle switch along its nullclines, to locate the steady states as extrema.

**Model.** Genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`.

**Method.** Integrate the potential by the trapezoidal rule along each nullcline and plot the accumulated value against X and against Y.

**Test.** The toggle-switch nullclines with `gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`.

**Show.** Plots of the accumulated potential along each nullcline, versus X and versus Y.

**Verification.** Each curve dips to a minimum at the stable states and rises to a maximum at the saddle; the two paths do not agree on the potential values, since the flow is not a true gradient.
:::

**Generated prompt (Python):** In Python, compute an effective potential for the two-variable toggle switch along its nullclines, to locate the steady states as extrema. The model is genes X and Y repress each other; each gene's transcription is a basal rate plus a repressive Hill function of the other, with linear degradation: `dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X`, `dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y`. Use integrate the potential by the trapezoidal rule along each nullcline and plot the accumulated value against X and against Y. Test it on The toggle-switch nullclines with `gX0 = 5, gX1 = 50, Yth = 100, nY = 4, kX = 0.1`; `gY0 = 4, gY1 = 40, Xth = 150, nX = 4, kY = 0.12`. Produce plots of the accumulated potential along each nullcline, versus X and versus Y. As a separate check, confirm that each curve dips to a minimum at the stable states and rises to a maximum at the saddle; the two paths do not agree on the potential values, since the flow is not a true gradient. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3H.1.1  (03-phase-plane/03h-multi-component.qmd)

::: {.callout-note title="Recipe 3H.1.1"}
**Objective.** Integrate a system with any number of components using a single RK4 routine.

**Model.** Any autonomous system of ODEs `dX/dt = f(X)`, where the state X and the derivative f(X) are vectors of the same length.

**Method.** The generic `RK4_generic` from Part 3A, unchanged; it works whenever the derivative function returns a vector the same length as the state.

**Test.** Apply it to the two-, three-, and many-gene systems of this chapter.

**Show.** The integrated trajectories for the two-, three-, and many-gene systems.

**Verification.** The same integrator handles each system correctly with no modification.
:::

**Generated prompt (Python):** In Python, integrate a system with any number of components using a single RK4 routine. The model is any autonomous system of ODEs `dX/dt = f(X)`, where the state X and the derivative f(X) are vectors of the same length. Use the generic `RK4_generic` from Part 3A, unchanged; it works whenever the derivative function returns a vector the same length as the state. Test it on Apply it to the two-, three-, and many-gene systems of this chapter. Produce the integrated trajectories for the two-, three-, and many-gene systems. As a separate check, confirm that the same integrator handles each system correctly with no modification. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3H.2.1  (03-phase-plane/03h-multi-component.qmd)

::: {.callout-note title="Recipe 3H.2.1"}
**Objective.** Simulate a two-gene negative-feedback loop (X activates Y, Y represses X).

**Model.** A two-gene negative-feedback loop in nondimensional form, with Hill coefficient 3 and unit degradation: `dx/dt = g/(1 + y^3) - x`, `dy/dt = h*x^3/(1 + x^3) - y`.

**Method.** Integrate the nondimensionalized system with the generic RK4.

**Test.** `g = 10`, `h = 10`.

**Show.** A phase-plane plot of trajectories spiraling into the single steady state.

**Verification.** Trajectories converge to a single stable steady state, spiraling in.
:::

**Generated prompt (Python):** In Python, simulate a two-gene negative-feedback loop (X activates Y, Y represses X). The model is a two-gene negative-feedback loop in nondimensional form, with Hill coefficient 3 and unit degradation: `dx/dt = g/(1 + y^3) - x`, `dy/dt = h*x^3/(1 + x^3) - y`. Use integrate the nondimensionalized system with the generic RK4. Test it on `g = 10`, `h = 10`. Produce a phase-plane plot of trajectories spiraling into the single steady state. As a separate check, confirm that trajectories converge to a single stable steady state, spiraling in. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3H.4.1  (03-phase-plane/03h-multi-component.qmd)

::: {.callout-note title="Recipe 3H.4.1"}
**Objective.** Simulate the toggle switch (X and Y repress each other) in nondimensionalized form.

**Model.** The toggle switch in nondimensional form, Hill coefficient 3 and unit degradation: `dx/dt = g/(1 + y^3) - x`, `dy/dt = h/(1 + x^3) - y`.

**Method.** Integrate the system with the generic RK4.

**Test.** `g = 5`, `h = 5`, from several initial conditions.

**Show.** A phase-plane plot of trajectories reaching one of the two stable steady states.

**Verification.** There are two stable steady states (x-high/y-low and x-low/y-high); which one a trajectory reaches depends on its initial condition.
:::

**Generated prompt (Python):** In Python, simulate the toggle switch (X and Y repress each other) in nondimensionalized form. The model is the toggle switch in nondimensional form, Hill coefficient 3 and unit degradation: `dx/dt = g/(1 + y^3) - x`, `dy/dt = h/(1 + x^3) - y`. Use integrate the system with the generic RK4. Test it on `g = 5`, `h = 5`, from several initial conditions. Produce a phase-plane plot of trajectories reaching one of the two stable steady states. As a separate check, confirm that there are two stable steady states (x-high/y-low and x-low/y-high); which one a trajectory reaches depends on its initial condition. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3H.5.1  (03-phase-plane/03h-multi-component.qmd)

::: {.callout-note title="Recipe 3H.5.1"}
**Objective.** Simulate a three-gene repressilator (a ring of three repressions) and see sustained oscillations.

**Model.** A ring of three mutual repressions in nondimensional form (X repressed by Z, Y by X, Z by Y), Hill coefficient 3 and unit degradation: `dx/dt = g/(1 + z^3) - x`, `dy/dt = h/(1 + x^3) - y`, `dz/dt = l/(1 + y^3) - z`.

**Method.** Integrate the three-variable ring with the generic RK4.

**Test.** `g = h = l = 5`.

**Show.** A time-series plot of the three genes oscillating in turn.

**Verification.** The three genes settle into a sustained limit-cycle oscillation, each peaking in turn in a fixed order.
:::

**Generated prompt (Python):** In Python, simulate a three-gene repressilator (a ring of three repressions) and see sustained oscillations. The model is a ring of three mutual repressions in nondimensional form (X repressed by Z, Y by X, Z by Y), Hill coefficient 3 and unit degradation: `dx/dt = g/(1 + z^3) - x`, `dy/dt = h/(1 + x^3) - y`, `dz/dt = l/(1 + y^3) - z`. Use integrate the three-variable ring with the generic RK4. Test it on `g = h = l = 5`. Produce a time-series plot of the three genes oscillating in turn. As a separate check, confirm that the three genes settle into a sustained limit-cycle oscillation, each peaking in turn in a fixed order. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 3H.6.1  (03-phase-plane/03h-multi-component.qmd)

::: {.callout-note title="Recipe 3H.6.1"}
**Objective.** Simulate a large system of many interacting species with the same integrator.

**Model.** A generalized Lotka-Volterra community of S species with a small immigration term, `dN_i/dt = N_i*(1 - sum_j a_ij*N_j) + D`, with self-interaction `a_ii = 1` and random off-diagonal interactions `a_ij`.

**Method.** Build the interaction matrix, then integrate with the generic RK4.

**Test.** `S = 50` species, dispersal `D = 1e-6`, off-diagonal `a_ij` drawn from Uniform(0, 2a) with mean interaction strength a = 0.08 (weak), 0.16, and 0.64 (strong).

**Show.** A time-series plot of the S species relaxing to steady state, some coexisting and some near zero.

**Verification.** The same RK4 routine scales to S variables, and the community relaxes to a steady state in which some species coexist while others fall to near-extinction.
:::

**Generated prompt (Python):** In Python, simulate a large system of many interacting species with the same integrator. The model is a generalized Lotka-Volterra community of S species with a small immigration term, `dN_i/dt = N_i*(1 - sum_j a_ij*N_j) + D`, with self-interaction `a_ii = 1` and random off-diagonal interactions `a_ij`. Use build the interaction matrix, then integrate with the generic RK4. Test it on `S = 50` species, dispersal `D = 1e-6`, off-diagonal `a_ij` drawn from Uniform(0, 2a) with mean interaction strength a = 0.08 (weak), 0.16, and 0.64 (strong). Produce a time-series plot of the S species relaxing to steady state, some coexisting and some near zero. As a separate check, confirm that the same RK4 routine scales to S variables, and the community relaxes to a steady state in which some species coexist while others fall to near-extinction. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 4A.1.1  (04-time-delays/04a-delayed-differential-equations.qmd)

::: {.callout-note title="Recipe 4A.1.1"}
**Objective.** Integrate a delayed exponential-growth model and see how the delay changes the approach to the steady state.

**Model.** Exponential growth whose rate is set by the population one delay time in the past, `dN/dt = r*N(t - tau)` (the steady state is N = 0).

**Method.** The Euler method for delay differential equations, `N_next = N + dt*r*N_delayed` with the delayed value taken `tau/dt` steps back, storing the whole trajectory including the history interval.

**Test.** Constant history `N(t) = 1` for `t <= 0`; `tau = 1`, `dt = 0.01`, integrated to t = 40; growth rates r = -0.3, -1.4, and -1.7.

**Show.** A plot of N(t) for r = -0.3, -1.4, and -1.7, showing monotonic, damped, and growing behavior.

**Verification.** As r drops past the critical value `-pi/2` (about -1.57), the decay of the N = 0 state changes from monotonic (r = -0.3) to a damped oscillation (r = -1.4) to an oscillation that grows without bound (r = -1.7).
:::

**Generated prompt (Python):** In Python, integrate a delayed exponential-growth model and see how the delay changes the approach to the steady state. The model is exponential growth whose rate is set by the population one delay time in the past, `dN/dt = r*N(t - tau)` (the steady state is N = 0). Use the Euler method for delay differential equations, `N_next = N + dt*r*N_delayed` with the delayed value taken `tau/dt` steps back, storing the whole trajectory including the history interval. Test it on Constant history `N(t) = 1` for `t <= 0`; `tau = 1`, `dt = 0.01`, integrated to t = 40; growth rates r = -0.3, -1.4, and -1.7. Produce a plot of N(t) for r = -0.3, -1.4, and -1.7, showing monotonic, damped, and growing behavior. As a separate check, confirm that as r drops past the critical value `-pi/2` (about -1.57), the decay of the N = 0 state changes from monotonic (r = -0.3) to a damped oscillation (r = -1.4) to an oscillation that grows without bound (r = -1.7). Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 4A.2.1  (04-time-delays/04a-delayed-differential-equations.qmd)

::: {.callout-note title="Recipe 4A.2.1"}
**Objective.** Integrate the same delayed model with a second-order method and compare it to the first-order result.

**Model.** Delayed exponential growth `dN/dt = r*N(t - tau)`.

**Method.** Heun's method for delay differential equations: an Euler predictor followed by a trapezoidal corrector, taken over whole time steps so the stored delayed values suffice.

**Test.** Constant history `N(t) = 1`; `tau = 1`, `dt = 0.01`, to t = 40; the unstable case r = -1.7, compared against the Euler result.

**Show.** A plot of the Euler and Heun trajectories for r = -1.7 overlaid.

**Verification.** The Euler and Heun trajectories agree early but drift apart later, with the first-order Euler curve overshooting the second-order one.
:::

**Generated prompt (Python):** In Python, integrate the same delayed model with a second-order method and compare it to the first-order result. The model is delayed exponential growth `dN/dt = r*N(t - tau)`. Use heun's method for delay differential equations: an Euler predictor followed by a trapezoidal corrector, taken over whole time steps so the stored delayed values suffice. Test it on Constant history `N(t) = 1`; `tau = 1`, `dt = 0.01`, to t = 40; the unstable case r = -1.7, compared against the Euler result. Produce a plot of the Euler and Heun trajectories for r = -1.7 overlaid. As a separate check, confirm that the Euler and Heun trajectories agree early but drift apart later, with the first-order Euler curve overshooting the second-order one. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 4B.1.1  (04-time-delays/04b-modeling-time-delays.qmd)

::: {.callout-note title="Recipe 4B.1.1"}
**Objective.** Simulate delay-driven logistic growth and locate the onset of sustained oscillation.

**Model.** The Hutchinson delay-logistic equation, logistic growth whose crowding term uses the population one delay in the past, `dN/dt = r*N(t)*(1 - N(t-tau)/B)`, with carrying capacity B.

**Method.** The second-order Heun integrator for delay differential equations.

**Test.** Constant history `N(t) = 1`; `tau = 1`, `B = 100`, `dt = 0.01`; growth rates r = 0.3, 1.5, pi/2 (about 1.57), and 1.7.

**Show.** A plot of N(t) for each r, from monotonic growth to a sustained oscillation.

**Verification.** Crossing the Hopf point `r*tau = pi/2`: r = 0.3 grows monotonically to B, r = 1.5 is a damped oscillation settling at B, r = pi/2 decays extremely slowly, and r = 1.7 becomes a sustained oscillation (a stable limit cycle).
:::

**Generated prompt (Python):** In Python, simulate delay-driven logistic growth and locate the onset of sustained oscillation. The model is the Hutchinson delay-logistic equation, logistic growth whose crowding term uses the population one delay in the past, `dN/dt = r*N(t)*(1 - N(t-tau)/B)`, with carrying capacity B. Use the second-order Heun integrator for delay differential equations. Test it on Constant history `N(t) = 1`; `tau = 1`, `B = 100`, `dt = 0.01`; growth rates r = 0.3, 1.5, pi/2 (about 1.57), and 1.7. Produce a plot of N(t) for each r, from monotonic growth to a sustained oscillation. As a separate check, confirm that crossing the Hopf point `r*tau = pi/2`: r = 0.3 grows monotonically to B, r = 1.5 is a damped oscillation settling at B, r = pi/2 decays extremely slowly, and r = 1.7 becomes a sustained oscillation (a stable limit cycle). Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 4B.2.1  (04-time-delays/04b-modeling-time-delays.qmd)

::: {.callout-note title="Recipe 4B.2.1"}
**Objective.** Simulate a self-repressing gene whose feedback acts after a delay, and watch oscillations emerge as the delay grows.

**Model.** Delayed negative autoregulation: basal plus repressive Hill production driven by the protein level a delay earlier, minus linear degradation, `dX/dt = g0 + g1/(1 + (X(t-tau)/Xth)^n) - k*X(t)`.

**Method.** The second-order Heun integrator for delay differential equations.

**Test.** Constant history `X(t) = 1`; `g0 = 10, g1 = 60, Xth = 200, n = 4, k = 0.1`, `dt = 0.01`, to t = 200; delay swept `tau = 5, 10, 15, 20`.

**Show.** A plot of X(t) for each delay tau, from a steady state to a sustained oscillation.

**Verification.** A short delay (tau = 5) settles at the steady state near 250; larger delays turn the response into a damped and then a sustained oscillation.
:::

**Generated prompt (Python):** In Python, simulate a self-repressing gene whose feedback acts after a delay, and watch oscillations emerge as the delay grows. The model is delayed negative autoregulation: basal plus repressive Hill production driven by the protein level a delay earlier, minus linear degradation, `dX/dt = g0 + g1/(1 + (X(t-tau)/Xth)^n) - k*X(t)`. Use the second-order Heun integrator for delay differential equations. Test it on Constant history `X(t) = 1`; `g0 = 10, g1 = 60, Xth = 200, n = 4, k = 0.1`, `dt = 0.01`, to t = 200; delay swept `tau = 5, 10, 15, 20`. Produce a plot of X(t) for each delay tau, from a steady state to a sustained oscillation. As a separate check, confirm that a short delay (tau = 5) settles at the steady state near 250; larger delays turn the response into a damped and then a sustained oscillation. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 4B.3.1  (04-time-delays/04b-modeling-time-delays.qmd)

::: {.callout-note title="Recipe 4B.3.1"}
**Objective.** Simulate a predator-prey model in which each species responds to the other after a lag.

**Model.** The Lotka-Volterra system with delayed cross-terms, `dN/dt = N(t)*(a - b*P(t-tau))`, `dP/dt = P(t)*(c*N(t-tau) - d)`.

**Method.** The generic multi-variable Heun integrator for delay differential equations, carrying a vector state and vector history.

**Test.** `a = 1, b = 0.03, c = 0.02, d = 1`; initial (N, P) = (30, 10); `dt = 0.01`, to t = 50; delay `tau = 0` and `tau = 0.01`.

**Show.** A phase-plane plot of the predator-prey orbit for tau = 0 and tau = 0.01.

**Verification.** Without delay (tau = 0) the orbit is a neutrally stable closed loop; even a one-step delay (tau = 0.01) makes it spiral outward.
:::

**Generated prompt (Python):** In Python, simulate a predator-prey model in which each species responds to the other after a lag. The model is the Lotka-Volterra system with delayed cross-terms, `dN/dt = N(t)*(a - b*P(t-tau))`, `dP/dt = P(t)*(c*N(t-tau) - d)`. Use the generic multi-variable Heun integrator for delay differential equations, carrying a vector state and vector history. Test it on `a = 1, b = 0.03, c = 0.02, d = 1`; initial (N, P) = (30, 10); `dt = 0.01`, to t = 50; delay `tau = 0` and `tau = 0.01`. Produce a phase-plane plot of the predator-prey orbit for tau = 0 and tau = 0.01. As a separate check, confirm that without delay (tau = 0) the orbit is a neutrally stable closed loop; even a one-step delay (tau = 0.01) makes it spiral outward. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 4C.1.1  (04-time-delays/04c-indirect-interactions.qmd)

::: {.callout-note title="Recipe 4C.1.1"}
**Objective.** Simulate a two-gene negative-feedback loop without any delay, as a baseline.

**Model.** A nondimensional two-node loop where X activates Y and Y represses X (Hill coefficient 3, unit degradation), `dx/dt = g/(1 + y^3) - x`, `dy/dt = h*x^3/(1 + x^3) - y`.

**Method.** The generic multi-variable Heun integrator (ordinary, no delay).

**Test.** `g = 10`, `h = 10`; initial (x, y) = (1, 1); `dt = 0.01`, to t = 10.

**Show.** A time-series plot of x and y relaxing to a steady state.

**Verification.** With no delay the loop simply relaxes to a stable steady state, with no oscillation.
:::

**Generated prompt (Python):** In Python, simulate a two-gene negative-feedback loop without any delay, as a baseline. The model is a nondimensional two-node loop where X activates Y and Y represses X (Hill coefficient 3, unit degradation), `dx/dt = g/(1 + y^3) - x`, `dy/dt = h*x^3/(1 + x^3) - y`. Use the generic multi-variable Heun integrator (ordinary, no delay). Test it on `g = 10`, `h = 10`; initial (x, y) = (1, 1); `dt = 0.01`, to t = 10. Produce a time-series plot of x and y relaxing to a steady state. As a separate check, confirm that with no delay the loop simply relaxes to a stable steady state, with no oscillation. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 4C.2.1  (04-time-delays/04c-indirect-interactions.qmd)

::: {.callout-note title="Recipe 4C.2.1"}
**Objective.** Add a delay to the two-node loop's repression and see it destabilize.

**Model.** The same two-node loop as 4C.1, but the repression of X by Y acts after a delay, `dx/dt = g/(1 + y(t-tau)^3) - x`, `dy/dt = h*x^3/(1 + x^3) - y`.

**Method.** The generic multi-variable Heun integrator for delay differential equations.

**Test.** `g = 10`, `h = 10`, `tau = 2`; constant history (x, y) = (1, 1); `dt = 0.01`, to t = 30.

**Show.** A time-series plot of x and y oscillating once the delay is added.

**Verification.** The delay destabilizes the steady state that was stable without it, giving a sustained oscillation.
:::

**Generated prompt (Python):** In Python, add a delay to the two-node loop's repression and see it destabilize. The model is the same two-node loop as 4C.1, but the repression of X by Y acts after a delay, `dx/dt = g/(1 + y(t-tau)^3) - x`, `dy/dt = h*x^3/(1 + x^3) - y`. Use the generic multi-variable Heun integrator for delay differential equations. Test it on `g = 10`, `h = 10`, `tau = 2`; constant history (x, y) = (1, 1); `dt = 0.01`, to t = 30. Produce a time-series plot of x and y oscillating once the delay is added. As a separate check, confirm that the delay destabilizes the steady state that was stable without it, giving a sustained oscillation. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 4C.3.1  (04-time-delays/04c-indirect-interactions.qmd)

::: {.callout-note title="Recipe 4C.3.1"}
**Objective.** Show that a chain of intermediate genes can stand in for a delay in producing oscillations.

**Model.** Undelayed repression rings that replace the delay with explicit intermediate genes. Three-node ring (X activates Y activates Z, Z represses X): `dx/dt = g/(1 + z^3) - x`, `dy/dt = h*x^3/(1 + x^3) - y`, `dz/dt = l*y^3/(1 + y^3) - z`. Four-node ring: the same chain with a fourth gene W between Z and X, `dx/dt = g/(1 + w^3) - x`, ..., `dw/dt = m*z^3/(1 + z^3) - w`. Every edge is a Hill interaction (coefficient 3) with unit degradation.

**Method.** The generic multi-variable Heun integrator (ordinary, no delay).

**Test.** All maximal rates equal to 10 (`g = h = l = 10`, and `m = 10` for four nodes); initial state all ones; `dt = 0.01`, to t = 30.

**Show.** Time-series plots of the three- and four-node rings, showing the oscillation and its lengthening period.

**Verification.** The three-node ring already oscillates without any delay; adding a fourth node lengthens the period, approaching the delayed two-node loop, so a delay and a chain of intermediate reactions capture the same lag.
:::

**Generated prompt (Python):** In Python, show that a chain of intermediate genes can stand in for a delay in producing oscillations. The model is undelayed repression rings that replace the delay with explicit intermediate genes. Three-node ring (X activates Y activates Z, Z represses X): `dx/dt = g/(1 + z^3) - x`, `dy/dt = h*x^3/(1 + x^3) - y`, `dz/dt = l*y^3/(1 + y^3) - z`. Four-node ring: the same chain with a fourth gene W between Z and X, `dx/dt = g/(1 + w^3) - x`, ..., `dw/dt = m*z^3/(1 + z^3) - w`. Every edge is a Hill interaction (coefficient 3) with unit degradation. Use the generic multi-variable Heun integrator (ordinary, no delay). Test it on All maximal rates equal to 10 (`g = h = l = 10`, and `m = 10` for four nodes); initial state all ones; `dt = 0.01`, to t = 30. Produce time-series plots of the three- and four-node rings, showing the oscillation and its lengthening period. As a separate check, confirm that the three-node ring already oscillates without any delay; adding a fourth node lengthens the period, approaching the delayed two-node loop, so a delay and a chain of intermediate reactions capture the same lag. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5A.1.1  (05-molecular-dynamics/05a-second-order-integrators.qmd)

::: {.callout-note title="Recipe 5A.1.1"}
**Objective.** Set up a harmonic oscillator as the test system for comparing integrators, and state what a correct simulation must preserve.

**Model.** A mass on a spring of stiffness k, `d^2x/dt^2 = -k*x`, written as the pair `dx/dt = v`, `dv/dt = -k*x`; the energy per unit mass `e = 0.5*k*x^2 + 0.5*v^2` is exactly constant for the true motion.

**Method.** Direct evaluation of the force `f(x) = -k*x`; the numerical integrators come in the following sections.

**Test.** `k = 0.1`; initial displacement `x0 = 1`, velocity `v0 = 2`.

**Show.** A plot of the exact oscillator motion x(t) and its constant energy, as the reference to match.

**Verification.** The exact motion is a pure sinusoid at constant energy, so any amplitude growth or energy drift in a simulation is a numerical artifact.
:::

**Generated prompt (Python):** In Python, set up a harmonic oscillator as the test system for comparing integrators, and state what a correct simulation must preserve. The model is a mass on a spring of stiffness k, `d^2x/dt^2 = -k*x`, written as the pair `dx/dt = v`, `dv/dt = -k*x`; the energy per unit mass `e = 0.5*k*x^2 + 0.5*v^2` is exactly constant for the true motion. Use direct evaluation of the force `f(x) = -k*x`; the numerical integrators come in the following sections. Test it on `k = 0.1`; initial displacement `x0 = 1`, velocity `v0 = 2`. Produce a plot of the exact oscillator motion x(t) and its constant energy, as the reference to match. As a separate check, confirm that the exact motion is a pure sinusoid at constant energy, so any amplitude growth or energy drift in a simulation is a numerical artifact. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5A.2.1  (05-molecular-dynamics/05a-second-order-integrators.qmd)

::: {.callout-note title="Recipe 5A.2.1"}
**Objective.** Integrate the harmonic oscillator with the Euler method and check whether it conserves energy.

**Model.** A mass on a spring of stiffness k, `d^2x/dt^2 = -k*x`, written as the pair `dx/dt = v`, `dv/dt = -k*x`; the energy per unit mass `e = 0.5*k*x^2 + 0.5*v^2` is exactly constant for the true motion.

**Method.** The explicit (forward) Euler integrator for Newton's equations, updating velocity then position, `v_next = v + dt*f`, `x_next = x + dt*v`.

**Test.** `k = 0.1`, initial `x0 = 1`, `v0 = 2`, to t = 100; step sizes `dt = 0.01` and `dt = 0.1`.

**Show.** Plots of x(t) and the energy e(t) at dt = 0.01 and dt = 0.1.

**Verification.** Euler does not conserve energy: at dt = 0.1 the amplitude and energy grow steadily, and even at dt = 0.01 the energy drifts slowly upward.
:::

**Generated prompt (Python):** In Python, integrate the harmonic oscillator with the Euler method and check whether it conserves energy. The model is a mass on a spring of stiffness k, `d^2x/dt^2 = -k*x`, written as the pair `dx/dt = v`, `dv/dt = -k*x`; the energy per unit mass `e = 0.5*k*x^2 + 0.5*v^2` is exactly constant for the true motion. Use the explicit (forward) Euler integrator for Newton's equations, updating velocity then position, `v_next = v + dt*f`, `x_next = x + dt*v`. Test it on `k = 0.1`, initial `x0 = 1`, `v0 = 2`, to t = 100; step sizes `dt = 0.01` and `dt = 0.1`. Produce plots of x(t) and the energy e(t) at dt = 0.01 and dt = 0.1. As a separate check, confirm that euler does not conserve energy: at dt = 0.1 the amplitude and energy grow steadily, and even at dt = 0.01 the energy drifts slowly upward. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5A.3.1  (05-molecular-dynamics/05a-second-order-integrators.qmd)

::: {.callout-note title="Recipe 5A.3.1"}
**Objective.** Integrate the harmonic oscillator with the leapfrog method and compare its energy behavior to Euler.

**Model.** A mass on a spring of stiffness k, `d^2x/dt^2 = -k*x`, written as the pair `dx/dt = v`, `dv/dt = -k*x`; the energy per unit mass `e = 0.5*k*x^2 + 0.5*v^2` is exactly constant for the true motion.

**Method.** The leapfrog integrator, carrying velocity at half-integer steps: a startup half-step `v_half = v0 + 0.5*dt*f0`, then `x_next = x + dt*v_half` and `v_next_half = v_half + dt*f_next`.

**Test.** `k = 0.1`, initial `x0 = 1`, `v0 = 2`, to t = 100; step sizes `dt = 0.01` and `dt = 0.1`.

**Show.** Plots of x(t) and the energy e(t) for leapfrog at both step sizes.

**Verification.** Leapfrog energy does not drift even at the large step (it is symplectic and time-reversible); a small energy ripple at dt = 0.1 comes only from x and v being stored at staggered times.
:::

**Generated prompt (Python):** In Python, integrate the harmonic oscillator with the leapfrog method and compare its energy behavior to Euler. The model is a mass on a spring of stiffness k, `d^2x/dt^2 = -k*x`, written as the pair `dx/dt = v`, `dv/dt = -k*x`; the energy per unit mass `e = 0.5*k*x^2 + 0.5*v^2` is exactly constant for the true motion. Use the leapfrog integrator, carrying velocity at half-integer steps: a startup half-step `v_half = v0 + 0.5*dt*f0`, then `x_next = x + dt*v_half` and `v_next_half = v_half + dt*f_next`. Test it on `k = 0.1`, initial `x0 = 1`, `v0 = 2`, to t = 100; step sizes `dt = 0.01` and `dt = 0.1`. Produce plots of x(t) and the energy e(t) for leapfrog at both step sizes. As a separate check, confirm that leapfrog energy does not drift even at the large step (it is symplectic and time-reversible); a small energy ripple at dt = 0.1 comes only from x and v being stored at staggered times. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5A.4.1  (05-molecular-dynamics/05a-second-order-integrators.qmd)

::: {.callout-note title="Recipe 5A.4.1"}
**Objective.** Integrate the harmonic oscillator with velocity Verlet, the standard molecular-dynamics integrator.

**Model.** A mass on a spring of stiffness k, `d^2x/dt^2 = -k*x`, written as the pair `dx/dt = v`, `dv/dt = -k*x`; the energy per unit mass `e = 0.5*k*x^2 + 0.5*v^2` is exactly constant for the true motion.

**Method.** The velocity-Verlet integrator, splitting the velocity update into two half-steps around the position update, `v_half = v + 0.5*dt*f`, `x_next = x + dt*v_half`, `v_next = v_half + 0.5*dt*f_next`.

**Test.** `k = 0.1`, initial `x0 = 1`, `v0 = 2`, to t = 100; step sizes `dt = 0.01` and `dt = 0.1`.

**Show.** Plots of x(t) and the energy e(t) for velocity Verlet at both step sizes.

**Verification.** Velocity Verlet conserves the energy at both step sizes and returns position and velocity at the same time points.
:::

**Generated prompt (Python):** In Python, integrate the harmonic oscillator with velocity Verlet, the standard molecular-dynamics integrator. The model is a mass on a spring of stiffness k, `d^2x/dt^2 = -k*x`, written as the pair `dx/dt = v`, `dv/dt = -k*x`; the energy per unit mass `e = 0.5*k*x^2 + 0.5*v^2` is exactly constant for the true motion. Use the velocity-Verlet integrator, splitting the velocity update into two half-steps around the position update, `v_half = v + 0.5*dt*f`, `x_next = x + dt*v_half`, `v_next = v_half + 0.5*dt*f_next`. Test it on `k = 0.1`, initial `x0 = 1`, `v0 = 2`, to t = 100; step sizes `dt = 0.01` and `dt = 0.1`. Produce plots of x(t) and the energy e(t) for velocity Verlet at both step sizes. As a separate check, confirm that velocity Verlet conserves the energy at both step sizes and returns position and velocity at the same time points. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5A.5.1  (05-molecular-dynamics/05a-second-order-integrators.qmd)

::: {.callout-note title="Recipe 5A.5.1"}
**Objective.** Integrate the harmonic oscillator with the original position-only Verlet method.

**Model.** A mass on a spring of stiffness k, `d^2x/dt^2 = -k*x`, written as the pair `dx/dt = v`, `dv/dt = -k*x`; the energy per unit mass `e = 0.5*k*x^2 + 0.5*v^2` is exactly constant for the true motion. This method does not store the velocity.

**Method.** The Verlet integrator, advancing position from the two previous positions, `x_{n+2} = 2*x_{n+1} - x_n + dt^2*f_{n+1}`, with a startup step `x_1 = x0 + dt*v0 + 0.5*dt^2*f0`.

**Test.** `k = 0.1`, initial `x0 = 1`, `v0 = 2`, to t = 100; step sizes `dt = 0.01` and `dt = 0.1`.

**Show.** A plot of x(t) for Verlet at both step sizes.

**Verification.** Verlet reproduces the oscillation at both step sizes while storing only positions, which is why it and velocity Verlet underlie most practical MD codes.
:::

**Generated prompt (Python):** In Python, integrate the harmonic oscillator with the original position-only Verlet method. The model is a mass on a spring of stiffness k, `d^2x/dt^2 = -k*x`, written as the pair `dx/dt = v`, `dv/dt = -k*x`; the energy per unit mass `e = 0.5*k*x^2 + 0.5*v^2` is exactly constant for the true motion. This method does not store the velocity. Use the Verlet integrator, advancing position from the two previous positions, `x_{n+2} = 2*x_{n+1} - x_n + dt^2*f_{n+1}`, with a startup step `x_1 = x0 + dt*v0 + 0.5*dt^2*f0`. Test it on `k = 0.1`, initial `x0 = 1`, `v0 = 2`, to t = 100; step sizes `dt = 0.01` and `dt = 0.1`. Produce a plot of x(t) for Verlet at both step sizes. As a separate check, confirm that verlet reproduces the oscillation at both step sizes while storing only positions, which is why it and velocity Verlet underlie most practical MD codes. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5B.1.1  (05-molecular-dynamics/05b-orbital-motion.qmd)

::: {.callout-note title="Recipe 5B.1.1"}
**Objective.** Simulate a single particle orbiting under a central attractive force.

**Model.** A central force of magnitude `G/r` pointing toward the origin, with Cartesian components `fx = -G*x/r^2`, `fy = -G*y/r^2`, and `r^2 = x^2 + y^2`.

**Method.** The generic vector velocity-Verlet integrator applied to one 2D particle.

**Test.** `G = 1`, `dt = 0.1`; two runs, initial position (4, 0) with velocity (0, 1) to t = 100, and (2, 2) with velocity (-1, 1) to t = 1000.

**Show.** Orbit plots in the plane for both initial conditions (the ellipse and the rosette).

**Verification.** The first initial condition gives a closed, nearly elliptical orbit; the second traces a precessing rosette that fills an annulus and never quite closes.
:::

**Generated prompt (Python):** In Python, simulate a single particle orbiting under a central attractive force. The model is a central force of magnitude `G/r` pointing toward the origin, with Cartesian components `fx = -G*x/r^2`, `fy = -G*y/r^2`, and `r^2 = x^2 + y^2`. Use the generic vector velocity-Verlet integrator applied to one 2D particle. Test it on `G = 1`, `dt = 0.1`; two runs, initial position (4, 0) with velocity (0, 1) to t = 100, and (2, 2) with velocity (-1, 1) to t = 1000. Produce orbit plots in the plane for both initial conditions (the ellipse and the rosette). As a separate check, confirm that the first initial condition gives a closed, nearly elliptical orbit; the second traces a precessing rosette that fills an annulus and never quite closes. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5B.2.1  (05-molecular-dynamics/05b-orbital-motion.qmd)

::: {.callout-note title="Recipe 5B.2.1"}
**Objective.** Simulate two equal masses orbiting under mutual attraction, and remove their overall drift.

**Model.** Two unit masses attracting each other with a `1/r`-type central force, `f1x = G*(x2-x1)/r12^2` and the equal-and-opposite `f2x = G*(x1-x2)/r12^2` (likewise in y), with `r12^2 = (x1-x2)^2 + (y1-y2)^2`.

**Method.** The generic vector velocity-Verlet integrator on the four coordinates, with a shift to the center-of-mass frame (subtract the mean velocity per axis).

**Test.** `G = 1`, `dt = 0.01`, to t = 100; initial positions (2, 0, -2, 0), velocities (0, 0.4, 0, -0.2); a second run centers the velocities to zero total momentum.

**Show.** Orbit plots before and after shifting to the center-of-mass frame.

**Verification.** With nonzero total momentum the pair orbits while drifting across the plane; subtracting the mean velocity yields a stationary pattern about the common center of mass.
:::

**Generated prompt (Python):** In Python, simulate two equal masses orbiting under mutual attraction, and remove their overall drift. The model is two unit masses attracting each other with a `1/r`-type central force, `f1x = G*(x2-x1)/r12^2` and the equal-and-opposite `f2x = G*(x1-x2)/r12^2` (likewise in y), with `r12^2 = (x1-x2)^2 + (y1-y2)^2`. Use the generic vector velocity-Verlet integrator on the four coordinates, with a shift to the center-of-mass frame (subtract the mean velocity per axis). Test it on `G = 1`, `dt = 0.01`, to t = 100; initial positions (2, 0, -2, 0), velocities (0, 0.4, 0, -0.2); a second run centers the velocities to zero total momentum. Produce orbit plots before and after shifting to the center-of-mass frame. As a separate check, confirm that with nonzero total momentum the pair orbits while drifting across the plane; subtracting the mean velocity yields a stationary pattern about the common center of mass. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5B.3.1  (05-molecular-dynamics/05b-orbital-motion.qmd)

::: {.callout-note title="Recipe 5B.3.1"}
**Objective.** Simulate the chaotic three-body problem.

**Model.** Three unit masses, each attracting every other by a `1/r^2`-scaled central force summed over partners, `f1x = G*((x2-x1)/r12^2 + (x3-x1)/r13^2)` (likewise for the other bodies and for y).

**Method.** The generic vector velocity-Verlet integrator on the six coordinates, started in the center-of-mass frame.

**Test.** `G = 1`, `dt = 0.01`, to t = 100; initial positions (-2, 0, 2, 0, 0, 3), velocities centered from (0.1, 0, -0.1, 0, 0, -0.05).

**Show.** An orbit plot of the three bodies, up to the escape.

**Verification.** The three bodies stay bound and weave for a while, then one gains enough energy to escape while the remaining pair recoils, the hallmark of three-body chaos.
:::

**Generated prompt (Python):** In Python, simulate the chaotic three-body problem. The model is three unit masses, each attracting every other by a `1/r^2`-scaled central force summed over partners, `f1x = G*((x2-x1)/r12^2 + (x3-x1)/r13^2)` (likewise for the other bodies and for y). Use the generic vector velocity-Verlet integrator on the six coordinates, started in the center-of-mass frame. Test it on `G = 1`, `dt = 0.01`, to t = 100; initial positions (-2, 0, 2, 0, 0, 3), velocities centered from (0.1, 0, -0.1, 0, 0, -0.05). Produce an orbit plot of the three bodies, up to the escape. As a separate check, confirm that the three bodies stay bound and weave for a while, then one gains enough energy to escape while the remaining pair recoils, the hallmark of three-body chaos. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5C.4.1  (05-molecular-dynamics/05c-box-of-particles.qmd)

::: {.callout-note title="Recipe 5C.4.1"}
**Objective.** Build a velocity-Verlet integrator that keeps particles inside a periodic box.

**Model.** Each unit-mass particle obeys `d^2x_i/dt^2 = sum over j of the pair force`, where the pair force is the Lennard-Jones force `F(r) = 1/r^13 - 1/r^7` projected along the minimum-image separation between particles i and j.

**Method.** The velocity-Verlet integrator with a periodic-boundary correction `x -> x - floor(x/a)*a`, applied to positions (not velocities) after each position update.

**Test.** Box side `a = 6.25` with 25 particles (the run in 5C.5); the folding rule is `bc_periodic(x, a) = x - floor(x/a)*a`.

**Show.** A short check that a particle stepping past a wall reappears on the opposite side.

**Verification.** A particle that leaves one side of the box reappears on the opposite side, while velocities are never wrapped.
:::

**Generated prompt (Python):** In Python, build a velocity-Verlet integrator that keeps particles inside a periodic box. The model is each unit-mass particle obeys `d^2x_i/dt^2 = sum over j of the pair force`, where the pair force is the Lennard-Jones force `F(r) = 1/r^13 - 1/r^7` projected along the minimum-image separation between particles i and j. Use the velocity-Verlet integrator with a periodic-boundary correction `x -> x - floor(x/a)*a`, applied to positions (not velocities) after each position update. Test it on Box side `a = 6.25` with 25 particles (the run in 5C.5); the folding rule is `bc_periodic(x, a) = x - floor(x/a)*a`. Produce a short check that a particle stepping past a wall reappears on the opposite side. As a separate check, confirm that a particle that leaves one side of the box reappears on the opposite side, while velocities are never wrapped. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5C.5.1  (05-molecular-dynamics/05c-box-of-particles.qmd)

::: {.callout-note title="Recipe 5C.5.1"}
**Objective.** Simulate a 2D box of Lennard-Jones particles from an ordered start into a disordered liquid.

**Model.** 25 unit-mass particles in a square box of side `a = 6.25` interacting through the Lennard-Jones force `F(r) = 1/r^13 - 1/r^7` with minimum-image periodic boundaries.

**Method.** The velocity-Verlet integrator with periodic boundaries (from 5C.4).

**Test.** 25 particles on a 5x5 grid at spacing 1.25; deterministic initial velocities spread over (-0.05, 0.05) by a golden-ratio sequence; a warm-up run to t = 10 then a production run to t = 100, both at `dt = 0.01`.

**Show.** Snapshots of the particle positions in the box, from the initial grid to the final disordered state.

**Verification.** Particles start on the regular grid and end in an irregular liquid-like arrangement; because the many-body dynamics are chaotic, the exact final configuration differs between R and Python while structural averages match.
:::

**Generated prompt (Python):** In Python, simulate a 2D box of Lennard-Jones particles from an ordered start into a disordered liquid. The model is 25 unit-mass particles in a square box of side `a = 6.25` interacting through the Lennard-Jones force `F(r) = 1/r^13 - 1/r^7` with minimum-image periodic boundaries. Use the velocity-Verlet integrator with periodic boundaries (from 5C.4). Test it on 25 particles on a 5x5 grid at spacing 1.25; deterministic initial velocities spread over (-0.05, 0.05) by a golden-ratio sequence; a warm-up run to t = 10 then a production run to t = 100, both at `dt = 0.01`. Produce snapshots of the particle positions in the box, from the initial grid to the final disordered state. As a separate check, confirm that particles start on the regular grid and end in an irregular liquid-like arrangement; because the many-body dynamics are chaotic, the exact final configuration differs between R and Python while structural averages match. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 5C.7.1  (05-molecular-dynamics/05c-box-of-particles.qmd)

::: {.callout-note title="Recipe 5C.7.1"}
**Objective.** Measure the radial distribution function of the simulated particle box.

**Model.** The radial distribution function `g(r) = <dN(r)> / (2*pi*r*dr*rho0)`, where `rho0 = N/a^2` is the mean density and `<dN(r)>` is the average number of particles in a shell of radius r and thickness dr around a reference particle.

**Method.** The radial distribution function g(r), computed by histogramming minimum-image pair distances and normalizing each bin by its shell area `2*pi*r*dr` and the mean density.

**Test.** The t = 100 trajectory from 5C.5 (25 particles, box side `a = 6.25`); bin width `dr = 0.05`, radii up to `a/2 = 3.125`.

**Show.** A plot of `g(r)` versus r, with its first and second peaks.

**Verification.** g(r) has the classic liquid shape: near zero below about r = 0.8, a sharp first peak near r = 1 (the Lennard-Jones equilibrium separation), a weaker second peak near r = 2, settling toward 1 at large r.
:::

**Generated prompt (Python):** In Python, measure the radial distribution function of the simulated particle box. The model is the radial distribution function `g(r) = <dN(r)> / (2*pi*r*dr*rho0)`, where `rho0 = N/a^2` is the mean density and `<dN(r)>` is the average number of particles in a shell of radius r and thickness dr around a reference particle. Use the radial distribution function g(r), computed by histogramming minimum-image pair distances and normalizing each bin by its shell area `2*pi*r*dr` and the mean density. Test it on The t = 100 trajectory from 5C.5 (25 particles, box side `a = 6.25`); bin width `dr = 0.05`, radii up to `a/2 = 3.125`. Produce a plot of `g(r)` versus r, with its first and second peaks. As a separate check, confirm that g(r) has the classic liquid shape: near zero below about r = 0.8, a sharp first peak near r = 1 (the Lennard-Jones equilibrium separation), a weaker second peak near r = 2, settling toward 1 at large r. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6A.3.1  (06-stochastic/06a-random-number-generators.qmd)

::: {.callout-note title="Recipe 6A.3.1"}
**Objective.** Generate exponentially distributed random numbers from uniform ones.

**Model.** The target is the exponential distribution on `y >= 0`.

**Method.** Inverse-transform sampling of the exponential distribution: for a uniform `x` on (0, 1), take `y = -ln(x)`.

**Test.** `n = 10000` uniform draws on (0, 1), transformed and shown as a density histogram.

**Show.** A density histogram of the samples against the exponential density.

**Verification.** The flat uniform input becomes the characteristic decaying exponential density after the transform.
:::

**Generated prompt (Python):** In Python, generate exponentially distributed random numbers from uniform ones. The model is the target is the exponential distribution on `y >= 0`. Use inverse-transform sampling of the exponential distribution: for a uniform `x` on (0, 1), take `y = -ln(x)`. Test it on `n = 10000` uniform draws on (0, 1), transformed and shown as a density histogram. Produce a density histogram of the samples against the exponential density. As a separate check, confirm that the flat uniform input becomes the characteristic decaying exponential density after the transform. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6A.4.1  (06-stochastic/06a-random-number-generators.qmd)

::: {.callout-note title="Recipe 6A.4.1"}
**Objective.** Generate standard-normal random numbers from uniform ones.

**Model.** The target is the standard normal distribution `N(0, 1)`.

**Method.** The polar (Marsaglia) Box-Muller method: draw a point uniform in the square [-1, 1]^2, keep it if `R2 = x^2 + y^2` lies in (0, 1), and return `x*sqrt(-2*ln(R2)/R2)` and `y*sqrt(-2*ln(R2)/R2)` as two independent normals.

**Test.** `n = 10000` accepted draws, shown as a density histogram.

**Show.** A density histogram of the samples against the standard normal density.

**Verification.** Both output streams follow the bell-shaped normal density.
:::

**Generated prompt (Python):** In Python, generate standard-normal random numbers from uniform ones. The model is the target is the standard normal distribution `N(0, 1)`. Use the polar (Marsaglia) Box-Muller method: draw a point uniform in the square [-1, 1]^2, keep it if `R2 = x^2 + y^2` lies in (0, 1), and return `x*sqrt(-2*ln(R2)/R2)` and `y*sqrt(-2*ln(R2)/R2)` as two independent normals. Test it on `n = 10000` accepted draws, shown as a density histogram. Produce a density histogram of the samples against the standard normal density. As a separate check, confirm that both output streams follow the bell-shaped normal density. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6B.1.1  (06-stochastic/06b-brownian-motion.qmd)

::: {.callout-note title="Recipe 6B.1.1"}
**Objective.** Simulate an ensemble of one-dimensional random walks and see how they spread.

**Model.** A walker starts at 0 and at each step moves `+dx` or `-dx` with equal probability, `x_next = x + dx*(+-1)`.

**Method.** A discrete symmetric 1D random walk, written as an explicit loop and as a cumulative sum of random steps.

**Test.** `dx = 1`, `dt = 1`, 1000 steps, 1000 independent walks; seed 12.

**Show.** A plot of several random-walk trajectories versus time.

**Verification.** Each walk wanders differently, but as a group they spread symmetrically about the origin.
:::

**Generated prompt (Python):** In Python, simulate an ensemble of one-dimensional random walks and see how they spread. The model is a walker starts at 0 and at each step moves `+dx` or `-dx` with equal probability, `x_next = x + dx*(+-1)`. Use a discrete symmetric 1D random walk, written as an explicit loop and as a cumulative sum of random steps. Test it on `dx = 1`, `dt = 1`, 1000 steps, 1000 independent walks; seed 12. Produce a plot of several random-walk trajectories versus time. As a separate check, confirm that each walk wanders differently, but as a group they spread symmetrically about the origin. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6B.3.1  (06-stochastic/06b-brownian-motion.qmd)

::: {.callout-note title="Recipe 6B.3.1"}
**Objective.** Repeat the random walk with continuous Gaussian steps and compare the spread to discrete steps.

**Model.** The same 1D walk, but each step is drawn from a Gaussian of mean 0 and standard deviation dx, `x_next = x + N(0, dx)`.

**Method.** A 1D random walk with Gaussian-distributed step sizes.

**Test.** `dx = 1`, `dt = 1`, 1000 steps, 1000 walks; seed 12; the spread compared at times 1, 101, ..., 901.

**Show.** A box plot of the walker positions at sampled times, next to the discrete-step walk.

**Verification.** The spread is statistically identical to the discrete-step walk, so the step distribution does not matter for the long-time behavior, only its mean and variance.
:::

**Generated prompt (Python):** In Python, repeat the random walk with continuous Gaussian steps and compare the spread to discrete steps. The model is the same 1D walk, but each step is drawn from a Gaussian of mean 0 and standard deviation dx, `x_next = x + N(0, dx)`. Use a 1D random walk with Gaussian-distributed step sizes. Test it on `dx = 1`, `dt = 1`, 1000 steps, 1000 walks; seed 12; the spread compared at times 1, 101, ..., 901. Produce a box plot of the walker positions at sampled times, next to the discrete-step walk. As a separate check, confirm that the spread is statistically identical to the discrete-step walk, so the step distribution does not matter for the long-time behavior, only its mean and variance. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6B.4.1  (06-stochastic/06b-brownian-motion.qmd)

::: {.callout-note title="Recipe 6B.4.1"}
**Objective.** Simulate Brownian motion in a plane.

**Model.** A walker at (0, 0) takes an independent Gaussian step in x and in y at each time, `x_next = x + N(0, 1)*dx`.

**Method.** Two-dimensional Brownian motion via independent per-axis Gaussian steps.

**Test.** `dx = (1, 1)`, `dt = 1`, 10000 steps; seed 12.

**Show.** A plot of the 2D walk path in the plane.

**Verification.** The path is a tangled 2D walk that drifts from the origin, with typical distance growing like sqrt(t).
:::

**Generated prompt (Python):** In Python, simulate Brownian motion in a plane. The model is a walker at (0, 0) takes an independent Gaussian step in x and in y at each time, `x_next = x + N(0, 1)*dx`. Use two-dimensional Brownian motion via independent per-axis Gaussian steps. Test it on `dx = (1, 1)`, `dt = 1`, 10000 steps; seed 12. Produce a plot of the 2D walk path in the plane. As a separate check, confirm that the path is a tangled 2D walk that drifts from the origin, with typical distance growing like sqrt(t). Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6B.5.1  (06-stochastic/06b-brownian-motion.qmd)

::: {.callout-note title="Recipe 6B.5.1"}
**Objective.** Simulate a Wiener process whose variance stays correct regardless of the step size.

**Model.** A Gaussian random walk `x_next = x + N(0, step_sd)` with the scaling `step_sd = sqrt(2*D*dt)`, so the variance grows as `2*D*t`.

**Method.** A Wiener-process simulation via Gaussian increments scaled by `sqrt(2*D*dt)`.

**Test.** `D = 0.5`, a deliberately large step `dt = 10`, to t = 1000, 1000 walks; seed 12.

**Show.** A box plot of the walker positions over time, showing variance growing as `2*D*t`.

**Verification.** Even at the large step the variance follows `2*D*t`, because the `sqrt(2*D*dt)` scaling keeps it on track.
:::

**Generated prompt (Python):** In Python, simulate a Wiener process whose variance stays correct regardless of the step size. The model is a Gaussian random walk `x_next = x + N(0, step_sd)` with the scaling `step_sd = sqrt(2*D*dt)`, so the variance grows as `2*D*t`. Use a Wiener-process simulation via Gaussian increments scaled by `sqrt(2*D*dt)`. Test it on `D = 0.5`, a deliberately large step `dt = 10`, to t = 1000, 1000 walks; seed 12. Produce a box plot of the walker positions over time, showing variance growing as `2*D*t`. As a separate check, confirm that even at the large step the variance follows `2*D*t`, because the `sqrt(2*D*dt)` scaling keeps it on track. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6C.1.1  (06-stochastic/06c-sde-integrators.qmd)

::: {.callout-note title="Recipe 6C.1.1"}
**Objective.** Integrate a stochastic differential equation and recover free Brownian motion as a test.

**Model.** A general SDE `dX = f(X)*dt + sqrt(2*D(X))*dW`, tested on pure diffusion (`f(X) = 0`, constant D).

**Method.** The Euler-Maruyama method, `X_next = X + f(X)*dt + sqrt(2*D(X))*dW` with `dW ~ N(0, dt)`.

**Test.** `f = 0`, `D = 10`, `X0 = 0`, to t = 200, `dt = 0.01`, 10 trajectories; seed 1.

**Show.** A plot of the ten Brownian trajectories versus time.

**Verification.** The trajectories reproduce the free Brownian motion of the previous chapter.
:::

**Generated prompt (Python):** In Python, integrate a stochastic differential equation and recover free Brownian motion as a test. The model is a general SDE `dX = f(X)*dt + sqrt(2*D(X))*dW`, tested on pure diffusion (`f(X) = 0`, constant D). Use the Euler-Maruyama method, `X_next = X + f(X)*dt + sqrt(2*D(X))*dW` with `dW ~ N(0, dt)`. Test it on `f = 0`, `D = 10`, `X0 = 0`, to t = 200, `dt = 0.01`, 10 trajectories; seed 1. Produce a plot of the ten Brownian trajectories versus time. As a separate check, confirm that the trajectories reproduce the free Brownian motion of the previous chapter. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6C.2.1  (06-stochastic/06c-sde-integrators.qmd)

::: {.callout-note title="Recipe 6C.2.1"}
**Objective.** Simulate the Ornstein-Uhlenbeck process and confirm its stationary variance.

**Model.** A restoring drift with constant noise, `dX = -k*X*dt + sqrt(2*D)*dW`, whose stationary variance is `<x^2> = D/k`.

**Method.** The Euler-Maruyama method applied to the Ornstein-Uhlenbeck process.

**Test.** `k = 1`, noise levels `D = 100, 25, 1`, `X0 = 0`, to t = 1000, `dt = 0.01`; seed 1.

**Show.** Trajectory plots for each D, and a plot of measured variance versus D against the line `<x^2> = D`.

**Verification.** Larger D gives noisier trajectories, and the measured variances fall on the line `<x^2> = D` (since k = 1), confirming `<x^2> = D/k`.
:::

**Generated prompt (Python):** In Python, simulate the Ornstein-Uhlenbeck process and confirm its stationary variance. The model is a restoring drift with constant noise, `dX = -k*X*dt + sqrt(2*D)*dW`, whose stationary variance is `<x^2> = D/k`. Use the Euler-Maruyama method applied to the Ornstein-Uhlenbeck process. Test it on `k = 1`, noise levels `D = 100, 25, 1`, `X0 = 0`, to t = 1000, `dt = 0.01`; seed 1. Produce trajectory plots for each D, and a plot of measured variance versus D against the line `<x^2> = D`. As a separate check, confirm that larger D gives noisier trajectories, and the measured variances fall on the line `<x^2> = D` (since k = 1), confirming `<x^2> = D/k`. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6C.3.1  (06-stochastic/06c-sde-integrators.qmd)

::: {.callout-note title="Recipe 6C.3.1"}
**Objective.** Implement the Milstein integrator, which improves on Euler-Maruyama when the noise depends on the state.

**Model.** An SDE `dX = f(X)*dt + s(X)*dW` with state-dependent noise `s(X) = sqrt(2*D(X))`.

**Method.** The Milstein method, adding the stochastic-Taylor correction `0.5*s(X)*s'(X)*(dW^2 - dt)` to the Euler-Maruyama step (this term vanishes when the noise is constant).

**Test.** Wiener increments `dW ~ N(0, dt)`; the integrator is applied to the gene-circuit model of 6C.4.

**Show.** The integrator applied to a test SDE, returning a trajectory (exercised in 6C.4).

**Verification.** The correction term improves convergence only when the noise depends on the state; for constant noise it reduces to Euler-Maruyama.
:::

**Generated prompt (Python):** In Python, implement the Milstein integrator, which improves on Euler-Maruyama when the noise depends on the state. The model is an SDE `dX = f(X)*dt + s(X)*dW` with state-dependent noise `s(X) = sqrt(2*D(X))`. Use the Milstein method, adding the stochastic-Taylor correction `0.5*s(X)*s'(X)*(dW^2 - dt)` to the Euler-Maruyama step (this term vanishes when the noise is constant). Test it on Wiener increments `dW ~ N(0, dt)`; the integrator is applied to the gene-circuit model of 6C.4. Produce the integrator applied to a test SDE, returning a trajectory (exercised in 6C.4). As a separate check, confirm that the correction term improves convergence only when the noise depends on the state; for constant noise it reduces to Euler-Maruyama. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6C.4.1  (06-stochastic/06c-sde-integrators.qmd)

::: {.callout-note title="Recipe 6C.4.1"}
**Objective.** Simulate expression noise in a bistable self-activating gene whose noise depends on the state.

**Model.** A self-activating gene SDE, Hill self-activation drift with square-root noise, `dX = [g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X]*dt + b*sqrt(X)*dW`.

**Method.** The Milstein method (the noise `b*sqrt(X)` depends on the state).

**Test.** `g0 = 10, g1 = 45, Xth = 200, n = 4, k = 0.15`; noise amplitudes `b = 5, 2, 0.5`; `X0 = 300`, to t = 1000, `dt = 0.01`; seed 1.

**Show.** Trajectory plots for each noise amplitude b, showing transitions between the two states.

**Verification.** The circuit is bistable (stable states near X = 100 and X = 300), and noise drives transitions between them whose frequency rises with b.
:::

**Generated prompt (Python):** In Python, simulate expression noise in a bistable self-activating gene whose noise depends on the state. The model is a self-activating gene SDE, Hill self-activation drift with square-root noise, `dX = [g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X]*dt + b*sqrt(X)*dW`. Use the Milstein method (the noise `b*sqrt(X)` depends on the state). Test it on `g0 = 10, g1 = 45, Xth = 200, n = 4, k = 0.15`; noise amplitudes `b = 5, 2, 0.5`; `X0 = 300`, to t = 1000, `dt = 0.01`; seed 1. Produce trajectory plots for each noise amplitude b, showing transitions between the two states. As a separate check, confirm that the circuit is bistable (stable states near X = 100 and X = 300), and noise drives transitions between them whose frequency rises with b. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6D.1.1  (06-stochastic/06d-stochastic-transitions.qmd)

::: {.callout-note title="Recipe 6D.1.1"}
**Objective.** Simulate a two-gene toggle switch with noise and watch it hop between states.

**Model.** A toggle-switch SDE with mutual repression and constant additive noise, `dX = [gX0 + gX1/(1 + (Y/Y0)^nY) - kX*X]*dt + b*dW_X`, `dY = [gY0 + gY1/(1 + (X/X0)^nX) - kY*Y]*dt + b*dW_Y`, with X and Y kept non-negative.

**Method.** The Euler-Maruyama method for a two-variable SDE.

**Test.** `g0 = 10, g1 = 40, X0 = Y0 = 100, n = 4, k = 0.1` for both genes; noise `b = 20`; initial (X, Y) = (50, 200), to t = 1000, `dt = 0.01`; seed 3.

**Show.** A time-series plot of X and Y and a phase-plane plot, showing hops between the two states.

**Verification.** The switch is bistable (low-X/high-Y and high-X/low-Y states), and noise drives transitions where the trajectory hops between them, swapping which gene is highly expressed.
:::

**Generated prompt (Python):** In Python, simulate a two-gene toggle switch with noise and watch it hop between states. The model is a toggle-switch SDE with mutual repression and constant additive noise, `dX = [gX0 + gX1/(1 + (Y/Y0)^nY) - kX*X]*dt + b*dW_X`, `dY = [gY0 + gY1/(1 + (X/X0)^nX) - kY*Y]*dt + b*dW_Y`, with X and Y kept non-negative. Use the Euler-Maruyama method for a two-variable SDE. Test it on `g0 = 10, g1 = 40, X0 = Y0 = 100, n = 4, k = 0.1` for both genes; noise `b = 20`; initial (X, Y) = (50, 200), to t = 1000, `dt = 0.01`; seed 3. Produce a time-series plot of X and Y and a phase-plane plot, showing hops between the two states. As a separate check, confirm that the switch is bistable (low-X/high-Y and high-X/low-Y states), and noise drives transitions where the trajectory hops between them, swapping which gene is highly expressed. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6D.2.1  (06-stochastic/06d-stochastic-transitions.qmd)

::: {.callout-note title="Recipe 6D.2.1"}
**Objective.** Count how often the noisy toggle switch jumps between its two states.

**Model.** The same toggle-switch SDE and trajectory as 6D.1 (seed 3, `b = 20`, `dt = 0.01`, to t = 1000).

**Method.** Counting transitions by fitting an ellipse to each state's cloud (from the covariance eigen-decomposition), assigning a point to a state when it lies inside that ellipse, and detecting a transition from the change in assigned state.

**Test.** The 6D.1 trajectory, with each state's ellipse axes enlarged by a factor of 1.5, splitting the plane along Y = X.

**Show.** The phase plane with the two fitted state ellipses, and the count and step of each detected transition.

**Verification.** The crude counter over-counts: a single boundary-jitter event registers as a cluster of transitions within a few steps, motivating a more robust statistic.
:::

**Generated prompt (Python):** In Python, count how often the noisy toggle switch jumps between its two states. The model is the same toggle-switch SDE and trajectory as 6D.1 (seed 3, `b = 20`, `dt = 0.01`, to t = 1000). Use counting transitions by fitting an ellipse to each state's cloud (from the covariance eigen-decomposition), assigning a point to a state when it lies inside that ellipse, and detecting a transition from the change in assigned state. Test it on The 6D.1 trajectory, with each state's ellipse axes enlarged by a factor of 1.5, splitting the plane along Y = X. Produce the phase plane with the two fitted state ellipses, and the count and step of each detected transition. As a separate check, confirm that the crude counter over-counts: a single boundary-jitter event registers as a cluster of transitions within a few steps, motivating a more robust statistic. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 6D.3.1  (06-stochastic/06d-stochastic-transitions.qmd)

::: {.callout-note title="Recipe 6D.3.1"}
**Objective.** Measure the transition rate of the noisy toggle switch from its mean first-passage time.

**Model.** The same two-variable toggle-switch SDE (`g0 = 10, g1 = 40`, threshold 100, `k = 0.1`, `n = 4` per gene, noise `b = 20`), implemented in compiled Fortran.

**Method.** The mean first-passage time: simulate until the trajectory crosses the separatrix, record the passage time, relax, and repeat, then take the rate `kappa = 1/(2*tau)` from the mean time tau.

**Test.** The compiled Fortran routine called from both R (`.Fortran`) and Python (`ctypes`); seed 11, total time `1e5`, `dt = 0.01`, relaxation time 100.

**Show.** The transition rates from R and from Python (identical for the same seed).

**Verification.** Both languages call the identical compiled routine with the same seed and return identical rates, showing the numerics live in the shared library; a longer total time gives more accurate rates.
:::

**Generated prompt (Python):** In Python, measure the transition rate of the noisy toggle switch from its mean first-passage time. The model is the same two-variable toggle-switch SDE (`g0 = 10, g1 = 40`, threshold 100, `k = 0.1`, `n = 4` per gene, noise `b = 20`), implemented in compiled Fortran. Use the mean first-passage time: simulate until the trajectory crosses the separatrix, record the passage time, relax, and repeat, then take the rate `kappa = 1/(2*tau)` from the mean time tau. Test it on The compiled Fortran routine called from both R (`.Fortran`) and Python (`ctypes`); seed 11, total time `1e5`, `dt = 0.01`, relaxation time 100. Produce the transition rates from R and from Python (identical for the same seed). As a separate check, confirm that both languages call the identical compiled routine with the same seed and return identical rates, showing the numerics live in the shared library; a longer total time gives more accurate rates. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7A.1.1  (07-pde/07a-modeling-diffusion.qmd)

::: {.callout-note title="Recipe 7A.1.1"}
**Objective.** Build the explicit finite-difference scheme for the 1D diffusion equation.

**Model.** The 1D diffusion equation `dP/dt = D*d2P/dX2` for a distribution `P(X, t)`, with no reaction term.

**Method.** The explicit forward-time centered-space finite-difference method, `P_i_next = P_i + D*(dt/dX^2)*(P_{i+1} + P_{i-1} - 2*P_i)`, with Dirichlet boundaries (P = 0 at both ends).

**Test.** The generic solver `pde_fd_diffusion(ngrid, X_all, dX, dt, D, ...)`; the stability factor `D*dt/dX^2` must stay small.

**Show.** The distribution P advanced several finite-difference steps (exercised in 7A.2).

**Verification.** The scheme is the diffusion analog of Euler integration for ODEs, advancing P one time step from its neighbors.
:::

**Generated prompt (Python):** In Python, build the explicit finite-difference scheme for the 1D diffusion equation. The model is the 1D diffusion equation `dP/dt = D*d2P/dX2` for a distribution `P(X, t)`, with no reaction term. Use the explicit forward-time centered-space finite-difference method, `P_i_next = P_i + D*(dt/dX^2)*(P_{i+1} + P_{i-1} - 2*P_i)`, with Dirichlet boundaries (P = 0 at both ends). Test it on The generic solver `pde_fd_diffusion(ngrid, X_all, dX, dt, D, ...)`; the stability factor `D*dt/dX^2` must stay small. Produce the distribution P advanced several finite-difference steps (exercised in 7A.2). As a separate check, confirm that the scheme is the diffusion analog of Euler integration for ODEs, advancing P one time step from its neighbors. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7A.2.1  (07-pde/07a-modeling-diffusion.qmd)

::: {.callout-note title="Recipe 7A.2.1"}
**Objective.** Diffuse an initially concentrated distribution and compare its spread to theory.

**Model.** The 1D diffusion equation `dP/dt = D*d2P/dX2`.

**Method.** The explicit finite-difference diffusion integrator (from 7A.1), tracking the mean, variance, and total probability.

**Test.** Domain length 100, `dX = 1`, `dt = 0.01`, `D = 1`; all probability initially at X = 0; compared with the theoretical variance `sigma^2 = 2*D*t`.

**Show.** A plot of P(X) spreading over time, and its variance versus t against `2*D*t`.

**Verification.** The point distribution spreads into a widening Gaussian whose variance grows as `2*D*t` at early times, while the total probability slowly leaks as the Dirichlet ends absorb it.
:::

**Generated prompt (Python):** In Python, diffuse an initially concentrated distribution and compare its spread to theory. The model is the 1D diffusion equation `dP/dt = D*d2P/dX2`. Use the explicit finite-difference diffusion integrator (from 7A.1), tracking the mean, variance, and total probability. Test it on Domain length 100, `dX = 1`, `dt = 0.01`, `D = 1`; all probability initially at X = 0; compared with the theoretical variance `sigma^2 = 2*D*t`. Produce a plot of P(X) spreading over time, and its variance versus t against `2*D*t`. As a separate check, confirm that the point distribution spreads into a widening Gaussian whose variance grows as `2*D*t` at early times, while the total probability slowly leaks as the Dirichlet ends absorb it. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7B.1.1  (07-pde/07b-reaction-diffusion.qmd)

::: {.callout-note title="Recipe 7B.1.1"}
**Objective.** Simulate Fisher's equation and see traveling waves of advance.

**Model.** Fisher's equation, logistic growth with diffusion, `du/dt = r*u*(1 - u/B) + D*d2u/dX2`.

**Method.** The explicit finite-difference reaction-diffusion integrator (a reaction term added to the diffusion step each iteration), with Dirichlet boundaries.

**Test.** Domain length 100, `dX = 1`, `dt = 0.01`, `D = 1`, r = 0.64, `B = 75`; a central patch, and separately a left-end patch, each `u = 50`.

**Show.** Plots of u(X) at successive times, showing traveling fronts from a central and a left-end patch.

**Verification.** A local patch grows to carrying capacity and emits two constant-speed fronts, while a left-end patch produces a single rightward traveling wave.
:::

**Generated prompt (Python):** In Python, simulate Fisher's equation and see traveling waves of advance. The model is fisher's equation, logistic growth with diffusion, `du/dt = r*u*(1 - u/B) + D*d2u/dX2`. Use the explicit finite-difference reaction-diffusion integrator (a reaction term added to the diffusion step each iteration), with Dirichlet boundaries. Test it on Domain length 100, `dX = 1`, `dt = 0.01`, `D = 1`, r = 0.64, `B = 75`; a central patch, and separately a left-end patch, each `u = 50`. Produce plots of u(X) at successive times, showing traveling fronts from a central and a left-end patch. As a separate check, confirm that a local patch grows to carrying capacity and emits two constant-speed fronts, while a left-end patch produces a single rightward traveling wave. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7B.2.1  (07-pde/07b-reaction-diffusion.qmd)

::: {.callout-note title="Recipe 7B.2.1"}
**Objective.** Solve the Fokker-Planck equation for the distribution of an SDE, rather than averaging many trajectories.

**Model.** The Fokker-Planck equation `dP/dt = -d(f*P)/dX + D*d2P/dX2`, tested on the Ornstein-Uhlenbeck well with drift `f = -k*X`, whose steady state is `P_ss(X) = sqrt(k/(2*pi*D))*exp(-k*X^2/(2*D))`.

**Method.** An explicit finite-difference Fokker-Planck integrator: a centered-difference drift term plus the finite-difference diffusion term, with Dirichlet boundaries.

**Test.** Domain length 100, `dX = 1`, `dt = 0.01`, `D = 1`; springs `k = 0.01` and `k = 0.03`; started from both a localized patch and a uniform distribution.

**Show.** Plots of P(X) relaxing to the steady-state Gaussian, for each spring stiffness k.

**Verification.** The distribution relaxes to the analytic steady-state Gaussian from any initial condition, and a stiffer spring confines it to a narrower peak.
:::

**Generated prompt (Python):** In Python, solve the Fokker-Planck equation for the distribution of an SDE, rather than averaging many trajectories. The model is the Fokker-Planck equation `dP/dt = -d(f*P)/dX + D*d2P/dX2`, tested on the Ornstein-Uhlenbeck well with drift `f = -k*X`, whose steady state is `P_ss(X) = sqrt(k/(2*pi*D))*exp(-k*X^2/(2*D))`. Use an explicit finite-difference Fokker-Planck integrator: a centered-difference drift term plus the finite-difference diffusion term, with Dirichlet boundaries. Test it on Domain length 100, `dX = 1`, `dt = 0.01`, `D = 1`; springs `k = 0.01` and `k = 0.03`; started from both a localized patch and a uniform distribution. Produce plots of P(X) relaxing to the steady-state Gaussian, for each spring stiffness k. As a separate check, confirm that the distribution relaxes to the analytic steady-state Gaussian from any initial condition, and a stiffer spring confines it to a narrower peak. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7C.1.1  (07-pde/07c-turing-instability.qmd)

::: {.callout-note title="Recipe 7C.1.1"}
**Objective.** Build a reaction-diffusion integrator for any number of coupled components.

**Model.** A generic n-component 1D reaction-diffusion system, `du_k/dt = f_k(u) + D_k*d2u_k/dX2`, with a reaction term per component and its own diffusion constant.

**Method.** A vectorized explicit forward-time centered-space integrator for multiple components, with periodic (wrap-around) boundaries.

**Test.** The generic solver `pde_fd_reaction_diffusion_multi(derivs, n, ngrid, dX, dt, D, ...)`, with D a length-n vector and the initial state an (n, ngrid) matrix.

**Show.** All components advanced together on the wrap-around grid (exercised in 7C.2).

**Verification.** The solver advances all components together on a wrap-around grid; from random initial perturbations, R and Python give patterns of the same wavelength at different positions.
:::

**Generated prompt (Python):** In Python, build a reaction-diffusion integrator for any number of coupled components. The model is a generic n-component 1D reaction-diffusion system, `du_k/dt = f_k(u) + D_k*d2u_k/dX2`, with a reaction term per component and its own diffusion constant. Use a vectorized explicit forward-time centered-space integrator for multiple components, with periodic (wrap-around) boundaries. Test it on The generic solver `pde_fd_reaction_diffusion_multi(derivs, n, ngrid, dX, dt, D, ...)`, with D a length-n vector and the initial state an (n, ngrid) matrix. Produce all components advanced together on the wrap-around grid (exercised in 7C.2). As a separate check, confirm that the solver advances all components together on a wrap-around grid; from random initial perturbations, R and Python give patterns of the same wavelength at different positions. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7C.2.1  (07-pde/07c-turing-instability.qmd)

::: {.callout-note title="Recipe 7C.2.1"}
**Objective.** Reproduce Turing pattern formation and its two failure modes in a substrate-depletion model.

**Model.** The Gierer-Meinhardt substrate-depletion model, `f(u, v) = u^2*v - u`, `g(u, v) = mu*(1 - u^2*v)`, with the substrate v diffusing fast (`Dv = 1`) and u slow (`Du = d`).

**Method.** The multi-component finite-difference reaction-diffusion integrator (from 7C.1), run in successive time blocks.

**Test.** Domain length 20, `dX = 0.2`, `dt = 0.01`; nearly uniform initial `u = v = 1` with +-0.1 noise (seed 10); a pattern-forming case (`d = 0.1`, `mu = 1.5`), a case with close diffusion constants (`d = 0.8`), and an oscillating case (`d = 0.3`, `mu = 0.9`, initial `u = 2`).

**Show.** Plots of u(X) at successive times for each case: a stationary pattern, a flat state, and a uniform oscillation.

**Verification.** The same model gives three outcomes set by the diffusion ratio and mu: a stationary periodic Turing pattern, a homogeneous steady state, and a spatially uniform temporal oscillation.
:::

**Generated prompt (Python):** In Python, reproduce Turing pattern formation and its two failure modes in a substrate-depletion model. The model is the Gierer-Meinhardt substrate-depletion model, `f(u, v) = u^2*v - u`, `g(u, v) = mu*(1 - u^2*v)`, with the substrate v diffusing fast (`Dv = 1`) and u slow (`Du = d`). Use the multi-component finite-difference reaction-diffusion integrator (from 7C.1), run in successive time blocks. Test it on Domain length 20, `dX = 0.2`, `dt = 0.01`; nearly uniform initial `u = v = 1` with +-0.1 noise (seed 10); a pattern-forming case (`d = 0.1`, `mu = 1.5`), a case with close diffusion constants (`d = 0.8`), and an oscillating case (`d = 0.3`, `mu = 0.9`, initial `u = 2`). Produce plots of u(X) at successive times for each case: a stationary pattern, a flat state, and a uniform oscillation. As a separate check, confirm that the same model gives three outcomes set by the diffusion ratio and mu: a stationary periodic Turing pattern, a homogeneous steady state, and a spatially uniform temporal oscillation. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7C.3.1  (07-pde/07c-turing-instability.qmd)

::: {.callout-note title="Recipe 7C.3.1"}
**Objective.** Form a Turing pattern from an activator-inhibitor model and contrast it with substrate depletion.

**Model.** The Gierer-Meinhardt activator-inhibitor model, `f(u, v) = u^2/v - u`, `g(u, v) = mu*(u^2 - v)`, with `Du = d` slow and `Dv = 1` fast (local self-enhancement, long-range inhibition).

**Method.** The multi-component finite-difference reaction-diffusion integrator (from 7C.1), run in successive time blocks.

**Test.** Domain length 20, `dX = 0.2`, `dt = 0.01`, `d = 0.1`, `mu = 1.5`; nearly uniform initial `u = v = 1` with +-0.1 noise (seed 10).

**Show.** A plot of u(X) and v(X) forming a stationary periodic pattern, peaking together.

**Verification.** It forms a stationary periodic Turing pattern, but with u and v peaking together rather than out of phase as in the substrate-depletion case.
:::

**Generated prompt (Python):** In Python, form a Turing pattern from an activator-inhibitor model and contrast it with substrate depletion. The model is the Gierer-Meinhardt activator-inhibitor model, `f(u, v) = u^2/v - u`, `g(u, v) = mu*(u^2 - v)`, with `Du = d` slow and `Dv = 1` fast (local self-enhancement, long-range inhibition). Use the multi-component finite-difference reaction-diffusion integrator (from 7C.1), run in successive time blocks. Test it on Domain length 20, `dX = 0.2`, `dt = 0.01`, `d = 0.1`, `mu = 1.5`; nearly uniform initial `u = v = 1` with +-0.1 noise (seed 10). Produce a plot of u(X) and v(X) forming a stationary periodic pattern, peaking together. As a separate check, confirm that it forms a stationary periodic Turing pattern, but with u and v peaking together rather than out of phase as in the substrate-depletion case. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7D.2.1  (07-pde/07d-pattern-formation-dictyostelium.qmd)

::: {.callout-note title="Recipe 7D.2.1"}
**Objective.** Simulate spiral cAMP waves in a Dictyostelium colony of excitable cells.

**Model.** The Kessler-Levine model: a cAMP field `dc/dt = a^2*(d2c/dX2 + d2c/dY2) - k*c + s`, coupled to discrete cells that cycle inactive -> excited (firing and secreting when local `c > c_T`) -> refractory -> inactive.

**Method.** A 2D finite-difference integrator with no-flux boundaries for the cAMP field, coupled each step to a vectorized cell state-machine update.

**Test.** A 101x101 grid, cell fraction 0.15, threshold `c_T = 1`, secretion `dc = 300` over `t_e = 2`, recovery `t_r = 20`, degradation `k = 0.5`, `dt = 0.01`, to t = 150.

**Show.** A snapshot of the 2D cAMP field with the rotating spiral wave and the cell states.

**Verification.** The coupled field and excitable cells self-organize into a rotating cAMP spiral wave, with excited cells riding the high-cAMP crest.
:::

**Generated prompt (Python):** In Python, simulate spiral cAMP waves in a Dictyostelium colony of excitable cells. The model is the Kessler-Levine model: a cAMP field `dc/dt = a^2*(d2c/dX2 + d2c/dY2) - k*c + s`, coupled to discrete cells that cycle inactive -> excited (firing and secreting when local `c > c_T`) -> refractory -> inactive. Use a 2D finite-difference integrator with no-flux boundaries for the cAMP field, coupled each step to a vectorized cell state-machine update. Test it on A 101x101 grid, cell fraction 0.15, threshold `c_T = 1`, secretion `dc = 300` over `t_e = 2`, recovery `t_r = 20`, degradation `k = 0.5`, `dt = 0.01`, to t = 150. Produce a snapshot of the 2D cAMP field with the rotating spiral wave and the cell states. As a separate check, confirm that the coupled field and excitable cells self-organize into a rotating cAMP spiral wave, with excited cells riding the high-cAMP crest. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7E.1.1  (07-pde/07e-2d-reaction-diffusion.qmd)

::: {.callout-note title="Recipe 7E.1.1"}
**Objective.** Build a two-component reaction-diffusion integrator on a 2D grid.

**Model.** A generic two-component 2D reaction-diffusion system, `du/dt = f(u, v) + Du*(d2u/dX2 + d2u/dY2)`, `dv/dt = g(u, v) + Dv*(d2v/dX2 + d2v/dY2)`.

**Method.** An explicit forward-time centered-space 2D integrator with periodic boundaries, the five-point Laplacian computed by whole-grid array shifts.

**Test.** The generic solver with a run-in-blocks driver; nearly uniform initial `u = v = 1` with +-0.1 noise (seed 10); 2D stability needs `D*dt/dX^2 < 1/4`, so `dt = 0.005`.

**Show.** u and v advanced on the 2D grid (exercised in 7E.2).

**Verification.** The vectorized 2D integrator evaluates the reaction on the whole grid at once, running faster than a point-by-point 1D approach.
:::

**Generated prompt (Python):** In Python, build a two-component reaction-diffusion integrator on a 2D grid. The model is a generic two-component 2D reaction-diffusion system, `du/dt = f(u, v) + Du*(d2u/dX2 + d2u/dY2)`, `dv/dt = g(u, v) + Dv*(d2v/dX2 + d2v/dY2)`. Use an explicit forward-time centered-space 2D integrator with periodic boundaries, the five-point Laplacian computed by whole-grid array shifts. Test it on The generic solver with a run-in-blocks driver; nearly uniform initial `u = v = 1` with +-0.1 noise (seed 10); 2D stability needs `D*dt/dX^2 < 1/4`, so `dt = 0.005`. Produce u and v advanced on the 2D grid (exercised in 7E.2). As a separate check, confirm that the vectorized 2D integrator evaluates the reaction on the whole grid at once, running faster than a point-by-point 1D approach. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7E.2.1  (07-pde/07e-2d-reaction-diffusion.qmd)

::: {.callout-note title="Recipe 7E.2.1"}
**Objective.** Grow a 2D Turing pattern from a substrate-depletion model.

**Model.** The 2D Gierer-Meinhardt substrate-depletion model, `f(u, v) = u^2*v - u`, `g(u, v) = mu*(1 - u^2*v)`, with `Du = d`, `Dv = 1`.

**Method.** The 2D finite-difference reaction-diffusion integrator (from 7E.1), run in successive blocks.

**Test.** A 101x101 grid, `dX = 0.2`, `dt = 0.005`, `d = 0.1`, `mu = 1.5`; nearly uniform initial `u = v = 1` with +-0.1 noise (seed 10).

**Show.** A 2D image of the u field forming a labyrinth of stripes.

**Verification.** The instability grows from near-uniform noise into a coarsening labyrinth of winding stripes.
:::

**Generated prompt (Python):** In Python, grow a 2D Turing pattern from a substrate-depletion model. The model is the 2D Gierer-Meinhardt substrate-depletion model, `f(u, v) = u^2*v - u`, `g(u, v) = mu*(1 - u^2*v)`, with `Du = d`, `Dv = 1`. Use the 2D finite-difference reaction-diffusion integrator (from 7E.1), run in successive blocks. Test it on A 101x101 grid, `dX = 0.2`, `dt = 0.005`, `d = 0.1`, `mu = 1.5`; nearly uniform initial `u = v = 1` with +-0.1 noise (seed 10). Produce a 2D image of the u field forming a labyrinth of stripes. As a separate check, confirm that the instability grows from near-uniform noise into a coarsening labyrinth of winding stripes. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 7E.3.1  (07-pde/07e-2d-reaction-diffusion.qmd)

::: {.callout-note title="Recipe 7E.3.1"}
**Objective.** Grow a 2D Turing pattern from an activator-inhibitor model and compare its morphology to stripes.

**Model.** The 2D Gierer-Meinhardt activator-inhibitor model, `f(u, v) = u^2/v - u`, `g(u, v) = mu*(u^2 - v)`, with `Du = d`, `Dv = 1`.

**Method.** The 2D finite-difference reaction-diffusion integrator (from 7E.1), run in successive blocks.

**Test.** A 101x101 grid, `dX = 0.2`, `dt = 0.005`, `d = 0.1`, `mu = 1.5`; nearly uniform initial `u = v = 1` with +-0.1 noise (seed 10), identical to 7E.2.

**Show.** A 2D image of the u field forming an array of spots.

**Verification.** From the same start and parameters as the substrate-depletion case, the activator-inhibitor kinetics instead settle into a regular array of spots, so the reaction kinetics select stripes versus spots.
:::

**Generated prompt (Python):** In Python, grow a 2D Turing pattern from an activator-inhibitor model and compare its morphology to stripes. The model is the 2D Gierer-Meinhardt activator-inhibitor model, `f(u, v) = u^2/v - u`, `g(u, v) = mu*(u^2 - v)`, with `Du = d`, `Dv = 1`. Use the 2D finite-difference reaction-diffusion integrator (from 7E.1), run in successive blocks. Test it on A 101x101 grid, `dX = 0.2`, `dt = 0.005`, `d = 0.1`, `mu = 1.5`; nearly uniform initial `u = v = 1` with +-0.1 noise (seed 10), identical to 7E.2. Produce a 2D image of the u field forming an array of spots. As a separate check, confirm that from the same start and parameters as the substrate-depletion case, the activator-inhibitor kinetics instead settle into a regular array of spots, so the reaction kinetics select stripes versus spots. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8A.1.1  (08-monte-carlo/08a-monte-carlo-methods.qmd)

::: {.callout-note title="Recipe 8A.1.1"}
**Objective.** Estimate pi by Monte Carlo sampling of a circle's area.

**Model.** The fraction of uniform points in the square [-1, 1]^2 that fall inside the unit circle is `pi/4`, so `pi = 4*P`.

**Method.** Monte Carlo estimation of pi by area sampling (the fraction of uniform points inside the unit circle, times four).

**Test.** Uniform points in [-1, 1]^2; a single estimate at `n = 1e4` and a convergence run to `n = 1e5`; seed 1.

**Show.** The running estimate of pi versus sample count, and a scatter of the sampled points in and out of the circle.

**Verification.** The estimate is noticeably off at n = 1e4 but settles near pi by n = 1e5, showing the slow 1/sqrt(n) Monte Carlo convergence.
:::

**Generated prompt (Python):** In Python, estimate pi by Monte Carlo sampling of a circle's area. The model is the fraction of uniform points in the square [-1, 1]^2 that fall inside the unit circle is `pi/4`, so `pi = 4*P`. Use monte Carlo estimation of pi by area sampling (the fraction of uniform points inside the unit circle, times four). Test it on Uniform points in [-1, 1]^2; a single estimate at `n = 1e4` and a convergence run to `n = 1e5`; seed 1. Produce the running estimate of pi versus sample count, and a scatter of the sampled points in and out of the circle. As a separate check, confirm that the estimate is noticeably off at n = 1e4 but settles near pi by n = 1e5, showing the slow 1/sqrt(n) Monte Carlo convergence. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8A.2.1  (08-monte-carlo/08a-monte-carlo-methods.qmd)

::: {.callout-note title="Recipe 8A.2.1"}
**Objective.** Estimate pi with Buffon's needle experiment.

**Model.** Dropping a needle of length l on a floor ruled with lines a distance d apart (l < d) crosses a line with probability `2*l/(pi*d)`, so `pi = 2*l/(P*d)`.

**Method.** Buffon's needle Monte Carlo estimation of pi, using rejection sampling for a uniform random needle direction.

**Test.** Line spacing `d = 1`, needle lengths `l = 0.2, 0.4, ..., 1.4`; `n = 1e4` then `n = 1e5` drops; seed 1.

**Show.** The estimate of pi for each needle length, versus sample count.

**Verification.** At n = 1e4 the estimate wanders around pi (and breaks for l > d), and n = 1e5 tightens it.
:::

**Generated prompt (Python):** In Python, estimate pi with Buffon's needle experiment. The model is dropping a needle of length l on a floor ruled with lines a distance d apart (l < d) crosses a line with probability `2*l/(pi*d)`, so `pi = 2*l/(P*d)`. Use buffon's needle Monte Carlo estimation of pi, using rejection sampling for a uniform random needle direction. Test it on Line spacing `d = 1`, needle lengths `l = 0.2, 0.4, ..., 1.4`; `n = 1e4` then `n = 1e5` drops; seed 1. Produce the estimate of pi for each needle length, versus sample count. As a separate check, confirm that at n = 1e4 the estimate wanders around pi (and breaks for l > d), and n = 1e5 tightens it. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8A.3.1  (08-monte-carlo/08a-monte-carlo-methods.qmd)

::: {.callout-note title="Recipe 8A.3.1"}
**Objective.** Integrate a function by Monte Carlo sampling and compare it to a deterministic rule.

**Model.** The integral of `f(x) = exp(-x^2)`; the reference is the midpoint rule, and the exact value over [0, inf) is `sqrt(pi)/2`.

**Method.** Monte Carlo integration by uniform sampling, `I = (x2 - x1)/n * sum f(x_i)`, alongside the deterministic midpoint rule.

**Test.** Interval [0, 3]; midpoint rule at `n = 1000`, Monte Carlo at `n = 1e4`; seed 1.

**Show.** The midpoint-rule and Monte Carlo estimates against `sqrt(pi)/2`.

**Verification.** Both land near `sqrt(pi)/2 = 0.886`; the midpoint rule is more accurate for this smooth 1D integrand, but Monte Carlo's error is dimension-independent.
:::

**Generated prompt (Python):** In Python, integrate a function by Monte Carlo sampling and compare it to a deterministic rule. The model is the integral of `f(x) = exp(-x^2)`; the reference is the midpoint rule, and the exact value over [0, inf) is `sqrt(pi)/2`. Use monte Carlo integration by uniform sampling, `I = (x2 - x1)/n * sum f(x_i)`, alongside the deterministic midpoint rule. Test it on Interval [0, 3]; midpoint rule at `n = 1000`, Monte Carlo at `n = 1e4`; seed 1. Produce the midpoint-rule and Monte Carlo estimates against `sqrt(pi)/2`. As a separate check, confirm that both land near `sqrt(pi)/2 = 0.886`; the midpoint rule is more accurate for this smooth 1D integrand, but Monte Carlo's error is dimension-independent. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8A.4.1  (08-monte-carlo/08a-monte-carlo-methods.qmd)

::: {.callout-note title="Recipe 8A.4.1"}
**Objective.** Reduce the variance of a Monte Carlo integral with importance sampling.

**Model.** The same integral of `f(x) = exp(-x^2)`, rewritten as `I = integral of f(x)/p(x) * p(x) dx` and estimated by the mean of `f(x_i)/p(x_i)` with `x_i` drawn from p.

**Method.** Monte Carlo importance sampling with an exponential density `p(x) = exp(-x)` (drawn by inverse transform), compared against uniform Monte Carlo.

**Test.** Range [0, 100], `n = 1e4`, 100 replicates each of importance and uniform sampling; seed 1.

**Show.** Box plots of the importance-sampling versus uniform estimates over many replicates.

**Verification.** Over the wide range, uniform sampling scatters badly by wasting points where f is negligible, while importance sampling concentrates points where f is large and converges with far less variance.
:::

**Generated prompt (Python):** In Python, reduce the variance of a Monte Carlo integral with importance sampling. The model is the same integral of `f(x) = exp(-x^2)`, rewritten as `I = integral of f(x)/p(x) * p(x) dx` and estimated by the mean of `f(x_i)/p(x_i)` with `x_i` drawn from p. Use monte Carlo importance sampling with an exponential density `p(x) = exp(-x)` (drawn by inverse transform), compared against uniform Monte Carlo. Test it on Range [0, 100], `n = 1e4`, 100 replicates each of importance and uniform sampling; seed 1. Produce box plots of the importance-sampling versus uniform estimates over many replicates. As a separate check, confirm that over the wide range, uniform sampling scatters badly by wasting points where f is negligible, while importance sampling concentrates points where f is large and converges with far less variance. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8B.1.1  (08-monte-carlo/08b-metropolis-algorithm.qmd)

::: {.callout-note title="Recipe 8B.1.1"}
**Objective.** Simulate a two-state promoter as a Markov chain and find its long-run state fractions.

**Model.** A two-state promoter Markov chain with transition matrix `T = [[0.7, 0.3], [0.1, 0.9]]` (rows from active, from inactive), so the active state leaves with probability 0.3 and the inactive with 0.1.

**Method.** Direct simulation of a two-state Markov chain (draw a uniform, pick the next state from the current row of T).

**Test.** The transition matrix above, starting active, `1e4` steps; seed 1.

**Show.** A plot of the promoter state over the first steps, and the running active/inactive fractions.

**Verification.** The trajectory bursts (short active spells, long silent stretches), and the running fractions converge to about 25% active and 75% inactive.
:::

**Generated prompt (Python):** In Python, simulate a two-state promoter as a Markov chain and find its long-run state fractions. The model is a two-state promoter Markov chain with transition matrix `T = [[0.7, 0.3], [0.1, 0.9]]` (rows from active, from inactive), so the active state leaves with probability 0.3 and the inactive with 0.1. Use direct simulation of a two-state Markov chain (draw a uniform, pick the next state from the current row of T). Test it on The transition matrix above, starting active, `1e4` steps; seed 1. Produce a plot of the promoter state over the first steps, and the running active/inactive fractions. As a separate check, confirm that the trajectory bursts (short active spells, long silent stretches), and the running fractions converge to about 25% active and 75% inactive. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8B.3.1  (08-monte-carlo/08b-metropolis-algorithm.qmd)

::: {.callout-note title="Recipe 8B.3.1"}
**Objective.** Sample a target distribution with Metropolis-Hastings and confirm it reproduces the target.

**Model.** The two-state bursting-gene promoter of 8B.1, now specified only by its stationary distribution `p_ss = (0.25, 0.75)` as the target `P(x)`.

**Method.** The Metropolis-Hastings algorithm with a symmetric two-state proposal (always propose the other state), so acceptance is `a = min(1, P(x')/P(x))`.

**Test.** Target `p_ss = (0.25, 0.75)`, starting active, `1e4` steps; seed 1.

**Show.** The sampled state fractions and the recovered transition matrix.

**Verification.** The sampler has a different transition matrix from the original chain yet reproduces the same stationary distribution, showing that many chains can share one stationary distribution.
:::

**Generated prompt (Python):** In Python, sample a target distribution with Metropolis-Hastings and confirm it reproduces the target. The model is the two-state bursting-gene promoter of 8B.1, now specified only by its stationary distribution `p_ss = (0.25, 0.75)` as the target `P(x)`. Use the Metropolis-Hastings algorithm with a symmetric two-state proposal (always propose the other state), so acceptance is `a = min(1, P(x')/P(x))`. Test it on Target `p_ss = (0.25, 0.75)`, starting active, `1e4` steps; seed 1. Produce the sampled state fractions and the recovered transition matrix. As a separate check, confirm that the sampler has a different transition matrix from the original chain yet reproduces the same stationary distribution, showing that many chains can share one stationary distribution. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8B.4.1  (08-monte-carlo/08b-metropolis-algorithm.qmd)

::: {.callout-note title="Recipe 8B.4.1"}
**Objective.** Sample a continuous Gaussian with Metropolis-Hastings and see how the step size affects the result.

**Model.** A continuous target `P(x) proportional to exp(-x^2)`, a Gaussian of variance 1/2.

**Method.** The Metropolis-Hastings algorithm with a uniform random-walk proposal `x' = x + Uniform(-dx_max, dx_max)` and acceptance `a = min(1, exp(-x'^2 + x^2))`.

**Test.** Start at x = 0, `1e4` steps, proposal widths `dx_max = 0.1, 0.5, 2, 10`; seed 1.

**Show.** Histograms of the samples at each proposal width against the true Gaussian.

**Verification.** The sampled histogram matches the true Gaussian best at an intermediate dx_max (about 50% acceptance); too-small steps barely explore and too-large steps mostly reject.
:::

**Generated prompt (Python):** In Python, sample a continuous Gaussian with Metropolis-Hastings and see how the step size affects the result. The model is a continuous target `P(x) proportional to exp(-x^2)`, a Gaussian of variance 1/2. Use the Metropolis-Hastings algorithm with a uniform random-walk proposal `x' = x + Uniform(-dx_max, dx_max)` and acceptance `a = min(1, exp(-x'^2 + x^2))`. Test it on Start at x = 0, `1e4` steps, proposal widths `dx_max = 0.1, 0.5, 2, 10`; seed 1. Produce histograms of the samples at each proposal width against the true Gaussian. As a separate check, confirm that the sampled histogram matches the true Gaussian best at an intermediate dx_max (about 50% acceptance); too-small steps barely explore and too-large steps mostly reject. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8B.5.1  (08-monte-carlo/08b-metropolis-algorithm.qmd)

::: {.callout-note title="Recipe 8B.5.1"}
**Objective.** Sample the 1D Ising model with Metropolis Monte Carlo.

**Model.** A 1D Ising chain of n spins `s_i = +-1` with open ends, energy `E = -J*sum s_i*s_{i+1}`, sampled from the Boltzmann distribution `p proportional to exp(-E/T)`.

**Method.** The Metropolis-Hastings algorithm with single-spin-flip proposals and acceptance `a = min(1, exp(-(E' - E)/T))`.

**Test.** `n = 10` spins, `J = 1`, temperature `T = 1`, 1000 steps, random initial spins; seed 1.

**Show.** A plot of the Ising energy over the Monte Carlo steps.

**Verification.** At T = 1 the sampler moves between the minimum-energy aligned configuration (`E = -9`) and higher-energy ones, with about 20% acceptance.
:::

**Generated prompt (Python):** In Python, sample the 1D Ising model with Metropolis Monte Carlo. The model is a 1D Ising chain of n spins `s_i = +-1` with open ends, energy `E = -J*sum s_i*s_{i+1}`, sampled from the Boltzmann distribution `p proportional to exp(-E/T)`. Use the Metropolis-Hastings algorithm with single-spin-flip proposals and acceptance `a = min(1, exp(-(E' - E)/T))`. Test it on `n = 10` spins, `J = 1`, temperature `T = 1`, 1000 steps, random initial spins; seed 1. Produce a plot of the Ising energy over the Monte Carlo steps. As a separate check, confirm that at T = 1 the sampler moves between the minimum-energy aligned configuration (`E = -9`) and higher-energy ones, with about 20% acceptance. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8C.3.1  (08-monte-carlo/08c-particles-in-a-box.qmd)

::: {.callout-note title="Recipe 8C.3.1"}
**Objective.** Sample equilibrium configurations of a Lennard-Jones particle box with Metropolis Monte Carlo.

**Model.** A 2D box of N particles interacting through the Lennard-Jones potential `U(r) = 1/(12*r^12) - 1/(6*r^6)` with minimum-image periodic boundaries, sampled from the Boltzmann distribution `proportional to exp(-U/T)` over positions.

**Method.** The Metropolis-Hastings algorithm with local single-particle moves and an incremental energy update.

**Test.** 25 particles on a 5x5 grid at spacing 1.25 (box side `a = 6.25`), move size `dxmax = 0.25`, temperature `T = 0.05`, 5000 equilibration then 10000 sampling steps; seed 10.

**Show.** A plot of the box energy over the equilibration and sampling phases.

**Verification.** During equilibration the energy falls as particles condense from the open grid into a denser cluster, then fluctuates around a low plateau.
:::

**Generated prompt (Python):** In Python, sample equilibrium configurations of a Lennard-Jones particle box with Metropolis Monte Carlo. The model is a 2D box of N particles interacting through the Lennard-Jones potential `U(r) = 1/(12*r^12) - 1/(6*r^6)` with minimum-image periodic boundaries, sampled from the Boltzmann distribution `proportional to exp(-U/T)` over positions. Use the Metropolis-Hastings algorithm with local single-particle moves and an incremental energy update. Test it on 25 particles on a 5x5 grid at spacing 1.25 (box side `a = 6.25`), move size `dxmax = 0.25`, temperature `T = 0.05`, 5000 equilibration then 10000 sampling steps; seed 10. Produce a plot of the box energy over the equilibration and sampling phases. As a separate check, confirm that during equilibration the energy falls as particles condense from the open grid into a denser cluster, then fluctuates around a low plateau. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8D.1.1  (08-monte-carlo/08d-gillespie-algorithm.qmd)

::: {.callout-note title="Recipe 8D.1.1"}
**Objective.** Build a Gillespie stochastic simulation for a general chemical reaction network.

**Model.** A generic reaction network: a state vector of molecule counts, each reaction j with a propensity `R_j` and a stoichiometry vector `gamma_j` that updates the state.

**Method.** The Gillespie stochastic simulation algorithm (the direct method): compute all propensities, draw a waiting time from an exponential with rate equal to their sum, pick reaction j with probability `R_j/R_tot`, then advance time and update the state.

**Test.** A generic driver taking user-supplied propensity and stoichiometry functions, with stop conditions `tmax` and a maximum iteration count.

**Show.** A time-and-state trajectory for a reaction network (exercised in 8D.2).

**Verification.** The routine generates exact sample trajectories of a reaction network from its master equation.
:::

**Generated prompt (Python):** In Python, build a Gillespie stochastic simulation for a general chemical reaction network. The model is a generic reaction network: a state vector of molecule counts, each reaction j with a propensity `R_j` and a stoichiometry vector `gamma_j` that updates the state. Use the Gillespie stochastic simulation algorithm (the direct method): compute all propensities, draw a waiting time from an exponential with rate equal to their sum, pick reaction j with probability `R_j/R_tot`, then advance time and update the state. Test it on A generic driver taking user-supplied propensity and stoichiometry functions, with stop conditions `tmax` and a maximum iteration count. Produce a time-and-state trajectory for a reaction network (exercised in 8D.2). As a separate check, confirm that the routine generates exact sample trajectories of a reaction network from its master equation. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8D.2.1  (08-monte-carlo/08d-gillespie-algorithm.qmd)

::: {.callout-note title="Recipe 8D.2.1"}
**Objective.** Simulate a constitutively transcribed gene with Gillespie and measure its intrinsic noise.

**Model.** A birth-death process `0 -> x` at rate g and `x -> 0` at rate `k*x`, with deterministic steady state `x_bar = g/k`.

**Method.** The Gillespie stochastic simulation algorithm, with time-weighted mean and standard deviation.

**Test.** `k = 0.1`, production `g = 100, 10, 1` (steady states 1000, 100, 10), `tmax = 100`; started at the steady state and separately from zero; seed 101.

**Show.** Trajectories of the molecule count for each g, and the standard deviation versus mean against `sqrt(x_bar)`.

**Verification.** The standard deviation follows `sqrt(x_bar)` (Poisson), so relative noise falls as `1/sqrt(x_bar)`, and from zero the trajectories climb along `x_bar*(1 - exp(-k*t))`.
:::

**Generated prompt (Python):** In Python, simulate a constitutively transcribed gene with Gillespie and measure its intrinsic noise. The model is a birth-death process `0 -> x` at rate g and `x -> 0` at rate `k*x`, with deterministic steady state `x_bar = g/k`. Use the Gillespie stochastic simulation algorithm, with time-weighted mean and standard deviation. Test it on `k = 0.1`, production `g = 100, 10, 1` (steady states 1000, 100, 10), `tmax = 100`; started at the steady state and separately from zero; seed 101. Produce trajectories of the molecule count for each g, and the standard deviation versus mean against `sqrt(x_bar)`. As a separate check, confirm that the standard deviation follows `sqrt(x_bar)` (Poisson), so relative noise falls as `1/sqrt(x_bar)`, and from zero the trajectories climb along `x_bar*(1 - exp(-k*t))`. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8D.3.1  (08-monte-carlo/08d-gillespie-algorithm.qmd)

::: {.callout-note title="Recipe 8D.3.1"}
**Objective.** Show that transcriptional bursting adds noise beyond the Poisson floor.

**Model.** A bursting birth-death process `0 -> n*x` at rate `g/n` and `x -> 0` at rate `k*x`, so each event adds n molecules while the mean production rate g stays fixed.

**Method.** The Gillespie stochastic simulation algorithm, varying the burst size at fixed mean rate.

**Test.** `g = 16`, `k = 0.1` (steady state 160), burst sizes `n = 1, 2, 4, 8`, started at steady state, `tmax = 2000`; seed 91.

**Show.** Trajectories for each burst size n, and the mean and standard deviation versus n.

**Verification.** The mean stays at 160 for every n, but the standard deviation grows with burst size, so larger rarer bursts inject more noise than single-molecule events.
:::

**Generated prompt (Python):** In Python, show that transcriptional bursting adds noise beyond the Poisson floor. The model is a bursting birth-death process `0 -> n*x` at rate `g/n` and `x -> 0` at rate `k*x`, so each event adds n molecules while the mean production rate g stays fixed. Use the Gillespie stochastic simulation algorithm, varying the burst size at fixed mean rate. Test it on `g = 16`, `k = 0.1` (steady state 160), burst sizes `n = 1, 2, 4, 8`, started at steady state, `tmax = 2000`; seed 91. Produce trajectories for each burst size n, and the mean and standard deviation versus n. As a separate check, confirm that the mean stays at 160 for every n, but the standard deviation grows with burst size, so larger rarer bursts inject more noise than single-molecule events. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 8D.4.1  (08-monte-carlo/08d-gillespie-algorithm.qmd)

::: {.callout-note title="Recipe 8D.4.1"}
**Objective.** Show that negative autoregulation suppresses gene-expression noise.

**Model.** A self-repressing gene with promoter states (unbound D, bound C) and molecules M: production from D at rate g0 and from C at the slower g1, degradation `M -> 0` at rate k, and dimer binding/unbinding `2M + D <-> C` at rates kon, koff.

**Method.** The Gillespie stochastic simulation algorithm for a gene circuit with promoter binding.

**Test.** `k = 0.1` throughout; self-inhibition (`g0 = 55, g1 = 5, kon = 0.002, koff = 90`), no feedback (`g0 = g1 = 30`), and low copy number (`g0 = 5.5, g1 = 0.5, kon = 0.2, koff = 90`); seed 101.

**Show.** Trajectories of the molecule count for each regime, with the mean and standard deviation.

**Verification.** The self-inhibiting gene has a lower standard deviation than the no-feedback case, and at low copy number the few-molecule binding makes the dynamics strongly bursty.
:::

**Generated prompt (Python):** In Python, show that negative autoregulation suppresses gene-expression noise. The model is a self-repressing gene with promoter states (unbound D, bound C) and molecules M: production from D at rate g0 and from C at the slower g1, degradation `M -> 0` at rate k, and dimer binding/unbinding `2M + D <-> C` at rates kon, koff. Use the Gillespie stochastic simulation algorithm for a gene circuit with promoter binding. Test it on `k = 0.1` throughout; self-inhibition (`g0 = 55, g1 = 5, kon = 0.002, koff = 90`), no feedback (`g0 = g1 = 30`), and low copy number (`g0 = 5.5, g1 = 0.5, kon = 0.2, koff = 90`); seed 101. Produce trajectories of the molecule count for each regime, with the mean and standard deviation. As a separate check, confirm that the self-inhibiting gene has a lower standard deviation than the no-feedback case, and at low copy number the few-molecule binding makes the dynamics strongly bursty. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 9A.2.1  (09-optimization/09a-mcmc-optimization.qmd)

::: {.callout-note title="Recipe 9A.2.1"}
**Objective.** Explore a multi-minimum landscape with Metropolis-Hastings at fixed temperature.

**Model.** Himmelblau's function `f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2`, which has four equal minima at f = 0.

**Method.** The Metropolis-Hastings algorithm at fixed temperature, with a uniform proposal displacement and acceptance `a = min(1, exp(-de/T))`.

**Test.** Start at (0, 0), `1e4` steps, step size 0.5, temperatures `T = 1, 10, 30, 50`; seed 1 per temperature.

**Show.** The sampled points at each temperature on the Himmelblau landscape.

**Verification.** No single temperature both explores widely and settles: T = 1 traps in one basin, while T = 50 roams but never pins a minimum.
:::

**Generated prompt (Python):** In Python, explore a multi-minimum landscape with Metropolis-Hastings at fixed temperature. The model is himmelblau's function `f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2`, which has four equal minima at f = 0. Use the Metropolis-Hastings algorithm at fixed temperature, with a uniform proposal displacement and acceptance `a = min(1, exp(-de/T))`. Test it on Start at (0, 0), `1e4` steps, step size 0.5, temperatures `T = 1, 10, 30, 50`; seed 1 per temperature. Produce the sampled points at each temperature on the Himmelblau landscape. As a separate check, confirm that no single temperature both explores widely and settles: T = 1 traps in one basin, while T = 50 roams but never pins a minimum. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 9A.3.1  (09-optimization/09a-mcmc-optimization.qmd)

::: {.callout-note title="Recipe 9A.3.1"}
**Objective.** Find a global minimum by cooling a Metropolis chain.

**Model.** Himmelblau's function `f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2`, which has four equal minima at f = 0.

**Method.** Simulated annealing, lowering the temperature during a single run on a geometric (or linear) cooling schedule with Metropolis acceptance.

**Test.** Start at (0, 0), `1e4` steps, step size 0.5, maximum temperature 50, geometric schedule (scaling 0.999); four runs with seeds 1 to 4.

**Show.** The path of each annealing run on the landscape, and the best value over steps.

**Verification.** Each run cools from a broad search into a single minimum, but which minimum depends on the random path, so different seeds land in different basins.
:::

**Generated prompt (Python):** In Python, find a global minimum by cooling a Metropolis chain. The model is himmelblau's function `f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2`, which has four equal minima at f = 0. Use simulated annealing, lowering the temperature during a single run on a geometric (or linear) cooling schedule with Metropolis acceptance. Test it on Start at (0, 0), `1e4` steps, step size 0.5, maximum temperature 50, geometric schedule (scaling 0.999); four runs with seeds 1 to 4. Produce the path of each annealing run on the landscape, and the best value over steps. As a separate check, confirm that each run cools from a broad search into a single minimum, but which minimum depends on the random path, so different seeds land in different basins. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 9A.4.1  (09-optimization/09a-mcmc-optimization.qmd)

::: {.callout-note title="Recipe 9A.4.1"}
**Objective.** Sample every basin in one run by exchanging replicas across temperatures.

**Model.** Himmelblau's function `f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2`, which has four equal minima at f = 0.

**Method.** Parallel tempering (replica exchange): run a Metropolis chain at each of several temperatures and swap adjacent-temperature replicas with acceptance `a = min(1, exp((f_i - f_j)*(1/T_i - 1/T_j)))`.

**Test.** 6 replicas at `T = 2.3, 5, 10, 20, 40, 80` with step sizes 0.3 to 2.5, all started at the origin, 100 swap rounds of 100 steps; seed 1.

**Show.** The sampled points of the coldest and hottest replicas on the landscape.

**Verification.** The coldest replica stays inside minima yet still hops between them via configurations passed down from hotter replicas, while the hottest roams the whole landscape.
:::

**Generated prompt (Python):** In Python, sample every basin in one run by exchanging replicas across temperatures. The model is himmelblau's function `f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2`, which has four equal minima at f = 0. Use parallel tempering (replica exchange): run a Metropolis chain at each of several temperatures and swap adjacent-temperature replicas with acceptance `a = min(1, exp((f_i - f_j)*(1/T_i - 1/T_j)))`. Test it on 6 replicas at `T = 2.3, 5, 10, 20, 40, 80` with step sizes 0.3 to 2.5, all started at the origin, 100 swap rounds of 100 steps; seed 1. Produce the sampled points of the coldest and hottest replicas on the landscape. As a separate check, confirm that the coldest replica stays inside minima yet still hops between them via configurations passed down from hotter replicas, while the hottest roams the whole landscape. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 9A.5.1  (09-optimization/09a-mcmc-optimization.qmd)

::: {.callout-note title="Recipe 9A.5.1"}
**Objective.** Cross barriers and settle into minima within one chain by letting the temperature wander.

**Model.** Himmelblau's function `f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2`, which has four equal minima at f = 0.

**Method.** Simulated tempering: a single chain whose temperature random-walks a discrete ladder, with temperature moves accepted as `a = min(1, (c'/c)*exp(-f*(1/T' - 1/T)))`.

**Test.** Start at (0, 0), `1e4` steps, mixing fraction 0.9, step size 1, a 16-rung ladder `T = 5, 10, ..., 80` with weights decreasing linearly from 5 to 1; seed 1.

**Show.** The single chain's path over the landscape and its temperature over steps.

**Verification.** Because the temperature wanders up and down, one trajectory both crosses barriers when hot and settles into minima when cold, visiting all four basins.
:::

**Generated prompt (Python):** In Python, cross barriers and settle into minima within one chain by letting the temperature wander. The model is himmelblau's function `f(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2`, which has four equal minima at f = 0. Use simulated tempering: a single chain whose temperature random-walks a discrete ladder, with temperature moves accepted as `a = min(1, (c'/c)*exp(-f*(1/T' - 1/T)))`. Test it on Start at (0, 0), `1e4` steps, mixing fraction 0.9, step size 1, a 16-rung ladder `T = 5, 10, ..., 80` with weights decreasing linearly from 5 to 1; seed 1. Produce the single chain's path over the landscape and its temperature over steps. As a separate check, confirm that because the temperature wanders up and down, one trajectory both crosses barriers when hot and settles into minima when cold, visiting all four basins. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 9B.1.1  (09-optimization/09b-dynamic-programming.qmd)

::: {.callout-note title="Recipe 9B.1.1"}
**Objective.** Compute the optimal end-to-end alignment of two sequences.

**Model.** Sequence-alignment scoring with a match score `s_match`, mismatch `s_mismatch`, and gap penalty g; the score matrix follows `F_{i,j} = max(F_{i-1,j-1} + S(A_i, B_j), F_{i,j-1} + g, F_{i-1,j} + g)`.

**Method.** The Needleman-Wunsch global alignment dynamic program (fill the score matrix, then trace back from the last cell).

**Test.** Sequences `A = GAATTCAGTTA`, `B = GGATCGA`; two schemes, match-only (`s_match = 1`, others 0) and (`s_match = 3, s_mismatch = -3, gap = -2`).

**Show.** The optimal alignment and its score for each scoring scheme.

**Verification.** Rewarding only matches maximizes matches regardless of gaps, while adding mismatch and gap penalties balances matches against gap cost, giving a different optimal alignment.
:::

**Generated prompt (Python):** In Python, compute the optimal end-to-end alignment of two sequences. The model is sequence-alignment scoring with a match score `s_match`, mismatch `s_mismatch`, and gap penalty g; the score matrix follows `F_{i,j} = max(F_{i-1,j-1} + S(A_i, B_j), F_{i,j-1} + g, F_{i-1,j} + g)`. Use the Needleman-Wunsch global alignment dynamic program (fill the score matrix, then trace back from the last cell). Test it on Sequences `A = GAATTCAGTTA`, `B = GGATCGA`; two schemes, match-only (`s_match = 1`, others 0) and (`s_match = 3, s_mismatch = -3, gap = -2`). Produce the optimal alignment and its score for each scoring scheme. As a separate check, confirm that rewarding only matches maximizes matches regardless of gaps, while adding mismatch and gap penalties balances matches against gap cost, giving a different optimal alignment. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 9B.2.1  (09-optimization/09b-dynamic-programming.qmd)

::: {.callout-note title="Recipe 9B.2.1"}
**Objective.** Find the best-matching subsegment of two sequences.

**Model.** The same match/mismatch/gap scoring, with negative scores clamped to zero, `H_{i,j} = max(H_{i-1,j-1} + S(A_i, B_j), H_{i,j-1} + g, H_{i-1,j} + g, 0)`.

**Method.** The Smith-Waterman local alignment dynamic program (trace back from the largest entry, stopping at a zero).

**Test.** Sequences `A = GAATTCAGTTA`, `B = GGATCGA`, scoring `s_match = 3, s_mismatch = -3, gap = -2`.

**Show.** The optimal local alignment, next to the global alignment for comparison.

**Verification.** The global alignment spans both sequences end to end, while the local alignment reports only the best-matching internal segment.
:::

**Generated prompt (Python):** In Python, find the best-matching subsegment of two sequences. The model is the same match/mismatch/gap scoring, with negative scores clamped to zero, `H_{i,j} = max(H_{i-1,j-1} + S(A_i, B_j), H_{i,j-1} + g, H_{i-1,j} + g, 0)`. Use the Smith-Waterman local alignment dynamic program (trace back from the largest entry, stopping at a zero). Test it on Sequences `A = GAATTCAGTTA`, `B = GGATCGA`, scoring `s_match = 3, s_mismatch = -3, gap = -2`. Produce the optimal local alignment, next to the global alignment for comparison. As a separate check, confirm that the global alignment spans both sequences end to end, while the local alignment reports only the best-matching internal segment. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 9B.3.1  (09-optimization/09b-dynamic-programming.qmd)

::: {.callout-note title="Recipe 9B.3.1"}
**Objective.** Find the exact shortest tour through a set of cities.

**Model.** The traveling salesman problem: the shortest closed route visiting each of n cities once, from a Euclidean distance matrix; the Held-Karp recurrence is `g(x, S) = min over i in S of (g(i, S without i) + d(i, x))`.

**Method.** Dynamic programming for the traveling salesman problem (the Held-Karp algorithm).

**Test.** 10 cities at fixed coordinates, the tour starting and ending at city 1.

**Show.** The shortest tour drawn through the ten cities, with its length.

**Verification.** It returns the exact shortest tour through all ten cities; the cost is `O(n^2 * 2^n)`, so it stays practical only for small n.
:::

**Generated prompt (Python):** In Python, find the exact shortest tour through a set of cities. The model is the traveling salesman problem: the shortest closed route visiting each of n cities once, from a Euclidean distance matrix; the Held-Karp recurrence is `g(x, S) = min over i in S of (g(i, S without i) + d(i, x))`. Use dynamic programming for the traveling salesman problem (the Held-Karp algorithm). Test it on 10 cities at fixed coordinates, the tour starting and ending at city 1. Produce the shortest tour drawn through the ten cities, with its length. As a separate check, confirm that it returns the exact shortest tour through all ten cities; the cost is `O(n^2 * 2^n)`, so it stays practical only for small n. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 9C.1.1  (09-optimization/09c-genetic-algorithm.qmd)

::: {.callout-note title="Recipe 9C.1.1"}
**Objective.** Minimize a rugged multimodal function with a genetic algorithm.

**Model.** The Rastrigin function `f(x1, x2) = 20 + sum(x_i^2 - 10*cos(2*pi*x_i))`, whose global minimum f = 0 sits at the origin amid many local minima.

**Method.** A genetic algorithm for continuous variables: crossover swaps whole variables between parents, mutation nudges one variable, and the fittest members are kept each generation (elitist selection).

**Test.** Two variables on `[-5.12, 5.12]^2`, population 50, crossover rate 0.7, 100 generations, mutation step 10% of the range; seed 1.

**Show.** The best score versus generation, converging toward zero.

**Verification.** The best score drops toward zero, converging on the global minimum at the origin despite the many local minima.
:::

**Generated prompt (Python):** In Python, minimize a rugged multimodal function with a genetic algorithm. The model is the Rastrigin function `f(x1, x2) = 20 + sum(x_i^2 - 10*cos(2*pi*x_i))`, whose global minimum f = 0 sits at the origin amid many local minima. Use a genetic algorithm for continuous variables: crossover swaps whole variables between parents, mutation nudges one variable, and the fittest members are kept each generation (elitist selection). Test it on Two variables on `[-5.12, 5.12]^2`, population 50, crossover rate 0.7, 100 generations, mutation step 10% of the range; seed 1. Produce the best score versus generation, converging toward zero. As a separate check, confirm that the best score drops toward zero, converging on the global minimum at the origin despite the many local minima. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 9C.2.1  (09-optimization/09c-genetic-algorithm.qmd)

::: {.callout-note title="Recipe 9C.2.1"}
**Objective.** Approximate the shortest tour of many cities with a genetic algorithm.

**Model.** The traveling salesman problem with the tour-length objective (the sum of consecutive city distances plus the return leg); a candidate solution is a permutation of the city indices.

**Method.** A genetic algorithm on permutations, using permutation-preserving mutation (swap two cities, or reverse a segment) and order-preserving crossover.

**Test.** A 10-city instance, population 10, crossover rate 0.8, 50 generations.

**Show.** The best tour and its length versus generation.

**Verification.** The genetic algorithm recovers the exact shortest tour found by dynamic programming on the small instance, and scales to larger instances where exact methods cannot.
:::

**Generated prompt (Python):** In Python, approximate the shortest tour of many cities with a genetic algorithm. The model is the traveling salesman problem with the tour-length objective (the sum of consecutive city distances plus the return leg); a candidate solution is a permutation of the city indices. Use a genetic algorithm on permutations, using permutation-preserving mutation (swap two cities, or reverse a segment) and order-preserving crossover. Test it on A 10-city instance, population 10, crossover rate 0.8, 50 generations. Produce the best tour and its length versus generation. As a separate check, confirm that the genetic algorithm recovers the exact shortest tour found by dynamic programming on the small instance, and scales to larger instances where exact methods cannot. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10A.2.1  (10-high-dim-data/10a-dimensionality-reduction.qmd)

::: {.callout-note title="Recipe 10A.2.1"}
**Objective.** Compute principal components from scratch and confirm two routes agree.

**Model.** The input is a data matrix X with samples in rows and features in columns.

**Method.** Principal component analysis two equivalent ways: the eigen-decomposition of the covariance matrix, and the singular value decomposition of the centered data.

**Test.** The 100x2 "strip" data from 10A.1 (points scattered along a line); centering only, no scaling.

**Show.** The component variances from both routes, and the PC1 direction drawn on the strip data.

**Verification.** The two routes return identical variances, and on the strip data PC1 carries almost all the variance and points along the strip, recovering the regression direction of 10A.1.
:::

**Generated prompt (Python):** In Python, compute principal components from scratch and confirm two routes agree. The model is the input is a data matrix X with samples in rows and features in columns. Use principal component analysis two equivalent ways: the eigen-decomposition of the covariance matrix, and the singular value decomposition of the centered data. Test it on The 100x2 "strip" data from 10A.1 (points scattered along a line); centering only, no scaling. Produce the component variances from both routes, and the PC1 direction drawn on the strip data. As a separate check, confirm that the two routes return identical variances, and on the strip data PC1 carries almost all the variance and points along the strip, recovering the regression direction of 10A.1. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10A.4.1  (10-high-dim-data/10a-dimensionality-reduction.qmd)

::: {.callout-note title="Recipe 10A.4.1"}
**Objective.** Project samples onto their leading principal components and summarize the trend with a principal curve.

**Model.** The input is the PCA scores of 26 samples by 500 genes, projected onto the first two components; each sample is labeled by its condition (A through I).

**Method.** A PCA scree plot (variance explained per component) with a 2D PC1-PC2 projection, and a fitted Hastie-Stuetzle principal curve.

**Test.** The gene-expression data (26 samples x 500 genes), colored by condition.

**Show.** A scree plot and the PC1-PC2 scatter colored by condition, with the fitted principal curve.

**Verification.** PC1 explains most of the variance, the samples fall in an ordered arc from condition A to I, and the principal curve summarizes that progression as a 1D coordinate per sample.
:::

**Generated prompt (Python):** In Python, project samples onto their leading principal components and summarize the trend with a principal curve. The model is the input is the PCA scores of 26 samples by 500 genes, projected onto the first two components; each sample is labeled by its condition (A through I). Use a PCA scree plot (variance explained per component) with a 2D PC1-PC2 projection, and a fitted Hastie-Stuetzle principal curve. Test it on The gene-expression data (26 samples x 500 genes), colored by condition. Produce a scree plot and the PC1-PC2 scatter colored by condition, with the fitted principal curve. As a separate check, confirm that pC1 explains most of the variance, the samples fall in an ordered arc from condition A to I, and the principal curve summarizes that progression as a 1D coordinate per sample. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10A.5.1  (10-high-dim-data/10a-dimensionality-reduction.qmd)

::: {.callout-note title="Recipe 10A.5.1"}
**Objective.** Embed high-dimensional samples in 2D three different ways and compare what each preserves.

**Model.** The input is the same 26 samples by 500 genes, colored by condition.

**Method.** Three 2D nonlinear embeddings: classical (Torgerson) multidimensional scaling, t-SNE, and UMAP.

**Test.** The gene-expression data; t-SNE perplexity 5, UMAP 5 neighbors; seed 1.

**Show.** The MDS, t-SNE, and UMAP 2D embeddings, points colored by condition.

**Verification.** All three separate the nine conditions; MDS lays them along the A-to-I progression (it preserves global distances), while t-SNE and UMAP emphasize tight local clusters whose between-cluster distances are not meaningful.
:::

**Generated prompt (Python):** In Python, embed high-dimensional samples in 2D three different ways and compare what each preserves. The model is the input is the same 26 samples by 500 genes, colored by condition. Use three 2D nonlinear embeddings: classical (Torgerson) multidimensional scaling, t-SNE, and UMAP. Test it on The gene-expression data; t-SNE perplexity 5, UMAP 5 neighbors; seed 1. Produce the MDS, t-SNE, and UMAP 2D embeddings, points colored by condition. As a separate check, confirm that all three separate the nine conditions; MDS lays them along the A-to-I progression (it preserves global distances), while t-SNE and UMAP emphasize tight local clusters whose between-cluster distances are not meaningful. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10B.1.1  (10-high-dim-data/10b-clustering.qmd)

::: {.callout-note title="Recipe 10B.1.1"}
**Objective.** Cluster points into groups with k-means implemented from scratch.

**Model.** Three Gaussian blobs in 2D (50 points each, at means 0, 1.5, and 3, standard deviation 0.5); the input is an n x 2 point matrix.

**Method.** k-means clustering (Lloyd's algorithm): random centroids, then alternating nearest-centroid assignment and centroid-mean updates, keeping the lowest within-cluster sum of squares over several restarts.

**Test.** `k = 3`, 10 restarts, up to 100 iterations; seed 123.

**Show.** A scatter of the points colored by cluster, with the three centroids marked.

**Verification.** k-means recovers the three blobs, placing a centroid at the middle of each, and the restarts guard against an unlucky initialization converging to a poor split.
:::

**Generated prompt (Python):** In Python, cluster points into groups with k-means implemented from scratch. The model is three Gaussian blobs in 2D (50 points each, at means 0, 1.5, and 3, standard deviation 0.5); the input is an n x 2 point matrix. Use k-means clustering (Lloyd's algorithm): random centroids, then alternating nearest-centroid assignment and centroid-mean updates, keeping the lowest within-cluster sum of squares over several restarts. Test it on `k = 3`, 10 restarts, up to 100 iterations; seed 123. Produce a scatter of the points colored by cluster, with the three centroids marked. As a separate check, confirm that k-means recovers the three blobs, placing a centroid at the middle of each, and the restarts guard against an unlucky initialization converging to a poor split. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10B.2.1  (10-high-dim-data/10b-clustering.qmd)

::: {.callout-note title="Recipe 10B.2.1"}
**Objective.** Cluster points by building a dendrogram instead of a flat partition.

**Model.** The input is the same three-blob data, via pairwise Euclidean distances.

**Method.** Agglomerative hierarchical clustering with Ward linkage, cut into a flat partition.

**Test.** The three-blob data, cut into `k = 3` clusters.

**Show.** The dendrogram, and the points colored by the three cut clusters.

**Verification.** The dendrogram shows three clear branches, and cutting into three recovers the blobs, agreeing with k-means.
:::

**Generated prompt (Python):** In Python, cluster points by building a dendrogram instead of a flat partition. The model is the input is the same three-blob data, via pairwise Euclidean distances. Use agglomerative hierarchical clustering with Ward linkage, cut into a flat partition. Test it on The three-blob data, cut into `k = 3` clusters. Produce the dendrogram, and the points colored by the three cut clusters. As a separate check, confirm that the dendrogram shows three clear branches, and cutting into three recovers the blobs, agreeing with k-means. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10B.3.1  (10-high-dim-data/10b-clustering.qmd)

::: {.callout-note title="Recipe 10B.3.1"}
**Objective.** Cluster genes by expression and display them as a heatmap.

**Model.** The input is the gene-expression data (500 genes x 26 samples), z-scored per gene.

**Method.** k-means clustering of the z-scored genes, shown as a cluster-ordered expression heatmap.

**Test.** `k = 5` clusters; heatmap color scale from -3 to 3; seed 1.

**Show.** A cluster-ordered heatmap of the gene-expression data.

**Verification.** The heatmap splits genes into co-expression modules that rise and fall together across the conditions, resolving the PCA progression of Part 10A into specific gene groups.
:::

**Generated prompt (Python):** In Python, cluster genes by expression and display them as a heatmap. The model is the input is the gene-expression data (500 genes x 26 samples), z-scored per gene. Use k-means clustering of the z-scored genes, shown as a cluster-ordered expression heatmap. Test it on `k = 5` clusters; heatmap color scale from -3 to 3; seed 1. Produce a cluster-ordered heatmap of the gene-expression data. As a separate check, confirm that the heatmap splits genes into co-expression modules that rise and fall together across the conditions, resolving the PCA progression of Part 10A into specific gene groups. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10B.4.1  (10-high-dim-data/10b-clustering.qmd)

::: {.callout-note title="Recipe 10B.4.1"}
**Objective.** Cluster points with a Gaussian mixture and pick the number of components automatically.

**Model.** A Gaussian mixture: the data modeled as a mixture of Gaussian components, each with its own mean, covariance, and weight, with soft assignments; the input is the three-blob data.

**Method.** Gaussian mixture models fit by expectation-maximization, with the number of components chosen by the Bayesian information criterion.

**Test.** The three-blob data, choosing among 1 to 6 components; seed 1.

**Show.** A scatter of the points colored by mixture component, and the BIC versus number of components.

**Verification.** The criterion selects three components matching the blobs, and the mixture recovers them with soft probabilities.
:::

**Generated prompt (Python):** In Python, cluster points with a Gaussian mixture and pick the number of components automatically. The model is a Gaussian mixture: the data modeled as a mixture of Gaussian components, each with its own mean, covariance, and weight, with soft assignments; the input is the three-blob data. Use gaussian mixture models fit by expectation-maximization, with the number of components chosen by the Bayesian information criterion. Test it on The three-blob data, choosing among 1 to 6 components; seed 1. Produce a scatter of the points colored by mixture component, and the BIC versus number of components. As a separate check, confirm that the criterion selects three components matching the blobs, and the mixture recovers them with soft probabilities. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10C.1.1  (10-high-dim-data/10c-network-algorithms.qmd)

::: {.callout-note title="Recipe 10C.1.1"}
**Objective.** Compute the structural properties of a network.

**Model.** The input is an undirected adjacency matrix; the example is Zachary's karate club (34 vertices).

**Method.** Network properties computed from scratch: degree distribution, clustering coefficient, shortest-path distance matrix, and diameter, plus betweenness centrality.

**Test.** Zachary's karate club, drawn with a Kamada-Kawai layout.

**Show.** The printed degree, clustering, diameter, and betweenness, and the network drawn with a Kamada-Kawai layout.

**Verification.** The club has a few high-degree hubs, mean degree about 4.6, diameter 5, and clustering about 0.57, and the two highest-betweenness nodes are the instructor and president who bridge the factions.
:::

**Generated prompt (Python):** In Python, compute the structural properties of a network. The model is the input is an undirected adjacency matrix; the example is Zachary's karate club (34 vertices). Use network properties computed from scratch: degree distribution, clustering coefficient, shortest-path distance matrix, and diameter, plus betweenness centrality. Test it on Zachary's karate club, drawn with a Kamada-Kawai layout. Produce the printed degree, clustering, diameter, and betweenness, and the network drawn with a Kamada-Kawai layout. As a separate check, confirm that the club has a few high-degree hubs, mean degree about 4.6, diameter 5, and clustering about 0.57, and the two highest-betweenness nodes are the instructor and president who bridge the factions. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10C.2.1  (10-high-dim-data/10c-network-algorithms.qmd)

::: {.callout-note title="Recipe 10C.2.1"}
**Objective.** Generate random networks from two classic models and compare their degree distributions.

**Model.** The Erdos-Renyi model `G(n, p)` connects each vertex pair independently with probability p; the Barabasi-Albert model grows by attaching each new vertex to m existing ones with probability proportional to their degree (preferential attachment).

**Method.** The Erdos-Renyi and Barabasi-Albert random-network models, implemented from scratch.

**Test.** `n = 1000`; Erdos-Renyi with `p = 0.006`, Barabasi-Albert with `m = 3`; seed 1.

**Show.** The degree distributions: a histogram for Erdos-Renyi and a log-log scatter for Barabasi-Albert.

**Verification.** Erdos-Renyi degrees cluster tightly around the mean, while Barabasi-Albert degrees follow a straight line on log-log axes (a scale-free power law), matching biological networks better.
:::

**Generated prompt (Python):** In Python, generate random networks from two classic models and compare their degree distributions. The model is the Erdos-Renyi model `G(n, p)` connects each vertex pair independently with probability p; the Barabasi-Albert model grows by attaching each new vertex to m existing ones with probability proportional to their degree (preferential attachment). Use the Erdos-Renyi and Barabasi-Albert random-network models, implemented from scratch. Test it on `n = 1000`; Erdos-Renyi with `p = 0.006`, Barabasi-Albert with `m = 3`; seed 1. Produce the degree distributions: a histogram for Erdos-Renyi and a log-log scatter for Barabasi-Albert. As a separate check, confirm that erdos-Renyi degrees cluster tightly around the mean, while Barabasi-Albert degrees follow a straight line on log-log axes (a scale-free power law), matching biological networks better. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10C.3.1  (10-high-dim-data/10c-network-algorithms.qmd)

::: {.callout-note title="Recipe 10C.3.1"}
**Objective.** Split a network into communities by maximizing modularity.

**Model.** The input is the karate-club adjacency matrix; the modularity objective is `Q = (1/2m)*sum_ij (A_ij - k_i*k_j/2m)*delta(c_i, c_j)`.

**Method.** Newman spectral community detection (bisection by the leading eigenvector of the modularity matrix), with the eigenvector sign fixed so the labeling is reproducible.

**Test.** Zachary's karate club, drawn with the same Kamada-Kawai layout as 10C.1.

**Show.** The network drawn with nodes colored by community, and the modularity Q.

**Verification.** The spectral split separates the club into communities of 18 and 16 with modularity about 0.37, closely matching the factional split Zachary recorded.
:::

**Generated prompt (Python):** In Python, split a network into communities by maximizing modularity. The model is the input is the karate-club adjacency matrix; the modularity objective is `Q = (1/2m)*sum_ij (A_ij - k_i*k_j/2m)*delta(c_i, c_j)`. Use newman spectral community detection (bisection by the leading eigenvector of the modularity matrix), with the eigenvector sign fixed so the labeling is reproducible. Test it on Zachary's karate club, drawn with the same Kamada-Kawai layout as 10C.1. Produce the network drawn with nodes colored by community, and the modularity Q. As a separate check, confirm that the spectral split separates the club into communities of 18 and 16 with modularity about 0.37, closely matching the factional split Zachary recorded. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10C.4.1  (10-high-dim-data/10c-network-algorithms.qmd)

::: {.callout-note title="Recipe 10C.4.1"}
**Objective.** Enumerate the attractors of a Boolean network.

**Model.** A Boolean network: each gene is on or off and updates synchronously by `x_i(t+1) = f_i(x(t))`; iterating from any start falls into a fixed-point or cyclic attractor.

**Method.** Boolean network dynamics with synchronous updating, enumerating attractors exhaustively from all initial states.

**Test.** The generic `boolean_attractors(update, n)` driver, which runs every one of the `2^n` start states to its first repeat.

**Show.** The attractors enumerated for a given update rule (applied in 10C.5 and 10C.6).

**Verification.** The routine returns every fixed-point and cyclic attractor, ready to apply to specific circuits.
:::

**Generated prompt (Python):** In Python, enumerate the attractors of a Boolean network. The model is a Boolean network: each gene is on or off and updates synchronously by `x_i(t+1) = f_i(x(t))`; iterating from any start falls into a fixed-point or cyclic attractor. Use boolean network dynamics with synchronous updating, enumerating attractors exhaustively from all initial states. Test it on The generic `boolean_attractors(update, n)` driver, which runs every one of the `2^n` start states to its first repeat. Produce the attractors enumerated for a given update rule (applied in 10C.5 and 10C.6). As a separate check, confirm that the routine returns every fixed-point and cyclic attractor, ready to apply to specific circuits. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10C.5.1  (10-high-dim-data/10c-network-algorithms.qmd)

::: {.callout-note title="Recipe 10C.5.1"}
**Objective.** Find the attractors of a Boolean toggle switch.

**Model.** A Boolean toggle switch of two mutually repressing genes, `X(t+1) = NOT Y(t)`, `Y(t+1) = NOT X(t)`.

**Method.** Boolean attractor enumeration (from 10C.4) over all four states.

**Test.** The toggle-switch update rule on `n = 2` genes.

**Show.** The toggle switch's attractors printed as arrow chains.

**Verification.** It finds two fixed points, `01` and `10` (the mutually exclusive on-states of a bistable switch), and an oscillation `00 -> 11 -> 00`.
:::

**Generated prompt (Python):** In Python, find the attractors of a Boolean toggle switch. The model is a Boolean toggle switch of two mutually repressing genes, `X(t+1) = NOT Y(t)`, `Y(t+1) = NOT X(t)`. Use boolean attractor enumeration (from 10C.4) over all four states. Test it on The toggle-switch update rule on `n = 2` genes. Produce the toggle switch's attractors printed as arrow chains. As a separate check, confirm that it finds two fixed points, `01` and `10` (the mutually exclusive on-states of a bistable switch), and an oscillation `00 -> 11 -> 00`. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

## Recipe 10C.6.1  (10-high-dim-data/10c-network-algorithms.qmd)

::: {.callout-note title="Recipe 10C.6.1"}
**Objective.** Find the attractors of a Boolean repressilator.

**Model.** A Boolean repressilator, a ring of three repressions, `X(t+1) = NOT Z(t)`, `Y(t+1) = NOT X(t)`, `Z(t+1) = NOT Y(t)`.

**Method.** Boolean attractor enumeration (from 10C.4) over all eight states.

**Test.** The repressilator update rule on `n = 3` genes.

**Show.** The repressilator's attractors printed as arrow chains.

**Verification.** It finds no fixed points, only two cyclic attractors: `000 -> 111 -> 000` and a six-state cycle, the discrete analog of the continuous limit cycle.
:::

**Generated prompt (Python):** In Python, find the attractors of a Boolean repressilator. The model is a Boolean repressilator, a ring of three repressions, `X(t+1) = NOT Z(t)`, `Y(t+1) = NOT X(t)`, `Z(t+1) = NOT Y(t)`. Use boolean attractor enumeration (from 10C.4) over all eight states. Test it on The repressilator update rule on `n = 3` genes. Produce the repressilator's attractors printed as arrow chains. As a separate check, confirm that it finds no fixed points, only two cyclic attractors: `000 -> 111 -> 000` and a six-state cycle, the discrete analog of the continuous limit cycle. Implement the method explicitly with short comments rather than calling a routine that does it in one step, and explain in one sentence why that check confirms the result.

---

