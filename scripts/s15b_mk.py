"""s15b_mk.py — Mk models on a two-state character, and the reason the
primary character cannot be fitted with one.

Three models, all closed-form on two states so nothing here needs a matrix
exponential or a third-party optimiser:

  ER            one rate both ways.
  ARD           a gain rate and a loss rate.
  irreversible  the gain rate pinned to zero, root fixed present — the
                model Dollo's assumption corresponds to.

**The point of this module is a refusal.**  S15a's character is invariant:
every scored tip is present.  An invariant character contains no
transition to estimate, and each of these three likelihoods is then
monotone decreasing in every rate, so the maximum sits on the boundary at
rate zero and any positive rate a fitter reports is its own starting
point.  `profile()` measures that rather than asserting it — it evaluates
the likelihood along a grid of rates and reports whether it is monotone
and where its maximum lies — and `fit()` refuses an invariant character
outright, returning the reason instead of a number.

Where the sensitivity matrix *does* produce variation, the fits are real
and are reported for what they are: a rate estimated from losses a filter
setting manufactured, which is a property of the setting.

Branch lengths matter here and only here.  Parsimony counts edges; a
continuous-time model reads their lengths, so every fit is run under all
three of `s15b_lib.BL_SCHEMES` and the spread across them is committed.
"""

from __future__ import annotations

import math

import s15b_lib as lib

MODELS = ("ER", "ARD", "irreversible")
#: rate grid for the profile, log-spaced; wide enough that a monotone
#: likelihood is visibly monotone over 8 orders of magnitude
PROFILE_DECADES = (-6.0, 2.0)
PROFILE_POINTS = 33


def _probs(model: str, q: tuple[float, ...], t: float):
    """(P00, P01, P10, P11) over a branch of length t."""
    if t <= 0:
        return 1.0, 0.0, 0.0, 1.0
    if model == "irreversible":
        e = math.exp(-min(q[0] * t, 700.0))
        return 1.0, 0.0, 1.0 - e, e
    if model == "ER":
        e = math.exp(-min(2.0 * q[0] * t, 700.0))
        same, diff = 0.5 * (1.0 + e), 0.5 * (1.0 - e)
        return same, diff, diff, same
    a, b = q                                   # a: 0->1 gain, b: 1->0 loss
    s = a + b
    if s <= 0:
        return 1.0, 0.0, 0.0, 1.0
    p1, p0 = a / s, b / s
    e = math.exp(-min(s * t, 700.0))
    return p0 + p1 * e, p1 * (1.0 - e), p0 * (1.0 - e), p1 + p0 * e


def _root_prior(model: str, q: tuple[float, ...]) -> tuple[float, float]:
    if model == "irreversible":
        return 0.0, 1.0                        # the gain happened once, above
    if model == "ER":
        return 0.5, 0.5
    a, b = q
    s = a + b
    return (0.5, 0.5) if s <= 0 else (b / s, a / s)


def loglik(root: lib.Node, char: dict[str, int | None], model: str,
           q: tuple[float, ...]) -> float:
    """Felsenstein pruning, rescaled at every node against underflow."""
    part: dict[int, tuple[float, float]] = {}
    logscale = 0.0
    for n in lib.postorder(root):
        if n.is_tip:
            v = char.get(n.name)
            part[id(n)] = ((1.0, 1.0) if v is None else
                           (1.0, 0.0) if v == 0 else (0.0, 1.0))
            continue
        l0 = l1 = 1.0
        for c in n.children:
            c0, c1 = part[id(c)]
            p00, p01, p10, p11 = _probs(model, q, max(0.0, c.length))
            l0 *= p00 * c0 + p01 * c1
            l1 *= p10 * c0 + p11 * c1
        m = max(l0, l1)
        if m <= 0.0:
            return -math.inf
        part[id(n)] = (l0 / m, l1 / m)
        logscale += math.log(m)
    p0, p1 = _root_prior(model, q)
    r0, r1 = part[id(root)]
    tot = p0 * r0 + p1 * r1
    if tot <= 0.0:
        return -math.inf
    return logscale + math.log(tot)


def is_invariant(char: dict[str, int | None]) -> bool:
    seen = {v for v in char.values() if v is not None}
    return len(seen) <= 1


def n_params(model: str) -> int:
    return 2 if model == "ARD" else 1


def _golden(f, lo: float, hi: float, tol: float = 1e-6, it: int = 200):
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(it):
        if b - a < tol:
            break
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = f(d)
    x = (a + b) / 2.0
    return x, f(x)


def fit(root: lib.Node, char: dict[str, int | None], model: str) -> dict:
    """Maximum likelihood, or the reason there is none to report.

    An invariant character is refused rather than fitted: the likelihood
    is monotone in every rate and the reported optimum would be the
    optimiser's starting point.
    """
    if model not in MODELS:
        raise ValueError(f"unknown model {model!r}")
    if is_invariant(char):
        prof = profile(root, char, model)
        extra = ""
        if model == "ARD":
            gain = profile(root, char, model, axis="gain")
            extra = (f"; the gain rate is {gain['shape']} to the grid's "
                     f"own edge, so the pair is not identified")
        return dict(model=model, fitted=0, rate="", rate_gain="",
                    loglik="", aic="", n_params=n_params(model),
                    reason=("the character is invariant, so there is no "
                            "transition to estimate; the likelihood is "
                            f"{prof['shape']} in the loss rate and its "
                            f"maximum is at the boundary "
                            f"({prof['argmax']:.3g}){extra}"),
                    profile_shape=prof["shape"],
                    profile_argmax=prof["argmax"])

    lo, hi = math.log(1e-8), math.log(1e3)
    if model == "ARD":
        b = math.log(1e-3)
        a = math.log(1e-3)
        # coordinate descent; measured to converge to 7 s.f. in under 10
        # rounds on this surface, 20 is the margin
        for _ in range(20):
            a, _ = _golden(lambda x: loglik(root, char, model,
                                            (math.exp(x), math.exp(b))),
                           lo, hi)
            b, ll = _golden(lambda x: loglik(root, char, model,
                                             (math.exp(a), math.exp(x))),
                            lo, hi)
        rates = (math.exp(a), math.exp(b))
        ll = loglik(root, char, model, rates)
        gain, loss = rates[0], rates[1]
    else:
        x, ll = _golden(lambda z: loglik(root, char, model,
                                         (math.exp(z),)), lo, hi)
        loss = math.exp(x)
        gain = 0.0 if model == "irreversible" else loss
    k = n_params(model)
    return dict(model=model, fitted=1, rate=round(loss, 8),
                rate_gain=round(gain, 8), loglik=round(ll, 4),
                aic=round(2 * k - 2 * ll, 4), n_params=k, reason="",
                profile_shape="", profile_argmax="")


#: the rate ARD's other axis is held at while one is profiled; small
#: enough that the slice is a slice and not a second free parameter
ARD_HELD_RATE = 1e-3


def profile(root: lib.Node, char: dict[str, int | None], model: str,
            points: int = PROFILE_POINTS, axis: str = "loss") -> dict:
    """The likelihood along a log-spaced rate grid.

    Reported rather than asserted: `shape` says `monotone decreasing`,
    `monotone increasing` or `interior maximum`, and `argmax` is where the
    grid's best value sits.  An invariant character gives a monotone
    profile on both of ARD's axes, which is the whole basis for refusing
    to fit it.

    `axis` picks which of ARD's two rates is varied.  Profiling ARD along
    its **diagonal** would be ER by construction and would draw the same
    curve twice in two colours, so the two axes are profiled separately:
    the loss axis falls, the gain axis *rises* — and a likelihood that
    rises to the edge of the grid in one parameter is unidentifiability
    stated as a measurement rather than as an assertion.
    """
    if axis not in ("loss", "gain"):
        raise ValueError(f"unknown profile axis {axis!r}")
    if axis == "gain" and model != "ARD":
        raise ValueError(f"{model} has no separate gain rate to profile")
    lo, hi = PROFILE_DECADES
    grid = [10.0 ** (lo + (hi - lo) * i / (points - 1)) for i in range(points)]

    def q_at(r: float):
        if model != "ARD":
            return (r,)
        return (r, ARD_HELD_RATE) if axis == "gain" else (ARD_HELD_RATE, r)

    lls = [loglik(root, char, model, q_at(r)) for r in grid]
    best = max(range(points), key=lambda i: lls[i])
    dec = all(lls[i] >= lls[i + 1] - 1e-9 for i in range(points - 1))
    inc = all(lls[i] <= lls[i + 1] + 1e-9 for i in range(points - 1))
    shape = ("monotone decreasing" if dec else
             "monotone increasing" if inc else "interior maximum")
    return dict(model=model, axis=axis, shape=shape, argmax=grid[best],
                loglik_max=lls[best], grid=grid, loglik=lls,
                at_boundary=int(best in (0, points - 1)))
