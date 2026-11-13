"""Equation balancing by null-space of the element/charge matrix."""

from __future__ import annotations

import itertools
from fractions import Fraction
from math import gcd

from .formula import Species


class BalanceError(ValueError):
    pass


def _fraction(value: float) -> Fraction:
    return Fraction(str(value)).limit_denominator(1000)


def _rref(matrix: list[list[Fraction]]) -> tuple[list[list[Fraction]], list[int]]:
    a = [row[:] for row in matrix]
    rows = len(a)
    cols = len(a[0])
    pivots: list[int] = []
    r = 0
    for c in range(cols):
        pivot = None
        for i in range(r, rows):
            if a[i][c] != 0:
                pivot = i
                break
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        value = a[r][c]
        a[r] = [x / value for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c] != 0:
                factor = a[i][c]
                a[i] = [x - factor * y for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return a, pivots


def _nullspace(matrix: list[list[Fraction]], cols: int) -> list[list[Fraction]]:
    reduced, pivots = _rref(matrix)
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for f in free:
        vector = [Fraction(0)] * cols
        vector[f] = Fraction(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][f]
        basis.append(vector)
    return basis


def _to_integers(vector: list[Fraction]) -> list[int]:
    denominator = 1
    for value in vector:
        denominator = denominator * value.denominator // gcd(denominator, value.denominator)
    ints = [int(value * denominator) for value in vector]
    common = 0
    for value in ints:
        common = gcd(common, abs(value))
    if common > 1:
        ints = [value // common for value in ints]
    return ints


def _positive_integer_solution(basis: list[list[Fraction]]) -> list[int] | None:
    if not basis:
        return None
    if len(basis) == 1:
        vector = basis[0]
        if all(v <= 0 for v in vector):
            vector = [-v for v in vector]
        if all(v >= 0 for v in vector) and any(v > 0 for v in vector):
            return _to_integers(vector)
        return None

    best: tuple[int, list[int]] | None = None
    for coeffs in itertools.product(range(0, 7), repeat=len(basis)):
        if not any(coeffs):
            continue
        combined = [
            sum((c * b[i] for c, b in zip(coeffs, basis)), Fraction(0))
            for i in range(len(basis[0]))
        ]
        if all(v > 0 for v in combined):
            ints = _to_integers(combined)
            score = sum(ints)
            if best is None or score < best[0]:
                best = (score, ints)
            if score == len(ints):
                break
        elif all(v < 0 for v in combined):
            combined = [-v for v in combined]
            ints = _to_integers(combined)
            score = sum(ints)
            if best is None or score < best[0]:
                best = (score, ints)
    return best[1] if best else None


def balance(
    reactants: list[Species], products: list[Species]
) -> tuple[list[int], list[int], list[str]]:
    """Return balanced coefficients for reactants and products.

    Raises :class:`BalanceError` when no positive integer solution exists.
    """
    species = reactants + products
    elements: list[str] = []
    for sp in species:
        for symbol in sp.counts:
            if symbol not in elements:
                elements.append(symbol)

    matrix: list[list[Fraction]] = []
    for symbol in elements:
        row = []
        for index, sp in enumerate(species):
            count = _fraction(sp.counts.get(symbol, 0.0))
            row.append(count if index < len(reactants) else -count)
        matrix.append(row)

    if any(sp.charge for sp in species):
        row = []
        for index, sp in enumerate(species):
            charge = Fraction(sp.charge)
            row.append(charge if index < len(reactants) else -charge)
        matrix.append(row)

    if not matrix:
        raise BalanceError("nothing to balance")

    basis = _nullspace(matrix, len(species))
    solution = _positive_integer_solution(basis)
    if solution is None:
        raise BalanceError(
            "no positive whole-number coefficients balance this reaction "
            "(check the formulas and that all products are listed)"
        )

    unused = [sp.display for sp, value in zip(species, solution) if value == 0]
    if unused:
        raise BalanceError(
            "no balanced equation uses every species listed; these would have "
            "to be left out entirely: " + ", ".join(unused)
        )

    r = solution[: len(reactants)]
    p = solution[len(reactants) :]
    notes: list[str] = []

    # Verify explicitly.
    for symbol in elements:
        left = sum(r[i] * reactants[i].counts.get(symbol, 0) for i in range(len(reactants)))
        right = sum(p[i] * products[i].counts.get(symbol, 0) for i in range(len(products)))
        if abs(left - right) > 1e-9:
            raise BalanceError(f"internal balance check failed for {symbol}")
    return r, p, notes


def equation_text(reactants, rcoeffs, products, pcoeffs) -> str:
    def side(species, coeffs):
        parts = []
        for sp, coeff in zip(species, coeffs):
            prefix = "" if coeff == 1 else f"{coeff} "
            parts.append(f"{prefix}{sp.display}")
        return " + ".join(parts)

    return f"{side(reactants, rcoeffs)} \u2192 {side(products, pcoeffs)}"
