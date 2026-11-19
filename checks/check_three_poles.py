"""Exact specialization at the first size with a nonconstant polynomial quotient.

Run with Python 3 and SymPy. The two computations use the pole recurrence and
polynomial division, respectively. This checks the normalization at both real
embeddings of Q(sqrt(5)). The accompanying manuscript gives the all-size proof.
"""

from math import prod

import sympy as sp


FIELD = sp.QQ.algebraic_field(sp.sqrt(5))
VARIABLE = sp.Symbol("Y")
ZERO = FIELD.zero
ONE = FIELD.one


def determinant(moments):
    a, b, c, d, e = moments
    return a * (c * e - d * d) - b * (b * e - c * d) + c * (b * d - c * c)


def normalized_coefficients(t):
    nodes = [t ** (2 * a) for a in (1, 2, 3)]
    weights = [
        ONE / prod(ONE - r / s for s in nodes if s != r)
        for r in nodes
    ]
    moments = [ZERO] + [j * t**j / (ONE + t ** (2 * j)) for j in range(1, 5)]
    tails = [
        sum(
            ((-1) ** ((m - 1) // 2) * t**m / (t**m - ONE) ** 2
             for m in range(1, 2 * a, 2)),
            ZERO,
        )
        for a in (1, 2, 3)
    ]
    # The nonzero E(3,d), in order d=1,...,5, are 6,2,6,1,2.
    clearer = (
        (t - ONE) ** 6 * (t + ONE) ** 2 * (t**2 + t + ONE) ** 6
        * (t**2 + ONE) * (t**4 + t**3 + t**2 + t + ONE) ** 2
    )
    samples = []
    for y in range(4):
        pole_values = [(-1) ** a * (y - tail) for a, tail in enumerate(tails, 1)]
        values = pole_values.copy()
        entries = []
        for k in range(5):
            via_recurrence = sum((w * v for w, v in zip(weights, values)), ZERO)
            via_residues = sum(
                (w * r**k * v for w, r, v in zip(weights, nodes, pole_values)),
                ZERO,
            )
            # Only x^4/D has a quotient of positive degree. The coefficient
            # of x there is the reciprocal of D's leading coefficient.
            if k == 4:
                via_residues -= prod(nodes) * moments[1]
            assert via_recurrence == via_residues
            entries.append(via_recurrence)
            values = [r * (v - moments[k]) for r, v in zip(nodes, values)]
        scaled = clearer * determinant(entries) / t**22
        samples.append((y, FIELD.to_sympy(scaled)))

    polynomial = sp.Poly(sp.interpolate(samples, VARIABLE), VARIABLE, extension=sp.sqrt(5))
    assert polynomial.degree() == 3
    assert FIELD.from_sympy(polynomial.LC()) == -t**2 * clearer
    return polynomial.all_coeffs()


def main():
    root = sp.sqrt(5)
    b = FIELD.from_sympy((3 + root) / 2)
    plus = normalized_coefficients(b)
    minus = normalized_coefficients(ONE / b)
    for a, conjugate in zip(plus, minus):
        assert sp.expand(a.xreplace({root: -root}) - conjugate) == 0
        # An element u+v*sqrt(5) belongs to Z[(1+sqrt(5))/2]
        # exactly when u-v and 2v are ordinary integers.
        v = a.coeff(root)
        u = a.subs(root, 0)
        assert (u - v).is_Integer and (2 * v).is_Integer
    print("PASS: all five moment entries at four probe values in both embeddings;")
    print("all four normalized coefficients, their conjugates and integrality;")
    print("and the signed leading coefficient at K=3.")


if __name__ == "__main__":
    main()
