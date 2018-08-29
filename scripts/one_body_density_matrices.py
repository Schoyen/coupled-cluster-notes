from sympy.physics.secondquant import (
    AntiSymmetricTensor,
    wicks,
    F,
    Fd,
    NO,
    evaluate_deltas,
    substitute_dummies,
)
from sympy import symbols, Dummy, latex, Rational

pretty_dummies_dict = {"above": "cdefgh", "below": "klmno", "general": "pqrstu"}

wicks_kwargs = {
    "simplify_dummies": True,
    "keep_only_fully_contracted": True,
    "simplify_kronecker_deltas": True,
}

sub_kwargs = {"new_indices": True, "pretty_indices": pretty_dummies_dict}


def get_ccsd_t_operators():
    i = symbols("i", below_fermi=True, cls=Dummy)
    a = symbols("a", above_fermi=True, cls=Dummy)

    t_ai = AntiSymmetricTensor("t", (a,), (i,))
    c_ai = NO(Fd(a) * F(i))

    i, j = symbols("i, j", below_fermi=True, cls=Dummy)
    a, b = symbols("a, b", above_fermi=True, cls=Dummy)

    t_abij = AntiSymmetricTensor("t", (a, b), (i, j))
    c_abij = NO(Fd(a) * Fd(b) * F(j) * F(i))

    T_1 = t_ai * c_ai
    T_2 = t_abij * c_abij

    return (T_1, T_2)


def get_ccsd_lambda_operators():
    i = symbols("i", below_fermi=True, cls=Dummy)
    a = symbols("a", above_fermi=True, cls=Dummy)

    l_ia = AntiSymmetricTensor("l", (i,), (a,))
    c_ia = NO(Fd(i) * F(a))

    i, j = symbols("i, j", below_fermi=True, cls=Dummy)
    a, b = symbols("a, b", above_fermi=True, cls=Dummy)

    l_ijab = AntiSymmetricTensor("t", (i, j), (a, b))
    c_ijab = NO(Fd(i) * Fd(j) * F(b) * F(a))

    L_1 = l_ia * c_ia
    L_2 = l_ijab * c_ijab

    return (L_1, L_2)


p, q = symbols("p, q", cls=Dummy)
c_pq = Fd(p) * F(q)

T = sum(get_ccsd_t_operators())
L = sum(get_ccsd_lambda_operators())

# Only keep non-zero terms
rho_eq = (1 - T) * c_pq * (1 + T)
rho_eq += (
    L
    * (1 - T + Rational(1, 2) * T ** 2)
    * c_pq
    * (1 + T + Rational(1, 2) * T ** 2)
)

rho = wicks(rho_eq, **wicks_kwargs).expand()
rho = evaluate_deltas(rho)
rho = substitute_dummies(rho, **sub_kwargs)

print(latex(rho))
