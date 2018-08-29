from sympy.physics.secondquant import (
    AntiSymmetricTensor,
    wicks,
    F,
    Fd,
    NO,
    evaluate_deltas,
    substitute_dummies,
    Commutator,
)
from sympy import symbols, Dummy, latex, Rational, Mul

pretty_dummies_dict = {
    "above": "abcdef",
    "below": "ijklmn",
    "general": "pqrstu",
}

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
    T_2 = Rational(1, 4) * t_abij * c_abij

    return (T_1, T_2)


def get_ccsd_lambda_operators():
    i = symbols("i", below_fermi=True, cls=Dummy)
    a = symbols("a", above_fermi=True, cls=Dummy)

    l_ia = AntiSymmetricTensor("l", (i,), (a,))
    c_ia = NO(Fd(i) * F(a))

    i, j = symbols("i, j", below_fermi=True, cls=Dummy)
    a, b = symbols("a, b", above_fermi=True, cls=Dummy)

    l_ijab = AntiSymmetricTensor("l", (i, j), (a, b))
    c_ijab = NO(Fd(i) * Fd(j) * F(b) * F(a))

    L_1 = l_ia * c_ia
    L_2 = Rational(1, 4) * l_ijab * c_ijab

    return (L_1, L_2)


p, q = symbols("p, q", cls=Dummy)
c_pq = Fd(p) * F(q)

T = sum(get_ccsd_t_operators())
L = sum(get_ccsd_lambda_operators())


def eval_equation(eq):
    eq = wicks(eq, **wicks_kwargs)
    eq = evaluate_deltas(eq.expand())
    eq = substitute_dummies(eq, **sub_kwargs)

    return eq


# Only keep non-zero terms

rho_eq = eval_equation(c_pq)
rho_eq += eval_equation(Commutator(c_pq, T))
rho_eq += eval_equation(L * c_pq)
comm = Commutator(c_pq, T)
rho_eq += eval_equation(L * comm)
comm = Commutator(comm, sum(get_ccsd_t_operators()))
rho_eq += Rational(1, 2) * eval_equation(L * comm)
# rho_eq = wicks((1 - T) * c_pq * (1 + T), **wicks_kwargs)
# rho_eq += wicks(
#    L
#    * (1 - T + Rational(1, 2) * T * T)
#    * c_pq
#    * (1 + T + Rational(1, 2) * T * T),
#    **wicks_kwargs
# )

rho = rho_eq.expand()
# rho = wicks(rho_eq, **wicks_kwargs)
for term in rho.args:
    print(evaluate_deltas(term))
rho = evaluate_deltas(rho)
rho = substitute_dummies(rho, **sub_kwargs)

print(latex(rho))