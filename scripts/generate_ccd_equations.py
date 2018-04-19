from sympy.physics.secondquant import (
        AntiSymmetricTensor, F, Fd, Commutator, wicks, evaluate_deltas,
        substitute_dummies, PermutationOperator, simplify_index_permutations
)
from sympy import symbols, Dummy, Rational, factorial, latex

pretty_dummies = {
    'above': 'cdefgh',
    'below': 'mno',
    'general': 'pqrstu'
}

wicks_kwargs = {
    "simplify_dummies": True,
    "keep_only_fully_contracted": True,
    "simplify_kronecker_deltas": True
}

sub_kwargs = {
    "new_indices": True,
    "pretty_indices": pretty_dummies
}

i, j, k, l = symbols("i, j, k, l", below_fermi=True)
a, b, c, d = symbols("a, b, c, d", above_fermi=True)

def get_hamiltonian():
    p, q, r, s = symbols("p, q, r, s", cls=Dummy)
    h = AntiSymmetricTensor("h", (p,), (q,))
    u = AntiSymmetricTensor("u", (p, r), (q, s))

    h = h * Fd(p) * F(q)
    u = u * Fd(p) * Fd(r) * F(s) * F(q)

    return h, Rational(1, 4) * u

def get_doubles_cluster_operator():
    i, j = symbols("i, j", below_fermi=True, cls=Dummy)
    a, b = symbols("a, b", above_fermi=True, cls=Dummy)

    t = AntiSymmetricTensor("t", (a, b), (i, j))
    t = t * Fd(a) * F(i) * Fd(b) * F(j)

    return [Rational(1, 4) * t]

def compute_hausdorff(h, cluster_func, num_terms=4):
    commutator = Commutator
    comm_term = h

    equation = comm_term

    for i in range(num_terms):
        t = sum(cluster_func())

        comm_term = wicks(commutator(comm_term, t))
        comm_term = substitute_dummies(evaluate_deltas(comm_term))

        equation += comm_term/factorial(i + 1)

    equation = equation.expand()
    equation = evaluate_deltas(equation)
    equation = substitute_dummies(
            equation, new_indices=True, pretty_indices=pretty_dummies)

    return equation


def get_energy_equation(equation_h, equation_u):
    energy = wicks(equation_h + equation_u, **wicks_kwargs)
    energy = substitute_dummies(energy, **sub_kwargs)

    return energy

def get_one_body_equation(equation_h, equation_u):
    one_body_eq = wicks(Fd(j) * F(b) *Fd(i) * F(a) * equation_h, **wicks_kwargs)

    p = PermutationOperator
    one_body_eq = simplify_index_permutations(one_body_eq, [p(a, b), p(i, j)])
    one_body_eq = substitute_dummies(one_body_eq, **sub_kwargs)

    return one_body_eq

def get_two_body_equation(equation_h, equation_u):
    two_body_eq = wicks(
            Fd(j) * F(b) * Fd(i) * F(a) * equation_u, **wicks_kwargs)

    p = PermutationOperator
    two_body_eq = simplify_index_permutations(two_body_eq, [p(i, j), p(a, b)])
    two_body_eq = substitute_dummies(two_body_eq, **sub_kwargs)

    return two_body_eq

def get_ccd_equations():
    h, u = get_hamiltonian()
    equation_h = compute_hausdorff(h, get_doubles_cluster_operator)
    equation_u = compute_hausdorff(u, get_doubles_cluster_operator)

    energy = get_energy_equation(equation_h, equation_u)
    one_body = get_one_body_equation(equation_h, equation_u)
    two_body = get_two_body_equation(equation_h, equation_u)

    return energy, one_body + two_body
