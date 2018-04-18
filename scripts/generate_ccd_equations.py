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
    t = t * Fd(a) * Fd(b) * F(i) * F(j)

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


i, j, k, l = symbols("i, j, k, l", below_fermi=True)
a, b, c, d = symbols("a, b, c, d", above_fermi=True)

#n_ai = AntiSymmetricTensor("n", (a,), (i,))
#n_ia = AntiSymmetricTensor("n", (i,), (a,))
#
#Dai = n_ai * Fd(a) * F(i)
#Dia = n_ia * Fd(i) * F(a)
#D0 = Dai + Dia

h, u = get_hamiltonian()

equation_h = compute_hausdorff(h, get_doubles_cluster_operator)
equation_u = compute_hausdorff(u, get_doubles_cluster_operator)

one_body_eq = wicks(
        Fd(j) * Fd(i) * F(b) * F(a) * equation_h, simplify_dummies=True,
        keep_only_fully_contracted=True, simplify_kronecker_deltas=True)

p = PermutationOperator
print (latex(simplify_index_permutations(one_body_eq, [p(a, b), p(i, j)])))

two_body_eq = wicks(
        Fd(j) * Fd(i) * F(b) * F(a) * equation_u, simplify_dummies=True,
        keep_only_fully_contracted=True, simplify_kronecker_deltas=True)

print ("\n")
p = PermutationOperator
print (latex(simplify_index_permutations(two_body_eq, [p(i, j), p(a, b)])))
