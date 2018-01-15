from sympy.physics.secondquant import (AntiSymmetricTensor, wicks,
        F, Fd, NO, evaluate_deltas, substitute_dummies, Commutator,
        simplify_index_permutations, PermutationOperator)
from sympy import (
    symbols, expand, pprint, Rational, latex, Dummy
)

pretty_dummies_dict = {
    'above': 'cdefgh',
    'below': 'mno',
    'general': 'pqrstu'
}

def getHamiltonian():
    p, q, r, s = symbols('p,q,r,s', cls=Dummy)
    h = AntiSymmetricTensor('h', (p,), (q,))  
    pq = Fd(p)*F(q)
    u = AntiSymmetricTensor('u', (p, r), (q, s))
    pqsr = Fd(p)*Fd(r)*F(s)*F(q)
    H1 = h*pq
    H2 = Rational(1, 4)*u*pqsr
    return (H1,H2)

def get_T_operators():
    """
    Returns a tuple (T1,T2) of unique operators.
    """
    i = symbols('i', below_fermi=True, cls=Dummy)
    a = symbols('a', above_fermi=True, cls=Dummy)
    t_ai = AntiSymmetricTensor('t', (a,), (i,))
    ai = Fd(a)*F(i)
    i, j = symbols('i,j', below_fermi=True, cls=Dummy)
    a, b = symbols('a,b', above_fermi=True, cls=Dummy)
    t_abij = AntiSymmetricTensor('t', (a, b), (i, j))
    abji = Fd(a)*F(i)*Fd(b)*F(j)

    T1 = t_ai*ai
    T2 = Rational(1, 4)*t_abij*abji
    return (T1, T2)

def computeHausdorff(H):

    print("Calculating 4 nested commutators")
    C = Commutator

    T1, T2 = get_T_operators()
    T = T2
    print("commutator 1...")
    comm1 = wicks(C(H, T))
    comm1 = evaluate_deltas(comm1)
    comm1 = substitute_dummies(comm1)

    T1, T2 = get_T_operators()
    T = T2
    print("commutator 2...")
    comm2 = wicks(C(comm1, T))
    comm2 = evaluate_deltas(comm2)
    comm2 = substitute_dummies(comm2)
    
    T1, T2 = get_T_operators()
    T = T2
    print("commutator 3...")
    comm3 = wicks(C(comm2, T))
    comm3 = evaluate_deltas(comm3)
    comm3 = substitute_dummies(comm3)

    T1, T2 = get_T_operators()
    T = T2
    print("commutator 4...")
    comm4 = wicks(C(comm3, T))
    comm4 = evaluate_deltas(comm4)
    comm4 = substitute_dummies(comm4)

    print("construct Hausdorff expansion...")
    eq = H + comm1 + comm2/2 + comm3/6 + comm4/24
    eq = eq.expand()
    eq = evaluate_deltas(eq)
    eq = substitute_dummies(eq, new_indices=True,
            pretty_indices=pretty_dummies_dict)
    print("*********************")
    print()
    return eq

def main():
    i, j, k, l = symbols('i,j,k,l', below_fermi=True)
    a, b, c, d = symbols('a,b,c,d', above_fermi=True)
    H1, H2 = getHamiltonian()
    eq_H1 = computeHausdorff(H1) #e^(-T)H1e^T
    eq_H2 = computeHausdorff(H2)
    print("CC energy:")
    Energy = wicks(eq_H1+eq_H2, simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
    Energy = substitute_dummies(Energy,new_indices=True, pretty_indices=pretty_dummies_dict)
    print(latex(Energy))
    print()

    """
    print("CC dE_H1/(dLia):")
    eqH1T1 = wicks(Fd(i)*F(a)*eq_H1, simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
    P = PermutationOperator
    eqH1T1 = simplify_index_permutations(eqH1T1, [P(a, b), P(i, j)])
    eqH1T1 = substitute_dummies(eqH1T1,new_indices=True, pretty_indices=pretty_dummies_dict)
    print(latex(eqH1T1))
    print() 

	print("CC dE_H2/(dLia):")
	eqH2T1 = wicks(Fd(i)*F(a)*eq_H2, simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
	P = PermutationOperator
	eqH2T1 = simplify_index_permutations(eqH2T1, [P(a, b), P(i, j)])
	eqH2T1 = substitute_dummies(eqH2T1,new_indices=True, pretty_indices=pretty_dummies_dict)
	print(latex(eqH2T1))
	print()
    """

    print("CC dE_H1/(dLijab):")
    eqH1T2 = wicks(Fd(j)*F(b)*Fd(i)*F(a)*eq_H1, simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
    P = PermutationOperator
    eqH1T2 = simplify_index_permutations(eqH1T2, [P(a, b), P(i, j)])
    eqH1T2 = substitute_dummies(eqH1T2,new_indices=True, pretty_indices=pretty_dummies_dict)
    print(latex(eqH1T2))
    print()

    print("CC dE_H2/(dLijab):")
    eqH2T2 = wicks(Fd(j)*F(b)*Fd(i)*F(a)*eq_H2, simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
    P = PermutationOperator
    eqH2T2 = simplify_index_permutations(eqH2T2, [P(i, j), P(a, b)])
    eqH2T2 = substitute_dummies(eqH2T2,new_indices=True, pretty_indices=pretty_dummies_dict)
    print(latex(eqH2T2))
    print()

if __name__ == "__main__":
    main()
