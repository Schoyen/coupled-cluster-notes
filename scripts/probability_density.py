from sympy.physics.secondquant import evaluate_deltas, substitute_dummies

from sympy import Rational, latex

from one_body_density_matrices import (
    get_ccsd_t_operators,
    get_ccsd_lambda_operators,
    eval_equation,
    sub_kwargs,
)

T = sum(get_ccsd_t_operators())
T_t = sum(get_ccsd_t_operators(ast_symb="t(t)"))
L = sum(get_ccsd_lambda_operators())
L_t = sum(get_ccsd_lambda_operators(ast_symb="l(t)"))

tilde_t_eq = Rational(1, 1)
tilde_t_eq += eval_equation(L_t)
tilde_t_eq += eval_equation(-L_t * T_t)
tilde_t_eq += eval_equation(L_t * T)
tilde_t_eq += eval_equation(-L_t * T_t * T)
tilde_t_eq += eval_equation(Rational(1, 2) * L_t * T_t ** 2)
tilde_t_eq += eval_equation(Rational(1, 2) * L_t * T ** 2)

tilde_t_eq = tilde_t_eq.expand()
tilde_t_eq = evaluate_deltas(tilde_t_eq)
tilde_t_eq = substitute_dummies(tilde_t_eq, **sub_kwargs)

print(latex(tilde_t_eq))
