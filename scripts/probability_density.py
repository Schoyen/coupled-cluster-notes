from sympy.physics.secondquant import evaluate_deltas, substitute_dummies

from sympy import Rational, latex

from one_body_density_matrices import (
    get_ccsd_t_operators,
    get_ccsd_lambda_operators,
    eval_equation,
    sub_kwargs,
)

T = sum(get_ccsd_t_operators())
T_2 = sum(get_ccsd_t_operators())
T_t = sum(get_ccsd_t_operators(ast_symb="t(t)"))
T_t_2 = sum(get_ccsd_t_operators(ast_symb="t(t)"))
L = sum(get_ccsd_lambda_operators())
L_t = sum(get_ccsd_lambda_operators(ast_symb="l(t)"))

tilde_t_eq = Rational(1, 1)
tilde_t_eq += eval_equation(-L_t * T_t)
tilde_t_eq += eval_equation(L_t * T)
tilde_t_eq += eval_equation(-L_t * T_t * T)
tilde_t_eq += eval_equation(Rational(1, 2) * L_t * T_t * T_t_2)
tilde_t_eq += eval_equation(Rational(1, 2) * L_t * T * T_2)

tilde_t_eq = tilde_t_eq.expand()
tilde_t_eq = evaluate_deltas(tilde_t_eq)
tilde_t_eq = substitute_dummies(tilde_t_eq, **sub_kwargs)

tilde_eq = Rational(1, 1)
tilde_eq += eval_equation(-L * T)
tilde_eq += eval_equation(L * T_t)
tilde_eq += eval_equation(-L * T * T_t)
tilde_eq += eval_equation(Rational(1, 2) * L * T * T_2)
tilde_eq += eval_equation(Rational(1, 2) * L * T_t * T_t_2)

tilde_eq = tilde_eq.expand()
tilde_eq = evaluate_deltas(tilde_eq)
tilde_eq = substitute_dummies(tilde_eq, **sub_kwargs)

print(latex(tilde_t_eq))
print(latex(tilde_eq))
