import pickle
import os
import re
from generate_ccd_equations import get_ccd_equations
from sympy import Mul, Rational
from sympy.physics.secondquant import AntiSymmetricTensor

filename = "ccd_equations.p"

if not os.path.isfile(filename):
    ccd_equations = get_ccd_equations()

    with open(filename, "wb") as f:
        pickle.dump(ccd_equations, f)
else:
    with open(filename, "rb") as f:
        ccd_equations = pickle.load(f)

energy_equation, amplitude_equation = ccd_equations


def generate_amplitude_function(equation):
    for arg in equation.args:
        if isinstance(arg, Mul):
            print (arg.args)
            print (type(arg.args[0]))
            print (arg.args[0])
        else:
            print ("YOLO")
            print (arg)

def _generate_term_function(term):
    if isinstance(term, Mul):
        print (term.args)
        for term in term.args:
            pass
    else:
        print ("YOLO")
        print (term)
        print (re.sub(r"_", "", str(term.upper[0])))

def generate_energy_function(equation):
    for arg in equation.args:
        _generate_term_function(arg)

generate_energy_function(energy_equation)
print ("\n")
generate_amplitude_function(amplitude_equation)
