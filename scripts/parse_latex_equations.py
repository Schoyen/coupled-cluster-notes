import re
import os
from generate_ccd_equations import get_ccd_equations
from sympy import latex

filename = "ccd_equations.dat"

if not os.path.isfile(filename):
    ccd_equations = [latex(eq) for eq in get_ccd_equations()]

    with open(filename, "w") as f:
        for eq in ccd_equations:
            f.write(eq + "\n")
else:
    with open(filename, "r") as f:
        ccd_equations = f.read().strip().split("\n")

energy_equation, amplitude_equation = ccd_equations
print (energy_equation)

def _remove_whitespace_and_empty_strings(string_list):
    return list(
            map(lambda x: x.strip(), filter(lambda x: x != "", string_list)))

def _split_fraction(string):
    if string.startswith(r"\frac{"):
        string = re.sub(r"\\frac\{(.*)\}\{(\d)\}", r"\1/\2", string)
    return string

def _split_fractions(string_list):
    return [_split_fraction(string) for string in string_list]

    return new_string_list
