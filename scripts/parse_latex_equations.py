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

_energy_equation, energy_equation = ccd_equations
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

def _create_einsum_call(string):
    if string == "+" or string == "-":
        return string

    print (string)
    divisor = re.findall(r"(?<=/)(\d)", string)

    assert len(divisor) < 2, \
            "The number of divisors is larger than 1 in expression: " \
            + "{0}".format(string)
    divisor = float(divisor[0]) if divisor else float(1)

    #matches = re.findall(r"(?<!\{)[a-zA-Z](?!\})", string)

    #for match in matches:
    #    print (match)

    #print (string)
    #indices = re.match(r".*\^\{(\w+)\}_\{(\w+)\}.*", string)
    #print (indices)
    #print (indices.group(1))
    #print (indices.group(2))
    print ("\n")

def _create_einsum_calls(string_list):
    return [_create_einsum_call(string) for string in string_list]


energy_equation = _remove_whitespace_and_empty_strings(
    re.split(r"([+-])", energy_equation))

energy_equation = _split_fractions(energy_equation)
print (energy_equation)
print ("\n")

print (_create_einsum_calls(energy_equation))
