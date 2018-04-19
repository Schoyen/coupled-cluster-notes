import pickle
import os
import re
from generate_ccd_equations import get_ccd_equations
from sympy import Mul, Number
from sympy.physics.secondquant import (
        AntiSymmetricTensor, PermutationOperator
)

filename = "ccd_equations.p"

if not os.path.isfile(filename):
    ccd_equations = get_ccd_equations()

    with open(filename, "wb") as f:
        pickle.dump(ccd_equations, f)
else:
    with open(filename, "rb") as f:
        ccd_equations = pickle.load(f)

energy_equation, amplitude_equation = ccd_equations


occupied = "ij"
virtual = "ab"

swapaxes = {"a": 0, "b": 1, "i": 2, "j": 3}



def _get_tensor_string(term):
    tensor = str(term.symbol)

    if tensor == "t":
        return tensor

    tensor += "["

    for ind in term.upper + term.lower:
        if ind.assumptions0.get("above_fermi"):
            tensor += "n:, "
        elif ind.assumptions0.get("below_fermi"):
            tensor += ":n, "
        else:
            raise NotImplementedError("Index must be above/below Fermi level")

    return tensor[:-2] + "]"

def _get_indices(term):
    indices = term.upper + term.lower

    return "".join(map(lambda x: re.sub(r"_", "", str(x)), indices))

class Wrapper:

    def __init__(self, term):
        self.args = [term]

    def __str__(self):
        return str(self.args[0])

def _generate_term_function(term, reduces_to):
    factor = "1"
    einsum_string = \
            "term = {factor} * np.einsum(\"{indices}\", {tensors})"
    rhs = "term"
    indices = ""
    tensors = ""

    if not isinstance(term, Mul):
        term = Wrapper(term)

    if sum([isinstance(_, AntiSymmetricTensor) for _ in term.args]) == 1:
        tensor = list(filter(
            lambda x: isinstance(x, AntiSymmetricTensor), term.args))[0]
        _indices = list(map(lambda x: str(x), tensor.upper + tensor.lower))

        if len(_indices) == len(set(_indices)):
            return "rhs += " + _get_tensor_string(tensor) + "\n"
            #for _term in term.args:
            #    if _term.is_number:
            #        factor += " * ({0})".format(str(float(_term)))

            #return "rhs += {0} * {1}\n".format(
            #        factor, _get_tensor_string(tensor))

    for _term in term.args:
        if _term.is_number:
            factor += " * ({0})".format(str(float(_term)))
        elif isinstance(_term, AntiSymmetricTensor):
            tensors += "{0}, ".format(_get_tensor_string(_term))
            indices += _get_indices(_term) + ", "
        elif isinstance(_term, PermutationOperator):
            rhs = "({0} - {0}.swapaxes({1}, {2}))".format(
                    rhs, *[swapaxes[str(arg)] for arg in _term.args])
        else:
            raise NotImplementedError("The current _term is not implemented")

    tensors = tensors[:-2]
    indices = indices[:-2] + " -> {0}".format(reduces_to)
    einsum_string = einsum_string.format(
            factor=factor, indices=indices, tensors=tensors)
    rhs = "{0}\nrhs += {1}\n".format(einsum_string, rhs)

    return rhs

def generate_terms(equation, reduces_to=""):
    code = ""
    for arg in equation.args:
        code += _generate_term_function(arg, reduces_to)

    return code

print (generate_terms(energy_equation))
print ("\n")
print (generate_terms(amplitude_equation, reduces_to="abij"))
