import pickle
import os
from generate_ccd_equations import get_ccd_equations

filename = "ccd_equations.p"

if not os.path.isfile(filename):
    ccd_equations = get_ccd_equations()

    with open(filename, "wb") as f:
        pickle.dump(ccd_equations, f)
else:
    with open(filename, "rb") as f:
        ccd_equations = pickle.load(f)

print (ccd_equations)
