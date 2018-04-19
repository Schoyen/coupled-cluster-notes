import os
from generate_ccd_equations import get_ccd_equations

filename = "ccd_equations.dat"

if not os.path.isfile(filename):
    ccd_equations = get_ccd_equations()

    with open(filename, "w") as f:
        for eq in ccd_equations:
            f.write(eq + "\n")
else:
    with open(filename, "r") as f:
        ccd_equations = f.read().split("\n")

print (ccd_equations)
