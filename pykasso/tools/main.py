"""
TODO
"""

import numpy as np

############
### GRID ###
############

pass

############
### MASK ###
############

pass

###############
### GEOLOGY ###
###############

pass

def read_gslib(file, nvar):
    with open(file) as f:
        print(f.readline(1))
    # data = np.genfromtxt(file, skip_header=3)
    # return GSLIB()
    return None

class GSLIB():
    pass

# lecteur de gslib 
# conversion en pickle_numpy