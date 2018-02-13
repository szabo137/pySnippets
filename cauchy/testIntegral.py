import numpy as np
from scipy import integrate, special

# P \int_{-1}^1 dx 1/(x - wvar) * (1 + sin(x))
print(integrate.quad(lambda x: np.sin(x), -10, 10, weight='cauchy', wvar=0))

# Check against known result
print(2*special.sici(1)[0])
