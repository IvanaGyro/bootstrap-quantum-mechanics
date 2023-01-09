from sage.rings.integer import Integer
from sage.rings.real_mpfr import RealField
# To enable the ability to find the 7th enery level at depth 60, we need at
# least 350 bits with using "hessenberg" algorithm to calculate the
# determinants.
# Here are some values of the required bits:
# - depth 60: 350 bits
# - depth 80: 450 bits
Real300 = RealField(Integer(500))
Real20 = RealField(Integer(40))
