from sympy import primerange
primes=list(primerange(2520,5000000))

print("primes")

from math import log,e


Euler_Mascheron=0.57721566490153286060651209008240243104215933593992

def not_RH_2(p):
    assert p>=2520
    return (3>=p*(e**Euler_Mascheron*2*log(log(p))-1))

for p in primes:
    if not_RH_2(p):
        print("RH not True")
    #print("dunno p:",p)

print("test2")
from sympy import prime
for n_p in range(369,348514):
    if not_RH_2(prime(n_p)):
        print("RH not True")