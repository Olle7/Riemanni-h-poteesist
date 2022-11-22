from math import *

M=0
x_2=1
x_1=0.5
while True:
    M+=1
    print(((sin(x_2*log(M+1)+2*(2*pi)/4)*(1-x_1)+sin(x_2*log(M+1)+1*(2*pi)/4)*x_2)*M**(1-x_1)*((1-x_1)**2+x_2**2)**(-1))-\
          ((sin(x_2*log(M)+2*(2*pi)/4)*(1-x_1)+sin(x_2*log(M)+1*(2*pi)/4)*x_2)*M**(1-x_1)*((1-x_1)**2+x_2**2)**(-1))
          ," ; M=",M)