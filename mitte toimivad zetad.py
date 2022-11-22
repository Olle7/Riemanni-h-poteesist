from math import *

def zeta6(*args,**kwargs):
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=300000

    if x_1<=0:
        raise NotImplementedError

    v_Re=(x_1**2-2*x_1+1+x_2**2)
    for n in range(2,iteratsioone+1):
        v_Re+=-sin(x_2*log(n)+1*(2*pi)/4)*(x_1**2-x_1-n*x_1+n+x_2**2)*n**-x_1

    v_Im=0#( sin(x_2*log(1)+1*(2*pi)/4)*x_2*(1-1) + sin(x_2*log(1)+2*(2*pi)/4)*(x_1**2-x_1-1*x_1+1+x_2**2) ) *1**(-x_1)
    for n in range(1,iteratsioone+1):
        v_Im+=( sin(x_2*log(n)+1*(2*pi)/4)*x_2*(n-1) + sin(x_2*log(n)+2*(2*pi)/4)*(x_1**2-x_1-n*x_1+n+x_2**2) ) *n**(-x_1)-\
        ( sin(x_2*log(n+1)+1*(2*pi)/4)*x_2 + sin(x_2*log(n+1)+2*(2*pi)/4)*(1-x_1) ) *n*(n+1)**(-x_1)

    return (x_1**2-2*x_1+1+x_2**2)**(-1)*(v_Re+1j*v_Im)

def zeta7(*args,**kwargs):
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=100

    if x_1<=0:
        raise NotImplementedError

    v_Re=(sin(x_2*log(1)+0*(2*pi)/4)*x_2*(1-1) + sin(x_2*log(1)+1*(2*pi)/4)*(x_1**2-x_1-1*x_1+1+x_2**2) ) *1**(-x_1)
    for n in range(1,iteratsioone+1):
        v_Re+=(sin(x_2*log(n+1)+0*(2*pi)/4)*x_2*(n+1-1) + sin(x_2*log(n+1)+1*(2*pi)/4)*(x_1**2-x_1-(n+1)*x_1+n+x_2**2) ) *(n+1)**(-x_1)-\
        ( sin(x_2*log(n+1)+0*(2*pi)/4)*x_2 + sin(x_2*log(n+1)+1*(2*pi)/4)*(1-x_1) ) *n*(n+1)**(-x_1)

    v_Im=0
    for n in range(1,iteratsioone+1):
        v_Im+=( sin(x_2*log(n)+1*(2*pi)/4)*x_2*(n-1) + sin(x_2*log(n)+2*(2*pi)/4)*(x_1**2-x_1-n*x_1+n+x_2**2) ) *n**(-x_1)-\
        ( sin(x_2*log(n+1)+1*(2*pi)/4)*x_2 + sin(x_2*log(n+1)+2*(2*pi)/4)*(1-x_1) ) *n*(n+1)**(-x_1)
    return (x_1**2-2*x_1+1+x_2**2)**(-1)*(v_Re+1j*v_Im)