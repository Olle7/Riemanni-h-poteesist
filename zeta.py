from math import sin,cos,log,pi,e,exp,floor,atan,factorial
def sgn(x):
    if x>0:
        return 1
    if x<0:
        return -1
    if x==0:
        return 0
from time import time
def zeta1_Dirichlet_series_Wikipediast(z, iteratsioone=10000):
    #avaldise_latex=r"$\forall_{x_1}(x_1 \in C \land Re(x_1)>0 \to \zeta(x_1)=(x_1-1)^{-1}*\sum_{n=1}^{\infty}(\frac{n}{(n+1)^{x_1}}-\frac{n-x_1}{n^{x_1}})$"
    if z.real<=0:
        raise NotImplementedError
    s=0
    for n in range(1,iteratsioone+1):
        s+=n/(n+1)**z-(n-z)/n**z
    return s*(z-1)**-1

#print(zeta1(3+5.5j))
default_iteratsioone=10**6*5



def zeta2_reaalosa_ja_imaginaarosa_eraldatud(*args, **kwargs):
#    latex=r"""$\forall_{x_1}(\forall_{x_2}(x_1 \in R \land x_2 \in R \land x_1>0 \to \\
#\zeta(x_1+i*x_2)=\\
#   ((x_1-1)^2+x_2^2)^{-1}*\sum_{n=1}^{\infty}(\\
#((x_1-1)*cos(x_2*ln(n+1))-x_2*sin(x_2*ln(n+1)))*n*(n+1)^{-x_1}-\\
#((n-x_1)*((x_1-1)*cos(x_2*ln(n))-x_2*sin(x_2*ln(n)))-x_2*((x_1-1)*sin(x_2*ln(n))+x_2*cos(x_2*ln(n))))*n^{-x_1}\\
#)\\
#-i*((x_1-1)^2+x_2^2)^{-1}*\sum_{n=1}^{\infty}(\\
#((x_1-1)*sin(x_2*ln(n+1))+x_2*cos(x_2*ln(n+1)))*n*(n+1)^{-x_1}-\\
#((n-x_1)*((x_1-1)*sin(x_2*ln(n))+x_2*cos(x_2*ln(n)))+x_2*((x_1-1)*cos(x_2*ln(n))-x_2*sin(x_2*ln(n))))*n^{-x_1}\\
#)=\\
#))$"""

    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=default_iteratsioone


    if x_1<=0:
        raise NotImplementedError

    v_Re=0
    for n in range(1,iteratsioone+1):
        v_Re+=((x_1-1)*cos(x_2*log(n+1))-x_2*sin(x_2*log(n+1)))*n*(n+1)**(-x_1)-\
        ((n-x_1)*((x_1-1)*cos(x_2*log(n))-x_2*sin(x_2*log(n)))-x_2*((x_1-1)*sin(x_2*log(n))+x_2*cos(x_2*log(n))))*n**(-x_1)
    v_Re*=((x_1-1)**2+x_2**2)**-1

    v_Im=0
    for n in range(1, iteratsioone+1):
        v_Im+=((x_1-1)*sin(x_2*log(n+1))+x_2*cos(x_2*log(n+1)))*n*(n+1)**(-x_1)-\
        ((n-x_1)*((x_1-1)*sin(x_2*log(n))+x_2*cos(x_2*log(n)))+x_2*((x_1-1)*cos(x_2*log(n))-x_2*sin(x_2*log(n))))*n**(-x_1)
    v_Im*=-((x_1-1)**2+x_2**2)**-1
    #return (v_Re,v_Im)
    return v_Re+1j*v_Im

#t??psus: zeta2,zeta3,zeta1
#print(zeta2(3,5.5))

def zeta3(*args,**kwargs):
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=default_iteratsioone

    if x_1<=0:
        raise NotImplementedError

    v_Re=0
    for n in range(1,iteratsioone+1):
        v_Re+=((x_1-1)*cos(x_2*log(n+1))-x_2*sin(x_2*log(n+1)))*n*(n+1)**(-x_1)-\
        ( (n*x_1-x_1**2+x_1-x_2**2-n)*cos(x_2*log(n)) - x_2*(n-1)*sin(x_2*log(n)) )*n**(-x_1)
    v_Re*=((x_1-1)**2+x_2**2)**-1

    v_Im=0
    for n in range(1, iteratsioone+1):
        v_Im+=((x_1-1)*sin(x_2*log(n+1))+x_2*cos(x_2*log(n+1)))*n*(n+1)**(-x_1)-\
        ( (n*x_1-x_1**2+x_1-x_2**2-n)*sin(x_2*log(n)) + x_2*(n-1)*cos(x_2*log(n)) )*n**(-x_1)
    v_Im*=-((x_1-1)**2+x_2**2)**-1
    #return (v_Re,v_Im)
    return v_Re+1j*v_Im

def zeta4(*args,**kwargs):
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=1000

    if x_1<=0:
        raise NotImplementedError

    v_Re=0
    for n in range(1,iteratsioone+1):
        v_Re+=( (x_1**2-x_1-n*x_1+n+x_2**2)*cos(x_2*log(n)) + x_2*(n-1)*sin(x_2*log(n)) )*n**(-x_1)+\
        (-(1-x_1)*cos(x_2*log(n+1))-x_2*sin(x_2*log(n+1)))*n*(n+1)**(-x_1)

    v_Im=0
    for n in range(1,iteratsioone+1):
        v_Im+=(-(x_1**2-x_1-n*x_1+n+x_2**2)*sin(x_2*log(n))+x_2*(n-1)*cos(x_2*log(n)))*n**(-x_1)+\
        ((1-x_1)*sin(x_2*log(n+1))-x_2*cos(x_2*log(n+1)))*n*(n+1)**(-x_1)
    return (x_1**2-2*x_1+1+x_2**2)**(-1)*(v_Re+1j*v_Im)

def zeta5(*args,**kwargs):
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=default_iteratsioone

    if x_1<=0:
        raise NotImplementedError

    v_Re=0
    v_Im=0
    for n in range(1,iteratsioone+1):
        v_Re+=(sin(x_2*log(n)+0*(2*pi)/4)*x_2*(n-1) + sin(x_2*log(n)+1*(2*pi)/4)*(x_1**2-x_1-n*x_1+n+x_2**2) ) *n**(-x_1)-\
        ( sin(x_2*log(n+1)+0*(2*pi)/4)*x_2 + sin(x_2*log(n+1)+1*(2*pi)/4)*(1-x_1) ) *n*(n+1)**(-x_1)

        v_Im+=( sin(x_2*log(n)+1*(2*pi)/4)*x_2*(n-1) + sin(x_2*log(n)+2*(2*pi)/4)*(x_1**2-x_1-n*x_1+n+x_2**2) ) *n**(-x_1)-\
        ( sin(x_2*log(n+1)+1*(2*pi)/4)*x_2 + sin(x_2*log(n+1)+2*(2*pi)/4)*(1-x_1) ) *n*(n+1)**(-x_1)

    return(x_1**2-2*x_1+1+x_2**2)**(-1)*(v_Re+1j*v_Im)

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
        iteratsioone=default_iteratsioone

    if x_1<=0:
        raise NotImplementedError

    v_Re=0
    v_Im=0
    for n in range(1,iteratsioone+1):
        v_Re+=sin(x_2*log(n)+1*(2*pi)/4)*(x_1**2-x_1-n*x_1+n+x_2**2)*n**(-x_1)\
             -sin(x_2*log(n+1)+1*(2*pi)/4)*(1-x_1)*n*(n+1)**(-x_1)

        v_Im+=sin(x_2*log(n)+2*(2*pi)/4)*(x_1**2-x_1-n*x_1+n+x_2**2)*n**(-x_1)\
             -sin(x_2*log(n+1)+2*(2*pi)/4)*(1-x_1)*n*(n+1)**(-x_1)

    v_Re-=sin(x_2*log(n+1)+0*(2*pi)/4)*x_2*n*(n+1)**(-x_1)
    v_Im-=sin(x_2*log(n+1)+1*(2*pi)/4)*x_2*n*(n+1)**(-x_1)
    return (v_Re+1j*v_Im)*(x_1**2-2*x_1+1+x_2**2)**(-1)

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
        iteratsioone=default_iteratsioone

    if x_1<=0:
        raise NotImplementedError

    v_Re=0
    v_Im=0
    for n in range(1,iteratsioone+1):
        v_Re+=sin(x_2*log(n)+1*(2*pi)/4)*n**(-x_1)
        v_Im+=sin(x_2*log(n)+2*(2*pi)/4)*n**(-x_1)

    v_Re-=(sin(x_2*log(n+1)+1*(2*pi)/4)*(1-x_1)+sin(x_2*log(n+1)+0*(2*pi)/4)*x_2)*n*(n+1)**(-x_1)*(x_1**2-2*x_1+x_2**2+1)**(-1)
    v_Im-=(sin(x_2*log(n+1)+2*(2*pi)/4)*(1-x_1)+sin(x_2*log(n+1)+1*(2*pi)/4)*x_2)*n*(n+1)**(-x_1)*(x_1**2-2*x_1+x_2**2+1)**(-1)
    return (v_Re+1j*v_Im)

def zeta8(*args,**kwargs):
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=default_iteratsioone

    if x_1<=0:
        raise NotImplementedError

    v_Re=0
    v_Im=0
    for n in range(1,iteratsioone+1):
        v_Re+=sin(x_2*log(n)+1*(2*pi)/4)*n**(-x_1)
        v_Im+=sin(x_2*log(n)+2*(2*pi)/4)*n**(-x_1)

    v_Re-=(sin(x_2*log(n+1)+1*(2*pi)/4)*(1-x_1)+sin(x_2*log(n+1)+0*(2*pi)/4)*x_2)*n**(1-x_1)*((1-x_1)**2+x_2**2)**(-1)
    v_Im-=(sin(x_2*log(n+1)+2*(2*pi)/4)*(1-x_1)+sin(x_2*log(n+1)+1*(2*pi)/4)*x_2)*n**(1-x_1)*((1-x_1)**2+x_2**2)**(-1)
    return (v_Re+1j*v_Im)

def zeta9_M_plus_1e_ei_ole(*args, **kwargs):
#    latex=r"""\lim_{M \to \infty}(\\
#\sum_{n=1}^{M}(sin(x_2*ln(n)+1*(2\pi)/4)*n^{-x_1})\\
#-(sin(x_2*ln(M)+1*(2\pi)/4)*(1-x_1)+sin(x_2*ln(M)+0*(2\pi)/4)*x_2)*M^{1-x_1}*((1-x_1)^2+x_2^2)^{-1}\\
#)\\
#+i*\lim_{M \to \infty}(\\
#\sum_{n=1}^{M}(sin(x_2*ln(n)+2*(2\pi)/4)*n^{-x_1})\\
#-(sin(x_2*ln(M)+2*(2\pi)/4)*(1-x_1)+sin(x_2*ln(M)+1*(2\pi)/4)*x_2)*M^{1-x_1}*((1-x_1)^2+x_2^2)^{-1}\\
#)=\\"""
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=default_iteratsioone

    if x_1<=0:
        raise NotImplementedError

    v_Re=0
    v_Im=0
    for n in range(1,iteratsioone+1):
        v_Re+=sin(x_2*log(n)+1*(2*pi)/4)*n**(-x_1)
        v_Im+=sin(x_2*log(n)+2*(2*pi)/4)*n**(-x_1)

    v_Re-=(sin(x_2*log(n)+1*(2*pi)/4)*(1-x_1)+sin(x_2*log(n)+0*(2*pi)/4)*x_2)*n**(1-x_1)*((1-x_1)**2+x_2**2)**(-1)
    v_Im-=(sin(x_2*log(n)+2*(2*pi)/4)*(1-x_1)+sin(x_2*log(n)+1*(2*pi)/4)*x_2)*n**(1-x_1)*((1-x_1)**2+x_2**2)**(-1)
    return (v_Re+1j*v_Im)

def zeta10(*args,**kwargs):
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if x_1<=0:
        raise NotImplementedError
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=default_iteratsioone
    M=log(iteratsioone)*abs(x_2)/(2*pi)
    iteratsioone=int(e**(int(M)*2*pi/abs(x_2)))
    print("M:",M,"; iteratsioone:",iteratsioone)

    v_Re=0
    v_Im=0
    for n in range(1,iteratsioone):
        v_Re+=sin(x_2*log(n)+1*(2*pi)/4)*n**(-x_1)
        v_Im+=sin(x_2*log(n)+2*(2*pi)/4)*n**(-x_1)
    v_Re-=e**(int(M)*2*pi*(1-x_1)/abs(x_2))*((1-x_1)**2+x_2**2)**(-1)*(1-x_1)
    v_Im-=e**(int(M)*2*pi*(1-x_1)/abs(x_2))*((1-x_1)**2+x_2**2)**(-1)*x_2

    return (v_Re+1j*v_Im)

def zeta11(*args,**kwargs):
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if x_1<=0:
        raise NotImplementedError
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=10**7
    M=log(iteratsioone)*abs(x_2)/(2*pi)
    iteratsioone=int(exp((int(M)*2*pi)/abs(x_2)))
    print("M:",M,"; iteratsioone:",iteratsioone)
    pi2=6.2831853071795864769252867665590057683943387987502116419498891846

    v_Re=0
    v_Im=0
    for n in range(1,iteratsioone):
        n_x1=n**x_1
        logx2=x_2*log(n)
        v_Re+= cos(logx2)/n_x1
        v_Im+=-sin(logx2)/n_x1

    v_Re-=exp(int(M)*pi2*(1-x_1)/abs(x_2))*((1-x_1)**2+x_2**2)**(-1)*(1-x_1)
    v_Im-=exp(int(M)*pi2*(1-x_1)/abs(x_2))*((1-x_1)**2+x_2**2)**(-1)*x_2

    return (v_Re+1j*v_Im)

def zeta12_q1_v??rdub_q2(*args, **kwargs):

    latex=r"""$\forall_{x_1}(\forall_{x_2}(x_1 \in R \land x_2 \in R \land x_1>0 \to \\
\zeta(x_1+i*x_2)=\\

\lim_{M \to \infty}(\\
\sum_{n=1}^{floor(e^{floor(M)*2\pi/abs(x_2)})}(sin(x_2*ln(n)+1*(2\pi)/4)*n^{-x_1})\\
-e^{floor(M)*2\pi*(1-x_1)/abs(x_2)}*((1-x_1)^2+x_2^2)^{-1}*(1-x_1)\\
)\\
+i*\lim_{M \to \infty}(\\
\sum_{n=1}^{floor(e^{floor(M)*2\pi/abs(x_2)})}(sin(x_2*ln(n)+2*(2\pi)/4)*n^{-x_1})\\
-e^{floor(M)*2\pi*(1-x_1)/abs(x_2)}*((1-x_1)^2+x_2^2)^{-1}*x_2\\
)=\\
))
$"""
    q=0.4
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if x_1<=0:
        raise NotImplementedError
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=default_iteratsioone

    M=(log(iteratsioone)/(2*pi)-q/x_2)*abs(x_2)

    iteratsioone=floor(e**((floor(M)/abs(x_2)+q/x_2)*2*pi))
    print("M:",M,"; iteratsioone:",iteratsioone)

    v_Re=0
    v_Im=0
    for n in range(1,iteratsioone):
        n_x1=n**x_1
        logx2=x_2*log(n)
        v_Re+=cos(logx2)/n_x1
        v_Im-=sin(logx2)/n_x1
    #for n in range(1,iteratsioone):
    #    v_Re+=sin(x_2*log(n)+1*(2*pi)/4)*n**(-x_1)
    #    v_Im+=sin(x_2*log(n)+2*(2*pi)/4)*n**(-x_1)
    v_Re-=(sin((q+1/4)*2*pi)*(1-x_1)+sin((q+0/4)*2*pi)*x_2)*e**((floor(M)/abs(x_2)+q/x_2)*2*pi*(1-x_1))*((1-x_1)**2+x_2**2)**(-1)
    v_Im-=(sin((q+2/4)*2*pi)*(1-x_1)+sin((q+1/4)*2*pi)*x_2)*e**((floor(M)/abs(x_2)+q/x_2)*2*pi*(1-x_1))*((1-x_1)**2+x_2**2)**(-1)

    return(v_Re+1j*v_Im)

def zeta12_1(*args, **kwargs):
    q=0.2451
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if x_1<=0:
        raise NotImplementedError
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=default_iteratsioone
    M=(log(iteratsioone)*x_2)/(2*pi)
    iteratsioone=floor(e**((floor(M))*2*pi/x_2))
    print("M:",M,"; iteratsioone:",iteratsioone)

    v_Re=0
    v_Im=0
    #for n in range(1,iteratsioone):
    #    n_x1=n**x_1
    #    logx2=x_2*log(n)
    #    v_Re+=cos(logx2)/n_x1
    #    v_Im-=sin(logx2)/n_x1
    for n in range(1,iteratsioone):
        v_Re+=sin(x_2*log(n)+(1+q)*(2*pi)/4)*n**(-x_1)
        v_Im+=sin(x_2*log(n)+(2+q)*(2*pi)/4)*n**(-x_1)
    v_Re-=(sin((q+1/4)*2*pi)*(1-x_1)+sin((q+0/4)*2*pi)*x_2)*e**((floor(M)+q)*2*pi*(1-x_1)/x_2)*((1-x_1)**2+x_2**2)**(-1)
    v_Im-=(sin((q+2/4)*2*pi)*(1-x_1)+sin((q+1/4)*2*pi)*x_2)*e**((floor(M)+q)*2*pi*(1-x_1)/x_2)*((1-x_1)**2+x_2**2)**(-1)

    return(v_Re+1j*v_Im)


def zeta13_tagumine_pool_nulliks(*args, **kwargs):

    latex=r"""$\forall_{x_1}(\forall_{x_2}(x_1 \in R \land x_2 \in R \land x_1>0 \to \\
\zeta(x_1+i*x_2)=\\
   \lim_{M \to \infty}( \sum_{n=1}^{floor(e^{(      floor(M)/2*2\pi/abs(x_2)-arctan((1-x_1)/x_2)/x_2)})}(sin(x_2*ln(n)+2\pi/4)*n^{-x_1}))\\
+i*\lim_{M \to \infty}(-\sum_{n=1}^{floor(e^{((floor(M)/2+1/4)*2\pi/abs(x_2)-arctan((1-x_1)/x_2)/x_2)})}(sin(x_2*ln(n))*n^{-x_1}))=\\
))
$"""

    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if x_1<=0:
        raise NotImplementedError
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=default_iteratsioone

    #M = (log(iteratsioone) * abs(x_2)) / (2 * pi) * 2 - max(q_1, q_2) / pi#vana

    M=max(
    (log(iteratsioone)*abs(x_2)+atan((1-x_1)/x_2)*sgn(x_2))*2/(2*pi),
    (log(iteratsioone)+atan((1-x_1)/x_2)/x_2)*2*abs(x_2)/(2*pi)-1/2
    )
    iteratsioone_1=floor(e**(floor(M)/2*2*pi/abs(x_2)-atan((1-x_1)/x_2)/x_2))
    iteratsioone_2=floor(e**((floor(M)/2+1/4)*2*pi/abs(x_2)-atan((1-x_1)/x_2)/x_2))

    print("M:",M,"; iteratsioone_1:",iteratsioone_1,"; iteratsioone_2:",iteratsioone_2)

    v_Re=0
    for n in range(1,iteratsioone_1+1):
        v_Re+=sin(x_2*log(n)+1*(2*pi)/4)*n**(-x_1)
    v_Im=0
    for n in range(1,iteratsioone_2+1):
        v_Im+=sin(x_2*log(n))*n**(-x_1)

    return(v_Re-1j*v_Im)

def zeta13_tagumine_pool_nulliks_optimeeritud(*args, **kwargs):

    latex=r"""\lim_{M \to \infty}( \sum_{n=1}^{floor(e^{(      floor(M)/2*2\pi-arctan((1-x_1)/x_2))/x_2})}(sin(x_2*ln(n)+2\pi/4)*n^{-x_1}))\\
+i* \lim_{M \to \infty}(-\sum_{n=1}^{floor(e^{((floor(M)/2+1/4)*2\pi-arctan((1-x_1)/x_2))/x_2})}(sin(x_2*ln(n))*n^{-x_1}))"""

    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if x_1<=0:
        raise NotImplementedError
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=default_iteratsioone

    q_1=-atan((1-x_1)/x_2)
    #q_2=+atan(x_2/(1-x_1))
    q_2=q_1+1/4*sgn((1-x_1)*x_2)*2*pi
    print("q_1-q_2=",q_1-q_2)

    M=(log(iteratsioone)*x_2)/(2*pi)*2-max(q_1,q_2)/pi

    iteratsioone_1=floor(e**((floor(M)*pi+q_1)/x_2))
    iteratsioone_2=floor(e**((floor(M)*pi+q_2)/x_2))

    print("M:",M,"; iteratsioone_1:",iteratsioone_1,"; iteratsioone_2:",iteratsioone_2)

    v_Re=0
    #for n in range(1,iteratsioone_1+1):
    #    v_Re+=sin(x_2*log(n)+1*(2*pi)/4)*n**(-x_1)
    v_Im=0
    #for n in range(1,iteratsioone_2+1):
    #    v_Im+=sin(x_2*log(n))*n**(-x_1)
    mini=min(iteratsioone_1,iteratsioone_2)
    for n in range(1,mini+1):
        n_x1=n**x_1
        logx2=x_2*log(n)
        v_Re+=cos(logx2)/n_x1
        v_Im-=sin(logx2)/n_x1
    if iteratsioone_1>iteratsioone_2:
        for n in range(mini+1,iteratsioone_1+1):
            v_Re+=cos(x_2*log(n))*n**-x_1
    else:
        for n in range(mini+1,iteratsioone_2+1):
            v_Im-=sin(x_2*log(n))*n**-x_1

    return(v_Re+1j*v_Im)

def zeta14_tagumine_pool_nulliks_tayloriga(*args, **kwargs):
    latex=r"""$\forall_{x_1}(\forall_{x_2}(x_1 \in R \land x_2 \in R \land x_1>0 \to \\
\zeta(x_1+i*x_2)=\\
\lim_{M \to \infty}( \sum_{n_0=0}^{floor(e^{floor(M)/2*2\pi/abs(x_2)-arctan((1-x_1)/x_2)/x_2})}(\sum_{n_1=0}^{(2*M)^M-1}(\\
(\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))\%2)*\\
(x_2*2)^{\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))}*\\
(-1)^{\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))/2}*\\
n_0^{\sum_{n_2=0}^{M-1}((2*n_2+1)*(floor(n_1/(2*M)^{n_2})\%(2*M)))}*\\
(n_0+1)^{-x_1}\\
(n_0+2)^{-\sum_{n_2=0}^{M-1}((2*n_2+1)*(floor(n_1/(2*M)^{n_2})\%(2*M)))}*\\
\Pi_{n_2=0}^{M-1}(factorial(floor(n_1/(2*M)^{n_2})\%(2*M))^{-1}*(2*n_2+1)^{-floor(n_1/(2*M)^{n_2})\%(2*M)} )\\
)))\\

+i*\lim_{M \to \infty}(\sum_{n_0=0}^{floor(e^{(floor(M)/2+1/4)*2\pi/abs(x_2)-arctan((1-x_1)/x_2)/x_2})}(\sum_{n_1=0}^{(2*M)^M-1}(\\
(\sum_{n_2=0}^{M-1}(1+floor(n_1/(2*M)^{n_2})\%(2*M))\%2)*\\
(x_2*2)^{\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))}*\\
(-1)^{(\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))+1)/2}*\\
n_0^{ \sum_{n_2=0}^{M-1}((2*n_2+1)*(floor(n_1/(2*M)^{n_2})\%(2*M)))}*\\
(n_0+1)^{-x_1}\\
(n_0+2)^{-\sum_{n_2=0}^{M-1}((2*n_2+1)*(floor(n_1/(2*M)^{n_2})\%(2*M)))}*\\
\Pi_{n_2=0}^{M-1}(factorial(floor(n_1/(2*M)^{n_2})\%(2*M))^{-1}*(2*n_2+1)^{-floor(n_1/(2*M)^{n_2})\%(2*M)} )\\
)))\\
))%forallide- lopusulud
$"""
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if x_1<=0:
        raise NotImplementedError
    if "iteratsioone" in kwargs:
        soovitud_iteratsioone=kwargs["iteratsioone"]
    else:
        soovitud_iteratsioone=20
        M=max(
    (log(soovitud_iteratsioone)*abs(x_2)+atan((1-x_1)/x_2)*sgn(x_2))*2/(2*pi),
    (log(soovitud_iteratsioone)+atan((1-x_1)/x_2)/x_2)*2*abs(x_2)/(2*pi)-1/2
    )
    iteratsioone_1=floor(e**(floor(M)/2*2*pi/abs(x_2)-atan((1-x_1)/x_2)/x_2))
    iteratsioone_2=floor(e**((floor(M)/2+1/4)*2*pi/abs(x_2)-atan((1-x_1)/x_2)/x_2))
    print("iteratsioone1=",iteratsioone_1)
    print("iteratsioone2=",iteratsioone_2)
    s0_real=0
    assert M>0
    M=int(M)
    print(M)
    for n0 in range(iteratsioone_1):
        print("n0=",n0)
        for n1 in range((2*M)**M):
            #if n1%100000==0:
            #    print(n1)
            p=1
            s2_0=0
            for n2 in range(M):
                s2_0+=(n1//(2*M)**n2)%(2*M)
            p*=(s2_0+1)%2
            if p==0:
                continue
            p*=(x_2*2)**s2_0
            p*=(-1)**(s2_0/2)
            s2_1=0
            for n2 in range(M):
                s2_1+=(2*n2+1)*(n1//(2*M)**n2)%(2*M)
            p*=n0**s2_1
            p*=(n0+1)**(-x_1)
            p*=(n0+2)**-s2_1
            for n2 in range(M):
                p*=factorial((n1//(2*M)**n2)%(2*M))**-1*(2*n2+1)**(-s2_0)
            s0_real+=p
    print("real=",s0_real)
    s0_imag=0
    for n0 in range(iteratsioone_2):
        print("n0=", n0)
        for n1 in range((2*M)**M):
            p=1
            s2_0=0
            for n2 in range(M):
                s2_0+=(n1//(2*M)**n2)%(2*M)
            p*=s2_0%2
            if p==0:
                continue
            p*=(x_2*2)**s2_0
            p*=(-1)**((s2_0+1)/2)
            s2_1=0
            for n2 in range(M):
                s2_1+=(2*n2+1)*(n1//(2*M)**n2)%(2*M)
            p*=n0**s2_1
            p*=(n0+1)**-x_1
            p*=(n0+2)**-s2_1
            for n2 in range(M):
                p*=factorial((n1//(2*M)**n2)%(2*M))**-1*(2*n2+1)**(-s2_0)
            s0_imag+=p
    return s0_real+1j*s0_imag

def zeta14_tagumine_pool_nulliks_tayloriga_optimeeritud(*args, **kwargs):
    latex=r"""$\forall_{x_1}(\forall_{x_2}(x_1 \in R \land x_2 \in R \land x_1>0 \to \\
\zeta(x_1+i*x_2)=\\
\lim_{M \to \infty}( \sum_{n_0=0}^{floor(e^{floor(M)/2*2\pi/abs(x_2)-arctan((1-x_1)/x_2)/x_2})}(\sum_{n_1=0}^{(2*M)^M-1}(\\
(\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))\%2)*\\
(x_2*2)^{\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))}*\\
(-1)^{\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))/2}*\\
n_0^{\sum_{n_2=0}^{M-1}((2*n_2+1)*(floor(n_1/(2*M)^{n_2})\%(2*M)))}*\\
(n_0+1)^{-x_1}\\
(n_0+2)^{-\sum_{n_2=0}^{M-1}((2*n_2+1)*(floor(n_1/(2*M)^{n_2})\%(2*M)))}*\\
\Pi_{n_2=0}^{M-1}(factorial(floor(n_1/(2*M)^{n_2})\%(2*M))^{-1}*(2*n_2+1)^{-floor(n_1/(2*M)^{n_2})\%(2*M)} )\\
)))\\

+i*\lim_{M \to \infty}(\sum_{n_0=0}^{floor(e^{(floor(M)/2+1/4)*2\pi/abs(x_2)-arctan((1-x_1)/x_2)/x_2})}(\sum_{n_1=0}^{(2*M)^M-1}(\\
(\sum_{n_2=0}^{M-1}(1+floor(n_1/(2*M)^{n_2})\%(2*M))\%2)*\\
(x_2*2)^{\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))}*\\
(-1)^{(\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))+1)/2}*\\
n_0^{ \sum_{n_2=0}^{M-1}((2*n_2+1)*(floor(n_1/(2*M)^{n_2})\%(2*M)))}*\\
(n_0+1)^{-x_1}\\
(n_0+2)^{-\sum_{n_2=0}^{M-1}((2*n_2+1)*(floor(n_1/(2*M)^{n_2})\%(2*M)))}*\\
\Pi_{n_2=0}^{M-1}(factorial(floor(n_1/(2*M)^{n_2})\%(2*M))^{-1}*(2*n_2+1)^{-floor(n_1/(2*M)^{n_2})\%(2*M)} )\\
)))\\
))%forallide- lopusulud
$"""
    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if x_1<=0:
        raise NotImplementedError
    if "iteratsioone" in kwargs:
        soovitud_iteratsioone=kwargs["iteratsioone"]
    else:
        soovitud_iteratsioone=25
        M=max(
    (log(soovitud_iteratsioone)*abs(x_2)+atan((1-x_1)/x_2)*sgn(x_2))*2/(2*pi),
    (log(soovitud_iteratsioone)+atan((1-x_1)/x_2)/x_2)*2*abs(x_2)/(2*pi)-1/2
    )

    print("M:",M)
    assert M > 0
    M = floor(M)

    iteratsioone_1=floor(e**(floor(M)/2*2*pi/abs(x_2)-atan((1-x_1)/x_2)/x_2))
    iteratsioone_2=floor(e**((floor(M)/2+1/4)*2*pi/abs(x_2)-atan((1-x_1)/x_2)/x_2))
    print("iteratsioone1=",iteratsioone_1)
    print("iteratsioone2=",iteratsioone_2)
    s0_real=0
    M2=2*M
    #print("n0=",n0)
    for n1 in range((2*M)**M):
        #if n1%100000==0:
        #    print(n1)
        p=1
        s2_0=0
        s2_1=0
        for n2 in range(M):
            k=(n1//M2**n2)%M2
            s2_0+=k
            s2_1+=(2*n2+1)*k
            p/=factorial(k)*(2*n2+1)**s2_0
        p*=(s2_0+1)%2
        if p==0:
            continue
        p*=(x_2*2)**s2_0
        p*=(-1)**(s2_0/2)
        sx=0
        for n0 in range(1,iteratsioone_1):
            sx+=(n0+1)**-x_1
            sx+=(n0/(n0+2))**s2_1
        p*=sx
        s0_real+=p
    print("real=",s0_real)
        
    s0_imag=0
    for n1 in range((2*M)**M):
        p=1
        s2_0=0
        s2_1=0
        for n2 in range(M):
            k=(n1//(2*M)**n2)%(2*M)
            s2_0+=k
            s2_1+=(2*n2+1)*k
            p*=factorial(k)**-1*(2*n2+1)**(-s2_0)
        p*=s2_0%2
        if p==0:
            continue
        p*=(x_2*2)**s2_0
        p*=(-1)**((s2_0+1)/2)
        sx=0
        for n0 in range(iteratsioone_2+1):
            sx+=(n0+1)**-x_1
            sx+=(n0/(n0+2))**s2_1
        p*=sx
        s0_imag+=p
    return s0_real+1j*s0_imag
def cos_y_ln_x_multinomial_5(x,y,M=6):
    latex=r"""$\forall_x(x>0 \to cos(y*ln(x))=\\
lim_{M \to \infty}(\sum_{n_1=0}^{(2*M)^M-1}(\\
(\sum_{n_2=0}^{M-1}(1+floor(n_1/(2*M)^{n_2})\%(2*M))\%2)*\\
(y*2)^{\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))}*\\
(-1)^{\sum_{n_2=0}^{M-1}(floor(n_1/(2*M)^{n_2})\%(2*M))/2}*\\
(x-1)^{ \sum_{n_2=0}^{M-1}((2*n_2+1)*(floor(n_1/(2*M)^{n_2})\%(2*M)))}*\\
(x+1)^{-\sum_{n_2=0}^{M-1}((2*n_2+1)*(floor(n_1/(2*M)^{n_2})\%(2*M)))}*\\
\Pi_{n_2=0}^{M-1}(factorial(floor(n_1/(2*M)^{n_2})\%(2*M))^{-1}*(2*n_2+1)^{-floor(n_1/(2*M)^{n_2})\%(2*M)} )\\
))\\%k2-summa ja lim'i lopusulg sel real.
)%forall lopusulg$"""
    s1=0
    for k1 in range((2*M)**M):
        p=1
        s2=0
        for k2 in range(M):
            s2+=(k1//(2*M)**k2)%(2*M)
        p*=(s2+1)%2
        #p*=s2<2*M
        if p==0:
            continue
        p*=(y*2)**s2
        #print("s2:",s2)
        p*=(-1)**(s2/2)


        s2=0
        for k2 in range(M):
            s2+=(2*k2+1)* ((k1//(2*M)**k2)%(2*M))
        p*=((x-1)/(x+1))**s2

        for k2 in range(M):
            p*=factorial((k1//(2*M)**k2)%(2*M))**(-1)*(2*k2+1)**(-((k1//(2*M)**k2)%(2*M)))
        s1+=p
    return s1
def sin_y_ln_x_multinomial_5(x,y,M=5):
    latex=r"""$\forall_x(x>0 \to sin(y*ln(x))=\\
lim_{M \to \infty}(\sum_{k_1=0}^{(2*M)^M-1}(\\
(\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))\%2)*\\
(y*2)^{\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))}*\\
(-1)^{(\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))+3)/2}*\\
(x-1)^{ \sum_{k_2=0}^{M-1}((2*k_2+1)*(floor(k_1/(2*M)^{k_2})\%(2*M)))}*\\
(x+1)^{-\sum_{k_2=0}^{M-1}((2*k_2+1)*(floor(k_1/(2*M)^{k_2})\%(2*M)))}*\\
\Pi_{k_2=0}^{M-1}(factorial(floor(k_1/(2*M)^{k_2})\%(2*M))^{-1}*(2*k_2+1)^{-floor(k_1/(2*M)^{k_2})\%(2*M)} )\\
))\\%k2-summa ja lim'i lopusulg sel real.
)%forall lopusulg$"""
    s1=0
    for k1 in range((2*M)**M):
        p=1
        s2=0
        for k2 in range(M):
            s2+=(k1//(2*M)**k2)%(2*M)
        p*=s2%2
        #p*=s2<2*M
        if p==0:
            continue
        p*=(y*2)**s2
        #print("s2:",s2)
        p*=(-1)**((s2+3)/2)

        s2=0
        for k2 in range(M):
            s2+=(2*k2+1)* ((k1//(2*M)**k2)%(2*M))
        p*=((x-1)/(x+1))**s2

        for k2 in range(M):
            p/=factorial((k1//(2*M)**k2)%(2*M))*(2*k2+1)**((k1//(2*M)**k2)%(2*M))
        s1+=p
    return s1


def zeta15(*args, **kwargs):

    if len(args)==1:
        x_1=args[0].real
        x_2=args[0].imag
    elif len(args)==2:
        x_1=args[0]
        x_2=args[1]
    if x_1<=0:
        raise NotImplementedError
    if "iteratsioone" in kwargs:
        iteratsioone=kwargs["iteratsioone"]
    else:
        iteratsioone=43

    #M = (log(iteratsioone) * abs(x_2)) / (2 * pi) * 2 - max(q_1, q_2) / pi#vana

    M=max(
    (log(iteratsioone)*abs(x_2)+atan((1-x_1)/x_2)*sgn(x_2))*2/(2*pi),
    (log(iteratsioone)+atan((1-x_1)/x_2)/x_2)*2*abs(x_2)/(2*pi)-1/2
    )
    iteratsioone_1=floor(e**(floor(M)/2*2*pi/abs(x_2)-atan((1-x_1)/x_2)/x_2))
    iteratsioone_2=floor(e**((floor(M)/2+1/4)*2*pi/abs(x_2)-atan((1-x_1)/x_2)/x_2))

    print("M:",M,"; iteratsioone_1:",iteratsioone_1,"; iteratsioone_2:",iteratsioone_2)

    v_Re=0
    for n in range(1,iteratsioone_1+1):
        print(cos_y_ln_x_multinomial_4_2(n, x_2), cos(x_2 * log(n)))
        #v_Re+=cos_y_ln_x_multinomial_5(n,x_2)*n**(-x_1)
        v_Re+=cos(x_2*log(n))*n**(-x_1)
    v_Im=0
    for n in range(1,iteratsioone_2+1):
        #print(sin_y_ln_x_multinomial_5(n,x_2),sin(x_2*log(n)))
        #v_Im+=sin_y_ln_x_multinomial_5(n,x_2)*n**(-x_1)
        v_Im+=sin(x_2*log(n))*n**(-x_1)

    return(v_Re-1j*v_Im)

def suhted(z):
    return zeta5(z.real,z.imag).real / zeta1_Dirichlet_series_Wikipediast(z).real, zeta5(z.real, z.imag).imag / zeta1_Dirichlet_series_Wikipediast(z).imag

#print(zeta1(1-2j,iteratsioone=10000))
#input()
#print(zeta6(3.54+5.21j))z


#print(suhted(4+1j))
#print(suhted(1-1j))
#print(suhted(2-2.5j))
zetad=[zeta1_Dirichlet_series_Wikipediast, zeta2_reaalosa_ja_imaginaarosa_eraldatud, zeta3, zeta4, zeta5, zeta6, zeta7, zeta8, zeta9_M_plus_1e_ei_ole, zeta10, zeta11, zeta12_q1_v??rdub_q2, zeta12_1, zeta13_tagumine_pool_nulliks_optimeeritud]
toimivad_zetad=[zeta1_Dirichlet_series_Wikipediast, zeta2_reaalosa_ja_imaginaarosa_eraldatud, zeta3, zeta4, zeta5, zeta6, zeta7, zeta8, zeta9_M_plus_1e_ei_ole,zeta13_tagumine_pool_nulliks]

def test(zeta=zeta2_reaalosa_ja_imaginaarosa_eraldatud):
    print("test:",zeta)
    zeta_v????rtused={
        0.5+6j:0.837223808066879519392085379114373422068052148992916728465961897+0.340218396943766415292940629924307409809030623790098941796105413j,
        3.54+5.21j:0.94145150136838204399961988000425280352177673217144192099708990+0.042261487481304266846570711131559236624545183101533092693667295j,
        4-5j:0.95121830700949569861514024576822624312072806034027557374540065-0.024932824507357722259855155244740571939365897230032890678838269j,
        2.93+40.22j:0.91927363005266867503823440502290173219111893589113817500208901-0.058873893494407299524527222529938166896010810611721303120020614j,
        0.7+11.33j:1.26076751803400086687692006648914300098585746351039557218708372-0.528294513391831425306205096828572148868207356151099491411787630j,
        10-6j:0.999502256328600532560326101538767375141432686780150206823351-0.000824796183729373433431091029946446889482127040642667084842333j,
        4+0.05j:1.082241936785460782743310833364791724793475643993049380159631-0.003444050373342241129776836665066108430086147497362107906399599j,
        0.1+12j:1.08189989714261188605953655665870352938936153821212593030587662-0.986030063319784609891282060634672603258161005766600811809208416j
    }
    for zeta_argument,zeta_v????rtus in zeta_v????rtused.items():
        vastus=zeta(zeta_argument)
        print("jagatised:",vastus.real/zeta_v????rtus.real,vastus.imag/zeta_v????rtus.imag,"vahed:",vastus.real-zeta_v????rtus.real,vastus.imag-zeta_v????rtus.imag)
print(77)
test(zeta13_tagumine_pool_nulliks)