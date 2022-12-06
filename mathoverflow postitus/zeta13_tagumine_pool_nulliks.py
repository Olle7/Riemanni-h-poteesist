from math import log,atan,floor,sin,pi,e

def sgn(x):
    if x>0:
        return 1
    if x<0:
        return -1
    if x==0:
        return 0

def zeta13(*args, **kwargs):
    #THE BIGGER ISNUMBER OF REQUESTED ITERATIONS THE MORE ACCURATE IS RESULT. KEYWORD "iterations" DETERMINES NUMBER OF REQUESTED ITERATIONS.
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
    if "iterations" in kwargs:
        iterations=kwargs["iterations"]
    else:
        iterations=10**6

    M=max(
    (log(iterations)*abs(x_2)+atan((1-x_1)/x_2)*sgn(x_2))*2/(2*pi),
    (log(iterations)+atan((1-x_1)/x_2)/x_2)*2*abs(x_2)/(2*pi)-1/2
    )
    iterations_1=floor(e**(floor(M)/2*2*pi/abs(x_2)-atan((1-x_1)/x_2)/x_2))
    iterations_2=floor(e**((floor(M)/2+1/4)*2*pi/abs(x_2)-atan((1-x_1)/x_2)/x_2))
    #print("M=",M)

    v_Re=0
    for n in range(1,iterations_1+1):
        v_Re+=sin(x_2*log(n)+1*(2*pi)/4)*n**(-x_1)
    v_Im=0
    for n in range(1,iterations_2+1):
        v_Im+=sin(x_2*log(n))*n**(-x_1)

    return(v_Re-1j*v_Im)


print(zeta13(0.734+1.13256j,iterations=10**7))