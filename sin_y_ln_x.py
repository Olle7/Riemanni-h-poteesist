from math import *
def sin_y_ln_x_maclaurin(x, y, M1=85,M2=300):
    latex=r"""$\forall_x(x>0 \to\\
sin(y*ln(x))=\sum_{k_1=0}^\infty((-1)^{k_1}*( 2*\sum_{{k_2}=0}^\infty((2*{k_2}+1)^{-1}*(x-1)^{2*{k_2}+1}*(x+1)^{-2*{k_2}-1}) )^{2*k_1+1}*y^{2*k_1+1}*\Gamma(2*k_1+2)^{-1})\\
)$"""
    s1=0
    for k1 in range(M1):
        s2=0
        for k2 in range(M2):
            s2+=(2*k2+1)**-1*(x-1)**(2*k2+1)*(x+1)**(-2*k2-1)
        s1+=(-1)**k1*2**(2*k1+1)*s2**(2*k1+1)*y**(2*k1+1)*gamma(2*k1+2)**-1
    return s1


def sin_x_maclaurin(x,M=85):
    latex=r"""$\forall_x(x \in R \to\\
sin(x)=\sum_{k=0}^\infty((-1)^k*x^{2*k+1}*\Gamma(2*k+2)^{-1})\\
)$"""
    s1=0
    for k1 in range(M):
        #s2=0
        #for k2 in range(M):
        #    s2+=(2*k2+1)**-1*(x-1)**(2*k2+1)*(x+1)**(-2*k2-1)
        s1+=(-1)**k1*x**(2*k1+1)*gamma(2*k1+2)**-1
    return s1

def sin_x_maclaurin_optimeeritud(x,M=85):
    s1=0
    for k1 in range(M):
        #s2=0
        #for k2 in range(M):
        #    s2+=(2*k2+1)**-1*(x-1)**(2*k2+1)*(x+1)**(-2*k2-1)
        s1+=(-1)**k1*x**(2*k1+1)*factorial(2*k1+1)**-1
    return s1

def sin_y_ln_x_maclaurin_optimeeritud(x, y, M1=10000,M2=100000):
    s = 0
    try:
        for k in range(M1):
            s+=(2*k+1)**-1*((x-1)/(x+1))**(2*k+1)
    except OverflowError:
        pass
    x=2*s*y

    s=0
    try:
        for k in range(M2):
            s+=(-1)**k*x**(2*k+1)*factorial(2*k+1)**-1
    except OverflowError:
        pass
    return s

def sin_y_ln_x_maclaurin_3(x,y,M1=10000,M2=100000):
    latex = r"""$\forall_x(x>0 \to\\
    sin(y*ln(x))=\sum_{k_1=0}^\infty(( \sum_{{k_2}=0}^\infty((2*k_2+1)^{-1}*(x-1)^{2*k_2+1}*(x+1)^{-2*k_2-1}) )^{2*k_1+1}*y^{2*k_1+1}*factorial(2*k_1+1)^{-1}*2^{2*k_1+1}(-1)^{k_1})\\
    )$"""

    s1=0
    try:
        for k1 in range(M2):
            s2=0
            try:
                for k2 in range(M1):
                    s2+=(2*k2+1)**-1*(x-1)**(2*k2+1)*(x+1)**(-2*k2-1)
            except OverflowError:
                pass
            ln_x=s2
            
            s1+=ln_x**(2*k1+1)*y**(2*k1+1)*factorial(2*k1+1)**-1*2**(2*k1+1)*(-1)**k1
    except OverflowError:
        pass
    return s1

def sin_y_ln_x_multinomial_1(x,y,M=6):##TOIMIB
    latex=r"""$\forall_x(x>0 \to\\
sin(y*ln(x))=lim_{M \to \infty}(\sum_{k_1=0}^{M-1}(\sum_{k_2=0}^{(2*k_1+2)^M-1}(\\
(\sum_{k_3=0}^{M-1}(\ floor(k_2/(2*k_1+2)^{k_3})\%(2*k_1+2)\ )=2*k_1+1)*\\%sulgudega 0
\Pi_{k_3=0}^{M-1}(factorial(\ floor(k_2/(2*k_1+2)^{k_3})\%(2*k_1+2)\ )^{-1}*((2*k_3+1)^{-1}*(x-1)^{2*k_3+1}*(x+1)^{-2*k_3-1})^{floor(k_2/(2*k_1+2)^{k_3})\%(2*k_1+2)} )\\%sulgudega 0
)*y^{2*k_1+1}*2^{2*k_1+1}*(-1)^{k_1}))\\%k1-summa, k2-summa ja lim'i lopusulg sel real.
)%forall lopusulg$"""

    s1=0
    for k1 in range(M):
        #print(k1)
        s2=0
        for k2 in range((2*k1+2)**M):
            s3=0
            for k3 in range(M):
                s3+=(k2//(2*k1+2)**(M-1-k3))%(2*k1+2)
            if s3!=2*k1+1:
                continue
            p3=1
            for k3 in range(M):
                p3*=factorial((k2//(2*k1+2)**(M-1-k3))%(2*k1+2))**-1* ( (2*k3+1)**-1*(x-1)**(2*k3+1)*(x+1)**(-2*k3-1)    )**((k2//(2*k1+2)**(M-1-k3))%(2*k1+2 ))
            s2+=p3
        s1+=s2*y**(2*k1+1)*2**(2*k1+1)*(-1)**k1
    return s1

def sin_y_ln_x_multinomial(x,y,M=6):
    latex=r"""$\forall_x(x>0 \to sin(y*ln(x))=\\
lim_{M \to \infty}(\sum_{k_1=0}^{M-1}(y^{2*k_1+1}*2^{2*k_1+1}*(-1)^{k_1}*\sum_{k_2=0}^{(2*k_1+2)^M-1}(\\
(\sum_{k_3=0}^{M-1}(\ floor(k_2/(2*k_1+2)^{k_3})\%(2*k_1+2)\ )=2*k_1+1)*\\%sulgudega 0
(x-1)^{ \sum_{k_3=0}^{M-1}((2*k_3+1)*(floor(k_2/(2*k_1+2)^{k_3})\%(2*k_1+2)))}*\\
(x+1)^{-\sum_{k_3=0}^{M-1}((2*k_3+1)*(floor(k_2/(2*k_1+2)^{k_3})\%(2*k_1+2)))}*\\
\Pi_{k_3=0}^{M-1}(factorial(\ floor(k_2/(2*k_1+2)^{k_3})\%(2*k_1+2)\ )^{-1}*(2*k_3+1)^{-floor(k_2/(2*k_1+2)^{k_3})\%(2*k_1+2)} )\\
)))\\%k1-summa, k2-summa ja lim'i lopusulg sel real.
)%forall lopusulg$"""

    s1=0
    for k1 in range(M):
        #print(k1)
        s2=0
        for k2 in range((2*k1+2)**M):
            p=1
            s3=0
            for k3 in range(M):
                s3+=(k2//(2*k1+2)**k3)%(2*k1+2)
            if s3!=2*k1+1:
                continue
            s3=0
            for k3 in range(M):
                s3+=(2*k3+1)*((k2//(2*k1+2)**k3)%(2*k1+2))
            p*=(x-1)**s3
            p*=(x+1)**-s3

            p3=1
            for k3 in range(M):
                p3*=factorial((k2//(2*k1+2)**k3)%(2*k1+2))**-1*(2*k3+1)**(-((k2//(2*k1+2)**k3)%(2*k1+2)))
            p*=p3
            s2+=p
        s1+=s2*y**(2*k1+1)*2**(2*k1+1)*(-1)**(k1)
    return s1

def sin_y_ln_x_multinomial_2(x,y,M=6):
    latex=r"""$\forall_x(x>0 \to sin(y*ln(x))=\\
-lim_{M \to \infty}(\sum_{k_1=0}^{(2*M)^M-1}(\\
(\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))\%2)*\\
(\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))<2*M)*\\%sulgudega 0
(y*2)^{\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))}*\\
(-1)^{(\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))-1)/2}*\\
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
        p*=s2<2*M
        if p==0:
            continue
        p*=(y*2)**s2
        #print("s2:",s2)
        p*=(-1)**((s2-1)/2)

        s2=0
        for k2 in range(M):
            s2+=(2*k2+1)* ((k1//(2*M)**k2)%(2*M))
        p*=(x-1)**s2
        p*=(x+1)**-s2

        for k2 in range(M):
            p*=factorial((k1//(2*M)**k2)%(2*M))**(-1)*(2*k2+1)**(-((k1//(2*M)**k2)%(2*M)))
        s1+=p
    return s1

def sin_y_ln_x_multinomial_4(x,y,M=6):
    latex=r"""$\forall_x(x>0 \to sin(y*ln(x))=\\
-lim_{M \to \infty}(\sum_{k_1=0}^{(2*M)^M-1}(\\
(\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))\%2)*\\
(y*2)^{\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))}*\\
(-1)^{(\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))-1)/2}*\\
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
        p*=(-1)**((s2-1)/2)

        s2=0
        for k2 in range(M):
            s2+=(2*k2+1)* ((k1//(2*M)**k2)%(2*M))
        p*=(x-1)**s2
        p*=(x+1)**-s2

        for k2 in range(M):
            p*=factorial((k1//(2*M)**k2)%(2*M))**(-1)*(2*k2+1)**(-((k1//(2*M)**k2)%(2*M)))
        s1+=p
    return s1

def sin_y_ln_x_multinomial_3(x,y,M=6):
    latex=r"""$\forall_x(x>0 \to sin(y*ln(x))=\\
-lim_{M \to \infty}(\sum_{k_1=0}^{(2*M)^M-1}(\\
(\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))\%2)*\\
(\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))<2*M)*\\%sulgudega 0
(y*2)^{\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))}*\\
(-1)^{(\sum_{k_2=0}^{M-1}(floor(k_1/(2*M)^{k_2})\%(2*M))-1)/2}*\\
(x-1)^{ \sum_{k_2=0}^{M-1}((2*k_2+1)*(floor(k_1/(2*M)^{k_2})\%(2*M)))}*\\
(x+1)^{-\sum_{k_2=0}^{M-1}((2*k_2+1)*(floor(k_1/(2*M)^{k_2})\%(2*M)))}*\\
\Pi_{k_2=0}^{M-1}(factorial(floor(k_1/(2*M)^{k_2})\%(2*M))^{-1}*(2*k_2+1)^{-floor(k_1/(2*M)^{k_2})\%(2*M)} )\\
))\\%k2-summa ja lim'i lopusulg sel real.
)%forall lopusulg$"""
    s1=0
    for k1 in range((M)**M):
        p=1
        s2=0
        for k2 in range(M):
            s2+=(k1//(M)**k2)%(M)
        p*=s2%2
        p*=s2<M-1
        if p==0:
            continue
        p*=(y)**s2
        #print("s2:",s2)
        p*=(-1)**((s2-1)/2)

        s2=0
        for k2 in range(M):
            s2+=(2*k2+1)* ((k1//(M)**k2)%(M))
        p*=(x-1)**s2
        p*=(x+1)**-s2

        for k2 in range(M):
            p*=factorial((k1//(M)**k2)%(M))**(-1)*(2*k2+1)**(-((k1//(M)**k2)%(M)))
        s1+=p
    return s1


def cos_y_ln_x_multinomial(x,y,M=6):
    latex=r"""$\forall_x(x>0 \to cos(y*ln(x))=\\
lim_{M \to \infty}(\sum_{k_1=0}^{M-1}(y^{2*k_1}*2^{2*k_1}*(-1)^{k_1}*\sum_{k_2=0}^{(2*k_1+1)^M-1}(\\
(\sum_{k_3=0}^{M-1}(\ floor(k_2/(2*k_1+1)^{k_3})\%(2*k_1+1)\ )=2*k_1)*\\%sulgudega 0
(x-1)^{ \sum_{k_3=0}^{M-1}((2*k_3+1)*(floor(k_2/(2*k_1+1)^{k_3})\%(2*k_1+1)))}*\\
(x+1)^{-\sum_{k_3=0}^{M-1}((2*k_3+1)*(floor(k_2/(2*k_1+1)^{k_3})\%(2*k_1+1)))}*\\
\Pi_{k_3=0}^{M-1}(factorial(\ floor(k_2/(2*k_1+1)^{k_3})\%(2*k_1+1)\ )^{-1}*(2*k_3+1)^{-floor(k_2/(2*k_1+1)^{k_3})\%(2*k_1+1)} )\\
)))\\%k1-summa, k2-summa ja lim'i lopusulg sel real.
)%forall lopusulg$"""

    s1=0
    for k1 in range(M):
        #print(k1)
        s2=0
        for k2 in range((2*k1+1)**M):
            p=1
            s3=0
            for k3 in range(M):
                s3+=(k2//(2*k1+1)**k3)%(2*k1+1)
            if s3!=2*k1:
                continue
            s3=0
            for k3 in range(M):
                s3+=(2*k3+1)*((k2//(2*k1+1)**k3)%(2*k1+1))
            p*=(x-1)**s3
            p*=(x+1)**-s3

            p3=1
            for k3 in range(M):
                p3*=factorial((k2//(2*k1+1)**k3)%(2*k1+1))**-1*(2*k3+1)**(-((k2//(2*k1+1)**k3)%(2*k1+1)))
            p*=p3
            s2+=p
        s1+=s2*y**(2*k1)*2**(2*k1)*(-1)**(k1)
    return s1

def sin_y_ln_x_multinomial_optimeeritud(x, y, M=7):
    latex=r"""$\forall_x(x>0 \to sin(y*ln(x))=\\
lim_{M \to \infty}(\sum_{k_1=0}^M(y^{2*k_1+1}*2^{2*k_1+1}*(-1)^{k_1}*\sum_{k_2=0}^{(2*k_1+1)^M}(\\
(\sum_{k_3=0}^{M}(\ floor(k_2/(2*k_1+1)^{k_3})\%(2*k_1+1)\ )=2*k_1+1)*\\%sulgudega 0
(x-1)^{ \sum_{k_3=0}^{M}((2*k_3+1)*(floor(k_2/(2*k_1+1)^{k_3})\%(2*k_1+1)))}*\\
(x+1)^{-\sum_{k_3=0}^{M}((2*k_3+1)*(floor(k_2/(2*k_1+1)^{k_3})\%(2*k_1+1)))}*\\
\Pi_{k_3=0}^{M}(factorial(\ floor(k_2/(2*k_1+1)^{k_3})\%(2*k_1+1)\ )^{-1}*(2*k_3+1)^{-floor(k_2/(2*k_1+1)^{k_3})\%(2*k_1+1)} )\\
)))\\%k1-summa, k2-summa ja lim'i lopusulg sel real.
)%forall lopusulg$"""

    s1=0
    for k1 in range(M):
        #print(k1)
        s2=0
        for k2 in range((2*k1+1)**M):
            p=1
            s3=0
            for k3 in range(M):
                s3+=(k2//(2*k1+1)**k3)%(2*k1+1)
            if s3!=2*k1+1:
                continue
            s3=0
            for k3 in range(M):
                s3+=(2*k3+1)*((k2//(2*k1+1)**k3)%(2*k1+1))
            p*=((x-1)/(x+1))**s3

            p3=1
            for k3 in range(M):
                p3*=factorial((k2//(2*k1+1)**k3)%(2*k1+1))**-1*(2*k3+1)**(-((k2//(2*k1+1)**k3)%(2*k1+1)))
            p*=p3
            s2+=p
        s1+=s2*y**(2*k1+1)*2**(2*k1+1)*(-1)**(k1)
    return s1

def kontrolli_sin_y_ln_x(sin_y_ln_x_maclaurin_funktsioon):
    for x_y in [(2.845902,1),(1,2),(4,2),(2,3),(0.234,1),(3.232,0.3),(e**pi,3),]:
        x,y=x_y[0],x_y[1]
        õige=sin(y*log(x))
        vastus=sin_y_ln_x_maclaurin_funktsioon(x,y)
        print(õige-vastus,"vastus:",vastus,"; õige:",õige,"; x:",x,"; y:",y)

def kontrolli_cos_y_ln_x(sin_y_ln_x_maclaurin_funktsioon):
    for x_y in [(2.845902,1),(1,2),(4,2),(2,3),(0.234,1),(3.232,0.3),(e**pi,3),]:
        x,y=x_y[0],x_y[1]
        õige=cos(y*log(x))
        vastus=sin_y_ln_x_maclaurin_funktsioon(x,y)
        print(õige-vastus,"vastus:",vastus,"; õige:",õige,"; x:",x,"; y:",y)

def kontrolli_sin_x(sin_x):
    for x in [1,2,3,pi,2*pi,3.1,-0.4,0.243]:
        õige=sin(x)
        vastus=sin_x(x)
        print(õige-vastus,"vastus:",vastus,"; õige:",õige,"; x:",x)

def ln_maclaurin(x,M=10000):#toimib
    latex=r"""$\forall_x(x>0 \to\\
ln(x)=2*\sum_{k=0}^\infty((2*k+1)^{-1}*(x-1)^{2*k+1}*(x+1)^{-2*k-1})\\
)$"""
    try:
        s=0
        for k in range(M):
            s+=(2*k+1)**-1*(x-1)**(2*k+1)*(x+1)**(-2*k-1)
    except OverflowError:
        pass
    return 2*s

#print(sin(log(3)))
#print(sin_y_ln_x_multinomial_optimeeritud(3,1,1))
#print(sin_y_ln_x_multinomial_optimeeritud(3,1,2))
#print(sin_y_ln_x_multinomial_optimeeritud(3,1,3))
#print(sin_y_ln_x_multinomial_optimeeritud(3,1,4))
#print(sin_y_ln_x_multinomial_optimeeritud(3,1,5))
#print(sin_y_ln_x_multinomial_optimeeritud(3,1,6))
#print(sin_y_ln_x_multinomial_optimeeritud(3,1,7))
#print(sin_y_ln_x_multinomial_optimeeritud(3,1,8))

#print(sin_y_ln_x_multinomial_1(2.345,1,6),sin_y_ln_x_multinomial(2.345,1,6))
#print(sin_y_ln_x_multinomial_1(2.345,2,6),sin_y_ln_x_multinomial(2.345,2,6))
#print(sin_y_ln_x_multinomial_1(3.5342,1.434,6),sin_y_ln_x_multinomial(3.5342,1.434,6))
#print(sin_y_ln_x_multinomial_1(1.3443,0.345,6),sin_y_ln_x_multinomial(1.3443,0.345,6))
#print(sin_y_ln_x_multinomial_1(0.5,0.5,6),sin_y_ln_x_multinomial(0.5,0.5,6))

#print(ln_maclaurin(e**4.2642))
kontrolli_sin_y_ln_x(sin_y_ln_x_multinomial_2)
#kontrolli_sin_x(sin_x_maclaurin_optimeeritud)