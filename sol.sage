import itertools
from subprocess import check_output
import fgb_sage
sys.set_int_max_str_digits(20000)

def flatter(M):
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    from re import findall
    return matrix(M.nrows(), M.ncols(), map(int, findall(b"-?\\d+", ret)))
N1 = 11195108435418195710792588075406654238662413452040893604269481198631380853864777816171135346615239847585274781942826320773907414919521767450698159152141823148113043170072260905000812966959448737906045653134710039763987977024660093279241536270954380974093998238962759705207561900626656220185252467266349413165950122829268464816965028949121220409917771866453266522778041491886000765870296070557269360794230165147639201703312790976341766891628037850902489808393224528144341687117276366107884626925409318998153959791998809250576701129098030933612584038842347204032289231076557168670724255156010233010888918002630018693299
N2 = 15100254697650550107773880032815145863356657719287915742500114525591753087962467826081728465512892164117836132237310655696249972190691781679185814089899954980129273157108546566607320409558512492474972517904901612694329245705071789171594962844907667870548108438624866788136327638175865250706483350097727472981522495023856155253124778291684107340441685908190131143526592231859940556416271923298043631447630144435617140894108480182678930181019645093766210388896642127572162172851596331016756329494450522133805279328640942549500999876562756779916153474958590607156569686953857510763692124165713467629066731049974996526071
alpha = .125
gamma = .42
b = 1 - alpha - gamma
b1 = 2^(1792-1775)
b2 = 2^(1792-1126)
bg1 = 2^(1792-905)
bg2 = 2^(1792-256)
bounds = [b2,2^905,2^256,2^1792]
PR.<x,y,z,w> = PolynomialRing(ZZ)
f = x*z + bg2*y*z + ZZ(N1)
K = (b2//b1)
G = Sequence([], f.parent())
m = 9
t = 4
s = 3
for i in range(m+1):
    for j in range(0,m-i+1):
        g = (y*z)^j*w^s*f^i*K^(m-i)*ZZ(N2)^(max(t-i,0))*ZZ(N1)^(-min(s,i+j))
        monomials = vector(g.monomials())
#         print(monomials)
        cof = vector(g.coefficients())
        for k in range(len(monomials)):
            monomial = monomials[k]
            while monomial % (z*w) == 0:
                monomial = N1*monomial//(z*w)
            monomials[k] = monomial
        g = cof*monomials
#         print(i,j,monomials)
#         assert g(xx,yy,zz,ww) % (K^m*q2^t) == 0
        G.append(g)
B, monomials = G.coefficient_matrix()
monomials = vector(monomials)
factors = [monomial(*bounds) for monomial in monomials]
for i, factor in enumerate(factors):
    B.rescale_col(i, factor)
print("dim:",B.nrows(),B.ncols())
B = flatter(B)
# except:
#     B = B.dense_matrix().LLL()
B = B.change_ring(QQ)
for i, factor in enumerate(factors):
    B.rescale_col(i, 1 / factor)
RR = B*monomials
# for i in range(len(RR)):
#     if RR[i](xx,yy,zz,ww)==0 and RR[i]!=0:
#         print(i)
RR[0] = z*w - ZZ(N1)
res = fgb_sage.groebner_basis(RR[0:20])
for r in res:
    for coef in r.coefficients():
        if N1%coef==0 and coef !=1 and abs(coef)!=N1:
            RR[1] = w - coef
            V = fgb_sage.groebner_basis(RR[0:20])
            print("x =",x - V[0])
            print("y =",y - V[1])
            print("z =",z - V[2])
            print("w =",w - V[3])
            break
