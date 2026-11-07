import numpy as np
import sympy as sp
from constantes import R
class PR: 
    def __init__(self,substancia):
        self.nome = "Peng-Robinson"  
        self.equacao = "P = RT/(V-b) - a/(V(V+b)+b(V-b))"
        self.substancia = substancia
        self.sigma = -(1+2**0.5)
        self.epsilon = -(1-2**0.5)
        self.c = 1/(2*(2**0.5)*self.b())

    def kappa(self):
        w = self.substancia.w
        kappa = 0.37464 + 1.54226*w - 0.26992*w**2
        return kappa

    def alpha(self,T):
        Tc = self.substancia.tc
        Tr = T/Tc
        kappa = self.kappa()
        alpha = (1 + kappa*(1-Tr**0.5))**2
        return alpha

    def ac(self):
        Tc = self.substancia.tc
        Pc = self.substancia.pc
        ac = (0.45724*R**2*Tc**2/Pc)
        return ac

    def a(self,T):
        ac = self.ac()
        alpha = self.alpha(T)
        a = ac*alpha
        return a

    def b(self):
        Tc = self.substancia.tc
        Pc = self.substancia.pc
        b = 0.07780*R*Tc/Pc
        return b

    def calcA(self,T,P):
        a = self.a(T)
        A = a*P/(R**2*T**2)
        return A

    def calcB(self,T,P):
        b = self.b()
        B = b*P/(R*T)
        return B

    def calcZ(self,T,P):
        A = self.calcA(T,P)
        B = self.calcB(T,P)
        # Coeficientes do polinômio cúbico em Z
        # Z^3 - (1-B)Z^2 + (A-3B^2-2B)Z - (AB-B^2-B^3) = 0
        raizes = np.roots([1, -(1-B), A-3*B**2-2*B, -(A*B-B**2-B**3)])
        raizes = self.filtrar_raizes(raizes)
        return raizes
    # Aqui estamos adotando a definição de que a raiz menor é a raiz líquida e a maior é a raiz de vapor. O mesmo vale para volume.

    def filtrar_raizes(self,raizes):
        raizes_reais = raizes[np.isreal(raizes)].real
        raizes_positivas = raizes_reais[raizes_reais > 0]
        if len(raizes_positivas) == 3:
            raiz_vapor = max(raizes_positivas)
            raiz_liquido = min(raizes_positivas)
            raizes = [raiz_liquido, raiz_vapor]
            return raizes
        elif len(raizes_positivas) == 2:
            raiz_vapor = max(raizes_positivas)
            raiz_liquido = min(raizes_positivas)
            raizes = [raiz_liquido, raiz_vapor]
            return raizes
        elif len(raizes_positivas) == 1:
            raiz_vapor = raizes_positivas[0]
            raizes = [raiz_vapor]
            return raizes
        else:
            return None

    def calcV(self,T,P):
        Z = self.calcZ(T,P)
        if Z is None:
            return None
        else:
            i = 0
            for z in Z:
                if len(Z) == 1:
                    V = z*R*T/P
                else:
                    if i == 0:
                        Vliq = z*R*T/P
                    else:
                        Vvap = z*R*T/P
                i += 1
            if len(Z) == 1:
                Vtotal = [V]
            else:
                Vtotal = [Vliq, Vvap]
            return Vtotal

    def dalpha_dT(self,T):
        Tc = self.substancia.tc
        Tr = T/Tc
        kappa = self.kappa()
        dalpha_dT = -kappa/((T*Tc)**0.5)*(1 + kappa*(1-Tr**0.5))
        return dalpha_dT
