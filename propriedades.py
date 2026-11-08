from eos import PR
from constantes import R
import numpy as np

class entropia: 
    def __init__(self):
        self.nome = "Entropia Residual"
    def residual(self, eos, T, P, n, fase): # Aqui vamos passar a instância da eos (equaação com substância), a temperatura, a pressão e a fase (líquido ou vapor) para calcular a entropia residual.
        # Cálculo da entropia residual usando a equação de estado de Peng-Robinson
        # Se houver valores de T e P, podemos calcular o residual direto.
        # n é a quantidade de matéria da substância.
        # Se n = 1, teremos base molar
        if T is not None and P is not None:
            #dS = Cv/T dT + (dP/dT)v dV
            # Nesta parte, calculamos a entropia residual, ou seja, a diferença entre a entropia real e a entropia ideal. 
            # Sres = S - Sideal
            # Primeiramente vamos achar o valor de de dalpha/dT, que é necessário para o cálculo da entropia residual.
            dalpha_dT = eos.dalpha_dT(T)
            #Agora vamos definir K como sendo o seguinte termo:
            ac = eos.ac()
            b = eos.b()
            # Podemos calcular o termo dependente de Z: 
            Z = eos.calcZ(T,P)
            if fase == "liquido":
                Z = min(Z)
            elif fase == "vapor":
                Z = max(Z)
            else:
                raise ValueError("Fase não reconhecida. Por favor, escolha entre 'liquido' ou 'vapor'.")
            # Pode-se calcular o volume molar a partir de Z:
            V = Z*R*T/P
            Zdependente = R*np.log(Z)
            # Agora podemos calcular a integral de infinito a V do residual: 
            # Fazendo no papel, chegamos na função analítica da entropia residual, que é a seguinte:
            integral_fparciais = eos.c*(np.log((V-eos.epsilon*b)/(V-eos.sigma*b)))
            Sres = R*np.log(V-b) - R*np.log(V) - dalpha_dT*ac*integral_fparciais + Zdependente
            Sres = Sres*n # Multiplicando pela quantidade de matéria da substância, caso seja necessário.
            return Sres
        else:
            raise ValueError("Valores de T e P não fornecidos. Por favor, forneça valores válidos para T e P.")
            
    def ideal (self, substancia, T1, T2, P1, P2, n): # Aqui vamos passar a instância da substância, a temperatura inicial e final e a pressão inicial e final para calcular a entropia ideal.
        # Aqui vamos calcular a entropia ideal, que é a entropia de uma substância em condições ideais.
        # Se n = 1, teremos base molar
        # A entropia ideal é calculada a partir da equação de estado ideal, que é a seguinte:
        Cp = substancia.coef_cp
        integral_Cp = 0
        if T1 is not None and T2 is not None and P1 is not None and P2 is not None:
            for i in range(len(Cp)):
                if i == 0:
                    integral_Cp += Cp[i]*np.log(T2/T1)
                else:
                    #Vamos fazer no loop for a integral dos coeficientes do CP: 
                    integral_Cp += Cp[i]/(i)*(T2**(i)-T1**(i))
            # Agora podemos calcular a entropia ideal:
            Sideal = integral_Cp - R*np.log(P2/P1)
            Sideal = Sideal*n # Multiplicando pela quantidade de matéria da substância, caso seja necessário.
        else:
            raise ValueError("Valores de T1, T2, P1 e P2 não fornecidos. Por favor, forneça valores válidos para T1, T2, P1 e P2.")
        return Sideal

class entalpia: 
    def __init__(self):
        self.nome = "Entalpia Residual"
    def residual(self, eos, T, P, fase): # Aqui vamos passar a instância da eos (equaação com substância), a temperatura, a pressão e a fase (líquido ou vapor) para calcular a entalpia residual.
        # Cálculo da entalpia residual usando a equação de estado de Peng-Robinson
        # Se houver valores de T e P, podemos calcular o residual direto.
        if T is not None and P is not None:
            #dH = Cp dT + (V-T(dP/dT)v)dP
            # Nesta parte, calculamos a entalpia residual, ou seja, a diferença entre a entalpia real e a entalpia ideal. 
            # Hres = H - Hideal
            # Primeiramente vamos achar o valor de de dalpha/dT, que é necessário para o cálculo da entalpia residual.
            dalpha_dT = eos.dalpha_dT(T)
            ac = eos.ac()
            b = eos.b()
            Z = eos.calcZ(T,P)
            if fase == "liquido":
                Z = min(Z)
            elif fase == "vapor":
                Z = max(Z)
            else:
                raise ValueError("Fase não reconhecida. Por favor, escolha entre 'liquido' ou 'vapor'.")
            V = Z*R*T/P
            Zdependente = R*(Z-1)
            integral_fparciais = eos.c*(np.log((V-eos.epsilon*b)/(V-eos.sigma*b)))
            Hres = -dalpha_dT*ac*integral_fparciais + Zdependente
            return Hres
        else:
            raise ValueError("Valores de T e P não fornecidos. Por favor, forneça valores válidos para T e P.")
    def ideal (self, substancia, T1, T2, P1, P2): # Aqui vamos passar a instância da substância, a temperatura inicial e final e a pressão inicial e final para calcular a entalpia ideal.
        # Cálculo da entalpia ideal usando a equação de estado de Peng-Robinson
        # Aqui vamos calcular a entalpia ideal, que é a entalpia de uma substância em condições ideais.
        # A entalpia ideal é calculada a partir da equação de estado ideal, que é a seguinte:
        Cp = substancia.coef_cp
        integral_Cp = 0
        if T1 is not None and T2 is not None and P1 is not None and P2 is not None:
            for i in range(len(Cp)):
                #Vamos fazer no loop for a integral dos coeficientes do CP: 
                integral_Cp += Cp[i]/(i+1)*(T2**(i+1)-T1**(i+1))
            # Agora podemos calcular a entalpia ideal:
            Hideal = integral_Cp*np.log(T2/T1) - R*np.log(P2/P1)
        else:
            raise ValueError("Valores de T1, T2, P1 e P2 não fornecidos. Por favor, forneça valores válidos para T1, T2, P1 e P2.")
        return Hideal

class energia_interna: 
    def __init__(self):
        self.nome = "Energia Interna Residual"
    def residual(self, eos, T, P,n, fase): # Aqui vamos passar a instância da eos (equaação com substância), a temperatura, a pressão e a fase (líquido ou vapor) para calcular a energia interna residual.
        # Cálculo da energia interna residual usando a equação de estado de Peng-Robinson
        # Se n = 1, teremos base molar
        # Se houver valores de T e P, podemos calcular o residual direto.
        if T is not None and P is not None:
            #dU = Cv dT + (dP/dT)v dV
            # Nesta parte, calculamos a energia interna residual, ou seja, a diferença entre a energia interna real e a energia interna ideal. 
            # Ures = U - Uideal
            # Primeiramente vamos achar o valor de de dalpha/dT, que é necessário para o cálculo da energia interna residual.
            dalpha_dT = eos.dalpha_dT(T)
            alpha = eos.alpha(T)
            ac = eos.ac()
            b = eos.b()
            Z = eos.calcZ(T,P)
            if fase == "liquido":
                Z = min(Z)
            elif fase == "vapor":
                Z = max(Z)
            else:
                raise ValueError("Fase não reconhecida. Por favor, escolha entre 'liquido' ou 'vapor'.")
            V = Z*R*T/P
            K = alpha - T*dalpha_dT
            integral_fparciais = eos.c*(np.log((V-eos.epsilon*b)/(V-eos.sigma*b)))
            Ures = K*ac*integral_fparciais
            Ures = Ures*n # Multiplicando pela quantidade de matéria da substância, caso seja necessário.
            return Ures
        else:
            raise ValueError("Valores de T e P não fornecidos. Por favor, forneça valores válidos para T e P.")
    def ideal (self, substancia, T1, T2, n): # Aqui vamos passar a instância da substância, a temperatura inicial e final e a pressão inicial e final para calcular a energia interna ideal.
        # Aqui vamos calcular a energia interna ideal, que é a energia interna de uma substância em condições ideais.
        # Se n = 1, teremos base molar
        # A energia interna ideal é calculada a partir da equação de estado ideal, que é a seguinte:
        Cp = substancia.coef_cp
        Cv = 0
        if T1 is not None and T2 is not None:
        # Vamos calcular o valor do CV
            for i in range(len(Cp)):
                if i == 0:
                    Cv = (Cp[i] - R) * (T2 - T1)
                else:
                    Cv += Cp[i]/(i+1)*(T2**(i+1)-T1**(i+1))
            Uideal = Cv 
            Uideal = Uideal*n # Multiplicando pela quantidade de matéria da substância, caso seja necessário.
        else:
            raise ValueError("Valores de T1, T2 não fornecidos. Por favor, forneça valores válidos para T1 e T2.")
        return Uideal