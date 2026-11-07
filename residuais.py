from eos import PR
from constantes import R
import numpy as np

class entropia: 
    def __init__(self):
        self.nome = "Entropia Residual"
    def residual(self, eos, T, P, fase): # Aqui vamos passar a instância da eos (equaação com substância), a temperatura, a pressão e a fase (líquido ou vapor) para calcular a entropia residual.
        # Cálculo da entropia residual usando a equação de estado de Peng-Robinson
        # Se houver valores de T e P, podemos calcular o residual direto.
        if T and P is not None:
            #dS = Cv/T dT - (dP/dT)v dV
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
            return Sres
        else:
            raise ValueError("Valores de T e P não fornecidos. Por favor, forneça valores válidos para T e P.")
            
    def ideal (self, substancia, T):
        # Cálculo da entropia ideal usando a equação de estado de Peng-Robinson
        eos = PR(substancia)
        Z = eos.calcZ(T, 101325)  # Pressão padrão de 1 atm
        A = eos.calcA(T, 101325)
        B = eos.calcB(T, 101325)
        # Cálculo da entropia ideal
        s_ideal = R * (Z - 1) - (A / (2 * np.sqrt(2) * B)) * np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B))
        return s_ideal