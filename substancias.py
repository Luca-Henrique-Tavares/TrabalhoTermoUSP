class Substancia: 
    def __init__(self,nome,tc,pc,w,coef_cp):
        self.nome = nome
        self.tc = tc # tc eh a temperatura critica do composto em K
        self.pc = pc #pc eh a pressao critica do composto em Pa
        self.w = w  # w eh o fator acentrico do composto
        self.coef_cp = coef_cp # coef_cp eh uma tupla com os coeficientes do cp

etileno = Substancia(nome="Etileno",tc=305.2,pc=4.85*10**6,w=0.098,coef_cp=[3.95,15.628*10**-2,-8.339*10**-5,17.657*10**-9])
npentano = Substancia(nome="N-pentano",tc=469.7,pc=3.36*10**6,w=0.251,coef_cp=[6.77,45.398*10**-2,-22.448*10**-5,42.259*10**-9])
