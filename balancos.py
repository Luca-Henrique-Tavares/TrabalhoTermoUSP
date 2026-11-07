def balanco_energia (dU,dH,Q,W,m): 
    # Se quiser base molar, basta passar m = 1. Se quiser base de massa, basta passar m = massa do sistema.
    # Usamos na disciplina que dU = mdH1-2 + Q + W
    # Se for transiente, vamos primeiro analisar qua(is) é(são) a(s) variável(is) que queremos calcular. Se for dU, dH, Q ou W, vamos isolar a variável desejada na equação de balanço de energia.
    # Vamos ver quais são none: 
    valores = { "dU": dU, "dH": dH, "Q": Q, "W": W, "m": m }
    incognita = []
    for variavel, valor in valores.items():
        if valor is None:
            incognita += [variavel]
        else:
            continue
    if len(incognita) > 1:
        raise ValueError(f"As variáveis {incognita} estão faltando. Por favor, forneça valores para todas as variáveis, exceto a que deseja calcular.")
    elif len(incognita) == 0:
        raise ValueError("Todas as variáveis fornecidas. Por favor, forneça valores para todas as variáveis, exceto a que deseja calcular.")
    else:
        incognita = incognita[0]
        if incognita == "dU":
            dU = m*dH + Q + W
            return dU
        elif incognita == "dH":
            if m == 0:
                raise ValueError("A massa é zero, portanto não é possível calcular dH. Por favor, forneça um valor de massa diferente de zero.")
            else:
                dH = (dU - Q - W)/m
                return dH
        elif incognita == "Q":
            Q = dU - m*dH - W
            return Q
        elif incognita == "W":
            W = dU - m*dH - Q
            return W
        else:
            raise ValueError("Variável desconhecida. Por favor, forneça valores para todas as variáveis, exceto a que deseja calcular.")

def balanco_entropia (dStt, Sgen, Q, T, m, dS):
    # Se quiser base molar, basta passar m = 1. Se quiser base de massa, basta passar m = massa do sistema.
    # Usamos na disciplina que dStt = Sgen + Q/T + mdS1-2
    # Se for transiente, vamos primeiro analisar qua(is) é(são) a(s) variável(is) que queremos calcular. Se for dStt, Sgen, Q ou T, vamos isolar a variável desejada na equação de balanço de entropia.
    # Vamos ver quais são none: 
    valores = { "dStt": dStt, "Sgen": Sgen, "Q": Q, "T": T, "m": m, "dS": dS }
    incognita = []
    for variavel, valor in valores.items():
        if valor is None:
            incognita += [variavel]
        else:
            continue
    if len(incognita) > 1:
        raise ValueError(f"As variáveis {incognita} estão faltando. Por favor, forneça valores para todas as variáveis, exceto a que deseja calcular.")
    elif len(incognita) == 0:
        raise ValueError("Todas as variáveis fornecidas. Por favor, forneça valores para todas as variáveis, exceto a que deseja calcular.")
    else:
        incognita = incognita[0]  
        if incognita == "dStt":
            dStt = Sgen + Q/T + m*dS
            return dStt
        elif incognita == "Sgen":
            Sgen = dStt - Q/T - m*dS
            return Sgen
        elif incognita == "Q":
            Q = (dStt - Sgen - m*dS)*T
            return Q
        elif incognita == "T":
            T = Q/(dStt - Sgen - m*dS)
            return T
        elif incognita == "dS":
            if m == 0:
                raise ValueError("A massa é zero, portanto não é possível calcular dS. Por favor, forneça um valor de massa diferente de zero.")
            else:
                dS = (dStt - Sgen - Q/T)/m
                return dS
        else:
            raise ValueError("Variável desconhecida. Por favor, forneça valores para todas as variáveis, exceto a que deseja calcular.")
