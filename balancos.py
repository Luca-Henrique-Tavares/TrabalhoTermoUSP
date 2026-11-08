def balanco_energia (dU,dH,Q,W): 
    # Usamos na disciplina que dU = mdH1-2 + Q + W
    # Vamos primeiro analisar qua(is) é(são) a(s) variável(is) que queremos calcular. Se for dU, dH, Q ou W, vamos isolar a variável desejada na equação de balanço de energia.
    # Vamos ver quais são none: 
    valores = { "dU": dU, "dH": dH, "Q": Q, "W": W}
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
            dU = dH + Q + W
            return dU
        elif incognita == "dH":
            dH = (dU - Q - W)
            return dH
        elif incognita == "Q":
            Q = dU - dH - W
            return Q
        elif incognita == "W":
            W = dU - dH - Q
            return W
        else:
            raise ValueError("Variável desconhecida. Por favor, forneça valores para todas as variáveis, exceto a que deseja calcular.")

def balanco_entropia (dStt, Sgen, Q, T, dS):
    # Usamos na disciplina que dStt = Sgen + Q/T + mdS1-2
    # Vamos primeiro analisar qua(is) é(são) a(s) variável(is) que queremos calcular. Se for dStt, Sgen, Q ou T, vamos isolar a variável desejada na equação de balanço de entropia.
    # Vamos ver quais são none: 
    valores = { "dStt": dStt, "Sgen": Sgen, "Q": Q, "T": T,"dS": dS }
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
            dStt = Sgen + Q/T + dS
            return dStt
        elif incognita == "Sgen":
            Sgen = dStt - Q/T - dS
            return Sgen
        elif incognita == "Q":
            Q = (dStt - Sgen - dS)*T
            return Q
        elif incognita == "T":
            T = Q/(dStt - Sgen - dS)
            return T
        elif incognita == "dS":
            dS = (dStt - Sgen - Q/T)
            return dS
        else:
            raise ValueError("Variável desconhecida. Por favor, forneça valores para todas as variáveis, exceto a que deseja calcular.")
