from propriedades import entalpia,entropia

def TporH(eos, T1, P1, P2, n, dH_alvo, T_min, T_max, fase):
    # Essa função será a responsável por calcular a função f(T) para entalpia
    # Essa função chamará o método da bisseção para encontrar a temperatura que corresponde a uma entalpia específica. 
    hres1 = entalpia().residual(eos, T1, P1, n, fase) # Calculando a entalpia residual na condição inicial
    def f(T):
        hideal = entalpia().ideal(eos.substancia, T1, T, P1, P2, n) # Calculando a entalpia ideal
        hres2 = entalpia().residual(eos, T, P2, n, fase) # Calculando a entalpia residual na condição final
        return hres2 + hideal - hres1 - dH_alvo 
    return bissecao(f, T_min, T_max) # Chamando o método da bisseção para encontrar a raiz da função f(T) = 0

def TporS(eos, T1, P1, P2, n, dS_alvo, T_min, T_max, fase):
    # Essa função será a responsável por calcular a função f(T) para entropia
    # Essa função chamará o método da bisseção para encontrar a temperatura que corresponde a uma entropia específica. 
    sres1 = entropia().residual(eos, T1, P1, n, fase) # Calculando a entropia residual na condição inicial
    def f(T):
        sideal = entropia().ideal(eos.substancia, T1, T, P1, P2, n) # Calculando a entropia ideal 
        sres2 = entropia().residual(eos, T, P2, n, fase) # Calculando a entropia residual na condição final
        return sres2 + sideal - sres1 - dS_alvo 
    return bissecao(f, T_min, T_max) # Chamando o método da bisseção para encontrar a raiz da função f(T) = 0
    
def bissecao(f, a, b, tol=1e-6):
    # Essa função será a responsável por implementar o método da bisseção para encontrar a raiz da função f(T) = 0.
    # A função biss é uma implementação do método da bisseção para encontrar a raiz de uma função contínua f(x) em um intervalo [a,b] onde f(a) e f(b) têm sinais opostos.
    if f(a) * f(b) >= 0:
        raise ValueError("A função deve ter sinais opostos em a e b.")
    fa = f(a)
    while (b - a) / 2 > tol:
        c = (a + b) / 2
        if f(c) == 0:
            return c
        elif fa * f(c) < 0:
            b = c
        else:
            a = c
            fa = f(c)
    return (a + b) / 2
