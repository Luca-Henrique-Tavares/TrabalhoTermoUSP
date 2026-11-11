# TrabalhoTermoUSP

Rotina computacional para cálculo de propriedades termodinâmicas de fluidos reais com a **equação de estado de Peng-Robinson**, desenvolvida para os exercícios computacionais de **PQI3420 – Termodinâmica I** (Escola Politécnica da USP, 2º semestre de 2026).

A rotina foi pensada para ser **genérica**: qualquer substância (com Tc, Pc, ω e cp de gás ideal) pode ser usada em quaisquer condições de T e P. As propriedades calculadas alimentam diretamente os balanços de energia e entropia de cada exercício.

---

## Sumário

- [Ideia central](#ideia-central)
- [Estrutura do projeto](#estrutura-do-projeto)
- [Módulos](#módulos)
  - [`constantes.py`](#constantespy)
  - [`substancias.py`](#substanciaspy)
  - [`eos.py`](#eospy)
  - [`propriedades.py`](#propriedadespy)
  - [`balancos.py`](#balancospy)
  - [`solvers.py`](#solverspy)
- [Notebooks dos exercícios](#notebooks-dos-exercícios)
- [Convenções](#convenções)
- [Como usar](#como-usar)
- [Validação](#validação)
- [Requisitos](#requisitos)

---

## Ideia central

Toda propriedade de um fluido real é escrita como **gás ideal + residual**:

$$
\Delta M = \Delta M^{ig} + M^R_2 - M^R_1 \qquad (M = U,\ H,\ S)
$$

- A **parte ideal** depende só da substância (cp) e dos estados (T, P).
- A **parte residual** mede o desvio do comportamento ideal e depende da equação de estado.

O código segue essa mesma separação: a substância, a equação de estado, as propriedades, os balanços e o solver ficam em módulos diferentes, e cada exercício só combina essas peças.

```
Substancia ──► PR (eos) ──► propriedades (ideal + residual) ──► balanços ──► resultado
                                     ▲
                                     └── solvers (quando uma T é incógnita)
```

---

## Estrutura do projeto

```
TrabalhoTermoUSP/
├── constantes.py        # constantes físicas e fatores de conversão
├── substancias.py       # classe Substancia + instâncias (etileno, n-pentano)
├── eos.py               # classe PR (Peng-Robinson): parâmetros, Z e V
├── propriedades.py      # entropia, entalpia e energia interna (ideal + residual)
├── balancos.py          # balanços de energia e de entropia
├── solvers.py           # bisseção e busca de T dado ΔH ou ΔS
├── ex_1.ipynb           # Questão 1
├── ex_2_config1.ipynb   # Questão 2 – Configuração I (válvula + trocador)
├── ex_2_config2.ipynb   # Questão 2 – Configuração II (turbina + trocador)
└── .gitignore
```

Dependências entre módulos (quem importa quem):

```
constantes ◄── eos ◄── propriedades ◄── solvers
     ▲                       ▲
     └──── propriedades      └── notebooks (importam todos)
substancias ◄── notebooks
balancos    ◄── notebooks
```

`constantes.py` não importa nada, o que evita importações circulares.

---

## Módulos

### `constantes.py`

Constantes usadas em todo o projeto, em SI.

| Nome | Valor | Uso |
|---|---|---|
| `R` | 8,3144 J/(mol·K) | constante dos gases |
| `BAR` | 10⁵ Pa | conversão de pressão: `P = 15 * BAR` |
| `kJ` | 10³ J | conversão de energia: `W = 10 * kJ` |

Conversões universais (que não dependem da substância) ficam aqui. Conversões que dependem da substância, como a de base molar para mássica, ficam na própria `Substancia`.

---

### `substancias.py`

Guarda os dados **fixos** de cada composto, que não mudam de um problema para outro.

**Classe `Substancia`**

| Atributo | Unidade | Descrição |
|---|---|---|
| `nome` | – | nome do composto |
| `tc` | K | temperatura crítica |
| `pc` | Pa | pressão crítica |
| `w` | – | fator acêntrico ω |
| `coef_cp` | – | coeficientes `[a, b, c, d]` de cp<sup>ig</sup> = a + bT + cT² + dT³, em J/(mol·K) |
| `mm` | g/mol | massa molar |

**Método**

- `Jmol_to_Jkg(grandeza)`: converte um valor por mol (J/mol) para por massa (J/kg), dividindo pela massa molar em kg/mol.

**Instâncias prontas:** `etileno` e `npentano`, criadas no próprio arquivo. Os outros módulos só importam:

```python
from substancias import etileno, npentano
```

T e P **não** ficam na substância: eles descrevem um ponto do processo, e a mesma substância aparece em vários pontos.

---

### `eos.py`

**Classe `PR`**: equação de estado de Peng-Robinson.

$$
P = \frac{RT}{V - b} - \frac{a(T)}{V(V+b) + b(V-b)}
$$

O objeto é criado uma vez por substância (`PR(etileno)`) e depois avaliado em quantos estados forem necessários. Por isso o construtor recebe só a substância; T e P entram como argumentos dos métodos.

**Atributos definidos no construtor**

| Atributo | Descrição |
|---|---|
| `substancia` | instância de `Substancia` |
| `sigma`, `epsilon` | constantes do denominador fatorado, −(1+√2) e −(1−√2) |
| `c` | fator 1/(2√2·b) da integral por frações parciais |

**Métodos**

| Método | Retorna | Expressão |
|---|---|---|
| `kappa()` | κ | 0,37464 + 1,54226ω − 0,26992ω² |
| `alpha(T)` | α(T) | [1 + κ(1 − √T<sub>r</sub>)]² |
| `dalpha_dT(T)` | dα/dT | −κ/√(T·Tc) · [1 + κ(1 − √T<sub>r</sub>)] (derivada feita à mão) |
| `ac()` | a<sub>c</sub> | 0,45724 R²Tc²/Pc (parte constante de a) |
| `a(T)` | a(T) | a<sub>c</sub> · α(T) |
| `b()` | b | 0,07780 RTc/Pc |
| `calcA(T, P)` | A | aP/(RT)² |
| `calcB(T, P)` | B | bP/(RT) |
| `calcZ(T, P)` | lista de Z | raízes reais da cúbica em Z |
| `filtrar_raizes(raizes)` | lista | mantém só raízes reais e positivas |
| `calcV(T, P)` | lista de V | V = ZRT/P para cada raiz |

**Como o Z é calculado.** Substituindo V = ZRT/P na PR, a equação vira uma cúbica em Z:

$$
Z^3 - (1-B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0
$$

As raízes são obtidas com `np.roots`. A cúbica é resolvida em Z (adimensional, da ordem de 1), e não em V, porque isso deixa o cálculo numericamente mais estável.

**Convenção das raízes:** `calcZ` sempre devolve uma **lista**:

- duas ou três raízes reais → `[Z_liquido, Z_vapor]` (menor e maior);
- uma raiz real → `[Z_unico]`.

Acima da temperatura crítica só existe uma fase, e a cúbica tem uma única raiz real.

---

### `propriedades.py`

Três classes com a mesma interface, uma por propriedade. Cada uma tem:

- `residual(eos, T, P, n, fase)`: propriedade residual M<sup>R</sup> = M − M<sup>ig</sup> no estado (T, P);
- `ideal(substancia, T1, T2, ...)`: variação de gás ideal entre dois estados.

O argumento `n` é a quantidade de matéria (use `n = 1` para base molar). O argumento `fase` é `"liquido"` ou `"vapor"` e escolhe a raiz de Z.

#### Dedução dos residuais

Os residuais partem das relações de Maxwell, integrando a T constante de V → ∞ (onde real e ideal coincidem) até o volume do estado:

$$
S^R = \int_{\infty}^{V} \left[\left(\frac{\partial P}{\partial T}\right)_V - \frac{R}{V}\right] dV + R\ln Z
\qquad
U^R = \int_{\infty}^{V} \left[T\left(\frac{\partial P}{\partial T}\right)_V - P\right] dV
$$

- O termo **R ln Z** na entropia corrige a referência: a integral compara real e ideal no **mesmo V**, mas o residual é definido na **mesma T e P**. Na energia interna esse termo não aparece, porque U<sup>ig</sup> não depende de V nem de P.
- O denominador da PR é fatorado (Bhaskara) em (V + (1−√2)b)(V + (1+√2)b), e a integral sai por **frações parciais**:

$$
I = \int_{\infty}^{V} \frac{dV}{V(V+b) + b(V-b)} = \frac{1}{2\sqrt{2}\,b}\,\ln\!\left[\frac{V + (1-\sqrt{2})b}{V + (1+\sqrt{2})b}\right]
$$

#### `entropia`

$$
S^R = R\ln\frac{V-b}{V} - a_c\,\frac{d\alpha}{dT}\, I + R\ln Z
$$

$$
\Delta S^{ig} = a\ln\frac{T_2}{T_1} + \sum_{i \ge 1} \frac{c_i}{i}\left(T_2^{\,i} - T_1^{\,i}\right) - R\ln\frac{P_2}{P_1}
$$

O termo i = 0 do cp vira logaritmo ao dividir por T, por isso é tratado separadamente no loop.

#### `energia_interna`

$$
U^R = a_c\left(\alpha - T\frac{d\alpha}{dT}\right) I
\qquad
\Delta U^{ig} = \int_{T_1}^{T_2} (c_p^{ig} - R)\, dT
$$

#### `entalpia`

$$
H^R = U^R + RT(Z - 1)
\qquad
\Delta H^{ig} = \sum_{i \ge 0} \frac{c_i}{i+1}\left(T_2^{\,i+1} - T_1^{\,i+1}\right)
$$

O `residual` da entalpia reaproveita o de energia interna, em vez de repetir a fórmula. Para gás ideal, U e H dependem só de T, então as partes ideais de U e H não têm termo de pressão.

---

### `balancos.py`

Duas funções que **isolam a incógnita** de um balanço. A incógnita é o argumento passado como `None`; termos ausentes são passados como `0`.

**`balanco_energia(dU, dH, Q, W, m)`**

$$
\Delta U = m\,\Delta H_{1\text{-}2} + Q + W
$$

**`balanco_entropia(dStt, Sgen, Q, T, m, dS)`**

$$
\Delta S_{tt} = m\,\Delta S_{1\text{-}2} + \frac{Q}{T} + S_{gen}
$$

Comportamento:

- **exatamente um `None`** → calcula e devolve essa variável;
- **nenhum ou mais de um `None`** → erro (um balanço é uma equação e resolve uma incógnita);
- divisão por `m = 0` é verificada antes de acontecer.

Como cada caso físico é escolhido pelos termos passados, não existe parâmetro de "tipo de problema":

- **sistema fechado** → `m = 0`;
- **regime permanente** → `dU = 0` (ou `dStt = 0`).

O `T` do balanço de entropia é a **temperatura da fronteira** por onde o calor passa.

---

### `solvers.py`

Usado quando uma temperatura de saída é desconhecida e determinada por uma condição sobre H ou S (válvula isentálpica, turbina isentrópica).

**`bissecao(f, a, b, tol=1e-6)`**

Método da bisseção genérico: acha x em [a, b] tal que f(x) = 0.

1. Verifica que f(a) e f(b) têm sinais opostos; senão, levanta erro.
2. Calcula o ponto médio c e descarta a metade do intervalo sem troca de sinal.
3. Repete enquanto (b − a)/2 > `tol`. Como a raiz está sempre dentro de [a, b], (b − a)/2 é um limite garantido para o erro do ponto médio.
4. Devolve o ponto médio do intervalo final.

**`TporH(eos, T1, P1, P2, n, dH_alvo, T_min, T_max, fase)`**

Acha T₂ tal que H(T₂, P₂) − H(T₁, P₁) = `dH_alvo`.

**`TporS(eos, T1, P1, P2, n, dS_alvo, T_min, T_max, fase)`**

Acha T₂ tal que S(T₂, P₂) − S(T₁, P₁) = `dS_alvo`.

As duas seguem a mesma estrutura:

- o residual da entrada, que não depende do T de tentativa, é calculado **uma vez**, fora da busca;
- uma função interna `f(T)` (closure) soma a parte ideal de T₁ até T e o residual da saída em T, subtrai o residual da entrada e o alvo;
- `f` é entregue à `bissecao`, que a avalia em vários T até achar o zero.

**Escolha do intervalo [T_min, T_max]:** o limite inferior é a temperatura crítica (garante fase única em todo o intervalo) e o superior é T₁ (válvula e turbina resfriam o gás). A checagem de sinais da bisseção confirma que a raiz está no intervalo.

---

## Notebooks dos exercícios

Cada notebook contém o enunciado, as hipóteses, os balanços em markdown e as chamadas aos módulos.

### `ex_1.ipynb` – Questão 1

Compressão isotérmica de 5 mol de n-pentano (500 K, 3 → 30 bar) em sistema fechado, com W = 10 kJ fornecido e vizinhança a 298 K.

1. Z nos dois estados (fase única: T > Tc).
2. ΔU = ΔU<sup>ig</sup> + U<sup>R</sup>₂ − U<sup>R</sup>₁.
3. Balanço de energia (sistema fechado) → Q.
4. ΔS = ΔS<sup>ig</sup> + S<sup>R</sup>₂ − S<sup>R</sup>₁.
5. Balanço de entropia com Q e T<sub>vizinhança</sub> → S<sub>gen</sub>.
6. Classificação pelo sinal de S<sub>gen</sub>.

### `ex_2_config1.ipynb` – Questão 2, Configuração I

Etileno de 480 K e 50 bar até 325 K e 15 bar, com válvula + trocador.

1. **Válvula** (regime permanente, Q = 0, W = 0) → ΔH = 0 → `TporH` para T₂.
2. **Trocador** (W = 0) → Q = ΔH de T₂ a 325 K.
3. Conversão para base mássica.

### `ex_2_config2.ipynb` – Questão 2, Configuração II

Mesmos estados de entrada e saída, com turbina + trocador.

1. **Turbina** (adiabática e reversível → S<sub>gen</sub> = 0, Q = 0) → ΔS = 0 → `TporS` para T₂.
2. Balanço de energia na turbina → W = ΔH.
3. **Trocador** → Q.
4. Conversão para base mássica.

Como os estados inicial e final são os mesmos nas duas configurações e H é função de estado, **Q + W da configuração II é igual ao Q da configuração I**. Na configuração II, parte da energia retirada sai como trabalho útil.

---

## Convenções

- **Unidades:** SI em todo o código (K, Pa, J, mol). Conversões só na entrada (`BAR`, `kJ`) e na saída (`Jmol_to_Jkg`, `/ kJ`).
- **Sinais:** Q e W positivos quando **entram** no sistema.
- **ΔH₁₋₂ no balanço de energia:** H<sub>entrada</sub> − H<sub>saída</sub>. Por isso, nos notebooks, a variação calculada pelas propriedades (saída − entrada) tem o sinal invertido antes de entrar no balanço.
- **Residual:** M<sup>R</sup> = M<sub>real</sub> − M<sup>ig</sup>, na mesma T e P.
- **Fase:** `"liquido"` usa a menor raiz de Z; `"vapor"`, a maior.

---

## Como usar

Os notebooks devem estar na mesma pasta dos módulos. Primeira célula recomendada:

```python
%load_ext autoreload
%autoreload 2
```

Exemplo de cálculo de ΔH de um estado a outro:

```python
from substancias import etileno
from eos import PR
from propriedades import entalpia
from constantes import BAR

eos = PR(etileno)

T1, P1 = 480, 50 * BAR
T2, P2 = 400, 15 * BAR

H_R1 = entalpia().residual(eos, T1, P1, 1, "vapor")
H_R2 = entalpia().residual(eos, T2, P2, 1, "vapor")
dH_ig = entalpia().ideal(etileno, T1, T2, P1, P2, 1)

dH = dH_ig + H_R2 - H_R1                     # J/mol
dH_massico = etileno.Jmol_to_Jkg(dH)         # J/kg
```

Exemplo de temperatura de saída de uma expansão isentrópica:

```python
from solvers import TporS

T2 = TporS(eos, T1, P1, P2, 1, 0, etileno.tc, T1, "vapor")
```

Para adicionar uma substância, basta criar uma nova instância de `Substancia` em `substancias.py` com Tc, Pc, ω, coeficientes de cp e massa molar.

---

## Validação

- **Derivada dα/dT:** comparada com diferença finita central sobre o próprio método `alpha` e com derivação simbólica (SymPy).
- **Residuais S<sup>R</sup> e U<sup>R</sup>:** as fórmulas analíticas foram comparadas com a integração numérica (`scipy.integrate.quad`, de ∞ até V) dos integrandos originais da PR. Os resultados coincidem.
- **Sinais esperados:** para gás, U<sup>R</sup> e H<sup>R</sup> são negativos (forças atrativas).
- **Solver:** testado com f(x) = x² − 2 (raiz √2) e com teste de ida e volta (o T recuperado coincide com o T usado para gerar o alvo).
- **Consistência física:** na Questão 2, Q + W da configuração II reproduz o Q da configuração I.

---

## Requisitos

- Python 3
- `numpy`
- `scipy` (só para validação com `quad`)
- `sympy` (só para conferir a derivada)
- Jupyter / `ipykernel` para os notebooks
