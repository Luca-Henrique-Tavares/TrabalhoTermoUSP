# Primeiro, vamos analizar o enunciado do problema:
# 1- Um conjunto cilindro-pistão contém inicialmente 5 mol de n-pentano a 500 K e 3 bar. O gás
# é comprimido isotermicamente até 30 bar, com um trabalho fornecido ao sistema de 10 kJ.
# Durante esse processo, ocorre transferência de calor com a vizinhança, que se encontra a uma
# temperatura de 298 K.
# Determine se esse processo é reversível, irreversível ou impossível. Admita que o n-pentano
# obedeça à equação de estado de Peng-Robinson.

from substancias import npentano
from eos import PR

# Primeira coisa: Temos T1, T2, P1, P2, n, W.
# Queremos saber se o processo é reverssível, irreversível ou impossível. (Sgen = ?)
# Vamos calcular o Z 