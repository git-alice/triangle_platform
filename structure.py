from sympy.abc import *
from sympy import *
from sympy.tensor import *

from variables import *

#СТРОЕНИЕ ВЕКТОРОВ
e['x'] = Matrix([1,0,0])
e['y'] = Matrix([0,1,0])
e['z'] = Matrix([0,0,1])
e['xi'] = Matrix([cos(alpha),sin(alpha),0])
e['eta'] = Matrix([-sin(alpha),cos(alpha),0])

n_wheel = lambda i: -sin(theta[i])*e['xi'] + cos(theta[i])*e['eta'] # нормаль к колесу
e_wheel = lambda i:  cos(theta[i])*e['xi'] + sin(theta[i])*e['eta'] # e_z и нормаль к колесу дополняет до правой тройки?


#СТРОЕНИЕ ПЛАТФОРМЫ НА КОЛЁСАХ
SP = lambda i: cos(beta[i]) * e['xi'] + sin(beta[i])*e['eta']
PC = lambda i: d*e_wheel(i) - h*e['z']
CD = lambda i: -r*e['z']
DC = lambda i: r*e['z']

#ТОЧКИ И ИХ ПРИНАДЛЕЖНОСТЬ
points = { # S, P, C, D - словари
	'S': S,
	'P': P,
	'C': C,
	'D': D
}
S['where'] = ['platform'];          S['name'] = 'S'
P['where'] = ['platform', 'fork'];  P['name'] = 'P'
C['where'] = ['fork', 'wheel'];     C['name'] = 'C'
D['where'] = ['wheel', 'floor'];    D['name'] = 'D'
