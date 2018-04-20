from sympy.abc import *
from sympy import *



# ИНИЦИАЛИЗАЦИЯ СИМВОЛОВ

# Время
t = symbols('t') 
# Индексы
i, P, C = symbols('i, P, C', cls=Idx) 
# (проекция вилки на e_z, проекция вилки на e_wheel, радиус, масса)
h, d, r = symbols('h,d,r') 
# Тензор инерции
# J = symbols('J') # J = IndexedBase('J') 
# Углы
psi=IndexedBase('psi');
beta = IndexedBase('beta');
alpha = symbols('alpha');
theta = IndexedBase('theta');
#моменты инерции
W = IndexedBase('W'); T = IndexedBase('T') 
# элементы матрицы моментов инерций для разных тел (ВРЕМЕННО)
a,b,c = symbols('a,b,c')


#СЛОВАРИ И ИХ СОДЕРЖИМОЕ  пример словаря: {'Аоексей': 21}

m = {}                                # словарь для масс
e, omega, v = {}, {}, {}              # вектора, омега, скорость
S, P, C, D  = {}, {}, {}, {}          # точки
eq = {}                               # словарь со всякими выражениями
delta = {};                           
nu = {}                               # псевдоскорости nu[1] и nu[2]
A = {}                                # части уравнения Д'Аламбера лагранжа для разных тел
velocity = {}                         # скорость
J = {}                                # моменты инерции для разных тел (осевые?)
K = {}                                # кинетический момент для тела K = I_z*omega
delta_r = {}                          # виртуальное перемещение для центра масс?
M = {}                                # момент сил
F = {}                                # силы, действующие на тело
omega_delta = {}                      # виртуальны поворот
coeff = {}                            # коэффиценты у Д'Аламюера-Лагранжа

# Виртуальные перемещения
delta['x'], delta['y'], delta['alpha'] = symbols('delta_x, delta_y, delta_alpha')
delta['theta'] = IndexedBase('delta_theta') 
delta['psi'] = IndexedBase('delta_psi');
# Псевдоскорости (взамен ẋ и ẏ)
nu[1], nu[2] = symbols('nu1, nu2') 
m['platform'], m['wheel'] = symbols('m1, m2')


#ЗАВИСИМОСТИ

x = x(t); y = y(t); alpha = alpha(t);
nu[1] = nu[1](t); nu[2] = nu[2](t);













