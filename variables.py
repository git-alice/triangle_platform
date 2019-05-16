import sympy.abc 
from sympy import symbols, IndexedBase, Function
from sympy.physics.vector import dynamicsymbols


# ИНИЦИАЛИЗАЦИЯ СИМВОЛОВ

# Время
t = symbols('t')
# Индексы
# i, P, C = symbols('i, P, C', cls=Idx) # i, P и С ???
# (проекция вилки на e_z, проекция вилки на e_wheel, радиус, масса)
h, d, r = symbols('h,d,r')
# Тензор инерции
# J = symbols('J') # J = IndexedBase('J')
# Углы
psi0, psi1, psi2 = dynamicsymbols('psi0, psi1, psi2')
beta0, beta1, beta2 = symbols('beta0, beta1, beta2')
theta0, theta1, theta2 = dynamicsymbols('theta0, theta1, theta2')
psi = [psi0, psi1, psi2]
beta = [beta0, beta1, beta2]
theta = [theta0, theta1, theta2]

#моменты инерции
W1, W2, W3 = dynamicsymbols('W1, W2, W3'); T1, T2, T3 = dynamicsymbols('T1, T2, T3')
T = {
    1: T1,
    2: T2,
    3: T3
}

W = {
    1: W1,
    2: W2,
    3: W3
}
# элементы матрицы моментов инерций для разных тел (ВРЕМЕННО)
a,b,c = symbols('a,b,c')


#СЛОВАРИ И ИХ СОДЕРЖИМОЕ  пример словаря: {'Алексей': 21}

m = {}                                # словарь для масс
e, omega, v = {}, {}, {}              # вектора, омега, скорость
S, P, C, D  = {}, {}, {}, {}          # точки
points = {}
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

# Думаю здесь проблем не должно возникнуть, НО МАЛО ЛИ
# Виртуальные перемещения
delta['x'], delta['y'], delta['alpha'] = symbols('delta_x, delta_y, delta_alpha')
delta_theta0, delta_theta1, delta_theta2 = dynamicsymbols('delta_theta0, delta_theta1, delta_theta2')
delta_psi0, delta_psi1, delta_psi2 = dynamicsymbols('delta_psi0, delta_psi1, delta_psi2')
delta['theta'] = [delta_theta0, delta_theta1, delta_theta2]
delta['psi'] = [delta_psi0, delta_psi1, delta_psi2];
# Псевдоскорости (взамен ẋ и ẏ)
m['platform'], m['wheel'] = symbols('m1, m2')


#ЗАВИСИМОСТИ
x, y, alpha = dynamicsymbols('x, y, alpha')
nu[1], nu[2] = dynamicsymbols('nu1, nu2')














