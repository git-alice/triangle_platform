#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 14:44:31 2018

@author: eval
"""

from sympy.abc import *
from sympy import *


# ИНИЦИАЛИЗАЦИЯ СИМВОЛОВ

# Время
t = symbols('t') 
# Индексы
# i, P, C = symbols('i, P, C', cls=Idx) 
# (проекция вилки на e_z, проекция вилки на e_wheel, радиус, масса)
h, d, r = symbols('h,d,r') 
# Тензор инерции
# J = symbols('J') # J = IndexedBase('J') 
# Углы
psi = symbols('psi');
beta = 0;
alpha = symbols('alpha');
theta = symbols('theta');
#управляющие моменты
W, T = symbols ('W, T');
# элементы матрицы моментов инерций для разных тел (ВРЕМЕННО)
a,b,c = symbols('a,b,c')
J0, J1, J3 = symbols('J0, J1, J3')

# Виртуальные перемещения
delta_x, delta_y, delta_alpha = symbols('delta_x, delta_y, delta_alpha')
delta_theta = symbols('delta_theta') 
delta_psi = symbols('delta_psi');
nu1, nu2 = symbols('nu1, nu2') 
dnu1, dnu2, ddalpha = symbols('dnu1, dnu2, ddalpha') 
m1, m2 = symbols('m1, m2')

#ЗАВИСИМОСТИ
x = x(t); y = y(t); alpha = alpha(t);theta = theta(t); psi = psi(t);
nu1 = nu1(t); nu2 = nu2(t);

#СТРОЕНИЕ ВЕКТОРОВ
e_x = Matrix([1,0,0])
e_y = Matrix([0,1,0])
e_z = Matrix([0,0,1])
e_xi = Matrix([cos(alpha),sin(alpha),0])
e_eta = Matrix([-sin(alpha),cos(alpha),0])

n_wheel = -sin(theta)*e_xi + cos(theta)*e_eta # нормаль к колесу
e_wheel = cos(theta)*e_xi + sin(theta)*e_eta # e_z и нормаль к колесу дополняет до правой тройки?


#СТРОЕНИЕ ПЛАТФОРМЫ НА КОЛЁСАХ
SP = a* (cos(beta)*e_xi + sin(beta)*e_eta)
PC = d*e_wheel - h*e_z
# check!!!
CP = -PC  
PS = -SP
DC = r*e_z

omega_platform= Derivative(alpha,t)*e_z
omega_fork    = omega_platform + Derivative(theta,t)*e_z
omega_wheel = omega_fork + Derivative(psi,t)*n_wheel

v_c = omega_wheel.cross(DC)
v_p = v_c + omega_fork.cross(CP)
v_s = v_p + omega_platform.cross(PS)

eq_deltax = delta_x - v_s[0].subs(Derivative(theta,t),delta_theta).subs(Derivative(alpha,t),delta_alpha).subs(Derivative(psi,t),delta_psi)
eq_deltay = delta_y - v_s[1].subs(Derivative(theta,t),delta_theta).subs(Derivative(alpha,t),delta_alpha).subs(Derivative(psi,t),delta_psi)
sols_delta = solve([eq_deltax,eq_deltay],[delta_theta,delta_psi], dict = True)
delta_theta = sols_delta[0][delta_theta];
delta_psi = sols_delta[0][delta_psi];

eq_nu1 = nu1 - v_s.dot(e_xi)
eq_nu2 = nu2 - v_s.dot(e_eta)
sols_nu = solve([eq_nu1,eq_nu2],[Derivative(theta,t),Derivative(psi,t)], dict = True)
dtheta = sols_nu[0][Derivative(theta,t)];
dpsi = sols_nu[0][Derivative(psi,t)];

#print(dtheta)
#print(dpsi)
# платформа
acc_s = Matrix([diff(nu1*e_xi[0] + nu2*e_eta[0],t),diff(nu1*e_xi[1] + nu2*e_eta[1],t),0])
delta_rs = delta_x * e_x + delta_y * e_y
dks = J0 * diff(alpha*e_z[2],t,2)
dp_term1 = m1 * acc_s.dot(delta_rs)
dp_term2 = (dks + T)*delta_alpha
dp = Poly(dp_term1 + dp_term2, [delta_x,delta_y,delta_alpha])

#вилка
df_term2 = T*(delta_alpha+delta_theta)
#print(delta_alpha+delta_theta);


# колесо
v_c1 = Matrix([0,0,0])
v_c1[0] = v_c[0].subs(Derivative(psi,t),dpsi).subs(Derivative(theta,t),dtheta);
v_c1[1] = v_c[1].subs(Derivative(psi,t),dpsi).subs(Derivative(theta,t),dtheta);
acc_c = Matrix([diff(v_c1[0],t),diff(v_c1[1],t),0])
delta_rc = Matrix([0,0,0])
delta_rc[0] = v_c[0].subs(Derivative(psi,t),delta_psi).subs(Derivative(theta,t),delta_theta).subs(Derivative(alpha,t),delta_alpha);
delta_rc[1] = v_c[1].subs(Derivative(psi,t),delta_psi).subs(Derivative(theta,t),delta_theta).subs(Derivative(alpha,t),delta_alpha);
acc_c[0] = acc_c[0].subs(Derivative(psi,t),dpsi).subs(Derivative(theta,t),dtheta);
acc_c[1] = acc_c[1].subs(Derivative(psi,t),dpsi).subs(Derivative(theta,t),dtheta);
dw_term1 = m2 * acc_c.dot(delta_rc)
#print(v_c[0]);
#print(diff(v_c1[0],t));
#print(acc_c[0]);
K = J1*omega_wheel.dot(e_wheel)*e_wheel + J3*omega_wheel.dot(n_wheel)*n_wheel + J1*omega_wheel.dot(e_z)*e_z
Ksubs = Matrix([0,0,0])
Ksubs[0] = K[0].subs(Derivative(psi,t),dpsi).subs(Derivative(theta,t),dtheta);
Ksubs[1] = K[1].subs(Derivative(psi,t),dpsi).subs(Derivative(theta,t),dtheta);
Ksubs[2] = K[2].subs(Derivative(psi,t),dpsi).subs(Derivative(theta,t),dtheta);
omega_delta_wheel = Matrix([0,0,0])
omega_delta_wheel[0] = omega_wheel[0].subs(Derivative(psi,t),delta_psi).subs(Derivative(theta,t),delta_theta).subs(Derivative(alpha,t),delta_alpha);
omega_delta_wheel[1] = omega_wheel[1].subs(Derivative(psi,t),delta_psi).subs(Derivative(theta,t),delta_theta).subs(Derivative(alpha,t),delta_alpha);
omega_delta_wheel[2] = omega_wheel[2].subs(Derivative(psi,t),delta_psi).subs(Derivative(theta,t),delta_theta).subs(Derivative(alpha,t),delta_alpha);
dK = trigsimp(diff(Ksubs,t) - W*e_z)
#print(dK)
#dw_term2 = trigsimp(dK.dot(omega_delta_wheel))
#print(dK)
#print('dw-term2 done')
#print(trigsimp(dw_term2));


#d = Poly(dp_term1 + dp_term2 + df_term2 + dw_term1, [delta_x,delta_y,delta_alpha])
d = Poly(trigsimp(dp_term1 + dp_term2 + df_term2 + dw_term1 + dw_term2), [delta_x,delta_y,delta_alpha])
#d = Poly((dp_term1 + dp_term2 + df_term2 + dw_term1 + dw_term2), [delta_x,delta_y,delta_alpha])
print('d  - done');



eq_dal_delta_x = d.coeff_monomial(delta_x).subs(Derivative(nu1,t),dnu1).subs(Derivative(nu2,t),dnu2).subs(Derivative(alpha,t,t),ddalpha);
eq_dal_delta_y = d.coeff_monomial(delta_y).subs(Derivative(nu1,t),dnu1).subs(Derivative(nu2,t),dnu2).subs(Derivative(alpha,t,t),ddalpha);
eq_dal_delta_alpha = d.coeff_monomial(delta_alpha).subs(Derivative(nu1,t),dnu1).subs(Derivative(nu2,t),dnu2).subs(Derivative(alpha,t,t),ddalpha);

print('Eq1')
print(eq_dal_delta_x.subs([(alpha,0), (theta,0), (nu1,0), (nu2,0), (Derivative(alpha,t),0)]));

print('Eq2')
print(eq_dal_delta_y.subs([(alpha,0), (theta,0), (nu1,0), (nu2,0), (Derivative(alpha,t),0)]));

print('Eq3')
print(eq_dal_delta_alpha.subs([(alpha,0), (theta,0), (nu1,0), (nu2,0), (Derivative(alpha,t),0)]));

EqMatrix, EqRHS = linear_eq_to_matrix([eq_dal_delta_x,eq_dal_delta_y,eq_dal_delta_alpha],[dnu1,dnu2,ddalpha]);
print('begin to invert')
EqS = EqMatrix.inv()
print('finished')
print(EqS[0,0]);

#dyn_eq = solve([eq_dal_delta_x,eq_dal_delta_y,eq_dal_delta_alpha],[Derivative(nu1,t), Derivative(nu2,t), Derivative(alpha,t,2)], dict = True)
#dyn_eq = solve(, dict = True)
#eq_nu1 = dyn_eq[0][Derivative(nu1,t)]
#eq_nu2 = dyn_eq[0][Derivative(nu2,t)]
#eq_ddalpha = dyn_eq[0][Derivative(alpha,t,2)]
#eq_nu1 = dyn_eq[0][dnu1]
#eq_nu2 = dyn_eq[0][dnu2]
#eq_ddalpha = dyn_eq[0][ddalpha]
#print(eq_nu1);
#print(trigsimp(eq_nu1));
print('end');
#print(eq_nu2);
#print(eq_ddalpha);
