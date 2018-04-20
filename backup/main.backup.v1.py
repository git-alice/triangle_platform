import latex 
import sympy.printing as printing
import sys
from sympy import *
from termcolor import colored, cprint
import numpy as np
import tqdm
import time


def scalar_product(a, b):
    return(np.dot(a,b))


def connections():





#init

ddalpha,dalpha, dnu1,dnu2,nu1,nu2,alpha,m,t,dx,dy,delta_alpha,J1,d,a,r = symbols('ddalpha,dalpha,dnu1,dnu2,nu1,nu2,alpha,m,t,dx,dy,delta_alpha,J1,d,a,r')
k = symbols('k', integer=True)
i = symbols('i', cls=Idx)
c = symbols('c', cls=Idx)
W = IndexedBase("W")
beta = IndexedBase('beta')
theta = IndexedBase('theta')
J2 = IndexedBase('J2')
T = IndexedBase('T')
dtheta = IndexedBase('dtheta')
dpsi=IndexedBase('dpsi')
ddpsi=IndexedBase('ddpsi')
F = IndexedBase('F')
a = symbols('a')
dr = symbols('dr')
dK = symbols('dK')
M = symbols('M')
n = symbols('n')
omega = symbols('omega')

    # delta_alpha - дельта альфа
    # dalpha - производная альфа по времени


# K = 2*c(t) + c(t)**2
# pprint(K)
# n = 3
pprint(Derivative(K,t).doit()) #.subs(Derivative(c(t),t),0)
Dalamber = scalar_product((m*a - summation(F[i],(i,1,n))), dr) + scalar_product(dK - M, omega)
pprint(Dalamber)


# for i in tqdm.tqdm(range(1000)):
    # time.sleep(0.01)
    # or other long operations





#СВЯЗИ v1
R_1_v1 = -nu1*cos(theta[i]) - nu2*sin(theta[i]) - a*sin(theta[i] - beta[i])*dalpha   -r*dpsi[i]
R_2_v1 = -dalpha*cos(beta[i] - theta[i]) - nu2*cos(theta[i]) + nu1*sin(theta[i])     -d*dtheta[i] 
#СВЯЗИ v2
R_1_v2 = -r*dpsi[i]*sin(theta[i]) - r*dtheta[i]*cos(theta[i]) - a*dalpha*cos(beta[i]) - nu2
R_2_v2 = -r*dpsi[i]*cos(theta[i]) - r*dtheta[i]*sin(theta[i]) - a*dalpha*sin(beta[i]) - nu1



A1 =  m*(
        (dnu1*cos(alpha) - nu1*sin(alpha)*dalpha - dnu2*sin(alpha) - nu2*cos(alpha)*dalpha)*dx + 
        (dnu1*sin(alpha) + nu1*cos(alpha)*dalpha + dnu2*cos(alpha) - nu2*sin(alpha)*dalpha)*dy  + 
        J1*ddalpha*delta_alpha
        );
A2 = summation(
              W[i]*(
                   delta_alpha + 
                   d**(-1)*(-delta_alpha*cos(beta[i] - theta[i]) - nu2*cos(theta[i]) + nu1*sin(theta[i]))
                   )
              ,(i,0,2));


# Dalamber.subs(a,-r*ddpsi[i]*())
A3_1 = summation(
            (-dnu1*cos(theta[i]))
            , (i,0,0));
A3_2 = summation((nu1*sin(theta[i])*d**(-1)*(-dalpha*cos(beta[i] - theta[i]) - nu2*cos(theta[i]) + nu1*sin(theta[i]))), (i,0,0));
A3_3 = summation((-diff(nu2,t)*sin(theta[i])), (i,0,0));
A3_4 = summation((nu2*cos(theta[i])*d**(-1)*(-dalpha*cos(beta[i] - theta[i]) - nu2*cos(theta[i]) + nu1*sin(theta[i]))), (i,0,0));
A3_5 = summation((-a*sin(theta[i] - beta[i])*ddalpha - a*cos(theta[i]- beta[i])*dalpha*dalpha), (i,0,0));
A3_6 = summation(((-dx*cos(alpha) + dy*sin(alpha))*cos(theta[i]) - (dy*cos(alpha) - dx*sin(alpha)) - a*sin(theta[i]-beta[i])*delta_alpha), (i,0,0));
A3_7 = summation((J2[i]*(dalpha + d**(-1)*(-dalpha*cos(beta[i] - theta[i]) - nu2*cos(theta[i]) + nu1*sin(theta[i])))*(delta_alpha + dtheta[i])), (i,0,0));
A3_8 = summation((r**(-1)*(-(dx*cos(alpha) + dy*sin(alpha))*cos(theta[i])) - (-dx*sin(alpha) + dy*cos(alpha)) - a*sin(theta[i]-beta[i])*delta_alpha), (i,0,0));
A3_9 = summation((J2[i]*r**(-1)*(-nu1*cos(theta[i]) - nu2*sin(theta[i]) - a*sin(theta[i] - beta[i])*dalpha ) - T[i]), (i,0,0));
A3 = m*(A3_1 + A3_2 + A3_3 + A3_4 + A3_5)*A3_6 + A3_7 + A3_8*A3_9;
A = A1+A2+A3

coeff_delta_alpha = (A1+A2).subs(dx,0).subs(dy,0).subs(delta_alpha,1).subs([(beta[j], pi/2+j*2*pi/3) for j in range(3)])
coeff_dy     = (A1+A2).subs(dx,0).subs(dy,1).subs(delta_alpha,0).subs([(beta[j], pi/2+j*2*pi/3) for j in range(3)])
coeff_dx     = (A1+A2).subs(dx,1).subs(dy,0).subs(delta_alpha,0).subs([(beta[j], pi/2+j*2*pi/3) for j in range(3)])
# cprint('\coeff alpha\n','magenta')
# print(coeff_delta_alpha)
# cprint('\coeff dy\n','magenta')
# print(coeff_dy)
# cprint('\coeff dx\n','magenta')
# print(coeff_dx)


print('\n\n')
system = [coeff_delta_alpha, coeff_dx, coeff_dy]
# print(system)
solve = solve(system, [dnu1, dnu2, ddalpha],dict=True)
rhs_dnu1, rhs_dnu2, rhs_ddalpha = tqdm.tqdm(solve[0][dnu1], solve[0][dnu2], solve[0][ddalpha])
pprint(rhs_dnu1)


# cprint('\ncoeffs\n','magenta')
# coeff_delta_alpha.coeff(nu1)

# print(around_delta_alpha)
# cprint('\ndiff(nu1(t),t)\n','magenta')
# print(around_delta_alpha.coeff(Derivative(nu1(t), t)))
# cprint('\ndnu2\n','magenta')
# print(around_delta_alpha.coeff(dnu2))
# cprint('\nalpha\n','magenta')
# print(around_delta_alpha.coeff(alpha))

# cprint('\ndx\n','magenta')
# cprint('\ndy\n','magenta')


# cprint('\nA1\n','magenta')
# pprint(A1)





# print(solve(A, delta_alpha))
# map(print(solve(A, dtheta[i])),[1,2,3])

# init_printing()
# with open('latex.tex','w') as f:
    # f.write(latex.print_tex(A))


