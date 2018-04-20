import latex 
import sympy.printing as printing
import sys
from sympy import *
from termcolor import colored, cprint


#init
nu1,nu2,alpha,m,t,dx,dy,dalpha,J1,d,a,r = symbols('nu1,nu2,alpha,m,t,dx,dy,dalpha,J1,d,a,r')
k = symbols('k',integer=True)
i = symbols('i',integer=True)
W = IndexedBase("W")
beta = IndexedBase('beta')
theta = IndexedBase('theta')
J2 = IndexedBase('J2')
T = IndexedBase('T')
dtheta = IndexedBase('dtheta')


A1 =  m*((diff(nu1(t),t)*cos(alpha) - nu1(t)*sin(alpha)*diff(alpha(t),t) - diff(nu2(t),t)*sin(alpha) -  nu2(t)*cos(alpha)*diff(alpha(t),t))*dx +
(diff(nu1(t),t)*sin(alpha) + nu1(t)*cos(alpha)*diff(alpha(t),t) + diff(nu2(t),t)*cos(alpha) -  nu2(t)*sin(alpha)*diff(alpha(t),t))*dy  + 
    J1*diff(alpha(t),t)*dalpha);
A2 = summation(W[i]*(dalpha + d**(-1)*(-dalpha*cos(beta[i] - theta[i]) - nu2*cos(theta[i]) + nu1*sin(theta[i]) ) ),(i,0,2));
A3_1 = summation((-diff(nu1(t),t)*cos(theta[i])), (i,0,2));
A3_2 = summation((nu1(t)*sin(theta[i])*d**(-1)*(-diff(alpha(t),t)*cos(beta[i] - theta[i]) - nu2(t)*cos(theta[i]) + nu1(t)*sin(theta[i]))), (i,0,2));
A3_3 = summation((-diff(nu2(t),t)*sin(theta[i])), (i,0,2));
A3_4 = summation((nu2(t)*cos(theta[i])*d**(-1)*(-diff(alpha(t),t)*cos(beta[i] - theta[i]) - nu2(t)*cos(theta[i]) + nu1(t)*sin(theta[i]))), (i,0,2));
A3_5 = summation((-a*sin(theta[i] - beta[i])*diff(alpha(t),t,2) - a*cos(theta[i]- beta[i])*diff(alpha(t),t)*diff(alpha(t),t)), (i,0,2));
A3_6 = summation(((-dx*cos(alpha) + dy*sin(alpha))*cos(theta[i]) - (dy*cos(alpha) - dx*sin(alpha)) - a*sin(theta[i]-beta[i])*dalpha), (i,0,2));
A3_7 = summation((J2[i]*(diff(alpha(t),t) + d**(-1)*(-diff(alpha(t),t)*cos(beta[i] - theta[i]) - nu2(t)*cos(theta[i]) + nu1(t)*sin(theta[i])))*(dalpha + dtheta[i])), (i,0,2));
A3_8 = summation((r**(-1)*(-(dx*cos(alpha) + dy*sin(alpha))*cos(theta[i])) - (-dx*sin(alpha) + dy*cos(alpha)) - a*sin(theta[i]-beta[i])*dalpha), (i,0,2));
A3_9 = summation((J2[i]*r**(-1)*(-nu1(t)*cos(theta[i]) - nu2(t)*sin(theta[i]) - a*sin(theta[i] - beta[i])*diff(alpha(t),t) ) - T[i]), (i,0,2));
A3 = m*(A3_1 + A3_2 + A3_3 + A3_4 + A3_5)*A3_6 + A3_7 + A3_8*A3_9;
A = A1+A2+A3

# 100
# 010
# 001

around_dalpha = A.subs(dx,0).subs(dy,0).subs(dalpha,1)
around_dy     = A.subs(dx,0).subs(dy,1).subs(dalpha,0)
around_dx     = A.subs(dx,1).subs(dy,0).subs(dalpha,0)

cprint('\naround alpha\n','magenta')
print(around_dalpha)
cprint('\ndiff(nu1(t),t)\n','magenta')
print(around_dalpha.coeff(Derivative(nu1(t), t)))
cprint('\ndiff(nu2(t),t)\n','magenta')
print(around_dalpha.coeff(diff(nu2(t),t)))
cprint('\nalpha\n','magenta')
print(around_dalpha.coeff(alpha))

# cprint('\ndx\n','magenta')
# cprint('\ndy\n','magenta')


# cprint('\nA1\n','magenta')
# pprint(A1)





# print(solve(A, dalpha))
# map(print(solve(A, dtheta[i])),[1,2,3])

# init_printing()
# with open('latex.tex','w') as f:
    # f.write(latex.print_tex(A))


