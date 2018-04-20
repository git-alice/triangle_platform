# import latex
# import sys
# import time
# import tqdm
# import numpy as np
# import inspect
# import types
# import sympy.printing as printing


from sympy import *
from sympy.tensor import *


from variables import *
from structure import *
from functions import *
from latex     import *


if __name__ == '__main__':

    k = 0 
    # угловые скорости (платформы, вилки, колеса) относительно пола
    omega['platform']= lambda i: Derivative(alpha,t)*e['z']
    omega['fork']    = lambda i: omega['platform'](i) + Derivative(theta[i],t)*e['z']
    omega['wheel']   = lambda i: omega['fork'](i) + Derivative(psi[i],t)*n_wheel(i)

    # УРАНЕНИЯ ЭЁЛЕРА И ОТСУТСТВИЕ ПРОСКАЛЬЗЫВАНИЯ 
    # (euler возвращает лямбда функцию)
    v[fkey(S)] = euler(S, P)
    v[fkey(P)] = euler(P, C)
    v[fkey(C)] = euler(C, D)
    v[fkey(D)] = lambda i: Matrix([0,0,0]) # проскальзывания нет


    #Полученные выражения из связей для nu1 nu2 delta_x delta_y
    eq[fkey(delta['x'])] = lambda i: scalar(v[fkey(S)](i), e['x'])
    eq[fkey(delta['y'])] = lambda i: scalar(v[fkey(S)](i), e['y'])
    eq[fkey(nu[1])]      = lambda i: scalar(v[fkey(S)](i), e['xi'])
    eq[fkey(nu[2])]      = lambda i: scalar(v[fkey(S)](i), e['eta'])

    eq['f(delta_x,delta_y)'] = lambda i: solve(
                          [Eq(eq[fkey(delta['x'])](i), delta['x']), Eq(eq[fkey(delta['y'])](i), delta['y'])],
                          [Derivative(psi[i],t), Derivative(theta[i],t)],
                          dict=True)[0]; # возвращает словарь с выражениями для diff(psi) и diff(theta)
    #.subs(Derivative(alpha,t), delta['alpha']) пдставить
    eq['f(nu1,nu2)'] = lambda i: solve(
                          [Eq(eq[fkey(nu[1])](i), nu[1]), Eq(eq[fkey(nu[2])](i), nu[2])],
                          [Derivative(psi[i],t), Derivative(theta[i],t)],
                          dict=True)[0]; # возвращает словарь с выражениями для diff(psi) и diff(theta)


    # cpprint(eq['f(delta_x,delta_y)'](i), pr_str = "eq['f(delta_x,delta_y)']           i=0")
    eq['diff(psi)']   = lambda i: eq['f(delta_x,delta_y)'](i)[Derivative(psi[i],t)].subs(Derivative(alpha,t), delta['alpha'])
    eq['diff(theta)']        = lambda i: eq['f(delta_x,delta_y)'](i)[Derivative(theta[i],t)].subs(Derivative(alpha,t), delta['alpha'])

    
    # cpprint(temp(0), pr_str="diff(psi)")
    # cpprint(eq['diff(theta)'](0), pr_str="diff(thtea) = g(delta[x], delta[y], delta[alpha])")


    #Д'Аламбер 
    # Dalamber = lambda i: scalar((m*velocity_s.diff(t) + F), delta_r) + scalar(K(i).diff(t) + M(i), omega_delta(i)) 



    #ТЕЛО 1 (платформа)
    velocity[fkey(S)]        = lambda i: nu[1]*e['xi'] + nu[2]*e['eta']  # !!!! (надо бы подставить, скорость которую уже считали выше) 
    F['platform']            = lambda i: zeros(3,1) # сил не действует
    delta_r[fkey(S)]         = subs_delta(lambda i: delta['x']*e['x'] + delta['y']*e['y'])
    J['platform']            = eye(3,3)*a  # !!!!! (симметричная) (временные коэффиценты)
    K['platform']            = lambda i: J['platform']*omega['platform'](i)
    omega_delta['platform']  = subs_delta(lambda i: omega['platform'](i).subs(Derivative(alpha,t),delta['alpha']))
    M['platform']            = lambda i: -W[0]*e['z']

    A['platform'] = dalamber(velocity[fkey(S)],
                             F['platform'],
                             delta_r[fkey(S)],
                             K['platform'],
                             M['platform'],
                             omega_delta['platform'])


    #ТЕЛО 2 (вилки)
    velocity[fkey(P)]    = lambda i: velocity[fkey(S)](i) + cross(omega['platform'](i), vec_by_2dots(S,P)(i))   # !!!!
    # velocity['fork']     = v[fkey(P)]
    F['fork']            = lambda i: zeros(3,1) # сил не действует
    delta_r[fkey(P)]     = subs_delta(lambda i: delta_r[fkey(S)](i) + cross(omega['platform'](i), vec_by_2dots(S,P)(i)))
    J['fork']            = zeros(3,3) # невесома
    K['fork']            = lambda i: J['fork']*omega['fork'](i)
    omega_delta['fork']  = subs_delta(lambda i: omega['fork'](i))
    M['fork']            = lambda i: -(-W[i]*e['z']) -(-T[i]*n_wheel(i)) # magic

    A['fork'] = dalamber(velocity[fkey(P)],
                         F['fork'],
                         delta_r[fkey(P)],
                         K['fork'],
                         M['fork'],
                         omega_delta['fork'])


    #ТЕЛО 3 (колёса)
    velocity[fkey(C)]     = lambda i: velocity[fkey(P)](i) + cross(omega['fork'](i), vec_by_2dots(P,C)(i))   # !!!!
    F['wheel']            = lambda i: zeros(3,1) # сил не действует
    delta_r[fkey(C)]      = subs_delta(lambda i: delta_r[fkey(P)](i) + cross(omega['fork'](i), vec_by_2dots(S,P)(i)))
    J['wheel']            = Matrix([[a,0,0],[0,b,0],[0,0,c]])
    K['wheel']            = lambda i: J['wheel']*omega['wheel'](i)
    omega_delta['wheel']  = subs_delta(lambda i: omega['wheel'](i))
    M['wheel']            = lambda i: -(-W[i]*e['z']) -(-T[i]*n_wheel(i)) # magic

    A['wheel'] = dalamber(velocity[fkey(C)],
                         F['wheel'],
                         delta_r[fkey(C)],
                         K['wheel'],
                         M['wheel'],
                         omega_delta['wheel'])




    

    # # Полный Д'Аламбер Лагранж
    A_full  = lambda i: A['platform'](i) + A['fork'](i) + A['wheel'](i)
    A_full_ = lambda i: Poly(A_full(i).subs(delta['psi'][i], eq['diff(psi)'](i)).subs(delta['theta'][i], eq['diff(theta)'](i)), 
        [delta['x'], delta['y'], delta['alpha']])

                                                         # ПРОВЕРИТЬ КОЭФФИЦЕННТЫ
    # print(type(A_full_(i).coeffs(dict=True)))
    print('COEFFS')
    A_coeffs = lambda i: A_full_(i).coeffs()

    A_coeffs = A_coeffs(0)

    coeff[delta['x']]     = lambda i: simplify(A_coeffs[0]).subs([alpha, beta[i], theta[i], nu[1], nu[2]], [0,0,0,0,0])
    coeff[delta['y']]     = lambda i: simplify(A_coeffs[1]).subs([alpha, beta[i], theta[i], nu[1], nu[2]], [0,0,0,0,0])
    coeff[delta['alpha']] = lambda i: simplify(A_coeffs[2]).subs([alpha, beta[i], theta[i], nu[1], nu[2]], [0,0,0,0,0])



    # coeff[delta['alpha']] = lambda i: (A_full_(i)).subs(delta['x'],0).subs(delta['y'],0).subs(delta['alpha'],1)#.subs([(beta[j], pi/2+j*2*pi/3) for j in range(3)])
    # coeff['0'] = lambda i: (A_full_(i)).subs(delta['x'],0).subs(delta['y'],0).subs(delta['alpha'],0)


    # ПЕЧАТЬ
    k = 0
    # cpprint(v[fkey(S)](k))
    # cpprint(v[fkey(P)](k))
    # cpprint(v[fkey(C)](k))
    # cpprint(v[fkey(D)](k))
    # cpprint(eq['diff(psi)'](0), pr_str="diff(psi)")


    # cpprint(A['wheel'](0), pr_str='ДЛЯ КОЛЕСА i=0')

    #  ПРОВЕРКА ВИРТУАЛЬНЫХ ПЕРЕМЕЩЕНИЙ И УГЛОВЫХ СКОРОСТЕЙ
    # cpprint(delta_r[fkey(S)](k), pr_str = 'delta_r [S]')
    # cpprint(delta_r[fkey(P)](k), pr_str = 'delta_r [P]')
    # cpprint(delta_r[fkey(C)](k), pr_str = 'delta_r [C]')
    # cpprint(omega_delta['platform'](k), pr_str = 'omega_delta [platform]')
    # cpprint(omega_delta['fork'](k), pr_str = 'omega_delta [fork]')
    # cpprint(omega_delta['wheel'](k), pr_str = 'omega_delta [wheel]')
    # cpprint(coeff[delta['alpha']](0), pr_str='coeff_delta_alpha')
    # print('collect')
    # print(A_full_(0).coeffs())
    # cpprint(coeff['0'](0), pr_str='NO')

    print(coeff[delta['x']](0))    
    print(coeff[delta['y']](0))
    print(coeff[delta['alpha']](0))

    cpprint(solve([coeff[delta['x']](0), coeff[delta['y']](0), coeff[delta['alpha']](0)], [Derivative(nu[1],t), Derivative(nu[2],t), Derivative(alpha,t,2)]), pr_str="EEEENDDDD")









    # pprint(e['xi'].diff(t))
    # pprint(simplify(m*velocity))
    # pprint(simplify(delta_r))

    # pprint(eqq + Matrix([Derivative(x,t),Derivative(y,t),0]))

    # pprint([eqq[i].coeff(Derivative(psi[1],t)) for i in range(3)]) #!!!

    # print(eqq.coeff(e['xi']))
    # print(eqq.coeff(e['eta']))
    # v_p = v_s + cross(omega_delta, e['z']) # v_s

    # pprint(psi.args[1])
    # print(omega[k])

    # nu1 = Derivative(x,t)*cos(alpha) + Derivative(y,t)*sin(alpha)
    # nu2 = -Derivative(x,t)*sin(alpha) + Derivative(y,t)*cos(alpha)





'''
    ###############
    # СЧИТАЯ \dot(x)e[x] +  \dot(y)e[y] = nu1*e[xi] + n2*e[eta]
    pseudo_vel_eq =  Derivative(x,t)*e['x'] +  Derivative(y,t)*e['y'] - nu1*e['xi'] - nu2*e['eta']

    # print(solve(eq.subs(Derivative(x,t),dx).subs(Derivative(x,t),dy),nu1))
    cprint('EQ','magenta')
    pprint(pseudo_vel_eq)
    cprint('EQ -> dx,dy','magenta')


    new_eq = eq.subs(Derivative(x,t),dx).subs(Derivative(y,t),dy)
    pprint(new_eq)
    print(solve((new_eq[0],new_eq[1]),(nu1,nu2))) 
'''