from sympy import *
from sympy.tensor import *
from termcolor import cprint

from variables import *
from structure import *

import structure    # для запуска функции по имени



def cpprint(obj, color='magenta', pr_str='***'):
    '''печатает объект, перед этим напечатав поясняющую цветнуюстроку'''
    cprint('\n' + ('[' + pr_str + ']').center(150,'.') + '\n', color)
    pprint(obj)


def fkey(dict_key):
    '''Нужна чтобы ключами словаря делать словарь, почему таак приходится делать не очень ясно
       https://stackoverflow.com/questions/13264511/typeerror-unhashable-type-dict'''
    if isinstance(dict_key, dict):
        return frozenset((key, fkey(value)) for key, value in dict_key.items())
    elif isinstance(dict_key, list):
        return tuple(fkey(value) for value in dict_key)
    return dict_key


def scalar(a, b):
    '''скалярное произведение'''
    return(a.dot(b))


def cross(a, b):
    '''векторное произведение'''
    return(a.cross(b))


def vec_invert(vec):
    '''из вектора AB получает BA
       (вектор это лямбда функция)'''
    return vec*(-1);

def intersection(dict_one, dict_two, key='where'):
    '''Находит пересечения значений по ключу в 2 словарях.
       Нужно для того чтобы писать формулу Эйлера и смотреть, принадлежат ли точки одному тело или нет'''
    return list(set(dict_one[key]) & set(dict_two[key]))


def vec_by_2dots(dot_one, dot_two):
    '''по точкам A, B  вернёт  AB или -BA'''
    vec_name = dot_one['name'] + dot_two['name']
    try:
        res_vec = getattr(structure, vec_name)
#         print(vec_name)
    except Exception as e:
        try:
            res_vec = lambda i: - getattr(structure, vec_name[::-1])(i)
#             print('-' + vec_name[::-1])
            # print('Обратный вектору {}?'.format(vec_name))
        except Exception as e:
            cprint('Такого вектора не задано','red')
            print(e)
#     print(res_vec(0))
    return res_vec


def subs_delta(eq):
    return lambda i: eq(i).subs(
        Derivative(alpha,t), delta['alpha']).subs(
        Derivative(psi[i],t), delta['psi'][i]).subs(
        Derivative(theta[i],t), delta['theta'][i])


def dalamber_subs(obj, eq):
    """ subs delta_psi delta_theta """
    return lambda i: obj(i).\
        subs(delta['psi'][i], eq['diff(psi)'](i)).\
        subs(delta['theta'][i], eq['diff(theta)'](i))


def dalamber_subs_nu(obj, eq):
    for i in range(3):
        obj = obj.\
        subs(Derivative(psi[i], t, 2), eq['diff(psi)_nu'](i).diff(t)).\
        subs(Derivative(psi[i], t, 2), eq['diff(psi)_nu'](i).diff(t)).\
        subs(Derivative(psi[i], t), eq['diff(psi)_nu'](i)).\
        subs(Derivative(theta[i], t), eq['diff(theta)_nu'](i))
    return obj


def euler(dot_one, dot_two):
    '''уравнение эйлера по 2 точкам для первой точки'''
    if dot_one == dot_two:
        cprint('Уравнениу Эйлера для одной точки', 'red')
#     print(vec_by_2dots(dot_two, dot_one)(0))
    return lambda i: v[dot_two['name']](i) + cross(omega[intersection(dot_one, dot_two)[0]](i),
                                                   vec_by_2dots(dot_two, dot_one)(i))


def dalamber(mass, velocity, F, delta_r, K, M, omega_delta):
    '''Д'Аламбера - Лагранжа'''
    dalamber = lambda i: scalar((mass*velocity(i).diff(t) - \
                                F(i)), delta_r(i)) + \
                         scalar(K(i).diff(t) - \
                                M(i), omega_delta(i)) 
    return dalamber