
import builtins
from sympy import *
# import dill as pickle  # dill позволяет записывать и функции в отличии от стандартного pickle
import cloudpickle as pickle # но даже дилл иногда не хочет работать

from variables import *
from structure import *


from termcolor import cprint
from pprint import pprint

# Для записи и чтения уравнений


def write_obj(obj, obj_name, comment='', folder='data/'):
    if comment:
        with open(folder + "comments.txt", "a") as fc:
            fc.write(obj_name + ': ' + comment + '\n')
    with open(folder + obj_name + '.pickle', 'wb') as fobj:
        pickle.dump(obj, fobj)


def read_obj(obj_name, folder='data/'):
    with open(folder + obj_name + '.pickle', 'rb') as fobj:
        return pickle.load(fobj)
    
def write_read_obg(obj, obj_name, calc):
    if calc:
        write_obj(obj, obj_name)
    else:
        read_obj(obj_name)
        return obj

def subs_init(eq):
    ''' '''
    eq = eq.subs({
            beta[0]: 0, beta[1]: 2*pi/3, beta[2]: 4*pi/3,  # раосторонний
            m['platform']: 2, 
            m['wheel']: 1,
            a: 3,
            b: 0.01,
            c: 0.02,
            r: 0.2
    })
    return eq


def subs_for_ode(eq):
    ''' getting rid by arrays | alpha'' presents like (alpha1, alpha1') '''
    print('Проверить чтобы это не заупскалось')
    nu1, nu2 = symbols('nu1, nu2')
    theta0, theta1, theta2 = symbols('theta0, theta1, theta2')
    psi0, psi1, psi2 = symbols('psi0, psi1, psi2')
    alpha1, alpha2 = symbols('alpha1, alpha2')
    
    eq = eq.subs({
        nu[1]: nu1,       nu[2]: nu2,
        alpha: alpha1,    Derivative(alpha, t): alpha2,
        psi[0]: psi0,     psi[1]: psi0,     psi[2]: psi2,
        theta[0]: theta0, theta[1]: theta1, theta[2]: theta2
    })
    return eq


def print_all_variables(obj):
    ''' Посмотрет какие переменные есть в выражении '''
    try:
        result = set(obj.free_symbols)
        if obj.find(alpha) or obj.find(Derivative(x, (t,1))) or obj.find(Derivative(x, (t,2))):
            result = result.union(obj.find(alpha))
        iter_set = result  
        for x in iter_set:
            if obj.find(Derivative(x, (t,1))):
                result = result.union(obj.find(Derivative(x, (t,1))))
            if obj.find(Derivative(x, (t,2))):
                result = result.union(obj.find(Derivative(x, (t,2))))
        return result
    except Exception as e:
        print(e)

def get_atoms(obj):
    return obj.atoms()


def c_print(text, obj, pretty=False):
    cprint('[' + text + '] ', 'red', end='');
    if pretty:
        pprint(obj)
    else:
        print(obj)