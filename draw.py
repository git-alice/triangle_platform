import numpy as np
import pandas as pd
from sympy import *
from sympy.tensor import *
import sympy.printing as printing

from variables import *
from structure import *
from functions import *
from latex     import *
from utils import read_obj, write_obj, subs_for_ode, subs_init, print_all_variables
import utils 
from ode_utils import numerical_solver
from pyodesys.symbolic import SymbolicSys
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.plotly as py

def data_for_plotly(*args):
    data = []
    for arg in args:
        data.append(go.Scatter(
            x = arg['x'],
            y = arg['y'],
            mode = 'lines',
            name = arg['name']))

    return data

def get_shapes(df, n):
    shapes = []
    for i in range(0, len(df.x), n):
        i = int(i)
#       треуголники
        shapes.append({
                'type': 'line',
                'xref': 'x',
                'yref': 'y',
                'x0': df.x.iloc[i] + df.SP0x.iloc[i],
                'y0': df.y.iloc[i] + df.SP0y.iloc[i],
                'x1': df.x.iloc[i] + df.SP1x.iloc[i],
                'y1': df.y.iloc[i] + df.SP1y.iloc[i],
                'line': {
                    'color': 'rgb(191, 63, 63)',
                    'width': 1,
                }})
        shapes.append({
                'type': 'line',
                'xref': 'x',
                'yref': 'y',
                'x0': df.x.iloc[i] + df.SP1x.iloc[i],
                'y0': df.y.iloc[i] + df.SP1y.iloc[i],
                'x1': df.x.iloc[i] + df.SP2x.iloc[i],
                'y1': df.y.iloc[i] + df.SP2y.iloc[i],
                'line': {
                    'color': 'rgb(63, 191, 191)',
                    'width': 1,
                }})
        shapes.append({
                'type': 'line',
                'xref': 'x',
                'yref': 'y',
                'x0': df.x.iloc[i] + df.SP2x.iloc[i],
                'y0': df.y.iloc[i] + df.SP2y.iloc[i],
                'x1': df.x.iloc[i] + df.SP0x.iloc[i],
                'y1': df.y.iloc[i] + df.SP0y.iloc[i],
                'line': {
                    'color': 'rgb(63, 191, 127)',
                    'width': 1,
                }})

# колеса
        shapes.append({
                'type': 'line',
                'xref': 'x',
                'yref': 'y',
                'x0': df.x.iloc[i] + df.SP0x.iloc[i],
                'y0': df.y.iloc[i] + df.SP0y.iloc[i],
                'x1': df.x.iloc[i] + df.SP0x.iloc[i] + df.PC0x.iloc[i],
                'y1': df.y.iloc[i] + df.SP0y.iloc[i] + df.PC0y.iloc[i],
                'line': {
                    'color': 'rgb(22, 100, 130)',
                    'width': 1,
                }})
        shapes.append({
                'type': 'line',
                'xref': 'x',
                'yref': 'y',
                'x0': df.x.iloc[i] + df.SP1x.iloc[i],
                'y0': df.y.iloc[i] + df.SP1y.iloc[i],
                'x1': df.x.iloc[i] + df.SP1x.iloc[i] + df.PC1x.iloc[i],
                'y1': df.y.iloc[i] + df.SP1y.iloc[i] + df.PC1y.iloc[i],
                'line': {
                    'color': 'rgb(22, 100, 130)',
                    'width': 1,
                }})
        shapes.append({
                'type': 'line',
                'xref': 'x',
                'yref': 'y',
                'x0': df.x.iloc[i] + df.SP2x.iloc[i],
                'y0': df.y.iloc[i] + df.SP2y.iloc[i],
                'x1': df.x.iloc[i] + df.SP2x.iloc[i] + df.PC2x.iloc[i],
                'y1': df.y.iloc[i] + df.SP2y.iloc[i] + df.PC2y.iloc[i],
                'line': {
                    'color': 'rgb(22, 100, 130)',
                    'width': 1,
                }})
    return shapes

def subs_symbols(obj):
    x = Function('x')(t)
    y = Function('y')(t)
    alpha = Function('alpha')(t)
    return obj.subs({x: symbols('x'),
                     y: symbols('y'),
                     alpha: symbols('alpha'),
                     psi[0]: symbols('psi_0'),
                     psi[1]: symbols('psi_1'),
                     psi[2]: symbols('psi_2'),
                     theta[0]: symbols('theta_0'),
                     theta[1]: symbols('theta_1'),
                     theta[2]: symbols('theta_2')
                    })

def subs_symbols_after(obj, yout, params=None):
    obj = obj.\
        subs(Function('alpha')(t), eq['alpha']).\
        subs(theta[0], symbols('theta_0')).\
        subs(theta[1], symbols('theta_1')).\
        subs(theta[2], symbols('theta_2')).\
        subs(psi[0], symbols('psi_0')).\
        subs(psi[1], symbols('psi_1')).\
        subs(psi[2], symbols('psi_2')).\
        subs({beta[0]: 0, beta[1]: 2*pi/3, beta[2]: 4*pi/3}).\
        subs({symbols('theta_0'): yout[3]}).\
        subs({symbols('theta_1'): yout[4]}).\
        subs({symbols('theta_2'): yout[5]})
    if params:
        obj = obj.subs(params)
    return obj 

def evaluate_two_vec(vec1, vec2, xout, yout, params):
    temp = [[], [], []]
    for i in range(3):
        for t_step, res_step in zip(xout, yout):
            temp[i].append(subs_symbols_after(vec1(i) + vec2(i), res_step, params).\
                           subs({t: t_step}).\
                           tolist())
    return temp

def evaluate_vec(vec, xout, yout, params):
    temp = [[], [], []]
    for i in range(3):
        print(vec(i))
        for t_step, res_step in zip(xout, yout):
            temp[i].append(subs_symbols_after(vec(i), res_step, params).\
                           subs({t: t_step}).\
                           tolist())
    return temp

def get_df(xout, yout, sc_coords, sp_coords, pc_coords):
    df = pd.DataFrame(yout, index=xout, columns=['psi0', 'psi1', 'psi2', 'theta0', 'theta1', 'theta2'])

    df['contact0x'] = pd.Series(to_np_float(np_concatenate(sc_coords[0][:,0])), index=df.index)
    df['contact0y'] = pd.Series(to_np_float(np_concatenate(sc_coords[0][:,1])), index=df.index)
    df['contact1x'] = pd.Series(to_np_float(np_concatenate(sc_coords[1][:,0])), index=df.index)
    df['contact1y'] = pd.Series(to_np_float(np_concatenate(sc_coords[1][:,1])), index=df.index)
    df['contact2x'] = pd.Series(to_np_float(np_concatenate(sc_coords[2][:,0])), index=df.index)
    df['contact2y'] = pd.Series(to_np_float(np_concatenate(sc_coords[2][:,1])), index=df.index)

    df['SP0x'] = pd.Series(to_np_float(np_concatenate(sp_coords[0][:,0])), index=df.index)
    df['SP0y'] = pd.Series(to_np_float(np_concatenate(sp_coords[0][:,1])), index=df.index)
    df['SP1x'] = pd.Series(to_np_float(np_concatenate(sp_coords[1][:,0])), index=df.index)
    df['SP1y'] = pd.Series(to_np_float(np_concatenate(sp_coords[1][:,1])), index=df.index)
    df['SP2x'] = pd.Series(to_np_float(np_concatenate(sp_coords[2][:,0])), index=df.index)
    df['SP2y'] = pd.Series(to_np_float(np_concatenate(sp_coords[2][:,1])), index=df.index)

    df['PC0x'] = pd.Series(to_np_float(np_concatenate(pc_coords[0][:,0])), index=df.index)
    df['PC0y'] = pd.Series(to_np_float(np_concatenate(pc_coords[0][:,1])), index=df.index)
    df['PC1x'] = pd.Series(to_np_float(np_concatenate(pc_coords[1][:,0])), index=df.index)
    df['PC1y'] = pd.Series(to_np_float(np_concatenate(pc_coords[1][:,1])), index=df.index)
    df['PC2x'] = pd.Series(to_np_float(np_concatenate(pc_coords[2][:,0])), index=df.index)
    df['PC2y'] = pd.Series(to_np_float(np_concatenate(pc_coords[2][:,1])), index=df.index)

    # calc x, y
    xy = evaluate_xy(eq['x'], eq['y'], xout)
    df['x'] = to_np_float(xy['x'])
    df['y'] = to_np_float(xy['y'])

    return df

def get_sc_coords(xout, yout, params):
    return np.array(evaluate_two_vec(SP, PC, xout, yout, params))

def get_sp_coords(xout, yout, params):
    return np.array(evaluate_vec(SP, xout, yout, params))

def get_pc_coords(xout, yout, params):
    return np.array(evaluate_vec(PC, xout, yout, params))

def evaluate_xy(eq_x, eq_y, xout):
    temp = {'x': [], 'y': []}
    for t_step in xout:
        if not type(eq_x) == float:
            temp['x'].append(eq_x.subs({t: t_step}))
        else:
            temp['x'].append(eq_x)
        if not type(eq_y) == float:
            temp['y'].append(eq_y.subs({t: t_step}))
        else:
            temp['y'].append(eq_y)
    return temp

def sum_two_list(a, b):
    return np.add(a, b)

def point_of_contact(x, y, params, ):
    pass

def np_concatenate(a):
    return list(np.concatenate(a))

def to_np_float(obj):
    obj = np.array(obj).astype(np.float64)
    return obj

def calc(eq_x, eq_y, eq_alpha, initial_conditions, params, t_end):
    eq['diff(psi)'] = read_obj('diff(psi)(nu)')
    eq['diff(theta)'] = read_obj('diff(theta)(nu)')
    print('diff(PSI) i = 0: ', eq['diff(psi)'](0))
    print('diff(THETA) i = 0: ', eq['diff(theta)'](0))

    # Представялем в виде массива
    eq['dot(psi)'] = [eq['diff(psi)'](0), eq['diff(psi)'](1), eq['diff(psi)'](2)]
    eq['dot(theta)'] = [eq['diff(theta)'](0), eq['diff(theta)'](1), eq['diff(theta)'](2)]

    # Подставляем геометрию платформы 
    for i in range(3):
        eq['dot(psi)'][i] = eq['dot(psi)'][i].subs({beta[0]: 0, beta[1]: 2*pi/3, beta[2]: 4*pi/3})
        eq['dot(theta)'][i] = eq['dot(theta)'][i].subs({beta[0]: 0, beta[1]: 2*pi/3, beta[2]: 4*pi/3})    

    x = Function('x')(t)
    y = Function('y')(t)
    alpha = Function('alpha')(t)
    # p = symbols('p')

    eq['x'] = eq_x
    eq['y'] = eq_y
    eq['alpha'] = eq_alpha
    eq['dot(alpha)'] = diff(eq_alpha).doit()

    # nu1, nu2
    eq['nu1'] = cos(alpha)*diff(x) + sin(alpha)*diff(y)
    eq['nu2'] = -sin(alpha)*diff(x) + cos(alpha)*diff(y)

    # nu1(alpha, x, y) -> nu1(t); nu2(alpha, x, y) -> nu2(t)
    for i in range(1,3):
        eq['nu' + str(i)] = eq['nu' + str(i)].subs({x: eq['x'],
                                                    y: eq['y'],
                                                    Derivative(alpha, t): eq['dot(alpha)'],
                                                    alpha: eq['alpha']})

    # 
    for i in range(3):
        eq['dot(psi)'][i] = eq['dot(psi)'][i].\
            subs({nu[1]: eq['nu1'],
                 nu[2]: eq['nu2'],
                 Derivative(alpha, t): eq['dot(alpha)'],
                 alpha: eq['alpha']}).\
            doit()
        eq['dot(theta)'][i] = eq['dot(theta)'][i].\
            subs({nu[1]: eq['nu1'],
                  nu[2]: eq['nu2'],
                  Derivative(alpha, t): eq['dot(alpha)'],
                  alpha: eq['alpha']}).\
            doit()  

    psi0, psi1, psi2, theta0, theta1, theta2, x, y, alpha = symbols('psi_0, psi_1, psi_2, theta_0, theta_1, theta_2, x, y, alpha')

    for i in range(3):
        eq['dot(psi)'][i] = subs_symbols(eq['dot(psi)'][i])
        eq['dot(theta)'][i] = subs_symbols(eq['dot(theta)'][i])

    left = [psi0, psi1, psi2, theta0, theta1, theta2]
    right = [eq['dot(psi)'][0], eq['dot(psi)'][1], eq['dot(psi)'][2], eq['dot(theta)'][0], eq['dot(theta)'][1], eq['dot(theta)'][2]]
    print(right)
    print('Количество совпадает:', len(left) == len(right))

    odesys = SymbolicSys(list(zip(left, right)), t, [p, r, d])
    xout, yout, info = odesys.integrate(t_end,
                                        initial_conditions,
                                        params,
                                        integrator='gsl',
                                        method='rk8pd',
                                        atol=1e-11,
                                        rtol=1e-12)

    return xout, yout, info, left, right