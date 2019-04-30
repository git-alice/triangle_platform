from sympy import *
from sympy.tensor import *
import sympy.printing as printing
from sympy.physics.vector import dynamicsymbols

from variables import *
from structure import *
from functions import *
from latex import *

import utils
from utils import write_obj, read_obj, c_print
from utils import print_all_variables as pa

class SymbolSolver(object):
	"""docstring for SymbolSolver"""
	def __init__(self):
		pass
		
	@staticmethod
	def it(f, symplify_it=True):
		if symplify_it:
			res = [simplify(f(i)) for i in range(3)]
		else:
			res = [f(i) for i in range(3)]
		return res

	def get_omegas(self, calc=True):
		""" угловые скорости (платформы, вилки, колеса) относительно пола """
		if calc:
			print(omega)
			omega['platform']= lambda i: Derivative(alpha,t)*e['z']
			omega['fork']    = lambda i: omega['platform'](i) + Derivative(theta[i],t)*e['z']
			omega['wheel']   = lambda i: omega['fork'](i) + Derivative(psi[i],t)*n_wheel(i)
			write_obj(omega, 'omega')
			# write_obj(omega['platform'], 'angular_velocity_platform')
			# write_obj(omega['fork'], 'angular_velocity_fork')
			# write_obj(omega['wheel'], 'angular_velocity_wheel')
		else:
			omega = read_obj('omega')
			# omega['platform'] = read_obj('angular_velocity_platform')
			# omega['fork'] = read_obj('angular_velocity_fork')
			# omega['wheel'] = read_obj('angular_velocity_wheel')
		return omega

	def get_euler_equations(self, calc=True):
		""" УРАНЕНИЯ ЭЁЛЕРА И ОТСУТСТВИЕ ПРОСКАЛЬЗЫВАНИЯ """
		if calc:
			v['S'] = euler(S, P)
			v['P'] = euler(P, C)
			v['C'] = euler(C, D)
			v['D'] = lambda i: Matrix([0,0,0]) # проскальзывания нет
			# write_obj(self.it(v['S']), 'point_S_velocity')
			# write_obj(self.it(v['P']), 'point_P_velocity')
			# write_obj(self.it(v['C']), 'point_C_velocity')
			# write_obj(self.it(v['D']), 'point_D_velocity')
			write_obj(v, 'v')
		else:
			v = read_obj('v')
			# v['S'] = read_obj('point_S_velocity')
			# v['P'] = read_obj('point_P_velocity')
			# v['C'] = read_obj('point_C_velocity')
			# v['D'] = read_obj('point_D_velocity')
		return v

	def get_constraint(self):
		if calc:
			#Полученные выражения из связей для nu1 nu2 delta_x delta_y
			eq['delta_x'] = lambda i: scalar(v['S'](i), e['x'])
			eq['delta_y'] = lambda i: scalar(v['S'](i), e['y'])
			eq['nu_1'] = lambda i: scalar(v['S'](i), e['xi'])
			eq['nu_2'] = lambda i: scalar(v['S'](i), e['eta'])
			write_obj(eq, 'eq')
		else:
			eq = read_obj('eq')
		return eq