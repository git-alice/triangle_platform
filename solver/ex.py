from pyodesys.symbolic import SymbolicSys
def f(t, y, p):
    return [y[1], -y[0] + p[0]*y[1]*(1 - y[0]**2)]

odesys = SymbolicSys.from_callback(f, 2, 1)
xout, yout, info = odesys.integrate(10, [1, 0], [1], integrator='odeint', nsteps=1000)
_ = odesys.plot_result()
import matplotlib.pyplot as plt; plt.show()  # doctest: +SKIP
