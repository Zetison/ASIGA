from sympy import *
from sympy.abc import r, z, theta, n

x_nm, x_nms, n3, R, L = symbols('x_nm x_nms n3 R L')
p = sin(n3*pi*(z+L/2)/L)*besselj(n, r*x_nm/R)*cos(n*theta)
k = sqrt((x_nm/R)**2 + (n3*pi/L)**2)
rdpdr = r*p.diff(r)
residual = rdpdr.diff(r)/r + p.diff(theta,2)/r**2 + p.diff(z,2) + k**2*p
print(simplify(residual))
print_latex(simplify(residual))

p = cos(n3*pi*(z+L/2)/L)*besselj(n, r*x_nms/R)*cos(n*theta)
k = sqrt((x_nms/R)**2 + (n3*pi/L)**2)
rdpdr = r*p.diff(r)
residual = rdpdr.diff(r)/r + p.diff(theta,2)/r**2 + p.diff(z,2) + k**2*p
print(simplify(residual))
print_latex(simplify(residual))

