import sympy as sym
import numpy as np
#symExpression = True
symExpression = False
if symExpression:
    x,y,c_xy,R1,R2,L,t = sym.symbols('x,y,c_xy,R1,R2,L,t')
else:
    R1 = 3
    R2 = 5
    L = 41
    t = 0.25*R2
    x,y,c_xy = sym.symbols('x,y,c_xy')

a = -(R2-R1)/L
#boundaryMethod = True
boundaryMethod = False
if boundaryMethod:
    c_z = (R1+R2+L)/2
    eq1 = sym.Eq((x/c_z)**2+(y/c_xy)**2,1)     # Equation for ellipse: (x/c_z)^2+(y/c_xy)^2=1
    eq2 = sym.Eq(y,a*(x-c_z+R1)+R1)           # Equation for straigth line: y = a*(x-c_z+R1)+R1
    eq3 = sym.Eq(2*x/c_z**2 + 2/c_xy**2*y*a,0) # Equation for derivative of ellipse at (x,y) which is equal to dy/dx = a
    result = sym.solve([eq1,eq2,eq3],(x,y,c_xy))        # Solve equations

    print(result[-1][2])
else:
    c_z = (R1+R2+L)/2 + t
    hyp = sym.sqrt(L**2+(R2-R1)**2)
    x1 = t*(R2-R1)/hyp
    y1 = t*L/hyp
    eq1 = sym.Eq((x/c_z)**2+(y/c_xy)**2,1)     # Equation for ellipse: (x/c_z)^2+(y/c_xy)^2=1
    eq2 = sym.Eq(y,a*(x-c_z+R1-x1)+y1+R1)           # Equation for straigth line: y = a*(x-x1)+y1
    eq3 = sym.Eq(2*x/c_z**2 + 2/c_xy**2*y*a,0) # Equation for derivative of ellipse at (x,y) which is equal to dy/dx = a
    result = sym.solve([eq1,eq2,eq3],(x,y,c_xy))        # Solve equations

    b = c_z
    print(a**4*b**4+1)
    print(result[-1][2])
