from sympy import *
from sympy.vector import CoordSys3D
from sympy.physics.vector import *


C = CoordSys3D('C')
t, a_1, a_2, a_3 = symbols('t,a_1,a_2,a_3', positive=True)
theta = symbols('theta', domain=Interval(0, pi))
phi = symbols('phi', domain=Interval(0, 2*pi))
if False:
    if False:
        r =  a_1*cos(theta)*C.i \
            +a_2*sin(theta)*C.j
        T = r.diff(theta).normalize()
        n_a = T.diff(theta).normalize()
        R = simplify(r + t*n_a)
        pprint(printing.octave.octave_code(simplify(n_a)))
    else:
        X_a =   a_1*sin(theta)*cos(phi)*C.i \
              + a_2*sin(theta)*sin(phi)*C.j \
              + a_3*cos(theta)*C.k
        dX_adtheta = X_a.diff(theta)
        dX_adphi = X_a.diff(phi)
        n_a = dX_adtheta.cross(dX_adphi).normalize()
        X_b = simplify(X_a + t*n_a)
        pprint(printing.octave.octave_code(simplify(X_b)))
        print_latex(simplify(X_b))
else:
    x_1, x_2, x_3 = symbols('x_1,x_2,x_3', real=True)
    q = sqrt(a_1**2*a_2**2*cos(theta)**2+sin(theta)**2*(a_1**2*a_3**2*sin(phi)**2+a_2**2*a_3**2*cos(phi)**2))
    expr1 = x_1 - (a_1+t*a_2*a_3/q)*sin(theta)*cos(phi)
    expr2 = x_2 - (a_2+t*a_1*a_3/q)*sin(theta)*sin(phi)
    expr3 = x_3 - (a_3+t*a_1*a_2/q)*cos(theta)
    Y = solve(expr1,y)
    print(Y)
    implicitForm = expr2.subs(y,Y)
    #print_latex(simplify(expr1))
    #print_latex(expr1)
    #print(' ')
    #print_latex(expr2)
    #print(' ')
    print_latex(implicitForm)
    #print_latex(collect(expand(expr1),y))
    #print(' ')
    #print_latex(collect(expand(expr2),y))
#pprint(R)
#octave_code(r+t*n_a)

#pprint(n_b.dot(n_a))
