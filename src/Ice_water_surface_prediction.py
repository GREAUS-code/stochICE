from gekko import GEKKO

"""
Input variables from HECRAS model
"""
n=0.03
n_ice=0.04
s=0.00038 # need this from hecras model (or ask the user to manually input it for each case)
t=0.2

"""
Input variables extracted from downstream cross-section of TAPE5.txt
"""
WSE_rivice=67.68
Q_rivice=50
A_rivice=61
P_rivice=38
B_rivice=36

"""
Preliminary calculations
"""

n_c=((1+(n_ice/n)**(3/2))/2)**(2/3)*n


"""
Rectangular cross-section approximation
"""

# Define initial guesses
D_i=10
A_glace=B_rivice*D_i
P_glace=2*B_rivice+2*D_i
R_glace=A_glace/P_glace
V_glace=1/(n_c)*(R_glace**(2/3))*(s**(1/2))

# setup up equation solver
m = GEKKO()  

D_i = m.Var(value = D_i, lb = 0)
A_glace = m.Var(value = A_glace, lb = 0)
P_glace = m.Var(value = P_glace)
R_glace = m.Var(value = R_glace)
V_glace = m.Var(value = V_glace)

# solve equations
m.Equation(B_rivice*D_i==A_glace)
m.Equation(2*B_rivice+2*D_i==P_glace)
m.Equation(A_glace/P_glace == R_glace)
m.Equation(1/(n_c)*(R_glace**(2/3))*(s**(1/2)) == V_glace)
m.Equation(V_glace*A_glace == Q_rivice)


m.solve(disp = True)

wse_ice=WSE_rivice-(A_rivice/B_rivice)+D_i.value[0]+t
print(D_i.value[0], V_glace.value[0]*A_glace.value[0], wse_ice)
