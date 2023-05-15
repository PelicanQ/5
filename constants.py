import calfem.core as cfc
import numpy as np
#Materialparametrar, koppar: 1, nylon: 2
E1 = 128e9    #Pa
E2 = 3e9      #Pa
ny1 = 0.36
ny2 = 0.39
alfa1 = 17.6e-6    #1/K
alfa2 = 80e-6      #1/K
rho1 = 8930     #kg/m^3   
rho2 = 1100     #kg/m^3
cp1 = 386       #J/kgK
cp2 = 1500      #J/kgK
k1 = 385        #W/m-K
k2 = 0.26       #W/m-K

a_c = 40        #W/m^2K
T_inf = 18       #grader C
T_0 = 18       #grader C
h_flow = -10**5       #W/m^2

#Sidor
L = 5e-3   #m
a = 0.1*L
b = 0.1*L
c = 0.3*L
d = 0.05*L
h = 0.15*L
t = 0.05*L
thickness = 5e-3 #m

ptype = 2 # plane strain
ep = [ptype, thickness]

D1 = np.eye(2) * k1
D2 = np.eye(2) * k2