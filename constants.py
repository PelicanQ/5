import calfem.core as cfc

#Materialparametrar, koppar: 1, nylon: 2
E1 = 128    #GPa
E2 = 3      #GPa
ny1 = 0.36
ny2 = 0.39
alfa1 = 17.6*10**-6    #1/K
alfa2 = 80*10**-6      #1/K
rho1 = 8930     #kg/m^3   
rho2 = 1100     #kg/m^3
cp1 = 386       #J/kgK
cp2 = 1500      #J/kgK
k1 = 385        #W/mK
k2 = 0.26       #W/mK

ac = 40         #W/m^2K
Tinf = 18       #grader C
h = 10**5       #W/m^2

#Sidor
L = 5    #mm
a = 0.1*L
b = 0.1*L
c = 0.3*L
d = 0.05*L
h = 0.15*L
t = 0.05*L

ptype = 1 # plane stress
ep = [ptype,t]

D1 = cfc.hooke(ptype, E1, ny1)
D2 = cfc.hooke(ptype, E2, ny2)