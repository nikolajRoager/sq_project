import numpy as np

m_p = 1.672621911E-27 #kg
q_p = 1.602176634E-19 #C
I   = 5 # A
N   = 1000
v   = 319515.54756336455 #m/s
mu0 = 4*np.pi*1E-7

B= mu0*I*N;
omega = q_p*B/m_p;
r =m_p*v/(q_p*B);

kg = 6.0221412901167394E26
s = 1e6
C = 1/(1.602176634E-19)
T = kg/(s*C)

print('1 kg =',kg,' u')
print('1 C  =',C,' e')
print('1 s  =',s,' μs')
print('1 T  =',kg/(s*C),' kg/(μs e)')

print('m_p  =',m_p,'kg ',m_p*kg,' u')
print('q_p  =',q_p,'C ' ,q_p*C,' e')
print('v  =',v,'m/s ',v/s,' m/(μs)')
print('ω  =',omega,' 1/s ',omega/s,' 1/(μs)')
print('T  =',1/omega,' s ',s/omega,' (μs)')
print('r  =',r,' m ',(m_p*kg*v/s)/(q_p*C*B*T),' m')
print('B  =',B,' T ',B*T,' u/(μs e)')
print('μ0  =',mu0,' kg m/C^2 ',mu0*kg/(C**2),'  u m/(e^2)')
