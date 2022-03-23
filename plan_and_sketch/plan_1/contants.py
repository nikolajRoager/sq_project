import numpy as np

m_p = 1.672621911E-27 #kg
q_p = 1.602176634E-19 #C
I   = 5 # A
N   = 1000
v   = 319515.54756336455 #m/s
mu0 = 4*np.pi*1E-7
eps0 =  8.8541878128e-12;

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
print('1 T  =',kg/(s*C),' u/(μs e)')
print('1 V/m  =',kg/((s**2)*C),' m u/(μs^2 e)')

print('m_p  =',m_p,'kg ',m_p*kg,' u')
print('q_p  =',q_p,'C ' ,q_p*C,' e')
print('v  =',v,'m/s ',v/s,' m/(μs)')
print('ω  =',omega,' 1/s ',omega/s,' 1/(μs)')
print('T  =',2*np.pi*1/omega,' s ',2*np.pi*s/omega,' (μs)')
print('r  =',r,' m ',(m_p*kg*v/s)/(q_p*C*B*T),' m')
print('B  =',B,' T ',B*T,' u/(μs e)')
print('μ0  =',mu0,' kg m/C^2 ',mu0*kg/(C**2),'  u m/(e^2)')
print('eps0  =',eps0,' S^2 C^2 /(m^3kg) ',eps0* (s**2) *(C**2)/(kg),' μs^2 e^2 /(m^3 u) ')

B=3.12e-5

print('B  =',B,' T ',B*T,' u/(μs e)')

omega = q_p*B/m_p;

print('ω  =',omega,' 1/s ',omega/s,' 1/(μs)')
print('T  =',2*np.pi/omega,' s ',2*np.pi*s/omega,' (μs)')

print((m_p*kg*v/s)/(q_p*C*B*T))

B=B/1000;

omega = q_p*B/(m_p);

print('ω  =',omega,' 1/s ',omega/s,' 1/(μs)')
print('T  =',2*np.pi/omega,' s ',2*np.pi*s/omega,' (μs)')


print((m_p*kg*v/s)/(q_p*C*B*T))

print(r*np.sin(0*np.pi/180));
print(r*np.sin(15*np.pi/180));
print(r*np.sin(30*np.pi/180));
print(r*np.sin(45*np.pi/180));
print(r*np.sin(60*np.pi/180));
print(r*np.sin(75*np.pi/180));
print(r*np.sin(90*np.pi/180));
