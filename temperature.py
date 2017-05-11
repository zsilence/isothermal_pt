def temperature(P1,P2,V1,V2,T1,T_f):
    import numpy as np
    Cp = 5193
    Cv = 3116
    pi = np.pi
    P_0 = 2*10**6
    R = 2078.5
    T_h = 300
    T_c = 77
    T_w = T_f
    f = 50
    nt = 360
    dt = 1/f/nt
    A = 0.25*pi*0.08**2
    alpha = 1000

    T2 = (Cv*(P2*V2/R-P1*V1/R)+P2*V2/R*Cp-alpha*A*dt*T_w+P1*(V2-V1))/(P1*V1/R/T1*Cp-alpha*A*dt)
    if P2*V2/R/T2 > P1*V1/R/T1:
        a = alpha*A*dt
        b = -(P1*V1/R/T1*Cp*T_f+alpha*A*dt*T_w-Cv*(P2*V2-P1*V1)/R-P1*(V2-V1))
        c = P2*V2/R*Cp*T_f
        T2 = (-b+np.sqrt(b**2-4*a*c)/2*a)
    return T2
