def tem(V,T_out):
    import numpy as np
    Cp = 5193
    Cv = 3116
    pi = np.pi
    P_0 = 2*10**6
    R = 2078.5
    T_h = 300
    T_c = 77

    V_c0 = 0.0001
    V_cs = 0.1*V_c0
    V_e0 = 0.00005
    V_es = 0.1*V_e0
    V_b0 = 0.00005
    V_bs = 0.1*V_b0
    f = 50
    nt = 360
    dt = 1/f/nt
    
