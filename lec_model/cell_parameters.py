"""
Cell numerical model
"""

# pylint: disable=R0903
# pylint: disable=R0902
# pylint: disable=R0913
# pylint: disable=R0914
# pylint: disable=C0103

class CellParameters:
    """
    Cell parameters
    """
    cell_height: float = 8.0 #cell height in um

    #Stretch parameters
    Astretch=1.2*0.3426 #amplitude of initial stretch
    stretch_exp=0.25#1.0*0.105 #stretch exponent (exp. decay)

    #Flux through stretch sensitive channels
    k_stretch=0.12 #0.08105
    kf=0.20 #0.1382
    kb=0.025 #0.04027
    #ksscc=1.15*1.025 #maksimal current through stretch sensitive channels

    ##Stretch induced phosphorylation
    Aalpha4=1.15*0.0282
    exp_alpha4=0.95*0.0138 #phosphorylation exponent (exp. decay)

    #Leakage current from ER and cytoplasm
    Jleak=0.1450

    #SARCO/Endoplasmic reticulum ATPase (SERCA) and plasma membrane Ca ATPase (PMCA) flux
    Vpump=5.341
    #Kpump=0.5030

    #Rynodine receptor dynamics
    Ka=0.37224
    Kb=0.63601
    Kc=0.05714
    #kryr=16.04 #maximal current through rynodine receptor sensitive channels

    ###IP3R3 receptor dynamics
    alpha1=40.0
    beta1=0.8
    beta4=0.01
    k2=0.5
    k3=0.5
    k1m=0.88
    k5=0.02
    #kip3r3=155.8 #maximal current through IP3 sensitive channels

    #PIP2 and IP3 dynamics and P2Y2 receptor dynamics
    PIP2tot=50000.0
    Rtot=20000.0
    kr=0.000175
    kp=0.03
    K1=6.5 #5.0
    ke=0.006
    K2=100.0
    ka=0.017
    sigma=0.001238
    Gtot=100000.0
    epsilon=0.85
    kd=0.15

    alpha=0.00002781
    rr=0.015
    K3=0.4
    #kdeg=1.25 #degradation rate of IP3

    Na=6.02252*10**23 #Avogadro's constant
