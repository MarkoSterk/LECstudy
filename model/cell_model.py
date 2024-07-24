"""
Cell numerical model
"""

# pylint: disable=R0903
# pylint: disable=R0902
# pylint: disable=R0913
# pylint: disable=R0914
# pylint: disable=C0103
import numpy as np

from .model_parameters import ModelParameters

class CellModel(ModelParameters):
    """
    Cell model parameters
    """
    cell_height: float = 8.0 #cell height in um

    #Stretch parameters
    Astretch=1.2*0.3426 #amplitude of initial stretch
    stretch_exp=1.0*0.105 #stretch exponent (exp. decay)

    #Flux through stretch sensitive channels
    k_stretch=0.12 #0.08105
    kf=0.20 #0.1382
    kb=0.025 #0.04027
    #ksscc=1.15*1.025 #maksimal current through stretch sensitive channels

    ##Stretch induced phosphorylation
    Aalpha4=1.15*0.0282
    exp_alpha4=0.95*0.0138 #phosphorylation exponent (exp. decay)

    ###Difusion of Ca2+ and IP3 through gap junctions
    Dca=512.7 #difusion constant of Ca2+ through GJ
    Dip3=913.9 #difusion constant for IP3 through GJ

    factor_diff=Dca/Dip3
    Cgjip3=10.0*0.025#0.025  #apparent constant for IP3 difusion
    Cgjca=factor_diff*Cgjip3 #apparent constant for Ca2+ difusion

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

    # Compartments

    def __init__(self, cell_num: int,
                volume: float,
                dist_from_stimulation: float,
                parameters: dict,
                neighbours: list = None):

        self.cell_num = cell_num
        self.volume = volume
        self.dist_from_stimulation = dist_from_stimulation

        if neighbours:
            self.neighbours = neighbours

        for key, value in parameters.items():
            setattr(self, key, value)

        self.time_of_activation = None
        self.simulation_time = []
        self.calcium_time_series = []
        self.ip3_time_series = []
        self.atp_time_series = []

        self.C: float = 0.125
        self.P: float = 0.01
        self.Osscc: float = 0.0
        self.stretch: float = 0.0
        self.Rs: float = 17000.0
        self.Rsp: float = 0.0
        self.G: float = 14.0
        self.PIP: float = 49997.0
        self.k_ip3_in: float = 0.0
        self.k_ca_out: float = 0.0

        self.CA_OUT: float = 0.0 #amount of added Ca2+ into stimulated cell(s)
        self.IP3_IN: float = 0.5771 #amount of added IP3 into stimulated cells

    def run_model_step(self, time: float,
                        step: int,
                        jgjca: float,
                        jgjip3: float,
                        L: float,
                        init_cells: list):
        """
        Calculates model step values
        """
        ###Saves binarized calcium time series
        ###Saves state of calcium to list
        if step%self.record_every == 0:
            self.simulation_time.append(time)
            self.calcium_time_series.append(self.C)
            self.ip3_time_series.append(self.P)
            self.atp_time_series.append(L)

        ##Saves activation time of cell
        time_step = self.record_every * self.dt ## time between two data points
        points_for_slope = int(self.time_interval_for_slope/time_step)
        if len(self.calcium_time_series)>points_for_slope+5 and not self.time_of_activation:
            slope = (self.calcium_time_series[-1] - self.calcium_time_series[-points_for_slope])/(points_for_slope*time_step)
            if slope>self.slopeTh:
                self.time_of_activation = time - int(points_for_slope/2)*time_step

        dstretch = None
        # Stimulation
        if self.time_of_stimulation <= time < (self.time_of_stimulation + self.stim_dur):
            self.stretch = self.Astretch * \
                np.e**(-self.stretch_exp*self.dist_from_stimulation)
            if self.cell_num in init_cells:
                self.k_ip3_in = self.IP3_IN/self.stim_dur
        else:
            dstretch = (-self.k_stretch*self.stretch)*self.dt
            self.k_ip3_in = 0.0

        alpha4 = self.Aalpha4 * np.e**(self.exp_alpha4*self.dist_from_stimulation)

        jpump = (self.Vpump*self.C**2) / (self.Kpump**2+self.C**2)

        w_inf = ((1.0+((self.Ka/self.C)**4)+((self.C/self.Kb)**3)) /
            (1.0+(1.0/self.Kc)+((self.Ka/self.C)**4)+((self.C/self.Kb)**3)))

        pryr = ((w_inf*(1.0+((self.C/self.Kb)**3))) /
            (1.0+((self.Ka/self.C)**4)+((self.C/self.Kb)**3)))

        jryr = self.kryr*pryr

        k1 = (self.alpha1*(self.C**3))/((self.beta1**3)+(self.C**3))

        k4 = (alpha4*self.P)/(self.beta4+self.P)

        phi = (1.0)/(1.0+((self.k2)/(self.k3+k4))*(1.0+k4/self.k5))

        O = (phi*self.P)/(((self.k1m+self.k2)/(k1))*phi+self.P)

        jip3r3 = self.kip3r3*O**4

        dOsscc = (self.stretch*self.kf-(self.stretch * self.kf+self.kb)*self.Osscc)*self.dt

        jsscc = self.ksscc*self.Osscc

        ca_noise =  self.ca_noise_amp*np.random.uniform(-1, 1)
        dC = (jgjca+self.Jleak-jpump+jryr +
                   jip3r3+jsscc+self.k_ca_out)*self.dt + ca_noise

        pr = (L*self.Rs)/(self.epsilon*self.Rtot*(self.K1+L))

        dG = (self.ka*(self.sigma+pr) * (self.Gtot-self.G)-self.kd*self.G)*self.dt

        dRsp = (L*(((self.kp*self.Rs)/(self.K1+L)) -((self.ke*self.Rsp)/(self.K2+L))))*self.dt

        dRs = (self.kr*self.Rtot-(self.kr+((self.kp*L) /
                    (self.K1+L)))*self.Rs-self.kr*self.Rsp)*self.dt

        rh = self.alpha*((self.C)/(self.K3+self.C))*self.G

        p_noise = self.ip3_noise_amp*np.random.uniform(-1, 1)
        dP = (rh*((1.0)/(self.Na*self.volume))*(10**(6)) *
                   self.PIP-self.kdeg*self.P+jgjip3+self.k_ip3_in)*self.dt + p_noise

        dPIP = (self.rr*self.PIP2tot-(rh+self.rr)*self.PIP -
                     self.rr*self.Na*self.volume*self.P*0.000001)*self.dt

        if dstretch:
            self.stretch+=dstretch
        self.Osscc += dOsscc
        self.C += dC
        self.G += dG
        self.Rsp += dRsp
        self.Rs += dRs
        self.P += dP
        self.PIP += dPIP
