"""
Initial conditions of the numerical model
"""
# pylint: disable=R0903



class InitialConditions:
    """
    Initial conditions of model
    """
    C: float = 0.125
    P: float = 0.01
    Osscc: float = 0.0
    stretch: float = 0.0
    Rs: float = 17000.0
    Rsp: float = 0.0
    G: float = 14.0
    PIP: float = 49997.0
    k_ip3_in: float = 0.0
    k_ca_out: float = 0.0

    CA_OUT: float = 0.0 #amount of added Ca2+ into stimulated cell(s)
    IP3_IN: float = 0.5771 #amount of added IP3 into stimulated cells
