"""
Experiments meta data
"""
####Lists of experiment names
from typing import Tuple


CONTROL_EXP_NAMES = [
    '02072019002002CN',
    '02072019007002CN',
    '03092019001002C',
    '03092019003002N',
    '03092019004002N',
    '8072019001001C',
    '08072019002001C',
    #'08072019002001CN',
    '09072019003001N',
    '09092019003002N',
    '10092019001002C',
    '10092019002001C',
    '10092019004001N',
    '15072019004001N',
    '16072019001002N',
    '16092019005002CN',
    '17092019002003CN',
    '17092019003001C',
    '17092019004001C',
    '22072019002002CN',
    '23092019001002N'
]

APIRAZA_EXP_NAMES = [
    '02072019002002apirazaCN',
    '08072019001005apirazaC',
    '10092019002001apirazaC',
    '10092019003001apirazaC',
    '15072019004004apirazaN',
    '17092019002003apirazaCN',
    '17092019003002apirazaC',#izloci?
    '17092019004001apiraza',#izloci?
    '22072019002002apirazaCN',
    '23092019001004apirazaN',
    '23092019002003apirazaN'
]

CBX_EXP_NAMES = [
    '03092019001001cbxC',
    '03092019003003cbxN',
    '09072019003001cbxN',
    '16072019001001cbxN'
]

###Analysis configurations
CONTROL_TS_CONFIGS = {
    'AMP_FACT': 0.5,
    'DISTANCE': 5,
    'WIDTH': 3,
    'PROMINENCE': 0.4,
    'REL_HEIGHT': 0.85
}

APIRAZA_TS_CONFIGS = {
    'AMP_FACT': 0.7,
    'DISTANCE': 5,
    'WIDTH': 5,
    'PROMINENCE': 0.3,
    'REL_HEIGHT': 0.75
}

CBX_TS_CONFIGS = {
    'AMP_FACT': 0.3,
    'DISTANCE': 15,
    'WIDTH': 10,
    'PROMINENCE': 0.6,
    'REL_HEIGHT': 0.75
}

###Method for selecting the appropriate experiment names (list)
###along with analysis configurations
def select_exp_names_and_configs(exp_type: str) -> Tuple[list, dict]:
    """
    Returns selected exp names list
    """
    if exp_type == 'control':
        return CONTROL_EXP_NAMES, CONTROL_TS_CONFIGS
    if exp_type == 'apiraza':
        return APIRAZA_EXP_NAMES, APIRAZA_TS_CONFIGS
    if exp_type == 'cbx':
        return CBX_EXP_NAMES, CBX_TS_CONFIGS

    raise Exception('This experiment type is not available')
