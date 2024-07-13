"""
Adds heterogeneity to cell parameters
"""
from scipy.stats import truncnorm
from model.cell_model import CellModel

#method for generating truncated normal distribution of random numbers
def get_truncated_normal(mean=0, sd=1, low=0.90, upp=1.10):
    """
    Generates numbers from a truncated normal distribution
    """
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

def cell_heterogeneity(cell: CellModel):
    """
    Sets cell heterogeneity
    """
    cell.Astretch = cell.Astretch*get_truncated_normal().rvs()
    cell.ksscc = cell.ksscc*get_truncated_normal().rvs()
    cell.Jleak = cell.Jleak*get_truncated_normal().rvs()
    cell.Gtot = cell.Gtot*get_truncated_normal().rvs()
    cell.kip3r3 = cell.kip3r3*get_truncated_normal().rvs()


