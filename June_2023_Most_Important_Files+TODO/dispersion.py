#Import Statements & Define Constants
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats as stats
from scipy.stats import lognorm
from scipy.stats import ks_2samp
import math
import rebound
import random
from numpy.random import seed, random
from scipy.stats import rayleigh
from scipy.stats import norm
import itertools
from spock import FeatureClassifier
from decimal import Decimal
import pandas as pd
from tqdm import tqdm
from matplotlib.ticker import EngFormatter
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import seaborn as sns
import statistics

def check_crossing(sim):
    ps = sim.particles
    for i1 in range(1,sim.N-1):
        i2 = i1+1 # next planet out 
        EMx = ps[i1].e*np.cos(ps[i1].pomega) - ps[i2].e*np.cos(ps[i2].pomega)
        EMy = ps[i1].e*np.sin(ps[i1].pomega) - ps[i2].e*np.sin(ps[i2].pomega)
        EM = np.sqrt(EMx**2 + EMy**2)
        EMcross = (ps[i2].a-ps[i1].a)/ps[i2].a
        if EM > EMcross:
            return True
    return False

# rori function
def p_ratios(P=None, logP=None, sim=None, log=True):
    """returns either the log or the normal period ratios of a system
       Param: a list of periods, log a, a, or a rebound sim
    """
    if logP != None:
        logpr = [logP[i+1]-logP[i] for i in range(len(logP)-1)]
        if log:
            return logpr
        else:
            return [10**x for x in logpr]
    if sim != None:
        ps = sim.particles[1:]
        P = [ps[i].P for i in range(len(ps))]   

    if log:
        return [np.log10(P[i+1]/P[i]) for i in range(len(P)-1)]
    else:
        return [P[i+1]/P[i] for i in range(len(P)-1)]
    
def dispersion(trios=None, systems=None, filter4=False, ddof=1): # input systems: orbital periods
    """returns the dispersion and error of a list of 
       trios as a tuple.
       For a single trio, only returns the dispersion
    """
    if systems != None:
        trios = get_trios(systems, filter4=filter4)
    
    if len(trios)==0:
        return 0, 0
    
    rel_var = []
    for t in trios:
        logpr = p_ratios(P=t, log=True)
        mu = np.mean(logpr)
        rel_var += [np.var(logpr)/mu**2]
    D = np.sqrt(np.mean(rel_var))
    if len(trios)>1: 
        error = np.std(rel_var, ddof=ddof)/np.sqrt(len(trios))#/(2*D)
    else:
        error = 10
    return D, error

if __name__ == '__main__':
    raise NotImplementedError
