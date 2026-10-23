# Internal inputs
from .EngineGeometry import EngineGeo, InjectorGeo, Tankage
from .cea import runCEA, MachNumber
from .Bisect import Bisect

# External imports
import math as m
from scipy.optimize import *

# This class constructes an engine object from key engine charactaristics like
# Fuel and Ox, engine geometery, and an injector object, or certan cirtical chamber dimensions and or thrust values and it will assume ideal geometery


class Engine:
    # Props will be taken from engine geo if pcritical geo is empty
    # This will have two execution modes, one from engine geomeyry, fuel and config, and one from ideal engine solving
    def __init__(
        self,
        tankage: Tankage,
        solveType: str,
        injector: InjectorGeo | None = None,
        geometry: EngineGeo | None = None,
        cstarEff: float = 1
    ):
        self.cstarEff = cstarEff
        self.tankage = tankage
        self.solveType = solveType

        # This is the solver for solving initial values from a given engine geometry, it calles a recompute function
        # this mean that it can be called from a provided new geometery if that happpens
        if self.solveType == "Geometry":
            self.injector: InjectorGeo = injector
            self.geometry: EngineGeo = geometry
            
            #solve with the geometry now
            self.GeoSolveSimple()

    def GeoSolveSimple(self):
        #we get the combustion results here [gas, m_dot] for a cstar of 1
        self.ceaResults = runCEA(self.tankage, self.geometry, self.injector)

        #if the cstar is not 1 we have more work to do
        if not (self.cstarEff == 1):
            #becasue we have the mass flow expected and ww know what the eff is we can compute the expected specific impulse
            #we can then solve for an of offset that will produce this value, and get the expected specific impulse.
            
            #we need to first compute the altered specific impulse, assuming the nozzle is the same\
            #first we need to compute the areas at the throat and exit, and gamma
            gamma = self.ceaResults[0].cp/self.ceaResults[0].cv
            R = self.ceaResults[0].cp
            
            Dt = self.geometry.Dt
            De = self.geometry.Re*2
            Ae = 0.25*m.pi*De**2

            Me = root_scalar(MachNumber, bracket=(1.001, 5), args=(Ae, Dt, gamma), xtol=(10**-10), method="brentq").root

            #We also need exit temperature 
            Tc = self.ceaResults[0].T
            Te = Tc*(1+(0.5*(gamma-1))*Me**2)**(-1)

            #final specifc impulse
            Isp = Me*m.sqrt(gamma*R*Te)

            #This is the targat spefic impulse from the cstar efficiency
            IspTrue = self.cstarEff*Isp
            
            #This will probably have a final value somewhere between 0-1, as this would alter OF to be lower
            OFoffset = root_scalar(IspD, bracket=(0.01, 1), args=(self.tankage, self.geometry, self.injector, IspTrue), xtol=(10**-6), method="brentq").root

            #now that we have the offset we can get the performance specifications using this value
            self.ceaResults = runCEA(self.tankage, self.geometry, self.injector, OFoffset)
            
        #now we get axial properties and assign them to variables
        



#function that uses the geometry cea call and returns specific impulse 
def IspD(OFoffset: float, tankage: Tankage, geometry: EngineGeo, injector: InjectorGeo, IspTargaet):
    #get the combustion props from the inputs
    ceaResults = runCEA(tankage, geometry, injector, OFoffset)
    
    #compute Isp as described in GeoSolveSimple
    gamma = ceaResults[0].cp/ceaResults[0].cv
    R = ceaResults[0].cp
            
    Dt = geometry.Dt
    De = geometry.Re*2
    Ae = 0.25*m.pi*De**2

    Me = root_scalar(MachNumber, bracket=(1.001, 5), args=(Ae, Dt, gamma), xtol=(10**-10), method="brentq").root

    Tc = ceaResults[0].T
    Te = Tc*(1+(0.5*(gamma-1))*Me**2)**(-1)

    Isp = Me*m.sqrt(gamma*R*Te)
    
    return (Isp-IspTargaet)