# Calculates the needed engine geometry paramaters based on the given infromation in the input file
"""
@author: phineas
"""

# External librarie
import math as m
import numpy as np
import CoolProp.CoolProp as CP

# Internal librarie
from ..Helpers.PropName import conCool

# conversion factors
d2r = m.pi / 180
psi2pa = 6894.71
in2m = 0.0254

# Engineometry object that can be constructed by stating a methog, ie, general, girssom, rao, and feeding it the needed inputs
# generalConstructorDemo
# ((0, 5, "5"), (), (), (), ())


# Grissom values must be input in the correct order
# (ExpRatio, ConRatio, Dt, Rcont, ThetaCont, Rexp, ThetaExp, Lchamber)
class EngineGeo:
    def __init__(self, method: str, inputs: tuple, cellnum: int):
        self.method = method
        self.inputs = inputs
        self.cellnum = cellnum

        # the first method here is "general", it uses an arbirary tuple of tuples of the form (x0, x1, func), to define a piecewise function
        # each list in the tuple should be a sequential piece of the piecewise function such that the x0 and x1 are the intervol and func is the function
        # the function is a string that is determined in python equation notation
        if self.method == "general":
            self.axval = np.linspace(0, inputs[-1][2], self.cellnum)
            fv = np.vectorize(self.generalfx)
            self.RatL = fv(self.axval)

            # create a higher resolution version of the point field for determining stuff if needed
            highdefL = np.linspace(0, inputs[-1][0], 10000)
            highdefR = fv(highdefL)

            # we now need some atributes from this
            self.Dt = min(highdefR)

        # grissom uses a method outlined in Grissom liquid film cooling.
        if self.method == "grissom":
            # we need LT to get declared and its over in the function for this, annoyingly, as I dont want to delcare all of those dummy variables
            self.grissomfx(0)

            # Produce the needed RatL
            self.axval = np.linspace(0, self.LT, self.cellnum)
            
            #Construct radius values with a factorized function
            fv = np.vectorize(self.grissomfx)
            self.RatL = fv(self.axval)
            
            #get various dimensions from that value that are important
            self.Rt = min(self.RatL)
            self.Dt = self.Rt*2
            self.Lt = self.axval[np.where(self.RatL == self.Rt)[0]]

            #exit radius
            self.Re = self.RatL[-1]
            
            
    # This is the function call that should be used externally the other ones will not work well if the other defentition is used
    def RatLf(self, x: float):
        if self.method == "genearl":
            return self.generalfx(x)

        if self.method == "grissomfx":
            return self.grissomfx(x)

    # this is the function that gets used to determine radius at some length for the general case
    def generalfx(self, x: float):
        # find what intervol of the piecewise function the current value is in
        i: int = 0
        while (
            not (self.inputs[i][0] <= x)
            and not (self.inputs[i][1] >= x)
            and i < len(self.inputs)
        ):
            i += 1
        section = self.inputs[i]

        # return the section evaluated at x
        return eval(section[2], {}, {"x": x})

    # This function uses the defined variables to compute the balues for grissom
    def grissomfx(self, x: float):
        # if this method is getting used we precompute the needed values and unpackage the input
        ExpRatio = self.inputs[0]
        ConRatio = self.inputs[1]
        Dt = self.inputs[2]
        Rcont = self.inputs[3]
        ThetaCont = self.inputs[4]
        Rexp = self.inputs[5]
        ThetaExp = self.inputs[6]
        Lchamber = self.inputs[7]

        # needed variables for the function to work must declare before we call the function
        Dc = Dt * m.sqrt(ConRatio)
        De = Dt * m.sqrt(ExpRatio)
        D2 = Dc - 2 * Rcont * (1 - m.cos(ThetaCont * d2r))
        D3 = Dt + 2 * Rexp * (1 - m.cos(ThetaCont * d2r))
        D5 = Dt + 2 * Rexp * (1 - m.cos(ThetaExp * d2r))
        L2 = Lchamber + Rcont * m.sin(ThetaCont * d2r)
        L3 = L2 + ((D2 - D3) / (2 * m.tan(ThetaCont * d2r)))
        Lt = L3 + Rexp * m.sin(ThetaCont * d2r)
        L5 = Lt + Rexp * m.sin(ThetaExp * d2r)
        self.LT = L5 + ((De - D5) / (2 * m.tan(ThetaExp * d2r)))

        # Start of contour calculation
        if x < Lchamber:
            return 0.5 * (Dc)
        if x < L2:
            return 0.5 * (Dc - 2 * (Rcont - m.sqrt(Rcont**2 - (x - Lchamber) ** 2)))
        if x < L3:
            return 0.5 * (D2 - 2 * (x - L2) * (m.tan(ThetaCont * d2r)))
        if x < L5:
            return 0.5 * (Dt + 2 * (Rexp - m.sqrt(Rexp**2 - (x - Lt) ** 2)))
        else:
            return 0.5 * (D5 + 2 * (x - L5) * (m.tan(ThetaExp * d2r)))


# inputs here should be tuples such that (fueltype, fuelTankTemp, fuelTankPressure)
class Tankage:
    def __init__(self, fuel: tuple, ox: tuple):
        # Set props from the fuel input tuplele
        self.fuel = fuel[0]
        self.fuelT = fuel[1]
        self.fuelP = fuel[2]

        # sets props from the ox input tuple
        self.ox = ox[0]
        self.oxT = ox[1]
        self.oxP = ox[2]


# Rocket engine injectorGeometry object, stores relevent info for injector
class InjectorGeo:
    # The tuples for fuel, ox, and blc holes, (HoleDiameter, HoleNumber, Cd, Angle)
    # The angle is measures such that 90 degres points directly at the engine wall
    def __init__(self, fuel: tuple, ox: tuple, BLC: tuple = (0, 0, 0, 0)):
        self.fuel = fuel
        self.ox = ox
        self.BLC = BLC

    # Gets the mass flow for the propellents should be upgraded to use the Dyer model for two phase flow so nutrus can be use
    # This functionr eturns a tuple of mass flows for (fuel, ox, BLC)
    def GetMassflow(self, Tankage: Tankage, pc: float):
        # Currently this just uses standard SPI for ease
        sizes = np.array([self.fuel[0], self.ox[0], self.BLC[0]])*in2m
        tankpress = np.array([Tankage.fuelP, Tankage.oxP, Tankage.fuelP])*psi2pa
        pcs = np.array([pc, pc, pc])
        numbers = np.array([self.fuel[1], self.ox[1], self.BLC[1]])
        Cds = np.array([self.fuel[2], self.ox[2], self.BLC[2]])

        
        # Gets the densities of the propellents in thanks
        FD = CP.PropsSI(
            "D", "P", Tankage.fuelP*psi2pa, "T", Tankage.fuelT, conCool(Tankage.fuel)
        )
        OxD = CP.PropsSI("D", "P", Tankage.oxP*psi2pa, "T", Tankage.oxT, conCool(Tankage.ox))
        densities = np.array([FD, OxD, FD])

        Areas = m.pi * (sizes / 2) ** 2
        MassflowRates = (
            Cds * Areas * np.sqrt(2 * densities * (tankpress - pcs))
        ) * numbers
        
        return MassflowRates


# Geo = EngineGeo("grissom", (8.4143, 20.25, 1, 0.75, 35, 0.5, 15, 5), 250)
# Geo = EngineGeo("grissom", (3, 9, 1, 0.5, 45, 0.5, 15, 3), 250)
# print(Geo.RatL)

# plt.plot(Geo.axval, Geo.RatL)
# plt.show()
