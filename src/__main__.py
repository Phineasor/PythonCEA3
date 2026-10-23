import math as m
import numpy as np
import cantera as ct
import CoolProp.CoolProp as CP

from .Utils.cea import runCEA, AxialValues, basicEngineChar
#from .Utils.EngineGeometry import R:wq

from .Inputs import InputValues as IV
from .Utils import Injector as Inj
from .Helpers import GraphGen
from .Helpers.PropName import conCool, conCant
from .Utils.EngineGeometry import EngineGeo, Tankage, InjectorGeo
from .Utils.Engine import Engine

# from .Utils.BLC2 import calcBLC2 #Temporarily Disabled
# from .Utils.Meshing import MeshGen as MG

in2m = 0.0254
psi2pa = 6894.71



# for i in range(IV.CellNum):
# print(AxialDistances[i])


if IV.Physical:
    #setup tankage and geometry
    tankage = Tankage((IV.Fuel, IV.FuelTankT, IV.FuelTankP), (IV.Ox, IV.OxTankT, IV.OxTankP))
    injector= InjectorGeo((IV.FuelOrificeDiameter, IV.FuelOrificeNum, IV.FuelOrificeCd, 0), 
                          (IV.OxOrificeDiameter, IV.OxOrificeNum, IV.FuelOrificeCd, 0), 
                          (IV.BLCOrificeDiameter, IV.BLCOrificeNum, IV.BLCOrificeCd))
                          
    contour = EngineGeo(method = "grissom", inputs = (IV.ExpRatio, IV.ConRatio, IV.Dt, IV.Rcont, IV.ThetaCont, IV.Rexp, IV.ThetaExp, IV.Lchamber), cellnum = IV.CellNum)

    Engine = Engine(tankage, "Geometry", injector, contour, IV.cStarEff)
    
    RadiusVal = contour.RatL
    
    print(
        "--------------------------------------------------------------------------------"
    )
    print("Physical driven rocket engine values")
    ceaOut = Engine.ceaResults
    print(ceaOut[0].report())
    γ = ceaOut[0].cp / ceaOut[0].cv
    R = ct.gas_constant / ceaOut[0].mean_molecular_weight
    val = AxialValues(ceaOut[0].T, ceaOut[0].P, ceaOut[0].density, ceaOut[0], contour)
    print("ExitPressure: " + str(val[1][IV.CellNum - 1]))

    v = (γ * R * val[0][-1]) ** 0.5 * val[4][-1]
    
    Mdots = injector.GetMassflow(tankage, ceaOut[0].P)
    FuelMdot = Mdots[0]
    OxMdot   = Mdots[1]
    BLCMdot  = Mdots[2]

    print("OF: " + str(OxMdot / FuelMdot))
    print("FuelMdot: " + str(FuelMdot))
    print("OxCMdot: " + str(OxMdot))
    print("FOMdot: " + str(OxMdot + FuelMdot))
    print("BLCMdot: " + str(BLCMdot))
    print("BLCMdot/BLCMdot+F: " + str(BLCMdot / (FuelMdot + BLCMdot)))
    print("pc: " + str(ceaOut[0].P * 0.000145038))
    print("Tc: " + str(ceaOut[0].T))
    print("Dc: " + str(ceaOut[0].density))
    print("Te: " + str(val[0][-1]))
    F = v * (FuelMdot + OxMdot) * 0.9 + (
        (m.pi * (RadiusVal[IV.CellNum - 1] * 0.0254) ** 2)
        * (val[1][IV.CellNum - 1] - 101000)
    )

    print("F : " + str(F))
    print("Ue: " + str(v))
    print("Cp: " + str(ceaOut[0].cp))
    print("Isp [S]: " + str(F / ((FuelMdot + OxMdot) * 9.81)))
    print("Viscosity: " + str(ceaOut[0].viscosity))
    print("thermal_conductivity: " + str(ceaOut[0].thermal_conductivity))
'''
print(
    "--------------------------------------------------------------------------------"
)
print("Standard Chamber Size, Pressure and Thrust driven engine paramaters")

List = basicEngineChar(
    IV.Pc, IV.AmbP, IV.OF, IV.F
)  # gets all of the known values and uses them to produce the needed engine params
print("Tc [K]: " + str(List[6]))
print("cStar: " + str(List[0]))
print("Thrust Coefficient: " + str(List[1]))
print("Total Mass flow rate [Kg/s]: " + str(List[2]))
MdotF = (1 / (1 + IV.OF)) * List[2]
MdotO = (IV.OF / (1 + IV.OF)) * List[2]
print("Fuel Mass flow [Kg/s]: " + str(MdotF))
print("Ox Mass flow rate [Kg/s]: " + str(MdotO))
print("Throat Area [m^2]: " + str(List[3]))
print("Throat diameter [in]: " + str(m.sqrt(List[3] / m.pi) * 2 / in2m))
print("ExpRatio: " + str(List[4]))
print("Exit diameter [in]: " + str(m.sqrt(List[5] / m.pi) * 2 / in2m))
Vc = (List[3] / (in2m**2)) * IV.LStar
print("Chamber Length [in]: " + str(Vc / (m.pi * (IV.Dc / 2) ** 2)))
F = IV.Pc * List[1] * List[3]
print("Thrust [N]: " + str(F))
print("Isp [s]: " + str(F / ((MdotF + MdotO) * 9.81)))

print(
    "--------------------------------------------------------------------------------"
)
# need oriface area
sqrtF = m.sqrt(
    2
    * CP.PropsSI("D", "T", IV.FuelTankT, "P", IV.Pc, conCool(IV.Fuel))
    * (IV.FuelTankP * psi2pa - IV.Pc)
)
sqrtO = m.sqrt(
    2
    * CP.PropsSI("D", "T", IV.OxTankT, "P", IV.Pc, conCool(IV.Ox))
    * (IV.OxTankP * psi2pa - IV.Pc)
)
AF = MdotF / (0.7 * sqrtF)
AO = MdotO / (0.7 * sqrtO)
print("Area of fuel oriface [in^2]: " + str(AF / (in2m**2)))
print("Area of ox   oriface [in^2]: " + str(AO / (in2m**2)))
# print(CP.PropsSI("D", "T", IV.FuelTankT, "P", IV.Pc, conCool(IV.Fuel)))
# print(CP.PropsSI("D", "T", IV.OxTankT, "P", IV.Pc, conCool(IV.Ox)))


# print(CP.PropsSI("T", "P", CombustionGas.P, "Q", 1, "Ethanol"))
# print(CP.PropsSI("T", "P", CombustionGas.P, "Q", 0.5, "Ethanol"))
# print(ceaOut[0]())

# print(AxialValues(CombustionGas.T, CombustionGas.P, CombustionGas.density, CombustionGas)[3])
# print(ceaOut[0].report())
# print(str(Pc/Inj.PSI2PA)+" : "+str(Pc))
# print("Pr = "+str(Pr))
# print("mu = "+str(CombustionGas.viscosity))
# print("Cp = "+str(CombustionGas.cp))
# print(CombustionGas.thermal_conductivity)
# print("Ox = "+str(264.172*10*OxMdot/CP.PropsSI("D", "T|liquid", IV.OxTankT, "P", IV.OxTankP*psi2pa, "O2")))
# print("Fuel = "+str(264.172*10*(FuelMdot+BLCMdot)/CP.PropsSI("D", "T|liquid", IV.FuelTankT, "P", IV.FuelTankP*psi2pa, "Ethanol")))
# print("Fuel gal/s = "+str(264.172*(FuelMdot+BLCMdot)/CP.PropsSI("D", "T|liquid", IV.FuelTankT, "P", IV.FuelTankP*psi2pa, "Ethanol")))
'''