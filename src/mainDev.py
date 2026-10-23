import cantera as ct


CombustionGas = ct.Solution("Dagaut_HP.yaml")
CombustionGas.TPY = 3000, (101000*30), "O2: 1"
print(ct.gas_constant/CombustionGas.mean_molecular_weight)

CombustionGas.equilibrate("TP")
print(ct.gas_constant/CombustionGas.mean_molecular_weight)
#print(CombustionGas.report())
