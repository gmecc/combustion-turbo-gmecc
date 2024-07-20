# combustion-turbo-gmecc
Program for calculating the fuel combustion process in the combustion chambers of gas turbine plants

The program is based on the liquid properties library `pyfluids`.

### Units
- temperature - degree Celsius _(°C)_;
- absolute pressure _(Pa)_

Parameter values ​​are specified as a tuple. For secondary parameters, the calculated range is determined in the tuple.

## Расчет калориметрической температуры горения топлива

```python
from combustion.combust import Combust
from scipy import constants as cst
from pyfluids import FluidsList

fuel = {
        'CH4': (FluidsList.Methane, .98, 35800000.),
        'C2H6': (FluidsList.Ethane, .0, 61300000.),
        'C3H8': (FluidsList.Propyne, .009, 91300000.),
        'C4H10': (FluidsList.Butene, .008, 118500000.),
        'C5H12': (FluidsList.CycloPentane, .0, 146500000.),
        'CO2': (FluidsList.CarbonDioxide, .0, .0),
        'N2': (FluidsList.Nitrogen, .003, .0),
        'O2': (FluidsList.Oxygen, .0, .0),
        'H2S': (FluidsList.SulfurDioxide, .0, 23400000.),
        'H2O': (FluidsList.Water, .0, .0),
        'H2': (FluidsList.Hydrogen, .0, 10800000.),
        }

# Расчет калориметрической температуры горения топлива
comb = Combust(fluid=fuel)
temp = comb._burnTemp(k_alpha=1.3, temp_air=300, pressure=cst.atm)
print(f"Калориметрическая температура горения топлива {temp:.0f} \u2103")
print("Состав газа:\n", comb.gas)
```

# Расчет коэффициента избытка воздуха
```python
k_alpha = comb.burnAlpha(temp_gas=950, temp_air=300, pressure=cst.atm)
print(f'Коэффициент избытка воздуха {k_alpha:0.3f}')
print("Состав газа:\n", comb.gas)
```

# Построение графика зависимости калориметрической температуры от коэффициента избытка воздуха
```python
comb.tempAlphaPl(k_alpha=(2, 5))
```

## About the author
Sergey Besedin,
dr. of sc., prof.

Konstantin Parfyonov, ingener
