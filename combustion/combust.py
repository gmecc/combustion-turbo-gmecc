# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 10:38:03 2024
# Расчет состава газа при сгорании природного газа
@author: User
"""

import pandas as pd
import numpy as np
import sys
from pyfluids import Fluid, FluidsList
from scipy import constants as cst
from scipy.optimize import root
import matplotlib.pyplot as plt
import combustion.spHvol

pd.set_option('display.float_format', '{:.4f}'.format)
pd.set_option('display.max_columns', None)

class Combust:
    def __init__(self, fluid):
        self.g = pd.Series([np.nan], index=['gas_vol'])
        self.f = pd.Series([np.nan], index=['mol_mass_fuel'])
        
        self.fuel = pd.DataFrame.from_dict(fluid, orient='index', 
                    columns=['fluid', 'mol', 'calorific'])
        self.fuel['mol_mass'] = np.nan
        self.fuel['mass'] = np.nan
        self.fuel['partial'] = np.nan
        self.fuel['density'] = np.nan
        self.fuel['enthalpy'] = np.nan
        
        if (self.fuel.mol[:].sum() != 1): sys.exit('НЕКОРРЕКТНЫЕ ИСХОДНЫЕ ДАННЫЕ!!!')
        
        for i in range(len(self.fuel)):
            self.fuel.mol_mass.iat[i] = Fluid(self.fuel.fluid.iat[i]).molar_mass
        
        self.f.mol_mass_fuel = (self.fuel.mol * self.fuel.mol_mass).sum()
        self.fuel.mass = self.fuel.mol * self.fuel.mol_mass / self.f.mol_mass_fuel # массовые доли
        
        self.gas = pd.DataFrame(columns=['fluid', 'mol_mass', 'mol', 'mass'],
                              index=['CO2', 'SO2', 'H2O', 'N2', 'O2'], dtype=float)
        
        self.gas['fluid'] = [FluidsList.CarbonDioxide,
                             FluidsList.SulfurDioxide,
                             FluidsList.Water,
                             FluidsList.Nitrogen,
                             FluidsList.Oxygen]
        
        for i in range(len(self.gas)):
            self.gas.mol_mass.iat[i] = Fluid(self.gas.fluid.iat[i]).molar_mass
        
        
        self.fluid_air = Fluid(FluidsList.Air)
        
        self.sph = pd.DataFrame(columns=['CO2', 'SO2', 'H2O', 'N2', 'O2'],
                index=range(5), dtype=float)
        
        self.sph.CO2 = spHvol.spHvol_CO2_coef
        self.sph.SO2 = spHvol.spHvol_SO2_coef
        self.sph.H2O = spHvol.spHvol_H2O_coef
        self.sph.N2 = spHvol.spHvol_N2_coef
        self.sph.O2 = spHvol.spHvol_O2_coef
        
        self.spHv_air_pol = np.polynomial.Polynomial(spHvol.spHvol_air_coef)
        
        
        # РАСЧЕТ СОСТАВА ПРОДУКТОВ СГОРАНИЯ
    def _massFraction(self, k_alpha):         
        self.k_alpha = k_alpha
        
        # Количество кислорода, необходимое для сжигания топлива м3/м3
        self.volume_O2 = (.5 * self.fuel.mol['H2'] + (2 * self.fuel.mol['CH4'] + 3.5 * 
                    self.fuel.mol['C2H6'] + 5 * self.fuel.mol['C3H8'] + 6.5 * self.fuel.mol['C4H10'] + 
                    8 * self.fuel.mol['C5H12']) + 1.5 * self.fuel.mol['H2S'] - self.fuel.mol['O2'])
        
        # Количество воздуха
        self.volume_air = (1 + 3.76) * self.volume_O2 * self.k_alpha
        
        # Объем продуктов сгорания при k_alpha=1
        self.volume_CO2 = (self.fuel.mol['CH4'] + 2 * self.fuel.mol['C2H6'] + 
                    3 * self.fuel.mol['C3H8'] + 4 * self.fuel.mol['C4H10'] + 
                    5 * self.fuel.mol['C5H12'] + self.fuel.mol['CO2'])
        self.volume_SO2 = self.fuel.mol['H2S']
        self.volume_H2O = (self.fuel.mol['H2O'] + self.fuel.mol['H2'] + self.fuel.mol['H2S'] + 
                    (2 * self.fuel.mol['CH4'] + 3 * self.fuel.mol['C2H6'] + 4 * self.fuel.mol['C3H8'] + 
                      5 * self.fuel.mol['C4H10'] + 6 * self.fuel.mol['C5H12']))
        self.volume_N2 = self.fuel.mol['N2'] + 3.76 * self.volume_O2 
        
        # Объем продуктов сгорания при k_alpha>1
        self.volume_N2_alpha = self.fuel.mol['N2'] + self.k_alpha * 3.76 * self.volume_O2
        self.volume_O2_alpha = (self.k_alpha - 1) * self.volume_O2
        
        # Объем продуктов сгорания
        self.g['gas_vol'] = (self.volume_CO2 + self.volume_SO2 + 
                    self.volume_H2O + self.volume_N2_alpha + 
                    self.volume_O2_alpha)
        # Состав продуктов сгорания - мольные доли
        self.gas.mol.at['CO2'] = self.volume_CO2 / self.g.gas_vol
        self.gas.mol.at['SO2'] = self.volume_SO2 / self.g.gas_vol
        self.gas.mol.at['H2O'] = self.volume_H2O / self.g.gas_vol
        self.gas.mol.at['N2'] = self.volume_N2_alpha / self.g.gas_vol
        self.gas.mol.at['O2'] = self.volume_O2_alpha / self.g.gas_vol
        self.mol_mass_gas = (self.gas.mol_mass * self.gas.mol).sum()
        
        # Массовые доли продуктов сгорания
        self.gas.mass = (self.gas.mol * self.gas.mol_mass / self.mol_mass_gas).round(4)
        self.gas.mass.at['N2'] = self.gas.mass.at['N2'] - self.gas.mass.sum() + 1
        self.g['gas_mass_sum'] = self.gas.mass.sum()
        self.f['calorific_vol'] = (self.fuel['mol'] * self.fuel['calorific']).sum()   # объемная теплота сгорания топлива
        self.f['calorific_mass'] = (self.fuel['mass'] * self.fuel['calorific']).sum() # массовая теплота сгорания топлива
        # self.gas.partial = self.gas.mol * self.press # расчет парциального давления газа
        
        # Удаление нулевых сторок и столбцов
        # nul_ind = self.gas[self.gas.mol == 0].index
        # self.gas = self.gas[self.gas.mol > 0] # удаление нулевых строк
        # self.fuel = self.fuel[self.fuel.mol > 0] # удаление нулевых строк
        # self.sph = self.sph.drop(columns=nul_ind) # удаление лишних столбцов 


        # ФУНКЦИЯ РАСЧЕТА ТЕМПЕРАТУРЫ ГОРЕНИЯ 
    def _tempFunc(self, temp):
        sph = np.zeros(len(self.gas))
        
        for i in range(len(self.gas)):
            spHeat_pol = np.polynomial.Polynomial(self.sph.iloc[:,i])
            sph[i] = spHeat_pol(temp).item()
        
        temp_resedual = abs((self.f.calorific_vol + self.spHv_air_pol(self.temp_air) * self.temp_air )/ 
                ((sph * self.gas.mol).sum() * self.g.gas_vol) - temp)
        return temp_resedual

        #РАСЧЕТ ТЕМПЕРАТУРЫ ГОРЕНИЯ
    def _burnTemp(self, k_alpha, temp_air=300, pressure=cst.atm): # расчет температуры (жаропроизводительность)
        self.temp_air = temp_air
        self.press = pressure
        self._massFraction(k_alpha=k_alpha) # расчет фракционного состава продуктов сгорания
        
        sol = root(self._tempFunc, 1000)
        # self.g['temp_gas'] = sol.x[0]
        return sol.x[0]

        # РАСЧЕТ К-АЛЬФА ПО ТЕПМЕРАТУРЕ ГОРЕНИЯ
    def burnAlpha(self, temp_gas, temp_air, pressure): # расчет k_alpha
        self.temperature = temp_gas
        self.pressure = pressure
    
        resedual = lambda k_alpha: abs(self._burnTemp(k_alpha=k_alpha, 
                temp_air=temp_air, pressure=pressure) - temp_gas)
    
        sol = root(resedual, 1)
        self.g['k_alpha'] = sol.x[0]
        return self.g.k_alpha
    
        # РАСЧЕТ ЗАВИСИМОСТИ И ПОСТОРЕНИЕ ГРАФИКА
    def tempAlphaPl(self, k_alpha, temp_air=300, pressure=cst.atm):
        n = 50
        tk = np.vectorize(self._burnTemp)
        ki = np.linspace(k_alpha[0], k_alpha[1], n)
        plt.plot(ki, tk(ki), color='red')
        plt.minorticks_on()
        plt.xlabel('k_alpha') 
        plt.ylabel('Temperature [\u2103]')
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        # plt.title('title', loc='left')
        plt.tight_layout() # оптимизируем поля и расположение объектов
        # plt.savefig('temp-k_alpha', dpi = 300)
        plt.show()



    
    
    
    













