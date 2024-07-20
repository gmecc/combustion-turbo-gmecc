# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 07:35:05 2024
Программа построения графиков средней объемной теплоемкости газов
@author: User
"""
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import spHvol

pd.set_option('display.float_format', '{:.3e}'.format)

class spHvol_plot:
    def __init__(self):
        self.sph = pd.DataFrame(columns=['CO2', 'SO2', 'H2O', 'N2', 'O2'],
                index=range(5), dtype=float)
        
        self.sph.CO2 = spHvol.spHvol_CO2_coef
        self.sph.SO2 = spHvol.spHvol_SO2_coef
        self.sph.H2O = spHvol.spHvol_H2O_coef
        self.sph.N2 = spHvol.spHvol_N2_coef
        self.sph.O2 = spHvol.spHvol_O2_coef
        
    def spHplot(self, fluid):
        temp = np.linspace(0, 2000)
        spHeat_pol = np.polynomial.Polynomial(self.sph[fluid])
        plt.plot(temp, spHeat_pol(temp), label=fluid)
        
    def Splot(self):
        
        for fluid in self.sph.columns:
            self.spHplot(fluid=fluid)
        
        plt.xlabel('$T, ^\circ C$') 
        plt.ylabel('$C_{p.mean}, Дж/(м3.К)$') 
        plt.minorticks_on()
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.title('Средняя удельная объемная теплоемкость', loc='left')
        plt.legend()
        plt.tight_layout() # оптимизируем поля и расположение объектов
        plt.show()
        
sph = spHvol_plot()
sph.Splot()
print(sph.sph)

