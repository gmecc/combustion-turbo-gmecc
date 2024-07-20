# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 07:35:05 2024
Программа построения графиков средней объемной теплоемкости газов
@author: User
"""
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
# import spHvol

spHvol_air_coef = [ 1.29916543e+03, 2.81610452e-02, 1.82765897e-04, -1.24381440e-07, 2.56717562e-11]
spHvol_CO2_coef = [ 1.63479959e+03,  9.75263813e-01, -5.45793612e-04,  1.83324681e-07, -2.67917924e-11]
spHvol_O2_coef = [ 1.30543577e+03,  1.74761949e-01,  6.42329236e-05, -8.82987243e-08,  2.28704857e-11]
spHvol_N2_coef = [ 1.30179986e+03, -1.08393937e-02,  2.16685585e-04, -1.35474341e-07,  2.67040689e-11]
spHvol_H2O_coef = [ 1.49735370e+03,  1.18320983e-01,  1.79154428e-04, -8.92543875e-08,  1.34614261e-11]
spHvol_SO2_coef = [ 1.87971903e+03,  5.92328524e-01, -9.65428295e-05, -9.09724113e-08, 3.29674615e-11]


pd.set_option('display.float_format', '{:.3e}'.format)

class spHvol_plot:
    def __init__(self):
        self.sph = pd.DataFrame(columns=['CO2', 'SO2', 'H2O', 'N2', 'O2'],
                index=range(5), dtype=float)
        
        self.sph.CO2 = spHvol_CO2_coef
        self.sph.SO2 = spHvol_SO2_coef
        self.sph.H2O = spHvol_H2O_coef
        self.sph.N2 = spHvol_N2_coef
        self.sph.O2 = spHvol_O2_coef
        
    def spHplot(self, fluid):
        temp = np.linspace(0, 2000)
        spHeat_pol = np.polynomial.Polynomial(self.sph[fluid])
        plt.plot(temp, spHeat_pol(temp), label=fluid)
        
    def Splot(self):
        
        for fluid in self.sph.columns:
            self.spHplot(fluid=fluid)
        
        plt.xlabel('$T, \u2103$')
        plt.ylabel('$C_{p.mean}, Дж/(м3.К)$')
        plt.minorticks_on()
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.title('Средняя удельная объемная теплоемкость', loc='left')
        plt.legend()
        plt.tight_layout() # оптимизируем поля и расположение объектов
        plt.show()
        
# sph = spHvol_plot()
# sph.Splot()
# print(sph.sph)

