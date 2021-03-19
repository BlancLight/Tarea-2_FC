# -*- coding: utf-8 -*-

#Estas son las librerias que se van a utilizar.
import numpy as np
import pandas as pds
import scipy.integrate as spint

M=28.9647/1000 #Masa Molar en kg/mol
R=8.314462 #Constante del gas ideal en J/(mol*K)
g=9.8 #Aceleracion de la gravedad en m*s^-2
yi=0 #Valor inicial de altura del mar en m
yf=3000 #Valor final de altura
n=31 #Cantidad de divisiones, se utiliza 31 para que divida el intervalo de 100 en 100
y=np.linspace(yi, yf, n) #Esta funcion nos genera arreglo con los valores espaciados de 100 en 100
P0=101325 #Valor inicial de la presion


def DifPresion(y,Presion_y):
    '''Esta funcion define la variacion de la temperatura con la altura T(y)

    Parámetros de la función
    ------------------------
    Presion_y : Variable dependiente que indica la presion con respecto a la altura.
    y : Variable independiente que indica la altura sobre el nivel del mar.

    Salida de la función
    --------------------
    DiferencialPresion : Funcion que es la ecuacion diferencial a resolver,
    esta contiene la variacion de la presion con la densidad,con la presion,
    la temperatura y la altura sobre el nivel del mar.
    '''
    DiferencialPresion=-(M*Presion_y*g)/(R*(293-(y/200)))
    return DiferencialPresion

#Se calculan los valores aproximados utilizando RK45 de Scipy
y_aprox_RK45 = spint.solve_ivp(DifPresion, [yi, yf], [P0],t_eval=y,method='RK45')


# Se genera la estructura de datos para presentar los resultados
conjuntodatos = pds.DataFrame({'Posicion':y,\
                'Aprox. RK45': y_aprox_RK45.y[0]},\
                columns=['Posicion','Aprox. RK45'])
print(" ")
print(conjuntodatos)