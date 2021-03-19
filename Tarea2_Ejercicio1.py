# coding=utf-8 
#Estas son las librerias que se van a utilizar.
import numpy as np
import pandas as pds

M = 28.9647/1000 #Masa Molar en kg/mol
R = 8.314462 #Constante del gas ideal en J/(mol*K)
g = 9.8 #Aceleracion de la gravedad en m*s^-2
yi = 0.0 #Valor inicial de altura del mar en m
yf = 3000.0 #Valor final de altura
h = 100 #Espaciado con el que se calcula cada aproximacion en m


y = np.arange(yi,yf+1,h) #Esta funcion permite crear un arreglo espaciado h entre cada valor.
P0 = 101325.0 #Valor inicial de presion

def Funcion_Temperatura(y):
    '''Esta funcion define la variacion de la temperatura con la altura T(y)

    Parámetros de la función
    ------------------------
    y : Variable independiente que indica la altura con respecto al nivel del mar.

    Salida de la función
    --------------------
    T: Funcion que contiene la variacion de la temperatura con la altitud.
    '''
    T = 293.0 - y/200.0
    return T

def Funcion_Densidad(Presion_y,y):
    '''Esta funcion define la variacion de la temperatura con la altura T(y)

    Parámetros de la función
    ------------------------
    Presion_y : Variable dependiente que indica la presion con respecto a la altura.
    y : Variable independiente que indica la altura sobre el nivel del mar.

    Salida de la función
    --------------------
    densidad : Funcion que contiene la variacion de la densidad con la presion,
    la temperatura y la altura sobre el nivel del mar.
    '''
    densidad = (M/(R*Funcion_Temperatura(y)))*Presion_y
    return densidad

def Ecuacion_Diferencial(Presion_y,y):
    '''Esta funcion define la variacion de la temperatura con la altura T(y)

    Parámetros de la función
    ------------------------
    Presion_y : Variable dependiente que indica la presion con respecto a la altura.
    y : Variable independiente que indica la altura sobre el nivel del mar.

    Salida de la función
    --------------------
    ec_diferencial : Funcion que es la ecuacion diferencial a resolver,
    esta contiene la variacion de la presion con la densidad,con la presion,
    la temperatura y la altura sobre el nivel del mar.
    '''
    ec_diferencial = -Funcion_Densidad(Presion_y,y)*g
    return ec_diferencial
    

def RK4(f, P0, y):
    '''Esta función realiza el cálculo de la solución para p(t) de la EDO de
    primer orden p'(y) = f(p(y),y)

    Parámetros de la función
    ------------------------
    f: función que corresponde al lado derecho de la EDO
    P0: valor de P(y) en P = P0 , lo que en conjunto representa la condición inicial
    y : el univeriso de valores de y en que se resolverá la EDO respectiva.

    Salida de la función
    --------------------
    P : arreglo que contiene los valores de la solución de P(y) en los puntos de interes.
    '''

    # Se calcula el valor del paso h
    h = (y[len(y)-1]-y[0])/len(y)

    # Se inicializa el arreglo con los valores aproximados utilizando la
    # cantidad de valores de y de interes
    P = np.zeros(len(y))

    # Se establece la condición inicial
    P[0] = P0

    # Se define el método de cálculo RK4
    for i in range(0, len(y)-1):
        k1 = h*f(P[i], y[i])
        k2 = h*f(P[i]+k1/2, y[i]+h/2)
        k3 = h*f(P[i]+k2/2, y[i]+h/2)
        k4 = h*f(P[i]+k3, y[i]+h)
        P[i+1] = P[i]+(k1+2*k2+2*k3+k4)/6
    return P


# Se calculan los valores aproximados de P(y) por RK4
p_aprox_RK4 = RK4(Ecuacion_Diferencial,P0,y)


# Se genera la estructura de datos para presentar los resultados
conjuntodatos = pds.DataFrame({'posicion': y, 'Aprox. RK4': p_aprox_RK4},\
                columns=['posicion', 'Aprox. RK4'])
print(" ")
print(conjuntodatos)
