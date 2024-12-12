#Impresión para verificar correcto funcionamiento entorno de programación
print("Inicializando entorno de programación Transformada de Fresnel...")

import numpy as np
import matplotlib.pyplot as plt

''' Definicion de las difefrentes funciones de transmitancia'''

def malla_Puntos_mascara_difractiva(resolucion, longitud_Arreglo):
    ''' CREACION DE LAS MALLAS DE PUNTOS Y LOS DELTAS PARA LOS PRODUCTOS ESPACIO FRECUENCIA '''
    x = np.linspace(-longitud_Arreglo / 2, longitud_Arreglo / 2, resolucion) #crea las mallas de puntos para el arreglo 
    y = np.linspace(-longitud_Arreglo / 2, longitud_Arreglo / 2, resolucion) 
    xx, yy = np.meshgrid(x, y) #crea una malla de puntos bidimensional 
    return xx, yy #retornamos la malla de puntos

def malla_Puntos_plano_medicion(resolucion, longitud_Arreglo):
    ''' CREACION DE LAS MALLAS DE PUNTOS Y LOS DELTAS PARA LOS PRODUCTOS ESPACIO FRECUENCIA '''
    x = np.linspace(-longitud_Arreglo / 2, longitud_Arreglo / 2, resolucion) #crea las mallas de puntos para el arreglo 
    y = np.linspace(-longitud_Arreglo / 2, longitud_Arreglo / 2, resolucion) 
    xx, yy = np.meshgrid(x, y) #crea una malla de puntos bidimensional 
    return xx, yy #retornamos la malla de puntos

#Se define una función con transmitancia asociada a una abertura circular
def funcion_Circulo(radio, centro, xx, yy): #definicion de la funcion para hacer circulo transparente
    ''' CREACCION DEL CONJUNTO DE PUNTOS DE LA MASCARA DE DFRACCION '''
    if centro is None: # que pasa si el centro no es definido en la funcion
        centro = (0, 0) #ubica el centro de la circunferencia en el origen
    distancia = (xx - centro[0])**2 + (yy - centro[1])**2 #calculamos la distancia desde el centro de la circunferencia a cada punto
    mascara = distancia <= radio**2 #los puntos de la mascara seran los puntos cuya distancia al centro es menor que el radio
    return mascara #devolvemos los puntos que cumplen la condicion para hacer parte de la mascara

def producto_espacio_frecuencia(longitud_onda,z,resolucion,longitud_Arreglo):

    delta_input = longitud_Arreglo/resolucion #Definiendo separación entre número total de muestras

    delta_output = (longitud_onda*z)/(resolucion*delta_input) #Definiendo separación entre muestras
                                                            #en plano de medición

    delta_muestreo = [delta_input,delta_output]

    return delta_muestreo



""" Definición de variables input para creación malla de puntos máscara difractiva"""
resolucion_input = 300
longitud_arreglo_input = 0.05 #UNIDADES: m

""" Definición de variables input para creación abertura circular"""
radio_input = 0.01 #UNIDADES: m
centro_input = None 

""" Creación de malla de puntos máscara difractiva"""
xx,yy = malla_Puntos_mascara_difractiva(resolucion_input,longitud_arreglo_input)

""" Creación de abertura circular """
mascara_circular = funcion_Circulo(radio_input,centro_input,xx,yy)

""" Parámetros de configuración del arreglo difractivo """
longitud_onda_input = 533E-9 #UNIDADES: m
z_input = 10 #UNIDADES: m
numero_im = 1j
numero_onda = ((2*np.pi)/longitud_onda_input)

""" Definición de variables input para creación malla de puntos plano medición"""

#Llamando a la función para determinar los deltas de input y output
delta_muestreo = producto_espacio_frecuencia(longitud_onda_input,z_input,resolucion_input,longitud_arreglo_input)

resolucion_medicion = resolucion_input 
longitud_arreglo_medicion = (delta_muestreo[1])*resolucion_medicion #UNIDADES: m

""" Creación de malla de puntos plano de medición"""
xx_medicion,yy_medicion = malla_Puntos_plano_medicion(resolucion_medicion,longitud_arreglo_input)

""" Definición de ecuación de Transformada de Fresnel"""

#Definiendo término de fase constante:
fase_constante = (np.exp(numero_im*numero_onda*z_input))/(numero_im*longitud_onda_input*z_input)

#Definiendo término de fase cuadrática en coordenadas de plano medicion:
fase_cuadrada_salida = np.exp((numero_im*numero_onda*(((xx_medicion)**2)+((yy_medicion)**2)))/(2*z_input))

#Definiendo término de transformada 

## Definiendo término correspondiente a "(campo entrada)x(factor de fase cuadrada) = campo mascara"

#### NOTA: (campo entrada) es resultado de la interacción del campo incidente con la transmitancia
#### de la máscara difractiva
campo_incidente = 1 #Se considera por facilidad esta característica en campo incidente 

campo_entrada = campo_incidente*mascara_circular

####Definiendo término de fase cuadrática en coordenadas de plano máscara difractiva:
fase_cuadrada_entrada = np.exp((numero_im*numero_onda*(((xx)**2)+((yy)**2)))/(2*z_input))

####Definiendo campo mascara
campo_mascara = campo_entrada*fase_cuadrada_entrada

## Definiendo termino correspondiente a transformada de Fourier del campo que pasa a través de la máscara difractiva

transformada_campo = np.fft.fft2(campo_mascara)
transformada_campo_centrada = np.fft.fftshift(transformada_campo)

## ---- ##
termino_transformada_campo = ((delta_muestreo[1])**2)*transformada_campo_centrada

#Definiendo resultado de campo difractado en plano de medición
campo_difractado = fase_constante*fase_cuadrada_salida*termino_transformada_campo

""" Calculando irradiancia del campo difractado """
irradiancia_campo_difractado = ((np.abs(campo_difractado))**2)

""" Graficando el resultado de difracción por transformada de Fresnel """

#GRAFICANDO MASCARA DIFRACTIVA DE ABERTURA CIRCULAR
plt.imshow(mascara_circular, extent=[-longitud_arreglo_input/2, longitud_arreglo_input/2, -longitud_arreglo_input/2, longitud_arreglo_input/2], cmap='gray')
plt.title("Máscara Circular")
plt.colorbar(label="Transmitancia")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()

#Patrón de difracción de Fresnel
plt.imshow(irradiancia_campo_difractado, extent=[-longitud_arreglo_medicion/2, longitud_arreglo_medicion/2, -longitud_arreglo_medicion/2, longitud_arreglo_medicion/2], cmap='gray', vmax= 1*(np.max(irradiancia_campo_difractado)))  # log1p para mejor visualización
plt.title("Patrón de Difracción")
plt.colorbar(label="Amplitud²")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()