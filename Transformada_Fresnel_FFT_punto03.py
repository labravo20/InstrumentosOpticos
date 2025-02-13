#Impresión para verificar correcto funcionamiento entorno de programación
print("Inicializando entorno de programación Transformada de Fresnel (Punto03, Práctica 01)...")

import numpy as np
import matplotlib.pyplot as plt
import LIBRERIA_Mascaras_Transmitancia as mask
import LIBRERIA_Funciones_importantes as function


""" Definición de variables input para creación malla de puntos máscara difractiva"""
resolucion_input = 3000
longitud_arreglo_input = 0.005 #UNIDADES: m


""" Creación de malla de puntos máscara difractiva"""
xx,yy = mask.malla_Puntos(resolucion_input,longitud_arreglo_input)

""" Creación de variables para definición de la transmitancia de amplitudd periódica """
m = 1 #Se asigna este valor para garantizar que la transmitancia va a estar en el intervalo [0,1]

L = 0.00005 #UNIDADES: m --> Se asigna un valor que pueda ser comparable con la longitud de 
           # la ventana del arreglo


""" Parámetros de configuración del arreglo difractivo """
longitud_onda_input = 533E-9 #UNIDADES: m

z_input = 0.05 #UNIDADES: m

N = 4

z_inputTransmitancia = (N*(L**2))/longitud_onda_input

numero_im = 1j
numero_onda = ((2*np.pi)/longitud_onda_input)

""" Definición de función de transmitancia de amplitud periódica """
transmitancia = (1/2)*(1 + (m*np.cos((2*np.pi*xx)/L)))

""" Definición de variables input para creación malla de puntos plano medición"""

#Llamando a la función para determinar los deltas de input y output
delta_muestreo = function.producto_espacio_frecuencia_TransformadaFresnel(longitud_onda_input,z_inputTransmitancia,resolucion_input,longitud_arreglo_input)

resolucion_medicion = resolucion_input 
longitud_arreglo_medicion = (delta_muestreo[1])*resolucion_medicion #UNIDADES: m

""" Creación de malla de puntos plano de medición"""
xx_medicion,yy_medicion = mask.malla_Puntos(resolucion_medicion,longitud_arreglo_medicion)

""" Definición de ecuación de Transformada de Fresnel"""

#Definiendo término de fase constante:
fase_constante = (np.exp(numero_im*numero_onda*z_inputTransmitancia))/(numero_im*longitud_onda_input*z_inputTransmitancia)

#Definiendo término de fase cuadrática en coordenadas de plano medicion:
fase_cuadrada_salida = np.exp((numero_im*numero_onda*(((xx_medicion)**2)+((yy_medicion)**2)))/(2*z_inputTransmitancia))

#Definiendo término de transformada 

## Definiendo término correspondiente a "(campo entrada)x(factor de fase cuadrada) = campo mascara"

#### NOTA: (campo entrada) es resultado de la interacción del campo incidente con la transmitancia
#### de la máscara difractiva
campo_incidente = 1 #Se considera por facilidad esta característica en campo incidente 

campo_entrada = campo_incidente*transmitancia

####Definiendo término de fase cuadrática en coordenadas de plano máscara difractiva:
fase_cuadrada_entrada = np.exp((numero_im*numero_onda*(((xx)**2)+((yy)**2)))/(2*z_inputTransmitancia))

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

#GRAFICANDO MASCARA DIFRACTIVA DE TRANSMITANCIA DE AMPLITUD PERIÓDICA

plt.imshow(transmitancia, extent=[-longitud_arreglo_input/2, longitud_arreglo_input/2, -longitud_arreglo_input/2, longitud_arreglo_input/2], cmap='gray')
plt.title("Máscara difractiva")
plt.colorbar(label="Transmitancia")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()

#PATRÓN DE DIFRACCIÓN USANDO TRANSFORMADA DE FRESNEL

plt.imshow(irradiancia_campo_difractado, extent=[-longitud_arreglo_medicion/2, longitud_arreglo_medicion/2, -longitud_arreglo_medicion/2, longitud_arreglo_medicion/2], cmap='gray')  
plt.title("Patrón de Difracción")
plt.colorbar(label="Amplitud²")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()