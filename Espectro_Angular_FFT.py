#Impresión para verificar correcto funcionamiento entorno de programación
print("Inicializando entorno de programación Espectro Angular...")

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

def funcion_Circulo(radio, centro, xx, yy): #definicion de la funcion para hacer circulo transparente
    ''' CREACCION DEL CONJUNTO DE PUNTOS DE LA MASCARA DE DFRACCION '''
    if centro is None: # que pasa si el centro no es definido en la funcion
        centro = (0, 0) #ubica el centro de la circunferencia en el origen
    distancia = (xx - centro[0])**2 + (yy - centro[1])**2 #calculamos la distancia desde el centro de la circunferencia a cada punto
    mascara = distancia <= radio**2 #los puntos de la mascara seran los puntos cuya distancia al centro es menor que el radio
    return mascara #devolvemos los puntos que cumplen la condicion para hacer parte de la mascara

def producto_espacio_frecuencia(longitud_onda,z,resolucion,longitud_Arreglo):

    delta_input = longitud_Arreglo/resolucion #Definiendo separación entre número total de muestras

    delta_output = (1)/(resolucion*delta_input) #Definiendo separación entre muestras
                                                            #en plano de medición

    delta_muestreo = [delta_input,delta_output]

    return delta_muestreo

""" Definición de variables input para creación malla de puntos máscara difractiva"""
resolucion_input = 1000
longitud_arreglo_input = 0.0025 #UNIDADES: m

""" Definición de variables input para creación abertura circular"""
radio_input = 0.0005 #UNIDADES: m
centro_input = None 

""" Creación de malla de puntos máscara difractiva"""
xx,yy = malla_Puntos_mascara_difractiva(resolucion_input,longitud_arreglo_input)

""" Creación de abertura circular """
mascara_circular = funcion_Circulo(radio_input,centro_input,xx,yy)

""" Parámetros de configuración del arreglo difractivo """
longitud_onda_input = 633E-9 #UNIDADES: m
z_input = 0.015 #UNIDADES: m
numero_im = 1j
numero_onda = ((2*np.pi)/longitud_onda_input)

""" Definición de variables input para creación malla de puntos plano medición"""

#Llamando a la función para determinar los deltas de input y output
delta_muestreo = producto_espacio_frecuencia(longitud_onda_input,z_input,resolucion_input,longitud_arreglo_input)

resolucion_medicion = resolucion_input 
longitud_arreglo_medicion = (delta_muestreo[1])*resolucion_medicion #UNIDADES: m

""" Creación de malla de puntos plano de medición"""
xx_espectro,yy_espectro = malla_Puntos_plano_medicion(resolucion_medicion,longitud_arreglo_medicion)

""" Definición de ecuación Espectro Angular"""

#Definiendo término de espectro angular en máscara difractiva

## Definiendo término correspondiente a "campo entrada"

#### NOTA: (campo entrada) es resultado de la interacción del campo incidente con la transmitancia
#### de la máscara difractiva
campo_incidente = 1 #Se considera por facilidad esta característica en campo incidente 

campo_entrada = campo_incidente*mascara_circular

## Definiendo termino correspondiente a transformada de Fourier del campo de entrada
transformada_campo = np.fft.fft2(campo_entrada)
transformada_campo_centrada = np.fft.fftshift(transformada_campo)

## ---- ##
espectro_angular_plano_mascara = ((delta_muestreo[0])**2)*transformada_campo_centrada

#Definiendo término de espectro angular en plano de medición

##Definiendo término de condición de evanescencia

### Definiendo término exponente que determina condición de evanescencia
cond_evanescencia = (np.sqrt(1-((longitud_onda_input**2)*(((xx_espectro)**2)+((yy_espectro)**2)))))

termino_cond_evanescencia = np.exp(numero_im*z_input*numero_onda*cond_evanescencia)

## ---- ##
espectro_angular_plano_medicion = termino_cond_evanescencia*espectro_angular_plano_mascara

#Definiendo término de transformada de Fourier inversa del espectro angular en el plano de medición
transformada_inversa_espectro_angular = np.fft.ifft2(espectro_angular_plano_medicion)

#Definiendo resultado de campo difractado en plano de medición
campo_difractado = ((delta_muestreo[1])**2)*transformada_inversa_espectro_angular

""" Calculando irradiancia del campo difractado """
irradiancia_campo_difractado = ((np.abs(campo_difractado))**2)

""" Graficando el resultado de difracción por espectro angular """

#GRAFICANDO MASCARA DIFRACTIVA DE ABERTURA CIRCULAR
plt.imshow(mascara_circular, extent=[-longitud_arreglo_input/2, longitud_arreglo_input/2, -longitud_arreglo_input/2, longitud_arreglo_input/2], cmap='gray')
plt.title("Máscara Circular")
plt.colorbar(label="Transmitancia")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()

#Patrón de difracción de Fresnel
plt.imshow(irradiancia_campo_difractado, extent=[-longitud_arreglo_input/2, longitud_arreglo_input/2, -longitud_arreglo_input/2, longitud_arreglo_input/2], cmap='gray', vmax= 1*(np.max(irradiancia_campo_difractado)))  # log1p para mejor visualización
plt.title("Patrón de Difracción")
plt.colorbar(label="Amplitud²")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()
