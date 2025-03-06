""" DESCRIPCIÓN 

Documento para trabajar la reconstrucción de hologramas. El proceso consiste en aplicar la transformada de Fourier 
al campo del holograma para filtrar, en el dominio de Fourier, las frecuencias asociadas a una de las imágenes GEMELAS
con información sobre el holograma. Posteriormente, el campo resultante de la filtración es multiplicado por el inverso de 
la onda plana de referencia usada para realizar la interferencia, pues de esta forma de descarta el aporte de la interferencia
antes mencionada. El resultado de las anteriores operaciones es el objeto original magnificado por el sistema de amplificación
del microscopio.

"""


""" Importando librerias y documentos necesarios """

import LIBRERIA_Mascaras_Transmitancia as mascaras
import LIBRERIA_Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import matplotlib.pyplot as plt
import LIBRERIA_Funciones_importantes as function
import LIBRERIA_Funciones_Graficacion as graph



' ################ INICIO SECCIÓN DE CARACTERIZACIÓN ARREGLO ################## '
##### LA SIGUIENTE SECCIÓN SE REFIERE A INPUTS QUE DEBEN SER ESPECIFICADOS POR EL USUARIO*

""" Definiendo parámetros de sensor """

#Numero de muestras/puntos distribuidos en el ancho del sensor
#resolucion_anchoSensorInput = 1280 
resolucion_anchoSensorInput = 1024  # Usando esta configuración se simplifica el cálculo para un sensor con distribución CUADRADA
                                    # NOTA: Se decide usar esta configuración porque ya se configuraron las características de 
                                    # reconstrucción adecuadas asociadas a este caso (específicamente las coordenadas del orden
                                    # +1 presente en la transformada de Fourier del holograma). 

#Numero de muestras/puntos distribuidos en el alto del sensor
resolucion_altoSensorInput = 1024

#Tamaño del pixel del sensor #UNIDADES: m
tamaño_PixelSensorInput = 5.2E-6



""" Definiendo parámetros de la iluminación coherente a utilizar  """

#Definiendo la longitud de onda asociada a la fuente en consideración
longitud_onda_input = 632.8E-9 #UNIDADES: m



""" Definiendo parámetros del OM usado para el sistema de magnificación """

#Magnificación objetivo de microscopio
magnificacion = 10

#Apertura numérica objetivo de microscopio
apertura_Numerica = 0.25 

#Distancia focal del objetivo de microscopio #UNIDADES: m
### NOTA: Si se desconoce la distancia focal del objetivo de microscopio escribir el número CERO (0) 
distancia_focalMO = 0

#Distancia focal asociada a la lente de tubo --> #UNIDADES: m
distancia_focalTL = 0.2   


' ################ FIN SECCIÓN DE CARACTERIZACIÓN ARREGLO ################## '


' ################ EMPIEZA SECCIÓN DE CÁLCULOS PARA CARACTERIZACIÓN ARREGLO ################## '

""" Calculando parámetros de muestreo del sensor utilizado """

#Se crea una lista a la cual se le asigna los valores de los delta de muestreo asociados al sensor
deltas_Sensor = {"deltaPlanoEntrada_X":tamaño_PixelSensorInput,"deltaPlanoEntrada_Y":tamaño_PixelSensorInput}

#Usando los datos anteriores se calcula el ancho y alto físicos del sensor
ancho_SensorInput = resolucion_anchoSensorInput*tamaño_PixelSensorInput
alto_SensorInput = resolucion_altoSensorInput*tamaño_PixelSensorInput



""" Definiendo parámetros de fuente """

#Calculando el número de onda
numero_onda_input = (2*np.pi)/longitud_onda_input


""" Definición de distancias del arreglo """

#Verficación de distancia focal asociada al objetivo de microscopio 
if (distancia_focalMO == 0):

    #En caso de distancia focal desconocida, esta se calcula haciendo uso de la relación de Magnificación
    distancia_focalMO = distancia_focalTL/magnificacion
else:

    distancia_focalMO = distancia_focalMO



""" Definiendo parámetro para el tamaño de la pupila """

#Calculando el ángulo de apertura
angulo_Apertura = np.arcsin(apertura_Numerica)

#Calculando el radio de la pupila asociada al objetivo de microscopio #UNIDADES: m
radio_pupilaInput =  distancia_focalMO*np.tan(angulo_Apertura) # Se define variable asociada al radio de la abertura que representará el
                                                               # radio del objetivo de microscopio --> Este se calcula haciendo uso de la
                                                               # abertura numérica y la longitud de tubo. 


' ################ FIN SECCIÓN DE CÁLCULOS PARA CARACTERIZACIÓN ARREGLO ################## '


' ################ EMPIEZA SECCIÓN DE RECONSTRUCCIÓN ################## '


""" Creación de mascaras de transmitancia """

# Cargar la imagen PNG como máscara de transmitancia
ruta_imagen_png = "/home/labravo/Desktop/Instrumentos ópticos/PROYECTO/HOLOGRAMAS/fibra_EnfoquePlanoEnfoqueLuz.tif" 
#ruta_imagen_png = "/home/labravo/Downloads/Holograma007.tif"  # Especifica la ruta de tu imagen
#ruta_imagen_png = "/home/labravo/Downloads/Mi_primer_holograma.PNG"  # Especifica la ruta de tu imagen
mascaraReconstruccion = function.cargar_imagen_png(ruta_imagen_png,resolucion_anchoSensorInput,resolucion_altoSensorInput)


# Crear la malla de puntos plano de entrada --> HOLOGRAMA A RECONSTRUIR 
xx_mascaraReconstruccion, yy_mascaraReconstruccion = mascaras.malla_Puntos(resolucion_anchoSensorInput, ancho_SensorInput,resolucion_altoSensorInput,
                                               alto_SensorInput)



""" Implementación de condición producto espacio frecuencia para cálculo de malla de puntos en dominio de Fourier """

# Crear malla de puntos de plano ESPECTRO DE FOURIER

delta_fx = (longitud_onda_input*distancia_focalMO)/(resolucion_anchoSensorInput * deltas_Sensor["deltaPlanoEntrada_X"]) #Tamaño de cada pixel en plano de Fourier (1/m)
delta_fy = (longitud_onda_input*distancia_focalMO)/(resolucion_altoSensorInput * deltas_Sensor["deltaPlanoEntrada_Y"]) #Tamaño de cada pixel en plano de Fourier (1/m)


#Calculando el ancho del arreglo en el plano de Fourier
anchoX_PlanoFourier = delta_fx*resolucion_anchoSensorInput
altoY_PlanoFourier = delta_fy*resolucion_altoSensorInput


# Crear la malla de puntos plano de entrada --> HOLOGRAMA A RECONSTRUIR 
xx_planoFourier, yy_planoFourier = mascaras.malla_Puntos(resolucion_anchoSensorInput, 
                                                         anchoX_PlanoFourier,
                                                         resolucion_altoSensorInput,
                                                         altoY_PlanoFourier)


""" Calculando la transformada de Fourier del holograma """
transformada_FourierReconstruccion = np.fft.fftshift(np.fft.fft2(mascaraReconstruccion))

#Calculando la intensidad de la transformada de Fourier
intensidad_TransformadaFourierRecosntruccion = (np.abs(transformada_FourierReconstruccion))**2



"""Se diseña e implementa una máscara de transmitancia para filtrar la imágen gemela de interés """
### Se diseñan diferentes condiciones de máscara de filtrado dependiendo del ángulo del haz de referencia usado para 
### la formación del holograma

#Definición del radio de la pupila en el plano de Fourier
#radio_pupilaEscaladaPlanoFourier= 0.0003 # --> USAF
radio_pupilaEscaladaPlanoFourier= 0.0004 # --> MUESTRAS CON ILUMINACIÓN 


# Se define el vector asociado a las coordenadas de la máscara de filtrado  --> CALIBRANDO COORDENADAS TRANSFORMADA DE FOURIER
#coordenadas_MascaraFiltrado = [-0.000258,0.000158] #coordenada TOMAS EXPERIMENTALES USAF
#coordenadas_MascaraFiltrado = [-0.00107,-0.0002] #coordenada TOMAS EXPERIMENTALES SENSOR CUADRADO
coordenadas_MascaraFiltrado = [-0.00081,-0.0002] #coordenada TOMAS EXPERIMENTALES SENSOR CUADRADO --> prueba DESPLAZAMIENTO
#coordenadas_MascaraFiltrado = [-0.000852,-0.0002] #coordenada TOMAS EXPERIMENTALES 

#Creación de la máscara de filtrado para el proceso de reconstrucción
mascaraReconstruccion_Filtrado = mascaras.funcion_Circulo(radio_pupilaEscaladaPlanoFourier,coordenadas_MascaraFiltrado,
                                                          xx_planoFourier,yy_planoFourier)

#Se implementa la filtración para obtener el campo filtrado
campo_Filtrado = mascaraReconstruccion_Filtrado*transformada_FourierReconstruccion

#Intensidad del campo resultante después del proceso de filtrado
intensidad_CampoFiltrado = (np.abs(campo_Filtrado))**2

""" Se vuelve al dominio espacial para determinar el campo complejo asociado al campo filtrado """

#Calculando la transformada de fourier inversa del campo shifteado
campo_FiltradoShifteado = np.fft.fftshift(campo_Filtrado)

#Calculando la transformada de Fourier inversa
campo_Reconstruccion = np.fft.ifft2(campo_FiltradoShifteado)

#Calculando intensidad del campo resultante en el dominio espacial
intensidad_CampoOpticoHolograma = (np.abs(campo_Reconstruccion))**2


""" Se diseña e implementa una onda plana, inversa al haz de referencia, para despreciar el aporte de la interferencia """

#Definición de onda plana para eliminar aporte de efectos de interferencia
#vector_cosenosDirectoresRef =  [-0.012656,0.0075936] #coordenada TOMAS EXPERIMENTALES USAF
vector_cosenosDirectoresRef =  [-0.0535,-0.01035] #coordenada TOMAS EXPERIMENTALES SENSOR CUADRADO
#vector_cosenosDirectoresRef =  [-0.0426,-0.01] #coordenada TOMAS EXPERIMENTALES

#Definición de vector de onda asociado al haz de referencia
vector_ondaRef = [numero_onda_input*vector_cosenosDirectoresRef[0], numero_onda_input*vector_cosenosDirectoresRef[1]]

# Definición de onda plana INVERSA al haz de referencia  
onda_PlanaRefInversa = np.exp(-1j*((vector_ondaRef[0]*xx_mascaraReconstruccion)  + (vector_ondaRef[1]*yy_mascaraReconstruccion)))


# Se multiplica el campo de interés (PLANO MEDICIÓN) con la onda plana generada
matriz_campoNOContribucionOndaPlana = campo_Reconstruccion*onda_PlanaRefInversa

#Intensidad del campo resultante de la reconstrucción
intensidad_matrizCampoNOcontribucionOndaPlana = (np.abs(matriz_campoNOContribucionOndaPlana))**2

' ################ FIN SECCIÓN DE RECONSTRUCCIÓN ################## '


' ################ EMPIEZA SECCIÓN DE GRAFICACIÓN ################## '

""" Graficando máscara de transmitancia"""
graph.graficar_transmitancia(mascaraReconstruccion,ancho_SensorInput,alto_SensorInput,"Imágen de entrada para reconstrucción")


#Graficando la transformada de Fourier del holograma
graph.graficar_intensidad(np.log(intensidad_TransformadaFourierRecosntruccion),anchoX_PlanoFourier,altoY_PlanoFourier,
                             "Transformada de Fourier del Holograma",1,1)

graph.graficar_transmitancia(mascaraReconstruccion_Filtrado,anchoX_PlanoFourier,altoY_PlanoFourier,"Mascara Filtrado")


graph.graficar_intensidad(intensidad_CampoFiltrado,anchoX_PlanoFourier,altoY_PlanoFourier,"Tranformada de Fourier filtrada")


graph.graficar_intensidad(intensidad_matrizCampoNOcontribucionOndaPlana,ancho_SensorInput,alto_SensorInput,
                          "Campo óptico del objeto")

graph.graficar_fase(np.angle(matriz_campoNOContribucionOndaPlana),ancho_SensorInput,alto_SensorInput,
                    "Distribución de fase Campo óptico objeto")


#MULTIPLICANDO POR TÉRMINO DE FASE 
# fase_adicionalPrueba = np.exp(1j*(1*(np.pi)/4))

# matriz_prueba = matriz_campoNOContribucionOndaPlana*fase_adicionalPrueba

# graph.graficar_fase(np.angle(matriz_prueba),ancho_SensorInput,alto_SensorInput,
#                     "Distribución de fase Campo óptico objeto + FASE PRUEBA")


' ################ EMPIEZA SECCIÓN DE GRAFICACIÓN ################## '