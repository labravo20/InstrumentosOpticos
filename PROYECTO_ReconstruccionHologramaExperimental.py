
""" Importando librerias y documentos necesarios """

import LIBRERIA_Mascaras_Transmitancia as mascaras
import LIBRERIA_Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import matplotlib.pyplot as plt
import LIBRERIA_Funciones_importantes as function
import LIBRERIA_Funciones_Graficacion as graph



######################## INICIA SECCIÓN DE RECONSTRUCCIÓN HOLOGRAMA ##############################

""" Definición de variables para cálculo de sistema 4F para reconstrucción del holograma """

#Definiendo distancia focal asociada a la primera lente
distancia_focalLente01 = 0.06

#Definiendo distancia focal asociada a la segunda lente
distancia_focalLente02 = 0.02

#Definiendo distancia de propagación arbitraria
distancia_propagacionAribitrariaReconstruccion = 0.02

' ################ INICIO SECCIÓN DE CARACTERIZACIÓN ARREGLO ################## '
##### LA SIGUIENTE SECCIÓN SE REFIERE A INPUTS QUE DEBEN SER ESPECIFICADOS POR EL USUARIO*

""" Definiendo parámetros de sensor """

#Numero de muestras/puntos distribuidos en el ancho del sensor
#resolucion_anchoSensorInput = 1280 
resolucion_anchoSensorInput = 1024

#Numero de muestras/puntos distribuidos en el alto del sensor
resolucion_altoSensorInput = 1024

#Tamaño del pixel del sensor #UNIDADES: m
tamaño_PixelSensorInput = 5.2E-6



""" Definiendo parámetros del montaje  """

#Se define una distancia de propagación arbitraria 'd' para simulación de sistema 4F --> #UNIDADES: m
distancia_propagacionAribitraria = 0.2



""" Definiendo parámetros de la iluminación coherente a utilizar  """

#Definiendo la longitud de onda asociada a la fuente en consideración
longitud_onda_input = 632.8E-9 #UNIDADES: m



""" Definiendo inclinación del ángulo de incidencia del haz de referencia para interferencia
en la formación del holograma"""

#Definiendo ángulo máximo --> LAS SIGUIENTES LINEAS DE CÓDIGO NO REQUIEREN INPUT DEL USUARIO
angulo_MaxHazReferencia = np.arcsin(longitud_onda_input/(2*tamaño_PixelSensorInput))
print("\nÁngulo máximo de inclinación del haz de referencia para garantizar buen funcionamiento sistema [GRADOS]:")
print(np.degrees(angulo_MaxHazReferencia))

#Definiendo ángulo de inclinación del haz de referencia respecto al haz objeto
# NOTA: Leer en terminal 'print' con información sobre el ángulo máximo 
cosenoDirector_HazReferenciaAlfa = np.cos(np.radians(0.25))

cosenoDirector_HazReferenciaBeta = np.cos(np.radians(0.25))


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
print(radio_pupilaInput)


""" Creación de mascara de transmitancia """

# Cargar la imagen PNG como máscara de transmitancia
#ruta_imagen_png = "/home/labravo/Desktop/Instrumentos ópticos/PROYECTO/HOLOGRAMAS/fibra_Enfoque_Luz1.tif"
ruta_imagen_png = "/home/labravo/Desktop/Instrumentos ópticos/PROYECTO/HOLOGRAMAS/fibra_EnfoquePlanoEnfoqueLuz.tif" 
#ruta_imagen_png = "/home/labravo/Downloads/Holograma007.tif"  # Especifica la ruta de tu imagen
#ruta_imagen_png = "/home/labravo/Downloads/Mi_primer_holograma.PNG"  # Especifica la ruta de tu imagen
mascaraReconstruccion = function.cargar_imagen_png(ruta_imagen_png,resolucion_anchoSensorInput,resolucion_altoSensorInput)



# Crear la malla de puntos plano de entrada --> HOLOGRAMA A RECONSTRUIR 
xx_mascaraReconstruccion, yy_mascaraReconstruccion = mascaras.malla_Puntos(resolucion_anchoSensorInput, ancho_SensorInput,resolucion_altoSensorInput,
                                               alto_SensorInput)


""" Graficando máscara de transmitancia"""
graph.graficar_transmitancia(mascaraReconstruccion,ancho_SensorInput,alto_SensorInput,"Imágen de entrada para reconstrucción")



# Crear malla de puntos de plano ESPECTRO DE FOURIER

delta_fx = (longitud_onda_input*distancia_focalMO)/(resolucion_anchoSensorInput * deltas_Sensor["deltaPlanoEntrada_X"]) #Tamaño de cada pixel en plano de Fourier (1/m)
delta_fy = (longitud_onda_input*distancia_focalMO)/(resolucion_altoSensorInput * deltas_Sensor["deltaPlanoEntrada_Y"]) #Tamaño de cada pixel en plano de Fourier (1/m)


#Calculando el ancho del arreglo en el plano de Fourier
anchoX_PlanoFourier = delta_fx*resolucion_anchoSensorInput
altoY_PlanoFourier = delta_fy*resolucion_altoSensorInput
print(delta_fy)

# Crear la malla de puntos plano de entrada --> HOLOGRAMA A RECONSTRUIR 
xx_planoFourier, yy_planoFourier = mascaras.malla_Puntos(resolucion_anchoSensorInput, 
                                                         anchoX_PlanoFourier,
                                                         resolucion_altoSensorInput,
                                                         altoY_PlanoFourier)


""" Calculando la transformada de Fourier del holograma """
transformada_FourierReconstruccion = np.fft.fftshift(np.fft.fft2(mascaraReconstruccion))

#Calculando la intensidad de la transformada de Fourier
intensidad_TransformadaFourierRecosntruccion = (np.abs(transformada_FourierReconstruccion))**2




"""Se diseña una máscara de transmitancia para filtrar la imágen gemela de interés """

#Definición del radio de la pupila en el plano de Fourier
#radio_pupilaEscaladaPlanoFourier= 0.0003 # --> USAF
radio_pupilaEscaladaPlanoFourier= 0.0005 # --> MUESTRAS CON ILUMINACIÓN 


# Se define el vector asociado a las coordenadas de la máscara de filtrado  --> CALIBRANDO COORDENADAS TRANSFORMADA DE FOURIER

#coordenadas_MascaraFiltrado = [-0.000258,0.000158] #coordenada TOMAS EXPERIMENTALES USAF
coordenadas_MascaraFiltrado = [-0.00107,-0.0002] #coordenada TOMAS EXPERIMENTALES CON ILUMINACIÓN
#coordenadas_MascaraFiltrado = [-0.000955,-0.00066] #coordenada TOMAS EXPERIMENTALES 

#Creación de la máscara de filtrado para el proceso de reconstrucción
mascaraReconstruccion_Filtrado = mascaras.funcion_Circulo(radio_pupilaEscaladaPlanoFourier,coordenadas_MascaraFiltrado,
                                                          xx_planoFourier,yy_planoFourier)


campo_Filtrado = mascaraReconstruccion_Filtrado*transformada_FourierReconstruccion

intensidad_CampoFiltrado = (np.abs(campo_Filtrado))**2


#Calculando la transformada de fourier inversa del campo shifteado
campo_FiltradoShifteado = np.fft.fftshift(campo_Filtrado)

#Calculando la transformada de Fourier inversa
campo_Reconstruccion = np.fft.ifft2(campo_FiltradoShifteado)

intensidad_CampoOpticoHolograma = (np.abs(campo_Reconstruccion))**2


#Definiceión de onda plana para eliminar aporte de efectos de interferencia
#vector_cosenosDirectoresRef =  [-0.012656,0.0075936] #coordenada TOMAS EXPERIMENTALES USAF
vector_cosenosDirectoresRef =  [-0.0535,-0.01035] #coordenada TOMAS EXPERIMENTALES CON ILUMINACIÓN
#vector_cosenosDirectoresRef =  [-0.04775,-0.033] #coordenada TOMAS EXPERIMENTALES

#Definición de vector de onda asociado al haz de referencia
vector_ondaRef = [numero_onda_input*vector_cosenosDirectoresRef[0], numero_onda_input*vector_cosenosDirectoresRef[1]]

# Definición de onda plana INVERSA al haz de referencia  
onda_PlanaRefInversa = np.exp(-1j*((vector_ondaRef[0]*xx_mascaraReconstruccion)  + (vector_ondaRef[1]*yy_mascaraReconstruccion)))


# Se multiplica el campo de interés (PLANO MEDICIÓN) con la onda plana generada
matriz_campoNOContribucionOndaPlana = campo_Reconstruccion*onda_PlanaRefInversa

intensidad_matrizCampoNOcontribucionOndaPlana = (np.abs(matriz_campoNOContribucionOndaPlana))**2

#Graficando la transformada de Fourier del holograma
graph.graficar_intensidad(np.log(intensidad_TransformadaFourierRecosntruccion),anchoX_PlanoFourier,altoY_PlanoFourier,
                             "Transformada de Fourier del Holograma",1,1)

graph.graficar_transmitancia(mascaraReconstruccion_Filtrado,anchoX_PlanoFourier,altoY_PlanoFourier,"Mascara Filtrado")


graph.graficar_intensidad(intensidad_CampoFiltrado,anchoX_PlanoFourier,altoY_PlanoFourier,"Tranformada de Fourier filtrada")

#graph.graficar_intensidad(intensidad_CampoOpticoHolograma,ancho_SensorInput,alto_SensorInput,"Campo Optico Filtrado")

#graph.graficar_fase(np.angle(campo_Reconstruccion),ancho_SensorInput,alto_SensorInput,"FASE")

graph.graficar_intensidad(intensidad_matrizCampoNOcontribucionOndaPlana,ancho_SensorInput,alto_SensorInput,
                          "Campo óptico del objeto")

graph.graficar_fase(np.angle(matriz_campoNOContribucionOndaPlana),ancho_SensorInput,alto_SensorInput,
                    "Distribución de fase Campo óptico objeto")


#MULTIPLICANDO POR TÉRMINO DE FASE 
fase_adicionalPrueba = np.exp(1j*(1*(np.pi)/4))

matriz_prueba = matriz_campoNOContribucionOndaPlana*fase_adicionalPrueba

graph.graficar_fase(np.angle(matriz_prueba),ancho_SensorInput,alto_SensorInput,
                    "Distribución de fase Campo óptico objeto + FASE PRUEBA")