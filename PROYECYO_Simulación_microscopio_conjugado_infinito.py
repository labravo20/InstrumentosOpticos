##TAREA --> CALCULAR RESOLUCIÓN DEL PROGRAMA TENIENDO EN CUENTA CARACTERÍSTICAS DEL MONTAJE REALIZADO

""" DESCRIPCIÓN 
--> Configuración de microscopio conjugado a infinito con iluminación coherente


En este documento se trabajará sobre un arreglo determinado por la siguiente estructura:
OBJETO --> *propagación(distancia focal MO)* --> OBJETIVO_MICROSCOPIO --> 
*propagación(distancia focal MO)* --> PUPILA (Determinada por la apertura del objetivo
de microscopio en consideración) --> *propagación(distancia arbitraria d)* --> LENTE_TL 
--> *propagación(distancia focal TL)* --> IMAGEN (Detección en cámara)

Debido a la condición de planos conjugados se debe sub dividir el proceso en dos tramos:
TRAMO 01: OBJETO --> *propagación(distancia focal MO)* --> OBJETIVO_MICROSCOPIO --> 
*propagación(distancia focal MO)* --> PUPILA (Determinada por la apertura del objetivo
de microscopio en consideración)

TRAMO 02: *propagación(distancia arbitraria d)* --> LENTE_TL 
--> *propagación(distancia focal TL)* --> IMAGEN (Detección en cámara)

"""

""" Se genera un print de mensaje inicial para verificar correcto funcionamiento del entorno """
print("Inicializando entorno de programación sistema de microscopio conjugado INFINITO...")



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
resolucion_anchoSensorInput = 1280 

#Numero de muestras/puntos distribuidos en el alto del sensor
resolucion_altoSensorInput = 1024

#Tamaño del pixel del sensor #UNIDADES: m
tamaño_PixelSensorInput = 5.2E-6


""" Definiendo parámetros objetivo microscopio """

#Magnificación objetivo de microscopio
magnificacion = 10

#Apertura numérica objetivo de microscopio
apertura_Numerica = 0.25 

#Distancia focal del objetivo de microscopio #UNIDADES: m
### NOTA: Si se desconoce la distancia focal del objetivo de microscopio escribir el número CERO (0) 
distancia_focalMO = 0

#Distancia focal asociada a la lente de tubo --> #UNIDADES: m
distancia_focalTL = 0.18   



""" Definiendo parámetros de la iluminación coherente a utilizar  """

#Definiendo la longitud de onda asociada a la fuente en consideración
longitud_onda_input = 632.8E-9 #UNIDADES: m



""" Definiendo parámetros del montaje  """

#Se define una distancia de propagación arbitraria 'd' para simulación de sistema 4F --> #UNIDADES: m
distancia_propagacionAribitraria = 0.18  #Se hace igual a la distancia focal TL para convertir el sistema en uno invariante al
                                         #desplazamiento --> ELIMINACIÓN DEL TÉRMINO DE FASE PARABÓLICO

' ################ FIN SECCIÓN DE CARACTERIZACIÓN ARREGLO ################## '


""" Calculando parámetros de muestreo del sensor utilizado """

#Se crea una lista a la cual se le asigna los valores de los delta de muestreo asociados al sensor
deltas_Sensor = {"deltaPlanoEntrada_X":tamaño_PixelSensorInput,"deltaPlanoEntrada_Y":tamaño_PixelSensorInput}

#Usando los datos anteriores se calcula el ancho y alto físicos del sensor
ancho_SensorInput = resolucion_anchoSensorInput*tamaño_PixelSensorInput
alto_SensorInput = resolucion_altoSensorInput*tamaño_PixelSensorInput



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


""" Definiendo parámetros de fuente """

#Calculando el número de onda
numero_onda_input = (2*np.pi)/longitud_onda_input



' ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL PRIMER TRAMO ------ '

""" Se calculan las matrices necesarias para estudiar el PRIMER TRAMO del arreglo difractivo
    OBJETO --> *propagación(distancia focal MO)* --> OBJETIVO_MICROSCOPIO --> 
    *propagación(distancia focal MO)* --> PUPILA (Determinada por la apertura del objetivo
    de microscopio en consideración) """

#Se calcula la matriz asociada al proceso de propagación desde el plano objeto
#hasta el plano de la lenteMO
matriz_propagacion01PrimerTramo = matriz.propagacion_MedioHomogeneo(distancia_focalMO)

#Se calcula la matriz asociada a la interacción con la lenteMO
matriz_lenteMO = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focalMO)

#Se calcula la matriz asociada al proceso de propagación desde el plano de la lente 01
#hasta el plano de la pupila
matriz_propagacion02PrimerTramo = matriz.propagacion_MedioHomogeneo(distancia_focalMO)



""" Calculando la matriz del sistema asociada al PRIMER TRAMO """

#NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
#en el arreglo.

#Se define la lista de matrices
lista_matricesPrimerTramoInvertida = [matriz_propagacion02PrimerTramo,matriz_lenteMO,matriz_propagacion01PrimerTramo]

#Se calcula la matriz del sistema teniendo en cuenta la lista de matrices definida anteriormente
matriz_SistemaPrimerTramo = matriz.matriz_Sistema(lista_matricesPrimerTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del PRIMER TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del sistema
camino_opticoCentralPrimerTramo = matriz.camino_Optico(lista_matricesPrimerTramoInvertida)


' ------ FIN SECCIÓN DE CÁLCULO MATRICIAL PRIMER TRAMO ------ '


' ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL SEGUNDO TRAMO ------ '

""" Se calculan las matrices necesarias para estudiar el SEGUNDO TRAMO del arreglo difractivo
    *propagación(distancia arbitraria d)* --> LENTE_TL --> *propagación(distancia focal TL)* 
    --> IMAGEN (Detección en cámara)"""

#Se calcula la matriz asociada al proceso de propagación desde el plano de la pupila 
# hasta el plano de la lente TL
matriz_propagacion01SegundoTramo = matriz.propagacion_MedioHomogeneo(distancia_propagacionAribitraria)

#Se calcula la matriz asociada a la interacción con la lente TL
matriz_lenteTL = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focalTL)

#Se calcula la matriz asociada al proceso de propagación desde el plano de la lente 02
#hasta el plano de medición
matriz_propagacion02SegundoTramo = matriz.propagacion_MedioHomogeneo(distancia_focalTL)



""" Calculando la matriz del sistema asociada al SEGUNDO TRAMO 
NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
en el arreglo."""

#Creando lista de matrices que describe el arreglo del SEGUNDO TRAMO
lista_matricesSegundoTramoInvertida = [matriz_propagacion02SegundoTramo,matriz_lenteTL,matriz_propagacion01SegundoTramo]

#Se calcula la matriz del sistema
matriz_SistemaSegundoTramo = matriz.matriz_Sistema(lista_matricesSegundoTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del SEGUNDO TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del Segundo tramo
camino_opticoCentralSegundoTramo = matriz.camino_Optico(lista_matricesSegundoTramoInvertida)

' ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL SEGUNDO TRAMO ------ '


' ----- EMPIEZA SECCIÓN DE CONFIGURACIÓN DE MALLAS DE PUNTOS PARA CADA PLANO ------ '

""" Creando mallas de puntos asociada a plano del sensor """
#Para poder conocer las condicones de creación de las mallas de puntos se debe considerar
#que los planos de interés son: PLANO OBJETO --> PLANO ANTERIOR A PUPILA --> PLANO SENSOR.
#Ahora bien, en este caso los datos asociados al sensor son aquellos que permitirán conocer
#las condiciones de todas las mallas de puntos, lo cual se logra usando las condiciones de
#producto espacio frecuencia en el caso específico de la transformada de Fresnel.
#Se calcula la malla de puntos asociada al plano de medición

#Creando malla de puntos plano del sensor/MEDICIÓN
xx_PlanoMedicion,yy_PlanoMedicion = mascaras.malla_Puntos(resolucion_anchoSensorInput,ancho_SensorInput,
                                            resolucion_altoSensorInput,alto_SensorInput)



""" Creando malla de puntos asociada a plano anterior a la pupila """

#Se llama función para determinar los deltas de muestreo del tramo SENSOR --> PUPILA
deltas_tramoPupilaMedicion = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_anchoSensorInput,
                                                                                             ancho_SensorInput,
                                                                                             matriz_SistemaSegundoTramo[0,1],
                                                                                             longitud_onda_input,
                                                                                             resolucion_altoSensorInput,
                                                                                             alto_SensorInput)

#Se calcula el ancho de la ventana del plano anterior a pupila
anchoX_VentanaPlanoPupila = resolucion_anchoSensorInput*deltas_tramoPupilaMedicion["deltaPlanoEntrada_X"]
altoY_VentanaPlanoPupila = resolucion_altoSensorInput*deltas_tramoPupilaMedicion["deltaPlanoEntrada_Y"]


#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoPupila, yy_PlanoPupila = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoPupila,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoPupila)



""" Creando malla de puntos asociada a plano de la máscara u objeto de entrada """

#Se llama función para determinar los deltas de muestreo del tramo PUPILA --> MÁSCARA
deltas_tramoMascaraPupila = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_anchoSensorInput,
                                                                                            anchoX_VentanaPlanoPupila,
                                                                                            matriz_SistemaPrimerTramo[0,1],
                                                                                            longitud_onda_input,
                                                                                            resolucion_altoSensorInput,
                                                                                            altoY_VentanaPlanoPupila)

#Se calcula el ancho de la ventana del plano de la máscara u objeto de entrada
anchoX_VentanaPlanoMascara = resolucion_anchoSensorInput*deltas_tramoMascaraPupila["deltaPlanoEntrada_X"]
altoY_VentanaPlanoMascara = resolucion_altoSensorInput*deltas_tramoMascaraPupila["deltaPlanoEntrada_Y"]


#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoMascara, yy_PlanoMascara = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoPupila,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoPupila)


' ----- FIN SECCIÓN DE CONFIGURACIÓN DE MALLAS DE PUNTOS PARA CADA PLANO ------ '


' ------ EMPIEZA SECCIÓN DE CÁLCULO RESULTADO DIFRACTIVO DE CADA TRAMO ------ '

""" Creando máscara de transmitancia asociada al objeto de estudio en el arreglo """

# Cargar la imagen PNG como máscara de transmitancia
ruta_imagen_png = "/home/labravo/Downloads/USAF_T-20.jpg"  # Especifica la ruta de imagen
mascara = function.cargar_imagen_png(ruta_imagen_png, resolucion_anchoSensorInput,resolucion_altoSensorInput)



""" Se calcula el resultado del proceso difractivo del PRIMER TRAMO"""

#Se calcula el campo de salida/en plano de la pupila --> Campo resultante del primer tramo 
campo_PlanoPupila = matriz.matriz_ABCD_Difraccion_Sensor(camino_opticoCentralPrimerTramo,mascara,
                                                 matriz_SistemaPrimerTramo[0,0],
                                                 matriz_SistemaPrimerTramo[0,1],
                                                 matriz_SistemaPrimerTramo[1,1],
                                                 xx_PlanoMascara,yy_PlanoMascara,
                                                 xx_PlanoPupila,yy_PlanoPupila,numero_onda_input,
                                                 deltas_tramoPupilaMedicion)



""" Se crea pupila para delimitar dominio del arreglo asociado a lentes FINITAS
NOTA: La pupila va a construirse a partir de la malla de puntos asociada al plano de la pupila."""

#Creación de máscara circular que representará el diafragma 
pupila = mascaras.funcion_Circulo(radio_pupilaInput, None, xx_PlanoPupila, yy_PlanoPupila)



""" Se calcula el campo de entrada para el SEGUNDO TRAMO del arreglo 
NOTA: El campo que llega a la pupila e interactua con la máscara asociada a la pupila
será el campo de entrada para el segundo tramo."""

#Calculando el campo de entrada al SEGUNDO TRAMO del arreglo
campo_entradaSegundoTramo = campo_PlanoPupila*pupila

#Calculando la intensidad del campo de entrada al SEGUNDO TRAMO del arreglo
intensidad_campoEntradaSegundoTramo = np.abs(campo_entradaSegundoTramo)**2



""" Se calcula el resultado del proceso difractivo del SEGUNDO TRAMO """

#Se calcula el campo de salida/en plano de medición --> Campo resultante de la difracción 
campo_PlanoMedicion = matriz.matriz_ABCD_Difraccion_Sensor_Shift(camino_opticoCentralSegundoTramo,
                                                    campo_entradaSegundoTramo,
                                                    matriz_SistemaSegundoTramo[0,0],
                                                    matriz_SistemaSegundoTramo[0,1],
                                                    matriz_SistemaSegundoTramo[1,1],xx_PlanoPupila,
                                                    yy_PlanoPupila,xx_PlanoMedicion,yy_PlanoMedicion,
                                                    numero_onda_input,deltas_Sensor)



#Se calcula la amplitud del campo de salida
amplitud_campoPlanoMedicion = np.abs(campo_PlanoMedicion)


#Se calcula la intensidad del campo de salida
intensidad_campoPlanoMedicion = amplitud_campoPlanoMedicion**2


' ------ FIN SECCIÓN DE CÁLCULO RESULTADO DIFRACTIVO DE CADA TRAMO ------ '


' ------ EMPIEZA SECCIÓN DE GRAFICACIÓN ------ '

""" Graficando máscara de transmitancia """
graph.graficar_transmitancia(mascara,anchoX_VentanaPlanoMascara,altoY_VentanaPlanoMascara,"Objeto de análisis")


""" Graficando intensidad del campo que entra al SEGUNDO TRAMO del arreglo"""
graph.graficar_intensidad(intensidad_campoEntradaSegundoTramo,anchoX_VentanaPlanoPupila,altoY_VentanaPlanoPupila,
                          "Intensidad del campo en plano de la PUPILA",1,0.6)


""" Graficando la intensidad del campo de salida del arreglo """
graph.graficar_intensidad(intensidad_campoPlanoMedicion,ancho_SensorInput,alto_SensorInput,"Intensidad del campo en plano de MEDICIÓN")


""" Graficar fase del campo de salida del arreglo """
graph.graficar_fase(np.angle(campo_PlanoMedicion),ancho_SensorInput,alto_SensorInput,"Fase del campo en el plano de medición")
