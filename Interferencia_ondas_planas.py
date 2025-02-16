
import numpy as np
import matplotlib.pyplot as plt
import LIBRERIA_Funciones_Graficacion as graph
import LIBRERIA_Funciones_importantes as function
import LIBRERIA_Mascaras_Transmitancia as mascaras
import LIBRERIA_Matrices_ABCD_Transferencia_rayos as matriz


""" Definiendo parámetros de sensor  DMM 37UX250-ML """

#Numero de muestras/puntos distribuidos en el ancho del sensor
resolucion_anchoSensorInput = 2048 

#Numero de muestras/puntos distribuidos en el alto del sensor
resolucion_altoSensorInput = 2048

#Tamaño del pixel del sensor
tamaño_PixelSensorInput = 3.45E-6

#Se crea una lista a la cual se le asigna los valores de los delta de muestreo asociados al sensor
deltas_Sensor = [tamaño_PixelSensorInput,tamaño_PixelSensorInput]

#Usando los datos anteriores se calcula el ancho y alto físicos del sensor
ancho_SensorInput = resolucion_anchoSensorInput*tamaño_PixelSensorInput
alto_SensorInput = resolucion_altoSensorInput*tamaño_PixelSensorInput



""" Definiendo parámetro para el tamaño de la pupila """

radio_pupilaInput = 0.07 #Se define variable asociada al radio de la abertura circular que representará
                           # el diafragma.



""" Definición de distancias del arreglo """

distancia_focal01 = 0.01  #Distancia focal asociada a la lente 01

distancia_focal02 = 0.05 #Distancia focal asociada a la lente 02

distancia_propagacionAribitraria = 0.01 #Se define una distancia de propagación arbitraria 



""" Definiendo parámetros de fuente """

#Definiendo la longitud de onda asociada a la fuente en consideración
longitud_onda_input = 632.8E-9 #UNIDADES: m

#Calculando el número de onda
numero_onda_input = (2*np.pi)/longitud_onda_input



""" ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL PRIMER TRAMO ------ """

""" Se calculan las matrices necesarias para estudiar el PRIMER TRAMO del arreglo difractivo
    OBJETO --> *propagación(distancia focal 01)* --> LENTE01 --> *propagación(distancia focal 01)* 
    --> PUPILA"""

#Se calcula la matriz asociada al proceso de propagación desde el plano objeto
#hasta el plano de la lente01
matriz_propagacion01PrimerTramo = matriz.propagacion_MedioHomogeneo(distancia_focal01)

#Se calcula la matriz asociada a la interacción con la lente 01
matriz_lente01 = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focal01)

#Se calcula la matriz asociada al proceso de propagación desde el plano de la lente 01
#hasta el plano de la pupila
matriz_propagacion02PrimerTramo = matriz.propagacion_MedioHomogeneo(distancia_focal01)



""" Calculando la matriz del sistema asociada al PRIMER TRAMO """

#NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
#en el arreglo.

#Se define la lista de matrices
lista_matricesPrimerTramoInvertida = [matriz_propagacion02PrimerTramo,matriz_lente01,matriz_propagacion01PrimerTramo]

#Se calcula la matriz del sistema teniendo en cuenta la lista de matrices definida anteriormente
matriz_SistemaPrimerTramo = matriz.matriz_Sistema(lista_matricesPrimerTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del PRIMER TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del sistema
camino_opticoCentralPrimerTramo = matriz.camino_Optico(lista_matricesPrimerTramoInvertida)



""" ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL SEGUNDO TRAMO ------ """

""" Se calculan las matrices necesarias para estudiar el SEGUNDO TRAMO del arreglo difractivo
    PUPILA--> *propagación(distancia arbitraria d)* --> LENTE02 --> *propagación(distancia focal 02)*
    --> IMAGEN """

#Se calcula la matriz asociada al proceso de propagación desde el plano de la pupila 
# hasta el plano de la lente 02
matriz_propagacion01SegundoTramo = matriz.propagacion_MedioHomogeneo(distancia_propagacionAribitraria)

#Se calcula la matriz asociada a la interacción con la lente 02
matriz_lente02 = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focal02)

#Se calcula la matriz asociada al proceso de propagación desde el plano de la lente 02
#hasta el plano de medición
matriz_propagacion02SegundoTramo = matriz.propagacion_MedioHomogeneo(distancia_focal02)



""" Calculando la matriz del sistema asociada al SEGUNDO TRAMO 
NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
en el arreglo."""

#Creando lista de matrices que describe el arreglo del SEGUNDO TRAMO
lista_matricesSegundoTramoInvertida = [matriz_propagacion02SegundoTramo,matriz_lente02,matriz_propagacion01SegundoTramo]

#Se calcula la matriz del sistema
matriz_SistemaSegundoTramo = matriz.matriz_Sistema(lista_matricesSegundoTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del SEGUNDO TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del Segundo tramo
camino_opticoCentralSegundoTramo = matriz.camino_Optico(lista_matricesSegundoTramoInvertida)



""" ----- EMPIEZA SECCIÓN DE CONFIGURACIÓN DE MALLAS DE PUNTOS PARA CADA PLANO ------ """

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
anchoX_VentanaPlanoPupila = resolucion_anchoSensorInput*deltas_tramoPupilaMedicion[0]
altoY_VentanaPlanoPupila = resolucion_altoSensorInput*deltas_tramoPupilaMedicion[1]

#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoPupila, yy_PlanoPupila = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoPupila,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoPupila)



""" Creando malla de puntos asociada a plano de la máscara u objeto de entrada """

#Se llama función para determinar los deltas de muestreo del tramo PUPILA --> MÁSCARA
deltas_tramoMascaraPupila = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_anchoSensorInput,
                                                                                            anchoX_VentanaPlanoPupila,
                                                                                            matriz_SistemaSegundoTramo[0,1],
                                                                                            longitud_onda_input,
                                                                                            resolucion_altoSensorInput,
                                                                                            altoY_VentanaPlanoPupila)

#Se calcula el ancho de la ventana del plano de la máscara u objeto de entrada
anchoX_VentanaPlanoMascara = resolucion_anchoSensorInput*deltas_tramoMascaraPupila[0]
altoY_VentanaPlanoMascara = resolucion_altoSensorInput*deltas_tramoMascaraPupila[1]

#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoMascara, yy_PlanoMascara = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoPupila,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoPupila)



#DEFINICIÓN DE VECTORES DE ONDA ASOCIADOS EN CADA CASO

k1 = (numero_onda_input,0)

k2 = (numero_onda_input,0)

#PARAMETRO DE PRUEBA: Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoMascara, yy_PlanoMascara = mascaras.malla_Puntos(3000,0.005, 3000,0.005)



# Generar el patrón de interferencia con diferentes fases
patron = function.interferencia_ondas_malla(xx_PlanoMascara,yy_PlanoMascara,k1,k2,fase1=0,fase2=np.pi/4)



# Mostrar el patrón de interferencia
graph.graficar_intensidad(patron,anchoX_VentanaPlanoMascara,altoY_VentanaPlanoMascara,
                          "Patron de interferencia",1,1)

