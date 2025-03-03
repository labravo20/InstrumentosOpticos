""" DESCRIPCIÓN 
--> Reconstrucción de Holograma generado a partir de la
interferencia entre la imagen resultante de una configuración de microscopio conjugado a infinito con iluminación coherente,  una
onda plana asociada al HAZ DE REFERENCIA.


En este documento se trabajará sobre un arreglo determinado por la siguiente estructura:
OBJETO --> *propagación(distancia focal MO)* --> OBJETIVO_MICROSCOPIO --> 
*propagación(distancia focal MO)* --> PUPILA (Determinada por la apertura del objetivo
de microscopio en consideración) --> *propagación(distancia arbitraria d)* --> LENTE_TL 
--> *propagación(distancia focal TL)* --> IMAGEN 

Debido a la condición de planos conjugados se debe sub dividir el proceso en dos tramos:
TRAMO 01: OBJETO --> *propagación(distancia focal MO)* --> OBJETIVO_MICROSCOPIO --> 
*propagación(distancia focal MO)* --> PUPILA (Determinada por la apertura del objetivo
de microscopio en consideración)

TRAMO 02: *propagación(distancia arbitraria d)* --> LENTE_TL 
--> *propagación(distancia focal TL)* --> IMAGEN

Posteriormente al tramo dos se procede a la formación del holograma, de forma que 
en el plano imagen se genera interferencia entre el campo de salida del anterior arreglo 
(haz objeto) y una onda plana de referencia (haz referencia) --> Esta interferencia se realiza
en el plano del sensor.

RECONSTRUCCIÓN DEL HOLOGRAMA: Se debe recuperar la información del campo complejo del objeto, para lograrlo se estudia
la transformada de Fourier del holograma resultante, se filtra el componente que contiene la información del campo, y finalmente
se contrarestan los efectos de interferencia con la onda plana usada antes para formar el holograma.

"""

""" Se genera un print de mensaje inicial para verificar correcto funcionamiento del entorno """
print("Inicializando entorno de programación sistema simulación holografía con MO CONJUGADO INFINITO...")


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
resolucion_anchoSensorInput = 1024

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
distancia_focalTL = 0.2   



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

' ------ FIN SECCIÓN DE CÁLCULO MATRICIAL SEGUNDO TRAMO ------ '


' ----- EMPIEZA SECCIÓN DE CONFIGURACIÓN DE MALLAS DE PUNTOS PARA CADA PLANO ------ '

""" Creando mallas de puntos asociada a plano del sensor """
#Para poder conocer las condicones de creación de las mallas de puntos se debe considerar
#que los planos de interés son: PLANO OBJETO --> PLANO ANTERIOR A PUPILA --> PLANO SENSOR.
#Ahora bien, en este caso los datos asociados al sensor son aquellos que permitirán conocer
#las condiciones de todas las mallas de puntos, lo cual se logra usando las condiciones de
#producto espacio frecuencia en el caso específico de la transformada de Fresnel.
#Se calcula la malla de puntos asociada al plano de medición

#Creando malla de puntos plano del sensor/MEDICIÓN
xx_PlanoImagen,yy_PlanoImagen = mascaras.malla_Puntos(resolucion_anchoSensorInput,ancho_SensorInput,
                                            resolucion_altoSensorInput,alto_SensorInput)



""" Creando malla de puntos asociada a plano anterior a la pupila """

#Se llama función para determinar los deltas de muestreo del tramo SENSOR --> PUPILA
deltas_tramoPupilaImagen = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_anchoSensorInput,
                                                                                             ancho_SensorInput,
                                                                                             matriz_SistemaSegundoTramo[0,1],
                                                                                             longitud_onda_input,
                                                                                             resolucion_altoSensorInput,
                                                                                             alto_SensorInput)

#Se calcula el ancho de la ventana del plano anterior a pupila
anchoX_VentanaPlanoPupila = resolucion_anchoSensorInput*deltas_tramoPupilaImagen["deltaPlanoEntrada_X"]
altoY_VentanaPlanoPupila = resolucion_altoSensorInput*deltas_tramoPupilaImagen["deltaPlanoEntrada_Y"]


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
                                                 deltas_tramoPupilaImagen)



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
campo_PlanoImagen = matriz.matriz_ABCD_Difraccion_Sensor_Shift(camino_opticoCentralSegundoTramo,
                                                    campo_entradaSegundoTramo,
                                                    matriz_SistemaSegundoTramo[0,0],
                                                    matriz_SistemaSegundoTramo[0,1],
                                                    matriz_SistemaSegundoTramo[1,1],xx_PlanoPupila,
                                                    yy_PlanoPupila,xx_PlanoImagen,yy_PlanoImagen,
                                                    numero_onda_input,deltas_Sensor)



#Se calcula la amplitud del campo de salida
amplitud_campoPlanoImagen = np.abs(campo_PlanoImagen)


#Se calcula la intensidad del campo de salida
intensidad_campoPlanoImagen = amplitud_campoPlanoImagen**2


' ------ FIN SECCIÓN DE CÁLCULO RESULTADO DIFRACTIVO DE CADA TRAMO ------ '


' ------ EMPIEZA SECCIÓN DE GRAFICACIÓN ------ '

""" Graficando máscara de transmitancia """
graph.graficar_transmitancia(mascara,anchoX_VentanaPlanoMascara,altoY_VentanaPlanoMascara,"Objeto de análisis")


""" Graficando la intensidad del campo de salida del arreglo """
graph.graficar_intensidad(intensidad_campoPlanoImagen,ancho_SensorInput,alto_SensorInput,"Intensidad del campo en plano IMÁGEN")


""" Graficando la información de fase de el campo en el plano imágen """
graph.graficar_fase(np.angle(campo_PlanoImagen),ancho_SensorInput,alto_SensorInput,"Fase de campo en Plano IMÁGEN")

' ------ FIN SECCIÓN DE GRAFICACIÓN ------ '



' ------ EMPIEZA SECCIÓN DE INTERFERENCIA PARA FORMACIÓN HOLOGRAMA ------ '

""" Definición del haz de referencia --> ONDA PLANA """

# Definición de onda plana INVERSA al haz de referencia  

### Se calcula cuál debe ser el orden de magnitud de la amplitud de la onda de referencia para que sea comparable con el campo
### asociado a la onda objeto
amplitudOndaReferencia = np.max(np.abs(campo_PlanoImagen))

onda_PlanaReferencia = (amplitudOndaReferencia)*np.exp(1j * numero_onda_input * 
                            ((xx_PlanoImagen * cosenoDirector_HazReferenciaAlfa) + 
                            (yy_PlanoImagen * cosenoDirector_HazReferenciaBeta)))



""" Se genera interferencia entre el haz de referencia (ONDA PLANA REF) y el haz objeto (CAMPO EN PLANO IMAGEN)"""

#Interferencia entre haz de referencia y haz objeto 
campo_interferenciaObjetoReferencia = onda_PlanaReferencia + campo_PlanoImagen
intensidad_interferenciaObjetoReferencia = (np.abs(campo_interferenciaObjetoReferencia))**2

' ------ FIN SECCIÓN DE INTERFERENCIA PARA FORMACIÓN HOLOGRAMA ------ '


' ------ EMPIEZA SECCIÓN DE GRAFICACIÓN ------ '

""" Graficando el campo resultante de la formación del holograma (interferencia entre haz de referencia y haz objeto)"""
graph.graficar_intensidad(intensidad_interferenciaObjetoReferencia,ancho_SensorInput,alto_SensorInput,"HOLOGRAMA")

' ------ FIN SECCIÓN DE GRAFICACIÓN ------ '


######################## INICIA SECCIÓN DE RECONSTRUCCIÓN HOLOGRAMA ##############################



' ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL PRIMER TRAMO RECONSTRUCCIÓN------ '


""" Se calculan las matrices necesarias para estudiar el PRIMER TRAMO del arreglo difractivo
    MASCARA --> *Propagación distancia focalMO* --> LENTE_MO  --> PUPILA """

#Se calcula la matriz asociada al proceso de propagación desde el plano objeto
#hasta el plano de la lente
matrizReconstruccion_propagacionPrimerTramo = matriz.propagacion_MedioHomogeneo(distancia_focalMO)

matrizReconstruccion_lente = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focalMO)

matrizReconstruccion_propagacion02 = matriz.propagacion_MedioHomogeneo(distancia_focalMO)



""" Calculando la matriz del sistema asociada al PRIMER TRAMO """

#NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
#en el arreglo.

#Se define la lista de matrices
lista_matricesReconstruccionPrimerTramoInvertida = [matrizReconstruccion_propagacion02,matrizReconstruccion_lente,matrizReconstruccion_propagacionPrimerTramo]

#En este caso la matriz del sistema es equivalente a la única matriz presente
matrizReconstruccion_SistemaPrimerTramo = matriz.matriz_Sistema(lista_matricesReconstruccionPrimerTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del PRIMER TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del sistema
caminoReconstruccion_opticoCentralPrimerTramo = matriz.camino_Optico(lista_matricesReconstruccionPrimerTramoInvertida)


' ------ FIN SECCIÓN DE CÁLCULO MATRICIAL PRIMER TRAMO RECONSTRUCCIÓN------- '


' ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL SEGUNDO TRAMO RECONSTRUCCIÓN------- '

""" Se calculan las matrices necesarias para estudiar el SEGUNDO TRAMO del arreglo difractivo
    PUPILA--> *propagación(distancia focal TL)* --> LENTE02 --> *propagación(distancia focal TL)*
    --> IMAGEN """

#Se calcula la matriz asociada al proceso de propagación desde el plano de la pupila 
# hasta el plano de la lente 02
matrizReconstruccion_propagacion01SegundoTramo = matriz.propagacion_MedioHomogeneo(distancia_propagacionAribitraria)

#Se calcula la matriz asociada a la interacción con la lente 02
matrizReconstruccion_lente02 = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focalTL)

#Se calcula la matriz asociada al proceso de propagación desde el plano de la lente 02
#hasta el plano de medición
matrizReconstruccion_propagacion02SegundoTramo = matriz.propagacion_MedioHomogeneo(distancia_focalTL)



""" Calculando la matriz del sistema asociada al SEGUNDO TRAMO 
NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
en el arreglo."""

#Creando lista de matrices que describe el arreglo del SEGUNDO TRAMO
lista_matricesReconstruccionSegundoTramoInvertida = [matrizReconstruccion_propagacion02SegundoTramo,matrizReconstruccion_lente02,matrizReconstruccion_propagacion01SegundoTramo]

#Se calcula la matriz del sistema
matrizReconstruccion_SistemaSegundoTramo = matriz.matriz_Sistema(lista_matricesReconstruccionSegundoTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del SEGUNDO TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del Segundo tramo
caminoReconstruccion_opticoCentralSegundoTramo = matriz.camino_Optico(lista_matricesReconstruccionSegundoTramoInvertida)

' ------ FIN SECCIÓN DE CÁLCULO MATRICIAL SEGUNDO TRAMO RECONSTRUCCIÓN- ------ '


' ----- EMPIEZA SECCIÓN DE CONFIGURACIÓN DE MALLAS DE PUNTOS PARA CADA PLANO RECONSTRUCCIÓN- ------ '

""" Creando máscara de plano objeto  """

# Crear la malla de puntos
xx_mascaraReconstruccion, yy_mascaraReconstruccion = mascaras.malla_Puntos(resolucion_anchoSensorInput, ancho_SensorInput,resolucion_altoSensorInput,
                                               alto_SensorInput)


""" Creando malla de puntos asociada a plano de la PUPILA """

#Se llama función para determinar los deltas de muestreo del tramo PUPILA --> MÁSCARA
deltasReconstruccion_tramoMascaraPupila = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_altoSensorInput,
                                                                                            ancho_SensorInput,
                                                                                            matrizReconstruccion_SistemaPrimerTramo[0,1],
                                                                                            longitud_onda_input,
                                                                                            resolucion_altoSensorInput,
                                                                                            alto_SensorInput)

#Se calcula el ancho de la ventana del plano de la máscara u objeto de entrada
anchoX_VentanaPlanoPupilaReconstruccion = resolucion_anchoSensorInput*deltasReconstruccion_tramoMascaraPupila["deltaPlanoEntrada_X"]
altoY_VentanaPlanoPupilaReconstruccion  = resolucion_altoSensorInput*deltasReconstruccion_tramoMascaraPupila["deltaPlanoEntrada_Y"]

#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoPupilaReconstruccion, yy_PlanoPupilaReconstruccion = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoPupilaReconstruccion,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoPupilaReconstruccion)


""" Creando malla de puntos asociada a plano Medición """

#Se llama función para determinar los deltas de muestreo del tramo SENSOR --> PUPILA
deltasReconstruccion_tramoPupilaMedicion = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_anchoSensorInput,
                                                                                             anchoX_VentanaPlanoPupilaReconstruccion,
                                                                                             matrizReconstruccion_SistemaSegundoTramo[0,1],
                                                                                             longitud_onda_input,
                                                                                             resolucion_altoSensorInput,
                                                                                             altoY_VentanaPlanoPupilaReconstruccion)

#Se calcula el ancho de la ventana del plano anterior a pupila
anchoX_VentanaPlanoMedicionReconstruccion = resolucion_anchoSensorInput*deltasReconstruccion_tramoPupilaMedicion["deltaPlanoEntrada_X"]
altoY_VentanaPlanoMedicionReconstruccion = resolucion_altoSensorInput*deltasReconstruccion_tramoPupilaMedicion["deltaPlanoEntrada_Y"]

#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoMedicionReconstruccion, yy_PlanoMedicionReconstruccion = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoMedicionReconstruccion,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoMedicionReconstruccion)

' ----- FIN SECCIÓN DE CONFIGURACIÓN DE MALLAS DE PUNTOS PARA CADA PLANO RECONSTRUCCIÓN------- '


' ------ EMPIEZA SECCIÓN DE CÁLCULO RESULTADO DIFRACTIVO DE CADA TRAMO RECONSTRUCCIÓN------- '

""" Creación de mascara de transmitancia """

# Cargar la imagen PNG como máscara de transmitancia
ruta_imagen_png = "/home/labravo/Downloads/Holograma007.tif"  # Especifica la ruta de tu imagen
#ruta_imagen_png = "/home/labravo/Downloads/Mi_primer_holograma.PNG"  # Especifica la ruta de tu imagen
mascaraReconstruccion = function.cargar_imagen_png(ruta_imagen_png,resolucion_anchoSensorInput,resolucion_altoSensorInput)
#mascaraReconstruccion = intensidad_interferenciaObjetoReferencia



""" Se calcula el resultado del proceso difractivo del PRIMER TRAMO"""

#Se calcula el campo de salida/en plano de la lente --> Campo resultante del primer tramo 
campoReconstruccion_PlanoPupila = matriz.matriz_ABCD_Difraccion_Sensor(caminoReconstruccion_opticoCentralPrimerTramo,mascaraReconstruccion,
                                                 matrizReconstruccion_SistemaPrimerTramo[0,0],
                                                 matrizReconstruccion_SistemaPrimerTramo[0,1],
                                                 matrizReconstruccion_SistemaPrimerTramo[1,1],xx_mascaraReconstruccion,
                                                 yy_mascaraReconstruccion,
                                                 xx_PlanoPupilaReconstruccion,yy_PlanoPupilaReconstruccion,numero_onda_input,
                                                 deltasReconstruccion_tramoMascaraPupila)



""" Se crea diafragma para delimitar dominio del arreglo asociado a lentes FINITAS
NOTA: El diafragma va a construirse a partir de la malla de puntos asociada al plano de la Lente."""

#Creación de máscara circular que representará el diafragma 
diafragmaReconstruccion = mascaras.funcion_Circulo(radio_pupilaInput, None, xx_PlanoPupilaReconstruccion, yy_PlanoPupilaReconstruccion)


""" Se calcula el campo de entrada para el SEGUNDO TRAMO del arreglo 
NOTA: El campo que llega a la lente e interactua con la máscara asociada al diafragma
será el campo de entrada para el segundo tramo."""

#Calculando el campo de entrada al SEGUNDO TRAMO del arreglo
campoReconstruccion_entradaSegundoTramo_sinFiltro = campoReconstruccion_PlanoPupila*diafragmaReconstruccion

#Calculando la intensidad del campo de entrada al SEGUNDO TRAMO del arreglo
intensidad_campoReconstruccionEntradaSegundoTramo_SinFiltro = np.abs(campoReconstruccion_entradaSegundoTramo_sinFiltro)**2



"""Se diseña una máscara de transmitancia para filtrar la imágen gemela de interés """

# Definición del Lado del rectangulo (dimensión de máscara de transmitancia)
lado_mascaraReconstruccion = 0.00033

# # Se calculan las coordenadas de la máscara de filtrado
# ## NOTA: Las coordenadas se calculan teniendo en cuenta la relación de cosenos directores con las coordenadas 
# ## frecuenciales.
# x_Filtrado = cosenoDirector_HazReferenciaAlfa*distancia_focalMO
# y_Filtrado = cosenoDirector_HazReferenciaBeta*distancia_focalMO

# # Se define el vector asociado a las coordenadas de la máscara de filtrado
# coordenadas_MascaraFiltrado = [-y_Filtrado,-x_Filtrado]
#coordenadas_MascaraFiltrado = [-0.00067,-0.00066]
coordenadas_MascaraFiltrado = [0.00029,-0.00018]

#Creación de la máscara de filtrado para el proceso de reconstrucción
mascaraReconstruccion_Filtrado = mascaras.funcion_Rectangulo(lado_mascaraReconstruccion,lado_mascaraReconstruccion,
                                                             coordenadas_MascaraFiltrado,xx_PlanoPupilaReconstruccion,
                                                             yy_PlanoPupilaReconstruccion)



""" Se calcula el campo de entrada para el SEGUNDO TRAMO del arreglo 
NOTA: El campo que llega a la pupila e interactua con la máscara asociada a la pupila
será el campo de entrada para el segundo tramo."""

#Calculando el campo de entrada al SEGUNDO TRAMO del arreglo
campoReconstruccion_entradaSegundoTramo = campoReconstruccion_PlanoPupila*diafragmaReconstruccion*mascaraReconstruccion_Filtrado

#Calculando la intensidad del campo de entrada al SEGUNDO TRAMO del arreglo
intensidad_campoReconstruccionEntradaSegundoTramo = np.abs(campoReconstruccion_entradaSegundoTramo)**2


""" Se calcula el resultado del proceso difractivo del SEGUNDO TRAMO """

#Se calcula el campo de salida/en plano de medición --> Campo resultante de la difracción 
campoReconstruccion_PlanoMedicion = matriz.matriz_ABCD_Difraccion_Sensor_Shift(caminoReconstruccion_opticoCentralSegundoTramo,
                                                    campoReconstruccion_entradaSegundoTramo,
                                                    matrizReconstruccion_SistemaSegundoTramo[0,0],
                                                    matrizReconstruccion_SistemaSegundoTramo[0,1],
                                                    matrizReconstruccion_SistemaSegundoTramo[1,1],xx_PlanoPupilaReconstruccion,
                                                    yy_PlanoPupilaReconstruccion,xx_PlanoMedicionReconstruccion,
                                                    yy_PlanoMedicionReconstruccion,
                                                    numero_onda_input,deltasReconstruccion_tramoPupilaMedicion)


#Se calcula la intensidad del campo de salida
intensidad_campoReconstruccionPlanoMedicion = np.abs(campoReconstruccion_PlanoMedicion)**2

' ------ FIN SECCIÓN DE CÁLCULO RESULTADO DIFRACTIVO DE CADA TRAMO RECONSTRUCCIÓN------- '

' ------ EMPIEZA SECCIÓN DE GRAFICACIÓN ------- '

""" Graficando máscara de transmitancia"""
graph.graficar_transmitancia(mascaraReconstruccion,ancho_SensorInput,alto_SensorInput,"Imágen de entrada para reconstrucción")


""" Graficando intensidad del campo que entra al SEGUNDO TRAMO del arreglo"""
graph.graficar_intensidad(intensidad_campoReconstruccionEntradaSegundoTramo_SinFiltro,anchoX_VentanaPlanoPupilaReconstruccion,
                              altoY_VentanaPlanoPupilaReconstruccion,"Transformada de Fourier del objeto",1,0.0001)


""" Graficando máscara de filtrado"""
graph.graficar_transmitancia(mascaraReconstruccion_Filtrado,anchoX_VentanaPlanoPupilaReconstruccion,
                             altoY_VentanaPlanoPupilaReconstruccion,"Máscara de filtrado")


""" Graficando campo asociado a transformada de Fourier filtrado para encontrar imágen real """
graph.graficar_intensidad(intensidad_campoReconstruccionEntradaSegundoTramo,anchoX_VentanaPlanoPupilaReconstruccion,
                             altoY_VentanaPlanoPupilaReconstruccion,"Transformada de Fourier del objeto",1,0.01)


""" Graficando la intensidad del campo de salida del arreglo """
graph.graficar_intensidad(intensidad_campoReconstruccionPlanoMedicion,anchoX_VentanaPlanoMedicionReconstruccion,
                          altoY_VentanaPlanoMedicionReconstruccion,"Intensidad campo de salida")

graph.graficar_fase(np.angle(campoReconstruccion_PlanoMedicion),anchoX_VentanaPlanoMedicionReconstruccion,altoY_VentanaPlanoMedicionReconstruccion,
                    "FASE CON CONTRIBUCION INTERFERENCIA")

camino_OpticoCampoObjeto = ((2*np.pi)/longitud_onda_input)*(np.angle(campoReconstruccion_PlanoMedicion))

graph.graficar_longitudCaminoOptico(camino_OpticoCampoObjeto,anchoX_VentanaPlanoMedicionReconstruccion,altoY_VentanaPlanoMedicionReconstruccion,
                          "Longitud de camino óptico")


' ------ FIN SECCIÓN DE GRAFICACIÓN ------- '


' ------ FIN SECCIÓN DE ELIMINACIÓN APORTE INTERFERENCIA ------- '

'Se realiza a continuación procedimiento para eliminar contribución de haz REFERENCIA'
# -----    ÁNGULO ENTRE HAZ DE REFERENCIA Y HAZ OBJETO --> 4.2255°  ------ #

#Definición de vector con cosenos directores asociado al haz referencia
#vector_cosenosDirectoresRef = [cosenoDirector_HazReferenciaAlfa,cosenoDirector_HazReferenciaBeta] 
#vector_cosenosDirectoresRef = [-0.0335,-0.033] 
vector_cosenosDirectoresRef = [0.0145,0.0019] 

#Definición de vector de onda asociado al haz de referencia
vector_ondaRef = [numero_onda_input*vector_cosenosDirectoresRef[0], numero_onda_input*vector_cosenosDirectoresRef[1]]

# Definición de onda plana INVERSA al haz de referencia  
onda_PlanaRefInversa = np.exp(-1j*((vector_ondaRef[0]*xx_PlanoMedicionReconstruccion)  + (vector_ondaRef[1]*yy_PlanoMedicionReconstruccion)))


# Se multiplica el campo de interés (PLANO MEDICIÓN) con la onda plana generada
matriz_campoNOContribucionOndaPlana = campoReconstruccion_PlanoMedicion*onda_PlanaRefInversa

intensidad_campoNOContribucionOndaPlana = (np.abs(matriz_campoNOContribucionOndaPlana))**2

graph.graficar_intensidad(intensidad_campoNOContribucionOndaPlana,anchoX_VentanaPlanoMedicionReconstruccion,altoY_VentanaPlanoMedicionReconstruccion,
                          "Imagen de objeto recuperado",1,1)


graph.graficar_fase(np.angle(matriz_campoNOContribucionOndaPlana),anchoX_VentanaPlanoMedicionReconstruccion,
                    altoY_VentanaPlanoMedicionReconstruccion,"Fase del campo complejo recuperado")


## FUNCIÓN PARA CÁLCULO DE CAMINO ÓPTICO ##

camino_OpticoCampoObjeto = ((2*np.pi)/longitud_onda_input)*(np.angle(matriz_campoNOContribucionOndaPlana))

graph.graficar_longitudCaminoOptico(camino_OpticoCampoObjeto,anchoX_VentanaPlanoMedicionReconstruccion,altoY_VentanaPlanoMedicionReconstruccion,
                          "Longitud de camino óptico")


