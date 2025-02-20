""" DESCRIPCIÓN 
--> Configuración de microscopio conjugado a finito con iluminación coherente


En este documento se trabajará sobre un arreglo determinado por la siguiente estructura:
OBJETO --> *propagación(distancia objeto)* --> OBJETIVO_MICROSCOPIO --> 
PUPILA (Determinada por la apertura del objetivo de microscopio en consideración) --> 
*propagación(distancia imagen)* --> IMAGEN 

Debido a la condición de planos conjugados se debe sub dividir el proceso en dos tramos:

TRAMO 01: 
OBJETO --> *propagación(distancia objeto)* --> OBJETIVO_MICROSCOPIO --> 
PUPILA (Determinada por la apertura del objetivo de microscopio en consideración)

TRAMO 02: 
PUPILA (Determinada por la apertura del objetivo de microscopio en consideración) --> 
*propagación(distancia imagen)* --> IMAGEN 

Posteriormente al tramo dos se procede a la formación del holograma, de forma que 
en el plano imagen se genera interferencia entre el campo de salida del anterior arreglo 
(haz objeto) y una onda plana de referencia (haz referencia). Después del paso anterior se 
procede a propagar el resultado de la interferencia una distancia arbitraria 'd'(implica presencia
de efectos difractivos), para llegar al plano donde se forma el holograma.

"""

""" Se genera un print de mensaje inicial para verificar correcto funcionamiento del entorno """
print("Inicializando entorno de programación sistema simulación holografía...")



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
resolucion_anchoSensorInput = 2848 

#Numero de muestras/puntos distribuidos en el alto del sensor
resolucion_altoSensorInput = 2848

#Tamaño del pixel del sensor #UNIDADES: m
tamaño_PixelSensorInput = 2.74E-6


""" Definiendo parámetros objetivo microscopio """

#Magnificación objetivo de microscopio
magnificacion = 60

#Apertura numérica objetivo de microscopio
apertura_Numerica = 0.85 

#Distancia focal del objetivo de microscopio #UNIDADES: m
distancia_focalMO = 0.2


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
angulo_HazReferencia = 6 #UNIDADES: Grados


""" Distancia de propagación para formación del holograma (PLANO IMAGEN --> PLANO HOLOGRAMA)"""
distancia_ImagenHolograma = 0.08 #UNIDADES: m

' ################ FIN SECCIÓN DE CARACTERIZACIÓN ARREGLO ################## '



""" Calculando parámetros de muestreo del sensor utilizado """

#Se crea una lista a la cual se le asigna los valores de los delta de muestreo asociados al sensor
deltas_Sensor = {"deltaPlanoEntrada_X":tamaño_PixelSensorInput,"deltaPlanoEntrada_Y":tamaño_PixelSensorInput}

#Usando los datos anteriores se calcula el ancho y alto físicos del sensor
ancho_SensorInput = resolucion_anchoSensorInput*tamaño_PixelSensorInput
alto_SensorInput = resolucion_altoSensorInput*tamaño_PixelSensorInput



""" Definición de distancias del arreglo """

#Distancia focal asociada a la lente de tubo --> #UNIDADES: m
distancia_Objeto = ((magnificacion + 1)*distancia_focalMO)/magnificacion   

#Se define una distancia de propagación arbitraria --> #UNIDADES: m
distancia_Imagen = (magnificacion + 1)*distancia_focalMO 



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
    OBJETO --> *propagación(distancia objeto)* --> OBJETIVO_MICROSCOPIO --> 
    PUPILA (Determinada por la apertura del objetivo de microscopio en consideración) """

#Se calcula la matriz asociada al proceso de propagación desde el plano objeto
#hasta el plano de la lenteMO
matriz_propagacionPrimerTramo = matriz.propagacion_MedioHomogeneo(distancia_Objeto)

#Se calcula la matriz asociada a la interacción con la lenteMO
matriz_lenteMO = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focalMO)



""" Calculando la matriz del sistema asociada al PRIMER TRAMO """

#NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
#en el arreglo.

#Se define la lista de matrices
lista_matricesPrimerTramoInvertida = [matriz_lenteMO,matriz_propagacionPrimerTramo]

#Se calcula la matriz del sistema teniendo en cuenta la lista de matrices definida anteriormente
matriz_SistemaPrimerTramo = matriz.matriz_Sistema(lista_matricesPrimerTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del PRIMER TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del sistema
camino_opticoCentralPrimerTramo = matriz.camino_Optico(lista_matricesPrimerTramoInvertida)

' ------ FIN SECCIÓN DE CÁLCULO MATRICIAL PRIMER TRAMO ------ '


' ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL SEGUNDO TRAMO ------ '

""" Se calculan las matrices necesarias para estudiar el SEGUNDO TRAMO del arreglo difractivo
    PUPILA (Determinada por la apertura del objetivo de microscopio en consideración) --> 
    *propagación(distancia imagen)* --> IMAGEN """

#Se calcula la matriz asociada al proceso de propagación desde el plano de la pupila 
# hasta el plano de la imagen
matriz_propagacionSegundoTramo = matriz.propagacion_MedioHomogeneo(distancia_Imagen)



""" Calculando la matriz del sistema asociada al SEGUNDO TRAMO 
NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
en el arreglo."""

#Creando lista de matrices que describe el arreglo del SEGUNDO TRAMO
lista_matricesSegundoTramoInvertida = [matriz_propagacionSegundoTramo]

#Se calcula la matriz del sistema
matriz_SistemaSegundoTramo = matriz.matriz_Sistema(lista_matricesSegundoTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del SEGUNDO TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del Segundo tramo
camino_opticoCentralSegundoTramo = matriz.camino_Optico(lista_matricesSegundoTramoInvertida)

' ------ FIN SECCIÓN DE CÁLCULO MATRICIAL SEGUNDO TRAMO ------ '


' ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL TRAMO FORMACIÓN HOLOGRAMA ------ '

""" Se calculan las matrices necesarias para estudiar el TRAMO DE FORMACIÓN DEL HOLOGRAMA
* Solo hay presencia de una componente de propagación entre el plano IMAGEN y el plano HOLOGRAMA
(separados una distancia 'distancia_ImagenHolograma') """

#Se calcula la matriz asociada al proceso de propagación desde el plano Imagen hasta el 
#plano Holograma
matriz_propagacionImagenHolograma = matriz.propagacion_MedioHomogeneo(distancia_ImagenHolograma)



""" Calculando la matriz del sistema asociada a la FORMACIÓN DEL HOLOGRAMA
NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
en el arreglo."""

#Creando lista de matrices que describe el arreglo del SEGUNDO TRAMO
lista_matricesFormacionHologramaInvertida = [matriz_propagacionImagenHolograma]

#Se calcula la matriz del sistema
matriz_SistemaFormacionHolograma = matriz.matriz_Sistema(lista_matricesFormacionHologramaInvertida)


""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del SEGUNDO TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del Segundo tramo
camino_opticoCentralTramoFormacionHolograma = matriz.camino_Optico(lista_matricesFormacionHologramaInvertida)


' ------ FIN SECCIÓN DE CÁLCULO MATRICIAL TRAMO FORMACIÓN HOLOGRAMA ------ '


' ----- EMPIEZA SECCIÓN DE CONFIGURACIÓN DE MALLAS DE PUNTOS PARA CADA PLANO ------ '

""" Creando mallas de puntos asociada a plano del sensor """
#Para poder conocer las condicones de creación de las mallas de puntos se debe considerar
#que los planos de interés son: PLANO OBJETO --> PLANO ANTERIOR A PUPILA --> PLANO IMAGEN -->
# PLANO HOLOGRAMA (Sensor).
#Ahora bien, en este caso los datos asociados al sensor son aquellos que permitirán conocer
#las condiciones de todas las mallas de puntos, lo cual se logra usando las condiciones de
#producto espacio frecuencia en el caso específico de la transformada de Fresnel.
#Se calcula la malla de puntos asociada al plano de medición

#Creando malla de puntos plano del sensor/MEDICIÓN
xx_PlanoHolograma,yy_PlanoHolograma = mascaras.malla_Puntos(resolucion_anchoSensorInput,ancho_SensorInput,
                                            resolucion_altoSensorInput,alto_SensorInput)



""" Creando malla de puntos asociada a plano IMAGEN """

#Se llama función para determinar los deltas de muestreo del tramo SENSOR --> IMÁGEN
deltas_tramoImagenHolograma = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_anchoSensorInput,
                                                                                              ancho_SensorInput,
                                                                                              matriz_SistemaFormacionHolograma[0,1],
                                                                                              longitud_onda_input,
                                                                                              resolucion_altoSensorInput,
                                                                                              alto_SensorInput)

#Se calcula el ancho de la ventana del plano imagen
anchoX_VentanaPlanoImagen = resolucion_anchoSensorInput*deltas_tramoImagenHolograma["deltaPlanoEntrada_X"]
altoY_VentanaPlanoImagen = resolucion_altoSensorInput*deltas_tramoImagenHolograma["deltaPlanoEntrada_Y"]


#Creando malla de puntos plano del imagen
xx_PlanoImagen,yy_PlanoImagen = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoImagen,
                                            resolucion_altoSensorInput,altoY_VentanaPlanoImagen)



""" Creando malla de puntos asociada a plano anterior a la pupila """

#Se llama función para determinar los deltas de muestreo del tramo SENSOR --> PUPILA
deltas_tramoPupilaImagen = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_anchoSensorInput,
                                                                                             anchoX_VentanaPlanoImagen,
                                                                                             matriz_SistemaSegundoTramo[0,1],
                                                                                             longitud_onda_input,
                                                                                             resolucion_altoSensorInput,
                                                                                             altoY_VentanaPlanoImagen)

#Se calcula el ancho de la ventana del plano anterior a pupila
anchoX_VentanaPlanoPupila = resolucion_anchoSensorInput*deltas_tramoPupilaImagen["deltaPlanoEntrada_X"]
altoY_VentanaPlanoPupila = resolucion_altoSensorInput*deltas_tramoPupilaImagen["deltaPlanoEntrada_Y"]


#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoPupila, yy_PlanoPupila = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoPupila,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoPupila)



""" Creando malla de puntos asociada a plano de la máscara u objeto de entrada """

#Se llama función para determinar los deltas de muestreo del tramo PUPILA --> MÁSCARA
deltas_tramoObjetoPupila = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_anchoSensorInput,
                                                                                            anchoX_VentanaPlanoPupila,
                                                                                            matriz_SistemaPrimerTramo[0,1],
                                                                                            longitud_onda_input,
                                                                                            resolucion_altoSensorInput,
                                                                                            altoY_VentanaPlanoPupila)

#Se calcula el ancho de la ventana del plano de la máscara u objeto de entrada
anchoX_VentanaPlanoObjeto = resolucion_anchoSensorInput*deltas_tramoObjetoPupila["deltaPlanoEntrada_X"]
altoY_VentanaPlanoObjeto = resolucion_altoSensorInput*deltas_tramoObjetoPupila["deltaPlanoEntrada_Y"]


#Se calcula la malla de puntos asociada al plano del objeto
xx_PlanoObjeto, yy_PlanoObjeto = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoPupila,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoPupila)


' ----- FIN SECCIÓN DE CONFIGURACIÓN DE MALLAS DE PUNTOS PARA CADA PLANO ------ '


' ------ EMPIEZA SECCIÓN DE CÁLCULO RESULTADO DIFRACTIVO DE CADA TRAMO ------ '

""" Creando máscara de transmitancia asociada al objeto de estudio en el arreglo """

# Cargar la imagen PNG como máscara de transmitancia
ruta_imagen_png = "/home/labravo/Downloads/USAF_T-20.jpg"  # Especifica la ruta de imagen
objeto = function.cargar_imagen_png(ruta_imagen_png, resolucion_anchoSensorInput,resolucion_altoSensorInput)



""" Se calcula el resultado del proceso difractivo del PRIMER TRAMO"""

#Se calcula el campo de salida/en plano de la pupila --> Campo resultante del primer tramo 
campo_PlanoPupila = matriz.matriz_ABCD_Difraccion_Sensor(camino_opticoCentralPrimerTramo,objeto,
                                                 matriz_SistemaPrimerTramo[0,0],
                                                 matriz_SistemaPrimerTramo[0,1],
                                                 matriz_SistemaPrimerTramo[1,1],
                                                 xx_PlanoObjeto,yy_PlanoObjeto,
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



""" Se calcula el resultado del proceso difractivo del SEGUNDO TRAMO """

#Se calcula el campo de salida/en plano de medición --> Campo resultante de la difracción 
campo_PlanoImagen = matriz.matriz_ABCD_Difraccion_Sensor_Shift(camino_opticoCentralSegundoTramo,
                                                    campo_entradaSegundoTramo,
                                                    matriz_SistemaSegundoTramo[0,0],
                                                    matriz_SistemaSegundoTramo[0,1],
                                                    matriz_SistemaSegundoTramo[1,1],xx_PlanoPupila,
                                                    yy_PlanoPupila,xx_PlanoImagen,yy_PlanoImagen,
                                                    numero_onda_input,deltas_tramoImagenHolograma)



#Se calcula la amplitud del campo de salida
amplitud_campoPlanoImagen = np.abs(campo_PlanoImagen)


#Se calcula la intensidad del campo de salida
intensidad_campoPlanoImagen = amplitud_campoPlanoImagen**2


' ------ FIN SECCIÓN DE CÁLCULO RESULTADO DIFRACTIVO DE CADA TRAMO ------ '


' ------ EMPIEZA SECCIÓN DE GRAFICACIÓN ------ '

""" Graficando máscara de transmitancia """

graph.graficar_transmitancia(objeto,anchoX_VentanaPlanoObjeto,altoY_VentanaPlanoObjeto,
                             "Muestra")


""" Graficando la intensidad del campo de salida del arreglo """

graph.graficar_intensidad(intensidad_campoPlanoImagen,ancho_SensorInput,alto_SensorInput,
                          "Intensidad recibida en el sensor",1,0.1)

' ------ FIN SECCIÓN DE GRAFICACIÓN ------ '


' ------ EMPIEZA SECCIÓN DE INTERFERENCIA ------ '

""" Definición del haz de referencia --> ONDA PLANA """

# Definición de onda plana INVERSA al haz de referencia  
onda_PlanaReferencia = np.exp(-1j*numero_onda_input*np.cos(np.radians(angulo_HazReferencia))
                              *np.sqrt((xx_PlanoImagen**2) + (yy_PlanoImagen**2)))

 

""" Se genera interferencia entre el haz de referencia (ONDA PLANA REF) y el haz objeto (CAMPO EN PLANO IMAGEN)"""

#Interferencia entre haz de referencia y haz objeto 
campo_interferenciaObjetoReferencia = onda_PlanaReferencia + campo_PlanoImagen

' ------ FIN SECCIÓN DE INTERFERENCIA ------ '


' ------ EMPIEZA SECCIÓN DE DIFRACCIÓN PARA FORMACIÓN HOLOGRAMA ------ '


""" Se calcula el resultado del proceso difractivo del TRAMO FORMACIÓN HOLOGRAMA """

#Se calcula el campo de salida/en plano de medición --> Campo resultante de la difracción 
campo_PlanoHolograma = matriz.matriz_ABCD_Difraccion_Sensor(camino_opticoCentralTramoFormacionHolograma,
                                                            campo_interferenciaObjetoReferencia,
                                                            matriz_SistemaFormacionHolograma[0,0],
                                                            matriz_SistemaFormacionHolograma[0,1],
                                                            matriz_SistemaFormacionHolograma[1,1],
                                                            xx_PlanoImagen,yy_PlanoImagen,xx_PlanoHolograma,
                                                            yy_PlanoHolograma,numero_onda_input,
                                                            deltas_Sensor)

#Se calcula la amplitud del campo de salida
amplitud_campoPlanoHolograma = np.abs(campo_PlanoHolograma)


#Se calcula la intensidad del campo de salida
intensidad_campoPlanoHolograma= amplitud_campoPlanoHolograma**2
 

' ------ FIN SECCIÓN DE DIFRACCIÓN PARA FORMACIÓN HOLOGRAMA ------ '


' ------ EMPIEZA SECCIÓN DE GRAFICACIÓN ------ '

""" Graficando la amplitud del campo asociada al HOLOGRAMA"""
graph.graficar_amplitud(amplitud_campoPlanoHolograma,ancho_SensorInput,alto_SensorInput,
                        "Amplitud del campo Holograma",1,0.3)

""" Graficando la intensidad del HOLOGRAMA """

graph.graficar_intensidad(intensidad_campoPlanoHolograma,ancho_SensorInput,alto_SensorInput,
                          "Intensidad recibida en el sensor",1,0.1)

' ------ FIN SECCIÓN DE GRAFICACIÓN ------ '