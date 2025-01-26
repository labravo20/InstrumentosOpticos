""" DESCRIPCIÓN...

En este documento se trabajará sobre un arreglo determinado por la siguiente estructura:
OBJETO --> *propagación(distancia focal 01)* --> LENTE01 --> *propagación(distancia focal 01)* -->
PUPILA --> *propagación(distancia arbitraria d)* --> LENTE02 --> *propagación(distancia focal 02)*

--> IMAGEN

Debido a la condición de planos conjugados se debe sub dividir el proceso en dos tramos:
TRAMO 01: OBJETO --> *propagación(distancia focal 01)* --> LENTE01 --> *propagación(distancia focal 01)* -->
PUPILA

TRAMO 02: *propagación(distancia arbitraria d)* --> LENTE02 --> *propagación(distancia focal 02)*
--> IMAGEN

"""

""" Se genera un print de mensaje inicial para verificar correcto funcionamiento del entorno """
print("Inicializando entorno de programación tercer punto SEGUNDA ENTREGA...")



""" Anotaciones importantes magnificación de imágenes """
#NOTA SOBRE CÁLCULO DISTANCIA IMAGEN Y OBJETO: (1/f) = (1/I) + (1/O)

#NOTA SOBRE MAGNIFICACIÓN (para calcular tamaño aproximado de imagen en función de tamaño objeto):

### M = I/O --> Tamaño de imagen será: (TAMAÑO OBJETO)*M



""" Importando librerias y documentos necesarios """

import Mascaras_Transmitancia as mascaras
import Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import Funciones_Graficacion as graficar
import Funciones_importantes as function


""" Definiendo parámetros de sensor  DMM 37UX250-ML """

#Numero de muestras/puntos distribuidos en el ancho del sensor
resolucion_anchoSensorInput = 2448 

#Numero de muestras/puntos distribuidos en el alto del sensor
resolucion_altoSensorInput = 2048

#NOTA IMPORTANTE --> En este caso para lograr una distribución de arreglo cuadrada (dada la 
# condición de caracterización donde se solicita un arreglo de entrada asociado a 125 um x 125 um)
# se usará unicamente el dato de menor distribución, es decir, se configurarán arreglos de 2048 x 2048, es decir:
resolucion_anchoSensorInput = 2048 

#Tamaño del pixel del sensor
tamaño_PixelSensorInput = 3.45E-6

#Se crea una lista a la cual se le asigna los valores de los delta de muestreo asociados al sensor
deltas_Sensor = [tamaño_PixelSensorInput,tamaño_PixelSensorInput]

#Usando los datos anteriores se calcula el ancho y alto físicos del sensor
ancho_SensorInput = resolucion_anchoSensorInput*tamaño_PixelSensorInput
alto_SensorInput = resolucion_altoSensorInput*tamaño_PixelSensorInput



""" Definiendo parámetros de máscara para procesamiento de imagen  """
#Radio del círculo asociado a la máscara circular
radio = 10E-5

#Se define radio interno del anillo
radio_internoAnillo = 1E-4

#Se define radio externo del anillo
radio_externoAnillo = 1.5E-4

# Se define el centro u origen para la configuración de la máscara 
centro = None  


""" Definiendo parámetro para el tamaño de la pupila """

radio_pupilaInput = 0.0035 #Se define variable asociada al radio de la abertura circular que representará
                           # el diafragma.



""" Definición de distancias del arreglo """

distancia_focal01 = 0.01  # Distancia focal asociada a la lente 01

#distancia_focal02 = 0.565 # Distancia focal asociada a la lente 02--> CALCULADA PARA LOGRAR 
#ANCHO DE VENTANA MÁSCARA DE 125E-6 m

distancia_focal02 = 0.1983 # Distancia focal asociada a la lente 02--> CALCULADA PARA LOGRAR 
#ANCHO DE VENTANA MÁSCARA DE 125E-6 m

distancia_propagacionAribitraria = 0.01 #Se define una distancia de propagación arbitraria 



""" Definiendo parámetros de fuente """

#Definiendo la longitud de onda asociada a la fuente en consideración
longitud_onda_input = 533E-9 #UNIDADES: m

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
                                                                                            matriz_SistemaPrimerTramo[0,1],
                                                                                            longitud_onda_input,
                                                                                            resolucion_altoSensorInput,
                                                                                            altoY_VentanaPlanoPupila)

#Se calcula el ancho de la ventana del plano de la máscara u objeto de entrada
anchoX_VentanaPlanoMascara = resolucion_anchoSensorInput*deltas_tramoMascaraPupila[0]
altoY_VentanaPlanoMascara = resolucion_altoSensorInput*deltas_tramoMascaraPupila[1]

#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoMascara, yy_PlanoMascara = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoMascara,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoMascara)



""" Creando máscara de transmitancia asociada al objeto de estudio en el arreglo """

# Cargar el archivo CSV
ruta_csv = "/home/labravo/Downloads/MuestraBio_E03.csv" # Reemplaza con la ruta de tu archivo
#mascara = function.cargar_documento_csv(ruta_csv,resolucion_anchoSensorInput,resolucion_altoSensorInput)
mascara = function.cargar_documento_csv_OPTION02(ruta_csv)

#Se calcula la intensidad asociada al campo de entrada
intensidad_mascara = np.abs(mascara)**2

#Se calcula la fase asociada al campo de entrada
mascara_fase = np.angle(mascara)

# graficar.graficar_fase(mascara_fase,anchoX_VentanaPlanoMascara,altoY_VentanaPlanoMascara,
#                        "Fase del campo de entrada")


""" ------ EMPIEZA SECCIÓN DE CÁLCULO RESULTADO DIFRACTIVO DE CADA TRAMO ------ """

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
pupila = mascaras.funcion_Circulo(radio_pupilaInput, centro, xx_PlanoPupila, yy_PlanoPupila)



""" Se calcula el campo de entrada para el SEGUNDO TRAMO del arreglo 
NOTA: El campo que llega a la pupila e interactua con la máscara asociada a la pupila
será el campo de entrada para el segundo tramo."""

#Calculando el campo que sale después de la interacción con la pupila
campo_salidaPupila = campo_PlanoPupila*pupila

#Calculando la intensidad del campo de entrada al SEGUNDO TRAMO del arreglo
intensidad_campoSalidaPupila = np.abs(campo_salidaPupila)**2



""" Se crea e implementa una 'máscara' adicional para procesar la imagen y eliminar el aporte 
proveniente de la fuente monocromática con la cual se eliminó la muestra """

# Crear la máscara del anillo de fase
# anillo_fase = mascaras.funcion_AnilloFase(radio_interno_anillo, radio_externo_anillo, fase_anillo, xx_PlanoPupila, yy_PlanoPupila)

# #Calculando el campo de entrada al SEGUNDO TRAMO del arreglo
# campo_entradaSegundoTramo = campo_salidaPupila*anillo_fase

# #Calculando la intensidad del campo de entrada al SEGUNDO TRAMO del arreglo
# intensidad_campoEntradaSegundoTramo = np.abs(campo_entradaSegundoTramo)**2

#Creación de una máscara rectangular de transmitancia para eliminar el aporte proviniente de la fuente monocromática
mascara_procesamiento = mascaras.funcion_CirculoInvertidoGaussian(radio,centro,xx_PlanoPupila,yy_PlanoPupila)

# #Se calcula el campo tras la interacción con la máscara de procesamiento
campo_entradaSegundoTramo = campo_salidaPupila*mascara_procesamiento

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

#Se calcula la fase asociada al campo en el plano de medición
campo_PlanoMedicionFase = np.angle(campo_PlanoMedicion)

#Se calcula la amplitud del campo de salida
amplitud_campoPlanoMedicion = np.abs(campo_PlanoMedicion)

#Se calcula la intensidad del campo de salida
intensidad_campoPlanoMedicion = amplitud_campoPlanoMedicion**2


#Versión NEGATIVA DE LA IMAGEN PLANO MEDICIÓN 

""" ------ EMPIEZA SECCIÓN DE GRAFICACIÓN ----- """

graficar.graficar_fase(np.angle(mascara),anchoX_VentanaPlanoMascara,altoY_VentanaPlanoMascara,"Fase muestra")

""" --- USANDO ANILLOS DE FASE ------ """

#Creación de una máscara rectangular de transmitancia para eliminar el aporte proviniente de la fuente monocromática
mascara_procesamiento = mascaras.funcion_AnilloFase(radio_internoAnillo,radio_externoAnillo,
                                                    -(np.pi)/2,xx_PlanoPupila,yy_PlanoPupila)

#Creación de máscara para agregar término de transparencia
mascara_transparencia = mascaras.funcion_Anillo(radio_internoAnillo,radio_externoAnillo,xx_PlanoPupila,
                                                yy_PlanoPupila,0.3)


#Se calcula el campo tras la interacción con la máscara de procesamiento
campo_entradaSegundoTramo = campo_salidaPupila*mascara_procesamiento*mascara_transparencia

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

#Se calcula la fase asociada al campo en el plano de medición
campo_PlanoMedicionFase = np.angle(campo_PlanoMedicion)

#Se calcula la amplitud del campo de salida
amplitud_campoPlanoMedicion = np.abs(campo_PlanoMedicion)

#Se calcula la intensidad del campo de salida
intensidad_campoPlanoMedicion = amplitud_campoPlanoMedicion**2



import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np


# Función para actualizar las máscaras, la intensidad y la transformada de Fourier
def actualizar_graficos(val):
    radioInterno = slider_radio1.val
    radioExterno = slider_radio2.val

    # Recalcular la máscara de anillos
    mascara_actualizada = mascaras.funcion_AnilloFase(radioInterno, radioExterno, -(np.pi)/2, 
                                                      anchoX_VentanaPlanoPupila, altoY_VentanaPlanoPupila)
    mascara_transparenciaActualizada = mascaras.funcion_Anillo(radioInterno, radioExterno, xx_PlanoPupila,
                                                                yy_PlanoPupila, 0.3)

    # Recalcular el campo de entrada al segundo tramo
    nuevo_campo_entradaSegundoTramo = campo_salidaPupila * mascara_actualizada * mascara_transparenciaActualizada

    # Calcular el campo de salida/en plano de medición
    nuevo_campo_PlanoMedicion = matriz.matriz_ABCD_Difraccion_Sensor_Shift(
        camino_opticoCentralSegundoTramo,
        nuevo_campo_entradaSegundoTramo,
        matriz_SistemaSegundoTramo[0, 0],
        matriz_SistemaSegundoTramo[0, 1],
        matriz_SistemaSegundoTramo[1, 1],
        xx_PlanoPupila, yy_PlanoPupila, xx_PlanoMedicion, yy_PlanoMedicion,
        numero_onda_input, deltas_Sensor
    )

    # Calcular la intensidad del campo de salida
    nuevo_intensidad_campoPlanoMedicion = (np.abs(nuevo_campo_PlanoMedicion))**2

    # Calcular la intensidad de la transformada de Fourier
    intensidad_transformada = (np.abs(nuevo_campo_entradaSegundoTramo))**2

    # Actualizar ambos gráficos
    img_salida.set_data(nuevo_intensidad_campoPlanoMedicion)
    img_transformada.set_data(intensidad_transformada)

    fig.canvas.draw_idle()


# Crear la figura con dos subgráficos (1 fila, 2 columnas)
fig, (ax_salida, ax_transformada) = plt.subplots(1, 2, figsize=(12, 6))
plt.subplots_adjust(left=0.1, bottom=0.25, wspace=0.4)  # Ajustar espaciado horizontal

# Primer gráfico: Intensidad del campo de salida
intensidad_inicial = intensidad_campoPlanoMedicion
img_salida = ax_salida.imshow(intensidad_inicial, extent=[-ancho_SensorInput / 2, ancho_SensorInput / 2,
                                                          -alto_SensorInput / 2, alto_SensorInput / 2],
                               cmap='gray', origin='lower',
                               vmax = 0.8*np.max(intensidad_inicial))
plt.colorbar(img_salida, ax=ax_salida, label="Intensidad")
ax_salida.set_title("Intensidad del campo de salida")
ax_salida.set_xlabel("x [m]")
ax_salida.set_ylabel("y [m]")

# Segundo gráfico: Transformada de Fourier
intensidad_transformada_inicial = np.abs(campo_entradaSegundoTramo)**2
img_transformada = ax_transformada.imshow(intensidad_transformada_inicial, extent=[-anchoX_VentanaPlanoPupila / 2, anchoX_VentanaPlanoPupila / 2,
                                                                                   -altoY_VentanaPlanoPupila / 2, altoY_VentanaPlanoPupila / 2],
                                          cmap='gray', origin='lower',
                                           vmax = 0.001*np.max(intensidad_transformada_inicial))
plt.colorbar(img_transformada, ax=ax_transformada, label="Intensidad")
ax_transformada.set_title("Intensidad de la Transformada de Fourier")
ax_transformada.set_xlabel("x [m]")
ax_transformada.set_ylabel("y [m]")

# Crear los sliders para modificar los radios
ax_radio1 = plt.axes([0.1, 0.15, 0.35, 0.03])  # Slider para el primer radio
slider_radio1 = Slider(ax_radio1, 'Radio Interno', 1E-9, 3E-6, valinit=radio_internoAnillo)

ax_radio2 = plt.axes([0.55, 0.15, 0.35, 0.03])  # Slider para el segundo radio
slider_radio2 = Slider(ax_radio2, 'Radio Externo', 1E-7, 3E-4, valinit=radio_externoAnillo)

# Conectar los sliders con la función de actualización
slider_radio1.on_changed(actualizar_graficos)
slider_radio2.on_changed(actualizar_graficos)

# Mostrar la ventana con los gráficos
plt.show()


""" --- USANDO ANILLOS SIN FASE ------ """

# #Creación de máscara para agregar término de transparencia
# mascara_transparencia = mascaras.funcion_Anillo(radio_internoAnillo,radio_externoAnillo,xx_PlanoPupila,
#                                                 yy_PlanoPupila,0)


# #Se calcula el campo tras la interacción con la máscara de procesamiento
# campo_entradaSegundoTramo = campo_salidaPupila*(1-mascara_transparencia)

# intensidad_campoEntradaSegundoTramo = np.abs(campo_entradaSegundoTramo)**2



# """ Se calcula el resultado del proceso difractivo del SEGUNDO TRAMO """

# #Se calcula el campo de salida/en plano de medición --> Campo resultante de la difracción 
# campo_PlanoMedicion = matriz.matriz_ABCD_Difraccion_Sensor_Shift(camino_opticoCentralSegundoTramo,
#                                                     campo_entradaSegundoTramo,
#                                                     matriz_SistemaSegundoTramo[0,0],
#                                                     matriz_SistemaSegundoTramo[0,1],
#                                                     matriz_SistemaSegundoTramo[1,1],xx_PlanoPupila,
#                                                     yy_PlanoPupila,xx_PlanoMedicion,yy_PlanoMedicion,
#                                                     numero_onda_input,deltas_Sensor)

# #Se calcula la fase asociada al campo en el plano de medición
# campo_PlanoMedicionFase = np.angle(campo_PlanoMedicion)

# #Se calcula la amplitud del campo de salida
# amplitud_campoPlanoMedicion = np.abs(campo_PlanoMedicion)

# #Se calcula la intensidad del campo de salida
# intensidad_campoPlanoMedicion = amplitud_campoPlanoMedicion**2



# import matplotlib.pyplot as plt
# from matplotlib.widgets import Slider
# import numpy as np


# # Función para actualizar las máscaras, la intensidad y la transformada de Fourier
# def actualizar_graficos(val):
#     radioInterno = slider_radio1.val
#     radioExterno = slider_radio2.val

#     # Recalcular la máscara de anillos
#     mascara_transparenciaActualizada = mascaras.funcion_Anillo(radioInterno, radioExterno, xx_PlanoPupila,
#                                                                 yy_PlanoPupila, 0)

#     # Recalcular el campo de entrada al segundo tramo
#     nuevo_campo_entradaSegundoTramo = campo_salidaPupila  * (1-mascara_transparenciaActualizada)

#     # Calcular el campo de salida/en plano de medición
#     nuevo_campo_PlanoMedicion = matriz.matriz_ABCD_Difraccion_Sensor_Shift(
#         camino_opticoCentralSegundoTramo,
#         nuevo_campo_entradaSegundoTramo,
#         matriz_SistemaSegundoTramo[0, 0],
#         matriz_SistemaSegundoTramo[0, 1],
#         matriz_SistemaSegundoTramo[1, 1],
#         xx_PlanoPupila, yy_PlanoPupila, xx_PlanoMedicion, yy_PlanoMedicion,
#         numero_onda_input, deltas_Sensor
#     )

#     # Calcular la intensidad del campo de salida
#     nuevo_intensidad_campoPlanoMedicion = (np.abs(nuevo_campo_PlanoMedicion))**2

#     # Calcular la intensidad de la transformada de Fourier
#     intensidad_transformada = (np.abs(nuevo_campo_entradaSegundoTramo))**2

#     # Actualizar ambos gráficos
#     img_salida.set_data(nuevo_intensidad_campoPlanoMedicion)
#     img_transformada.set_data(intensidad_transformada)

#     fig.canvas.draw_idle()


# # Crear la figura con dos subgráficos (1 fila, 2 columnas)
# fig, (ax_salida, ax_transformada) = plt.subplots(1, 2, figsize=(12, 6))
# plt.subplots_adjust(left=0.1, bottom=0.25, wspace=0.4)  # Ajustar espaciado horizontal

# # Primer gráfico: Intensidad del campo de salida
# intensidad_inicial = intensidad_campoPlanoMedicion
# img_salida = ax_salida.imshow(intensidad_inicial, extent=[-ancho_SensorInput / 2, ancho_SensorInput / 2,
#                                                           -alto_SensorInput / 2, alto_SensorInput / 2],
#                                cmap='gray', origin='lower',
#                                vmax = 0.8*np.max(intensidad_inicial))
# plt.colorbar(img_salida, ax=ax_salida, label="Intensidad")
# ax_salida.set_title("Intensidad del campo de salida")
# ax_salida.set_xlabel("x [m]")
# ax_salida.set_ylabel("y [m]")

# # Segundo gráfico: Transformada de Fourier
# intensidad_transformada_inicial = np.abs(campo_entradaSegundoTramo)**2
# img_transformada = ax_transformada.imshow(intensidad_transformada_inicial, extent=[-anchoX_VentanaPlanoPupila / 2, anchoX_VentanaPlanoPupila / 2,
#                                                                                    -altoY_VentanaPlanoPupila / 2, altoY_VentanaPlanoPupila / 2],
#                                           cmap='gray', origin='lower',
#                                            vmax = 0.01*np.max(intensidad_transformada_inicial))
# plt.colorbar(img_transformada, ax=ax_transformada, label="Intensidad")
# ax_transformada.set_title("Intensidad de la Transformada de Fourier")
# ax_transformada.set_xlabel("x [m]")
# ax_transformada.set_ylabel("y [m]")

# # Crear los sliders para modificar los radios
# ax_radio1 = plt.axes([0.1, 0.15, 0.35, 0.03])  # Slider para el primer radio
# slider_radio1 = Slider(ax_radio1, 'Radio Interno', 1E-9, 3E-2, valinit=radio_internoAnillo)

# ax_radio2 = plt.axes([0.55, 0.15, 0.35, 0.03])  # Slider para el segundo radio
# slider_radio2 = Slider(ax_radio2, 'Radio Externo', 1E-4, 3E-2, valinit=radio_externoAnillo)

# # Conectar los sliders con la función de actualización
# slider_radio1.on_changed(actualizar_graficos)
# slider_radio2.on_changed(actualizar_graficos)

# # Mostrar la ventana con los gráficos
# plt.show()