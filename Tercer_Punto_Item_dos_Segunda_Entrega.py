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
import matplotlib.pyplot as plt
import Funciones_importantes as function


""" Definiendo parámetros de sensor  DMM 37UX250-ML """

#Numero de muestras/puntos distribuidos en el ancho del sensor
resolucion_anchoSensorInput = 2448 

#Numero de muestras/puntos distribuidos en el alto del sensor
resolucion_altoSensorInput = 2048

#Tamaño del pixel del sensor
tamaño_PixelSensorInput = 3.45E-6

#Se crea una lista a la cual se le asigna los valores de los delta de muestreo asociados al sensor
deltas_Sensor = [tamaño_PixelSensorInput,tamaño_PixelSensorInput]

#Usando los datos anteriores se calcula el ancho y alto físicos del sensor
ancho_SensorInput = resolucion_anchoSensorInput*tamaño_PixelSensorInput
alto_SensorInput = resolucion_altoSensorInput*tamaño_PixelSensorInput



""" Definiendo parámetros de máscara difractiva """

#Radio del círculo asociado a la máscara circular
radio = 0.0005 

#Lados asociados a la máscara rectangular 
lado_Rectangulo01 = 5E-5
lado_Rectangulo02 = 5E-5

# Se define el centro u origen para la configuración de la máscara 
centro = None  



""" Definiendo parámetro para el tamaño de la pupila """

radio_pupilaInput = 0.035 #Se define variable asociada al radio de la abertura circular que representará
                           # el diafragma.



""" Definición de distancias del arreglo """

distancia_focal01 = 0.01  #Distancia focal asociada a la lente 01

distancia_focal02 = 0.565 #Distancia focal asociada a la lente 02--> CALCULADA PARA LOGRAR 
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



""" Creando máscara de transmitancia asociada al objeto de estudio en el arreglo """

#Creación de una máscara circular de transmitancia
#mascara = mascaras.funcion_Circulo(radio, centro, xx_PlanoMascara, yy_PlanoMascara)

#Creación de una máscara rectangular de transmitancia
#mascara = mascaras.funcion_Rectangulo(lado_Rectangulo01,lado_Rectangulo02,centro,xx_PlanoMascara,yy_PlanoMascara)

#Creación de una máscara con un corazón de transmitancia
#mascara = mascaras.funcion_Corazon(centro,xx_PlanoMascara,yy_PlanoMascara,radio)

# Cargar la imagen PNG como máscara de transmitancia
ruta_imagen_png = "/home/labravo/Downloads/Ruido_E03.png"  # Especifica la ruta de tu imagen
#mascara = function.cargar_imagen_png(ruta_imagen_png, resolucion_anchoSensorInput,resolucion_altoSensorInput)

from scipy.ndimage import zoom

# Cargar el archivo CSV
ruta_csv = "/home/labravo/Downloads/MuestraBio_E03.csv"  # Reemplaza con la ruta de tu archivo

# Leer el archivo y reemplazar 'i' con 'j'
with open("/home/labravo/Downloads/MuestraBio_E03.csv", "r") as file:
    contenido = file.read().replace("i", "j")

# Guardar el archivo actualizado
with open("/home/labravo/Downloads/MuestraBio_E03.csv", "w") as file:
    file.write(contenido)

datos_csv = np.genfromtxt(ruta_csv, delimiter=',', dtype=complex)
resolucion_x_csv, resolucion_y_csv = datos_csv.shape

# Factores de escala para ajustar a las dimensiones del plano de la máscara
factor_x = resolucion_anchoSensorInput / resolucion_x_csv
factor_y = resolucion_altoSensorInput / resolucion_y_csv

# Interpolar el CSV para ajustar su tamaño
datos_csv_ajustados = zoom(datos_csv, (factor_y, factor_x), order=1)  #


# Asignar la máscara ajustada
mascara = datos_csv_ajustados

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



""" Graficando máscara de transmitancia asignada a abertura circular"""

# plt.imshow(mascara, extent=[-anchoX_VentanaPlanoMascara/2, anchoX_VentanaPlanoMascara/2,
#                              -altoY_VentanaPlanoMascara/2, altoY_VentanaPlanoMascara/2], 
#                              cmap='gray')
# plt.title("Máscara")
# plt.colorbar(label="Transmitancia")
# plt.xlabel("X (m)")
# plt.ylabel("Y (m)")
# plt.show()

# Visualizar la máscara ajustada
plt.imshow(
    (np.abs(datos_csv_ajustados)**2),
    extent=[
        -anchoX_VentanaPlanoMascara / 2, anchoX_VentanaPlanoMascara / 2,
        -altoY_VentanaPlanoMascara / 2, altoY_VentanaPlanoMascara / 2
    ],
    cmap="gray",
    vmin = 1.8*(np.min((np.abs(datos_csv_ajustados))**2)))

plt.title("Máscara ajustada desde CSV")
plt.colorbar(label="Transmitancia")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()

""" Graficando intensidad del campo que entra al SEGUNDO TRAMO del arreglo"""

plt.imshow(intensidad_campoEntradaSegundoTramo, 
           extent=[-anchoX_VentanaPlanoPupila/2, anchoX_VentanaPlanoPupila/2,
                -altoY_VentanaPlanoPupila/2, altoY_VentanaPlanoPupila/2], 
            cmap='gray',
            vmax= 0.001*np.max(intensidad_campoEntradaSegundoTramo))
plt.title("Campo PUPILA")
plt.colorbar(label="Amplitud")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()


""" Graficando la intensidad del campo de salida del arreglo """

plt.imshow(intensidad_campoPlanoMedicion, extent=[-ancho_SensorInput/2, 
                                                  ancho_SensorInput/2, 
                                                  -alto_SensorInput/2, 
                                                  alto_SensorInput/2], 
           cmap='gray',
           vmax = 1*(np.max(intensidad_campoPlanoMedicion)),
           vmin = 1.8*(np.min(intensidad_campoPlanoMedicion)))
plt.title("Intensidad")
plt.colorbar(label="Intensidad")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()



