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
print("Inicializando entorno de programación segundo punto TECERA ENTREGA...")



""" Importando librerias y documentos necesarios """

import LIBRERIA_Mascaras_Transmitancia as mascaras
import LIBRERIA_Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import matplotlib.pyplot as plt
import LIBRERIA_Funciones_importantes as function
import LIBRERIA_Funciones_Graficacion as graph


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
                                                                                            matriz_SistemaPrimerTramo[0,1],
                                                                                            longitud_onda_input,
                                                                                            resolucion_altoSensorInput,
                                                                                            altoY_VentanaPlanoPupila)
print(matriz_SistemaSegundoTramo[0,1])
#Se calcula el ancho de la ventana del plano de la máscara u objeto de entrada
anchoX_VentanaPlanoMascara = resolucion_anchoSensorInput*deltas_tramoMascaraPupila[0]
altoY_VentanaPlanoMascara = resolucion_altoSensorInput*deltas_tramoMascaraPupila[1]

#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoMascara, yy_PlanoMascara = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoPupila,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoPupila)



""" Creando máscara de transmitancia asociada al objeto de estudio en el arreglo --> HOLOGRAMA """

# Cargar la imagen PNG como máscara de transmitancia
ruta_imagen_png = "/home/labravo/Downloads/Hologram.tiff"  # Especifica la ruta de tu imagen
mascara = function.cargar_imagen_png(ruta_imagen_png, resolucion_anchoSensorInput,resolucion_altoSensorInput)




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
pupila = mascaras.funcion_Circulo(radio_pupilaInput, None, xx_PlanoPupila, yy_PlanoPupila)

#Se crea máscara para filtrar el espectro de Fourier y obtener imágen real
#mascara_Filtrado = mascaras.funcion_Circulo(0.0005,[-0.0027,0.002],xx_PlanoMascara,yy_PlanoMascara)
#mascara_Filtrado = mascaras.funcion_Circulo(0.001,[-0.0026,0.002],xx_PlanoMascara,yy_PlanoMascara)
mascara_Filtrado = mascaras.funcion_Rectangulo(0.0015,0.0015,[-0.0031,0.002],xx_PlanoMascara,yy_PlanoMascara)

""" Se calcula el campo de entrada para el SEGUNDO TRAMO del arreglo 
NOTA: El campo que llega a la pupila e interactua con la máscara asociada a la pupila
será el campo de entrada para el segundo tramo."""

#Calculando el campo de entrada al SEGUNDO TRAMO del arreglo antes de filtrar
campo_entradaSegundoTramo_SinFiltro = campo_PlanoPupila*pupila

#Calculando la intensidad del campo de entrada al SEGUNDO TRAMO del arreglo
intensidad_campoEntradaSegundoTramo_SinFiltro = np.abs(campo_entradaSegundoTramo_SinFiltro)**2

#Calculando el campo de entrada al SEGUNDO TRAMO del arreglo
campo_entradaSegundoTramo = campo_PlanoPupila*pupila*mascara_Filtrado

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
graph.graficar_transmitancia(mascara,anchoX_VentanaPlanoMascara,altoY_VentanaPlanoMascara,"Holograma")


"""Graficando amplitud del campo que en plano de Fourier"""
graph.graficar_amplitud(np.abs(campo_PlanoPupila),anchoX_VentanaPlanoPupila,
                              altoY_VentanaPlanoPupila,"Transformada de Fourier del objeto",1,0.001)

""" Graficando intensidad del campo que entra al SEGUNDO TRAMO del arreglo"""
graph.graficar_intensidad(intensidad_campoEntradaSegundoTramo_SinFiltro,anchoX_VentanaPlanoPupila,
                              altoY_VentanaPlanoPupila,"Transformada de Fourier del objeto",1,0.00001)


""" Graficando máscara de filtrado"""
graph.graficar_transmitancia(mascara_Filtrado,anchoX_VentanaPlanoPupila,altoY_VentanaPlanoPupila,"Máscara de filtrado")


""" Graficando campo asociado a transformada de Fourier filtrado para encontrar imágen real """
graph.graficar_intensidad(intensidad_campoEntradaSegundoTramo,anchoX_VentanaPlanoPupila,
                             altoY_VentanaPlanoPupila,"Transformada de Fourier del objeto",1,0.001)


""" Graficando la intensidad del campo de salida del arreglo """
graph.graficar_intensidad(intensidad_campoPlanoMedicion,ancho_SensorInput,alto_SensorInput,"Intensidad campo de salida")




""" ###############Guardando matriz resultante como archivo txt######################### """
#print(campo_PlanoMedicion) #Para corroborar la información de salida
#print(type(campo_PlanoMedicion))



# Guardar la parte real e imaginaria en archivos separados
# np.savetxt("matriz_real.txt", campo_PlanoMedicion.real, fmt="%.18e")
# np.savetxt("matriz_imaginaria.txt", campo_PlanoMedicion.imag, fmt="%.18e")

# print("Matrices guardadas correctamente en 'matriz_real.txt' y 'matriz_imaginaria.txt'.")

# # import os
# # print("Archivos guardados en:", os.path.abspath("matriz_real.txt"))

# # Cargar las partes real e imaginaria
# real_cargado = np.loadtxt("matriz_real.txt")
# imaginario_cargado = np.loadtxt("matriz_imaginaria.txt")

# # Reconstruir la matriz compleja
# matriz_recuperada = real_cargado + 1j * imaginario_cargado

#print("Matriz reconstruida:")
#print(matriz_recuperada)


import matplotlib.image as mpimg
from scipy.ndimage import zoom



################ Parametros discretización ##########



num_pixels_x = 2048
num_pixels_y = 2048

delta_x = 3.45E-6 # Tamaño de cada pixel (m) 
delta_y = 3.45E-6 # Tamaño de cada pixel (m) 

delta_fx = 1/(num_pixels_x * delta_x) #Tamaño de cada pixel en plano de Fourier (1/m)
delta_fy = 1/(num_pixels_y * delta_y) #Tamaño de cada pixel en plano de Fourier (1/m)


# Malla de coordenadas espaciales (plano de entrada) ---> Plano de salida sistema 4f
x=np.arange(-num_pixels_x//2,num_pixels_x//2) 
y=np.arange(-num_pixels_y//2,num_pixels_y//2)
x,y=x*delta_x,y*delta_y
X, Y = np.meshgrid(x, y)

# malla de coordenadas espectrales
f_x=np.arange(-num_pixels_x//2,num_pixels_x//2)
f_y=np.arange(-num_pixels_y//2,num_pixels_y//2)
f_x,f_y=f_x*delta_fx,f_y*delta_fy
F_X, F_Y = np.meshgrid(f_x, f_y)



#Paso 1: Recibimos el patron de difracción del objeto: es una dsitribución de intensidad }



'Tomemos la matriz que respresenta la imagen y luego vamos a sacar su raíz para hallar el campo óptico a la salida'

#Se define la matriz de puntos para el análisis
matriz_campo = campo_PlanoMedicion

#Verificación del tamaño de la matriz de estudio
print("\n Tamaño  matriz campo óptico de salida:",matriz_campo.shape)



# 1. Aplicamos FFT para hallar A[p,q,z]
'NOTA: nuestro espectro angular de salida sale con un shift a causa de la fft, por lo que antes de dividir la matriz'
' punto a punto debemos shiftear el resultado de la fft'

espectro_angular_salida = np.fft.fftshift(np.fft.fft2(matriz_campo)) 


for distancia_propagacionAribitraria in np.arange(0.01, 0.68, 0.0001):


    # 2. Dividimos(matrices) punto a punto por la función de propagación para hallar A[p,q,0]

    funcion_transferencia = np.exp((1j*distancia_propagacionAribitraria*numero_onda_input)*(np.sqrt(1-(longitud_onda_input**2)*(F_X**2 + F_Y**2))))


    #Graficando la fase de nuestra función de transferencia
    plt.imshow(np.angle(funcion_transferencia),cmap='gray') # Graficamos la fase para comprobar que 
                                                            # efectivamente obtenemos una fase circular 
                                                            # sin réplicas o distorsiones
    plt.title("Fase función de transferencia")
    plt.show()


    #Calculando en espectro angular en el plano de entrada al sistema
    espectro_angular_entrada = espectro_angular_salida / funcion_transferencia
    print("\n Matriz Espectro angular entrada : \n", espectro_angular_entrada)

    # 3. Aplicamos IFFT para hallar U[x,y,0]

    "NOTA: la ifft recibe una matriz sin shiftear (desordenada) por lo que debemos desordenar nuestra matriz A[p,q,0] para"
    "aplicarle ifft. Recordar que al final de la ifft no se hace shift pues este algoritmo automáticamente reordena"

    campo_optico_entrada = np.fft.ifft2(np.fft.fftshift(espectro_angular_entrada))
    print("\n Matriz campo óptico entrada: \n", campo_optico_entrada)


    # Graficar la magnitud campo optico de entrada 
    'NOTA: En este caso sólo nos interesa la amplitud, por lo que tomamos la magnitud'
    plt.figure(figsize=(6, 6))
    plt.imshow(np.abs(campo_optico_entrada), cmap='gray', extent=[x[0], x[-1], y[0], y[-1]]) #extend se utiliza para escalar el rango de x tomado (-51.2 -51.2) en la grafica a realizar- x[0] es el valor de mas a la iz mientras que x[-1] es el ultimo valor del intervalo x
    plt.title("campo optico de entrada")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.colorbar()
    plt.show()
