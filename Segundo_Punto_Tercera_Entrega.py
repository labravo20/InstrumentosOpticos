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
print("Inicializando entorno de programación Matrices ABCD transferencia de rayos...")


""" Anotaciones importantes magnificación de imágenes """
#NOTA SOBRE CÁLCULO DISTANCIA IMAGEN Y OBJETO: (1/f) = (1/I) + (1/O)

#NOTA SOBRE MAGNIFICACIÓN (para calcular tamaño aproximado de imagen en función de tamaño objeto):

### M = I/O --> Tamaño de imagen será: (TAMAÑO OBJETO)*M



""" Importando librerias y documentos necesarios """

import LIBRERIA_Mascaras_Transmitancia as mascaras
import LIBRERIA_Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import matplotlib.pyplot as plt
import LIBRERIA_Funciones_importantes as function
from PIL import Image
import LIBRERIA_Funciones_Graficacion as graph


""" Definiendo parámetros de máscara difractiva """

resolucion_Input = 2048  # Número de puntos en la malla --> Asociado a comparación con una referencia de cámara
longitud_ArregloInput = (3.45E-6)*resolucion_Input  #Tamaño físico del área de la ventana
radio = 0.05  # Radio del círculo 
centro = None  # El centro será el origen si es None



""" Definiendo parámetro para el tamaño del diafragma """

radio_diafragmaInput = 0.07 #Se define variable asociada al radio de la abertura circular que representará
                           # el diafragma.



""" Definición de distancias del arreglo """

distancia_focal = 0.01

distancia_focal02 = 0.01

distancia_propagacionAribitraria = 0.01

""" Definiendo parámetros de fuente """

#Definiendo la longitud de onda asociada a la fuente en consideración
longitud_onda_input = 632.8E-9 #UNIDADES: m

#Calculando el número de onda
numero_onda_input = (2*np.pi)/longitud_onda_input



""" Creando máscara de plano objeto  """

# Crear la malla de puntos
xx_mascara, yy_mascara = mascaras.malla_Puntos(resolucion_Input, longitud_ArregloInput)

# Crear la máscara circular
#mascara = mascaras.funcion_Circulo(radio, centro, xx_mascara, yy_mascara)
#mascara = mascaras.funcion_Rectangulo(radio,radio,centro,xx_mascara,yy_mascara)
#mascara = mascaras.funcion_Corazon(centro,xx_mascara,yy_mascara,radio)

# Cargar la imagen PNG como máscara de transmitancia
ruta_imagen_png = "/home/labravo/Downloads/Hologram.tiff"  # Especifica la ruta de tu imagen
mascara = function.cargar_imagen_png(ruta_imagen_png,resolucion_Input,resolucion_Input)



""" Se calculan las matrices necesarias para estudiar el PRIMER TRAMO del arreglo difractivo
    OBJETO --> *Propagación* --> LENTE  --> PUPILA """

#Se calcula la matriz asociada al proceso de propagación desde el plano objeto
#hasta el plano de la lente
matriz_propagacionPrimerTramo = matriz.propagacion_MedioHomogeneo(distancia_focal)

matriz_lente = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focal)

matriz_propagacion02 = matriz.propagacion_MedioHomogeneo(distancia_focal)



""" Calculando la matriz del sistema asociada al PRIMER TRAMO """

#NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
#en el arreglo.

#Se define la lista de matrices
lista_matricesPrimerTramoInvertida = [matriz_propagacion02,matriz_lente,matriz_propagacionPrimerTramo]

#En este caso la matriz del sistema es equivalente a la única matriz presente
matriz_SistemaPrimerTramo = matriz.matriz_Sistema(lista_matricesPrimerTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del PRIMER TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del sistema
camino_opticoCentralPrimerTramo = matriz.camino_Optico(lista_matricesPrimerTramoInvertida)



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



""" ########## ANÁLISIS DE MALLAS DE PUNTOS ########### """

""" Creando malla de puntos asociada a plano de la PUPILA """

#Se llama función para determinar los deltas de muestreo del tramo PUPILA --> MÁSCARA
deltas_tramoMascaraPupila = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_Input,
                                                                                            longitud_ArregloInput,
                                                                                            matriz_SistemaPrimerTramo[0,1],
                                                                                            longitud_onda_input,
                                                                                            resolucion_Input,
                                                                                            longitud_ArregloInput)

#Se calcula el ancho de la ventana del plano de la máscara u objeto de entrada
anchoX_VentanaPlanoPupila = resolucion_Input*deltas_tramoMascaraPupila[0]
altoY_VentanaPlanoPupila = resolucion_Input*deltas_tramoMascaraPupila[1]

#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoPupila, yy_PlanoPupila = mascaras.malla_Puntos(resolucion_Input,anchoX_VentanaPlanoPupila,
                                                       resolucion_Input,altoY_VentanaPlanoPupila)


""" Creando malla de puntos asociada a plano Medición """

#Se llama función para determinar los deltas de muestreo del tramo SENSOR --> PUPILA
deltas_tramoPupilaMedicion = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_Input,
                                                                                             anchoX_VentanaPlanoPupila,
                                                                                             matriz_SistemaSegundoTramo[0,1],
                                                                                             longitud_onda_input,
                                                                                             resolucion_Input,
                                                                                             altoY_VentanaPlanoPupila)

#Se calcula el ancho de la ventana del plano anterior a pupila
anchoX_VentanaPlanoMedicion = resolucion_Input*deltas_tramoPupilaMedicion[0]
altoY_VentanaPlanoMedicion = resolucion_Input*deltas_tramoPupilaMedicion[1]

#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoMedicion, yy_PlanoMedicion = mascaras.malla_Puntos(resolucion_Input,anchoX_VentanaPlanoPupila,
                                                       resolucion_Input,altoY_VentanaPlanoPupila)




#####################################################

""" Se calcula el resultado del proceso difractivo del PRIMER TRAMO"""

#Se calcula el campo de salida/en plano de la lente --> Campo resultante del primer tramo 
campo_PlanoPupila = matriz.matriz_ABCD_Difraccion(camino_opticoCentralPrimerTramo,mascara,
                                                 matriz_SistemaPrimerTramo[0,0],
                                                 matriz_SistemaPrimerTramo[0,1],
                                                 matriz_SistemaPrimerTramo[1,1],xx_mascara,yy_mascara,
                                                 xx_PlanoPupila,yy_PlanoPupila,numero_onda_input,
                                                 deltas_tramoMascaraPupila)




""" Se crea diafragma para delimitar dominio del arreglo asociado a lentes FINITAS
NOTA: El diafragma va a construirse a partir de la malla de puntos asociada al plano de la Lente."""

#Creación de máscara circular que representará el diafragma 
diafragma = mascaras.funcion_Circulo(radio_diafragmaInput, centro, xx_PlanoPupila, yy_PlanoPupila)


""" Se calcula el campo de entrada para el SEGUNDO TRAMO del arreglo 
NOTA: El campo que llega a la lente e interactua con la máscara asociada al diafragma
será el campo de entrada para el segundo tramo."""

#Calculando el campo de entrada al SEGUNDO TRAMO del arreglo
campo_entradaSegundoTramo_sinFiltro = campo_PlanoPupila*diafragma

#Calculando la intensidad del campo de entrada al SEGUNDO TRAMO del arreglo
intensidad_campoEntradaSegundoTramo_SinFiltro = np.abs(campo_entradaSegundoTramo_sinFiltro)**2


mascara_Filtrado = mascaras.funcion_Rectangulo(0.0002,0.0002,[-0.000622,0.000395],xx_PlanoPupila,yy_PlanoPupila)


""" Se calcula el campo de entrada para el SEGUNDO TRAMO del arreglo 
NOTA: El campo que llega a la pupila e interactua con la máscara asociada a la pupila
será el campo de entrada para el segundo tramo."""


#Calculando el campo de entrada al SEGUNDO TRAMO del arreglo
campo_entradaSegundoTramo = campo_PlanoPupila*diafragma*mascara_Filtrado

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
                                                    numero_onda_input,deltas_tramoPupilaMedicion)


#Se calcula la amplitud del campo de salida
amplitud_campoPlanoMedicion = np.abs(campo_PlanoMedicion)

#Se calcula la intensidad del campo de salida
intensidad_campoPlanoMedicion = amplitud_campoPlanoMedicion**2


""" Graficando máscara de transmitancia asignada a abertura circular"""
#graph.graficar_transmitancia(mascara,longitud_ArregloInput,longitud_ArregloInput,"Holograma")


"""Graficando amplitud del campo que en plano de Fourier"""
#graph.graficar_amplitud(np.abs(campo_PlanoPupila),anchoX_VentanaPlanoPupila,
#                              altoY_VentanaPlanoPupila,"Transformada de Fourier del objeto",1,0.001)

""" Graficando intensidad del campo que entra al SEGUNDO TRAMO del arreglo"""
#graph.graficar_intensidad(intensidad_campoEntradaSegundoTramo_SinFiltro,anchoX_VentanaPlanoPupila,
#                              altoY_VentanaPlanoPupila,"Transformada de Fourier del objeto",1,0.00001)


""" Graficando máscara de filtrado"""
#graph.graficar_transmitancia(mascara_Filtrado,anchoX_VentanaPlanoPupila,altoY_VentanaPlanoPupila,"Máscara de filtrado")


""" Graficando campo asociado a transformada de Fourier filtrado para encontrar imágen real """
#graph.graficar_intensidad(intensidad_campoEntradaSegundoTramo,anchoX_VentanaPlanoPupila,
#                             altoY_VentanaPlanoPupila,"Transformada de Fourier del objeto",1,0.001)


""" Graficando la intensidad del campo de salida del arreglo """
#graph.graficar_intensidad(intensidad_campoPlanoMedicion,anchoX_VentanaPlanoMedicion,altoY_VentanaPlanoMedicion,"Intensidad campo de salida")


#### SECCIÓN DE RECONSTRUCCIÓN DEL HOLOGRAMA ####

import matplotlib.image as mpimg
from scipy.ndimage import zoom



################ Parametros discretización ##########


num_pixels_x = 2048
num_pixels_y = 2048

delta_x = 3.45e-6 # Tamaño de cada pixel (m) 
delta_y = 3.45e-6 # Tamaño de cada pixel (m) 

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



'Tomar la matriz que respresenta la imagen y luego vamos a sacar su raíz para hallar el campo óptico a la salida'

#Se define la matriz de puntos para el análisis
matriz_campo = campo_PlanoMedicion

#Verificación del tamaño de la matriz de estudio
print("\n Tamaño  matriz campo óptico de salida:",matriz_campo.shape)


'Se realiza a continuación procedimiento para eliminar contribución de haz REFERENCIA'
# -----    ÁNGULO ENTRE HAZ DE REFERENCIA Y HAZ OBJETO --> 4.2255°  ------ #

#Definición de vector con cosenos directores asociado al haz referencia
#vector_cosenosDirectoresRef = [0.0622,0.0395] # NOTA--> Valores obtenidos en cálculo analítico
vector_cosenosDirectoresRef = [0.06,0.0388] 

#Definición de vector de onda asociado al haz de referencia
vector_ondaRef = [numero_onda_input*vector_cosenosDirectoresRef[0], numero_onda_input*vector_cosenosDirectoresRef[1]]

# Definición de onda plana INVERSA al haz de referencia  
onda_PlanaRefInversa = np.exp(-1j*((vector_ondaRef[0]*F_X)  + (vector_ondaRef[1]*F_Y)))


# Se multiplica el campo de interés (PLANO MEDICIÓN) con la onda plana generada
matriz_campoNOContribucionOndaPlana = matriz_campo*onda_PlanaRefInversa


# 1. Aplicamos FFT para hallar A[p,q,z]
'NOTA: nuestro espectro angular de salida sale con un shift a causa de la fft, por lo que antes de dividir la matriz'
' punto a punto debemos shiftear el resultado de la fft'

espectro_angular_salida = np.fft.fftshift(np.fft.fft2(matriz_campoNOContribucionOndaPlana)) 



for distancia_propagacionAribitraria in np.arange(0.0845, 0.087, 0.00025):


    # 2. Dividimos(matrices) punto a punto por la función de propagación para hallar A[p,q,0]

    funcion_transferencia = np.exp((1j*distancia_propagacionAribitraria*numero_onda_input)*(np.sqrt(1-(longitud_onda_input**2)*(F_X**2 + F_Y**2))))



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