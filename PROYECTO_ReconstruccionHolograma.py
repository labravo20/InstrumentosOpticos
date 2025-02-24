""" Importando librerias y documentos necesarios """

import LIBRERIA_Mascaras_Transmitancia as mascaras
import LIBRERIA_Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import matplotlib.pyplot as plt
import LIBRERIA_Funciones_importantes as function
import LIBRERIA_Funciones_Graficacion as graph
import matplotlib.image as mpimg
from scipy.ndimage import zoom

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
magnificacion = 20

#Apertura numérica objetivo de microscopio
apertura_Numerica = 0.4 

#Distancia focal del objetivo de microscopio #UNIDADES: m
### NOTA: Si se desconoce la distancia focal del objetivo de microscopio escribir el número CERO (0) 
distancia_focalMO = 0

#Distancia focal asociada a la lente de tubo --> #UNIDADES: m
distancia_focalTL = 0.18   



""" Definiendo parámetros del montaje  """

#Se define una distancia de propagación arbitraria 'd' para simulación de sistema 4F --> #UNIDADES: m
distancia_propagacionAribitraria = 0.01 

#Se define distancia de propagación para formación del holograma (PLANO IMAGEN --> PLANO HOLOGRAMA)
distancia_ImagenHolograma = 0.0001 #UNIDADES: m



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
angulo_HazReferencia = 1 #UNIDADES: Grados

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


################### INICIO SECCIÓN PARA RECONSTRUCCIÓN DEL HOLOGRAMA #########################

' ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL PRIMER TRAMO ------ '

'NOTA --> Para la construcción del sistema 4f se usarán distancias arbitrarias que no afectan los resultados del proceso de reconstrucción.'

""" Se calculan las matrices necesarias para estudiar el PRIMER TRAMO del arreglo difractivo
    OBJETO --> *Propagación* --> LENTE  --> PUPILA """

#Se calcula la matriz asociada al proceso de propagación desde el plano objeto
#hasta el plano de la lente
matriz_propagacionPrimerTramo = matriz.propagacion_MedioHomogeneo(distancia_focalMO)

matriz_lente = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focalMO)

matriz_propagacion02 = matriz.propagacion_MedioHomogeneo(distancia_focalMO)



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

' ------ FIN SECCIÓN DE CÁLCULO MATRICIAL PRIMER TRAMO ------ '


' ------ EMPIEZA SECCIÓN DE CÁLCULO MATRICIAL SEGUNDO TRAMO ------ '

""" Se calculan las matrices necesarias para estudiar el SEGUNDO TRAMO del arreglo difractivo
    PUPILA--> *propagación(distancia arbitraria d)* --> LENTE02 --> *propagación(distancia focal 02)*
    --> IMAGEN """

#Se calcula la matriz asociada al proceso de propagación desde el plano de la pupila 
# hasta el plano de la lente 02
matriz_propagacion01SegundoTramo = matriz.propagacion_MedioHomogeneo(distancia_propagacionAribitraria)

#Se calcula la matriz asociada a la interacción con la lente 02
matriz_lente02 = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focalTL)

#Se calcula la matriz asociada al proceso de propagación desde el plano de la lente 02
#hasta el plano de medición
matriz_propagacion02SegundoTramo = matriz.propagacion_MedioHomogeneo(distancia_focalTL)



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


' ------ FIN SECCIÓN DE CÁLCULO MATRICIAL SEGUNDO TRAMO ------ '


' ----- EMPIEZA SECCIÓN DE CONFIGURACIÓN DE MALLAS DE PUNTOS PARA CADA PLANO ------ '

""" Creando máscara de plano objeto  """

# Crear la malla de puntos
xx_mascaraReconstruccion, yy_mascaraReconstruccion = mascaras.malla_Puntos(resolucion_anchoSensorInput, ancho_SensorInput,
                                               resolucion_altoSensorInput,alto_SensorInput)



""" Creando malla de puntos asociada a plano de la PUPILA """

#Se llama función para determinar los deltas de muestreo del tramo PUPILA --> MÁSCARA
deltas_tramoMascaraPupilaReconstruccion = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_anchoSensorInput,
                                                                                            ancho_SensorInput,
                                                                                            matriz_SistemaPrimerTramo[0,1],
                                                                                            longitud_onda_input,
                                                                                            resolucion_altoSensorInput,
                                                                                            alto_SensorInput)

#Se calcula el ancho de la ventana del plano de la máscara u objeto de entrada
anchoX_VentanaPlanoPupilaReconstruccion = resolucion_anchoSensorInput*deltas_tramoMascaraPupilaReconstruccion["deltaPlanoEntrada_X"]
altoY_VentanaPlanoPupilaReconstruccion = resolucion_altoSensorInput*deltas_tramoMascaraPupilaReconstruccion["deltaPlanoEntrada_Y"]

#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoPupilaReconstruccion, yy_PlanoPupilaReconstruccion = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoPupilaReconstruccion,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoPupilaReconstruccion)


""" Creando malla de puntos asociada a plano Medición """

#Se llama función para determinar los deltas de muestreo del tramo SENSOR --> PUPILA
deltas_tramoPupilaMedicionReconstruccion = function.producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_anchoSensorInput,
                                                                                             anchoX_VentanaPlanoPupilaReconstruccion,
                                                                                             matriz_SistemaSegundoTramo[0,1],
                                                                                             longitud_onda_input,
                                                                                             resolucion_altoSensorInput,
                                                                                             altoY_VentanaPlanoPupilaReconstruccion)

#Se calcula el ancho de la ventana del plano anterior a pupila
anchoX_VentanaPlanoMedicionReconstruccion = resolucion_anchoSensorInput*deltas_tramoPupilaMedicionReconstruccion["deltaPlanoEntrada_X"]
altoY_VentanaPlanoMedicionReconstruccion = resolucion_altoSensorInput*deltas_tramoPupilaMedicionReconstruccion["deltaPlanoEntrada_Y"]

#Se calcula la malla de puntos asociada al plano de la pupila
xx_PlanoMedicionReconstruccion, yy_PlanoMedicionReconstruccion = mascaras.malla_Puntos(resolucion_anchoSensorInput,anchoX_VentanaPlanoMedicionReconstruccion,
                                                       resolucion_altoSensorInput,altoY_VentanaPlanoMedicionReconstruccion)

' ----- FIN SECCIÓN DE CONFIGURACIÓN DE MALLAS DE PUNTOS PARA CADA PLANO ------ '


' ------ EMPIEZA SECCIÓN DE CÁLCULO RESULTADO DIFRACTIVO DE CADA TRAMO ------ '

""" Creando máscara de transmitancia asociada al objeto de estudio en el arreglo """

# Cargar la imagen PNG como máscara de transmitancia
ruta_imagen_png = "/home/labravo/Downloads/Holograma.png"  # Especifica la ruta de imagen
mascara = function.cargar_imagen_png(ruta_imagen_png, resolucion_anchoSensorInput,resolucion_altoSensorInput)


""" Se calcula el resultado del proceso difractivo del PRIMER TRAMO"""

#Se calcula el campo de salida/en plano de la lente --> Campo resultante del primer tramo 
campo_PlanoPupilaReconstruccion = matriz.matriz_ABCD_Difraccion_Sensor(camino_opticoCentralPrimerTramo, mascara,
                                                 matriz_SistemaPrimerTramo[0,0],
                                                 matriz_SistemaPrimerTramo[0,1],
                                                 matriz_SistemaPrimerTramo[1,1],xx_mascaraReconstruccion,yy_mascaraReconstruccion,
                                                 xx_PlanoPupilaReconstruccion,yy_PlanoPupilaReconstruccion,numero_onda_input,
                                                 deltas_tramoMascaraPupilaReconstruccion)


intensidad_campo_PlanoPupilaReconstruccion = (np.abs(campo_PlanoPupilaReconstruccion))**2

""" Se crea diafragma para delimitar dominio del arreglo asociado a lentes FINITAS
NOTA: El diafragma va a construirse a partir de la malla de puntos asociada al plano de la Lente."""

#Creación de máscara circular que representará el diafragma 
pupilaReconstruccion = mascaras.funcion_Circulo(radio_pupilaInput, None, xx_PlanoPupilaReconstruccion, yy_PlanoPupilaReconstruccion)


""" Se calcula el campo de entrada para el SEGUNDO TRAMO del arreglo 
NOTA: El campo que llega a la lente e interactua con la máscara asociada al diafragma
será el campo de entrada para el segundo tramo."""

#Calculando el campo de entrada al SEGUNDO TRAMO del arreglo
campo_entradaSegundoTramoReconstruccion_sinFiltro = campo_PlanoPupilaReconstruccion*pupilaReconstruccion

#Calculando la intensidad del campo de entrada al SEGUNDO TRAMO del arreglo
intensidad_campoEntradaSegundoTramoReconstruccion_SinFiltro = np.abs(campo_entradaSegundoTramoReconstruccion_sinFiltro)**2


mascara_Filtrado = mascaras.funcion_Rectangulo(0.0002,0.0002,[-0.000622,0.000395],xx_PlanoPupilaReconstruccion,yy_PlanoPupilaReconstruccion)


""" Se calcula el campo de entrada para el SEGUNDO TRAMO del arreglo 
NOTA: El campo que llega a la pupila e interactua con la máscara asociada a la pupila
será el campo de entrada para el segundo tramo."""


#Calculando el campo de entrada al SEGUNDO TRAMO del arreglo
campo_entradaSegundoTramoReconstruccion = campo_PlanoPupilaReconstruccion*pupilaReconstruccion*mascara_Filtrado

#Calculando la intensidad del campo de entrada al SEGUNDO TRAMO del arreglo
intensidad_campoEntradaSegundoTramoReconstruccion = np.abs(campo_entradaSegundoTramoReconstruccion)**2


""" Se calcula el resultado del proceso difractivo del SEGUNDO TRAMO """

#Se calcula el campo de salida/en plano de medición --> Campo resultante de la difracción 
campo_PlanoMedicionReconstruccion = matriz.matriz_ABCD_Difraccion_Sensor_Shift(camino_opticoCentralSegundoTramo,
                                                    campo_entradaSegundoTramoReconstruccion,
                                                    matriz_SistemaSegundoTramo[0,0],
                                                    matriz_SistemaSegundoTramo[0,1],
                                                    matriz_SistemaSegundoTramo[1,1],xx_PlanoPupilaReconstruccion,
                                                    yy_PlanoPupilaReconstruccion,xx_PlanoMedicionReconstruccion,yy_PlanoMedicionReconstruccion,
                                                    numero_onda_input,deltas_tramoPupilaMedicionReconstruccion)


#Se calcula la amplitud del campo de salida
amplitud_campoPlanoMedicionReconstruccion = np.abs(campo_PlanoMedicionReconstruccion)

#Se calcula la intensidad del campo de salida
intensidad_campoPlanoMedicionReconstruccion = amplitud_campoPlanoMedicionReconstruccion**2

' ------ FIN SECCIÓN DE CÁLCULO RESULTADO DIFRACTIVO DE CADA TRAMO ------ '


' ------ EMPIEZA SECCIÓN DE GRAFICACIÓN ------ '

""" Graficando mascara de entrada al sistema 4f """
graph.graficar_transmitancia(mascara,ancho_SensorInput,alto_SensorInput,"HOLOGRAMA")

""" Graficando intensidad del campo que entra al SEGUNDO TRAMO del arreglo"""
graph.graficar_intensidad(intensidad_campoEntradaSegundoTramoReconstruccion_SinFiltro,anchoX_VentanaPlanoPupilaReconstruccion,
                              altoY_VentanaPlanoPupilaReconstruccion,"Transformada de Fourier del objeto",1,0.0001)

' ------ FIN SECCIÓN DE GRAFICACIÓN ------ '



################ Parametros discretización ##########


num_pixels_x = 2848
num_pixels_y = 2848

delta_x = 2.74E-6 # Tamaño de cada pixel (m) 
delta_y = 2.74E-6 # Tamaño de cada pixel (m) 

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
matriz_campo = campo_PlanoMedicionReconstruccion

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