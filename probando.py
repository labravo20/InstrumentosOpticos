""" DESCRIPCIÓN 

Documento para trabajar la reconstrucción de hologramas. El proceso consiste en aplicar la transformada de Fourier 
al campo del holograma para filtrar, en el dominio de Fourier, las frecuencias asociadas a una de las imágenes GEMELAS
con información sobre el holograma. Posteriormente, el campo resultante de la filtración es multiplicado por el inverso de 
la onda plana de referencia usada para realizar la interferencia, pues de esta forma de descarta el aporte de la interferencia
antes mencionada. El resultado de las anteriores operaciones es el objeto original magnificado por el sistema de amplificación
del microscopio.

"""


""" Importando librerias y documentos necesarios """

import LIBRERIA_Mascaras_Transmitancia as mascaras
import LIBRERIA_Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import matplotlib.pyplot as plt
import LIBRERIA_Funciones_importantes as function
import LIBRERIA_Funciones_Graficacion as graph
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches

' ################ INICIO SECCIÓN DE CARACTERIZACIÓN ARREGLO ################## '
##### LA SIGUIENTE SECCIÓN SE REFIERE A INPUTS QUE DEBEN SER ESPECIFICADOS POR EL USUARIO*

""" Definiendo parámetros de sensor """

#Numero de muestras/puntos distribuidos en el ancho del sensor
#resolucion_anchoSensorInput = 1280 
resolucion_anchoSensorInput = 1024  # Usando esta configuración se simplifica el cálculo para un sensor con distribución CUADRADA
                                    # NOTA: Se decide usar esta configuración porque ya se configuraron las características de 
                                    # reconstrucción adecuadas asociadas a este caso (específicamente las coordenadas del orden
                                    # +1 presente en la transformada de Fourier del holograma). 

#Numero de muestras/puntos distribuidos en el alto del sensor
resolucion_altoSensorInput = 1024

#Tamaño del pixel del sensor #UNIDADES: m
tamaño_PixelSensorInput = 5.2E-6



""" Definiendo parámetros de la iluminación coherente a utilizar  """

#Definiendo la longitud de onda asociada a la fuente en consideración
longitud_onda_input = 632.8E-9 #UNIDADES: m



""" Definiendo parámetros del OM usado para el sistema de magnificación """

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


' ################ EMPIEZA SECCIÓN DE CÁLCULOS PARA CARACTERIZACIÓN ARREGLO ################## '

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


' ################ FIN SECCIÓN DE CÁLCULOS PARA CARACTERIZACIÓN ARREGLO ################## '


' ################ EMPIEZA SECCIÓN DE IMPORTACIÓN Y PREPARACIÓN DE HOLOGRAMAS ################## '

""" Creación de mascaras de transmitancia """

# Cargar la imagen como máscara de transmitancia --> USAF PRUEBA
#ruta_imagen = "/home/labravo/Desktop/Instrumentos ópticos/PROYECTO/CARACTERIZACION_FASE/HOLOGRAMAS/test_USAF250nm.tif" #--> AJUSTAR ÁNGULO
#ruta_imagen = "/home/labravo/Desktop/Instrumentos ópticos/PROYECTO/CARACTERIZACION_FASE/HOLOGRAMAS/test_USAF300nm.tif" #--> FUNCIONA 
#ruta_imagen = "/home/labravo/Desktop/Instrumentos ópticos/PROYECTO/CARACTERIZACION_FASE/HOLOGRAMAS/test_USAF350nm.tif" #--> FUNCIONA    
#ruta_imagenReferencia = "/home/labravo/Desktop/Instrumentos ópticos/PROYECTO/CARACTERIZACION_FASE/HOLOGRAMAS/referenciaUSAF.tif" 


# Cargar la imagen como máscara de transmitancia --> USAF PRUEBA 
ruta_imagen = "/home/labravo/Desktop/Instrumentos ópticos/PROYECTO/HOLOGRAMAS/tejido_epitelialBucalJOSE4.0(HOLOGRAMA).tif" 
#ruta_imagen = "/home/labravo/Desktop/Instrumentos ópticos/PROYECTO/HOLOGRAMAS/tejido_epitelialBucalJOSE3.0(HOLOGRAMA).tif" 


""" Implementación de condición producto espacio frecuencia para cálculo de malla de puntos en dominio de Fourier """

# Crear malla de puntos de plano ESPECTRO DE FOURIER

delta_fx = (longitud_onda_input*distancia_focalMO)/(resolucion_anchoSensorInput * deltas_Sensor["deltaPlanoEntrada_X"]) #Tamaño de cada pixel en plano de Fourier (1/m)
delta_fy = (longitud_onda_input*distancia_focalMO)/(resolucion_altoSensorInput * deltas_Sensor["deltaPlanoEntrada_Y"]) #Tamaño de cada pixel en plano de Fourier (1/m)


#Calculando el ancho del arreglo en el plano de Fourier
anchoX_PlanoFourier = delta_fx*resolucion_anchoSensorInput
altoY_PlanoFourier = delta_fy*resolucion_altoSensorInput


# Crear la malla de puntos plano de entrada --> HOLOGRAMA A RECONSTRUIR 
xx_planoFourier, yy_planoFourier = mascaras.malla_Puntos(resolucion_anchoSensorInput, 
                                                        anchoX_PlanoFourier,
                                                        resolucion_altoSensorInput,
                                                        altoY_PlanoFourier)


#CREACIÓN DE FUNCIÓN PARA RECONSTRUCCIÓN 
def reconstruccion_Holograma(ruta_imagenHolograma):

    ' ################ EMPIEZA SECCIÓN DE RECONSTRUCCIÓN ################## '

    mascaraReconstruccion = function.cargar_imagen_png(ruta_imagenHolograma,resolucion_anchoSensorInput,resolucion_altoSensorInput)


    # Crear la malla de puntos plano de entrada --> HOLOGRAMA A RECONSTRUIR 
    xx_mascaraReconstruccion, yy_mascaraReconstruccion = mascaras.malla_Puntos(resolucion_anchoSensorInput, ancho_SensorInput,resolucion_altoSensorInput,
                                                alto_SensorInput)



    """ Calculando la transformada de Fourier del holograma """
    transformada_FourierReconstruccion = np.fft.fftshift(np.fft.fft2(mascaraReconstruccion))

    #Calculando la intensidad de la transformada de Fourier
    intensidad_TransformadaFourierRecosntruccion = (np.abs(transformada_FourierReconstruccion))**2

    #Normalizando la intensidad de la transformada de Fourier
    intensidad_TransformadaFourierRecosntruccionNormalizada = intensidad_TransformadaFourierRecosntruccion / np.max(intensidad_TransformadaFourierRecosntruccion)



    """Se diseña e implementa una máscara de transmitancia para filtrar la imágen gemela de interés """
    ### Se diseñan diferentes condiciones de máscara de filtrado dependiendo del ángulo del haz de referencia usado para 
    ### la formación del holograma

    #Definición del radio de la pupila en el plano de Fourier
    #radio_pupilaEscaladaPlanoFourier= 0.0003 # --> USAF
    radio_pupilaEscaladaPlanoFourier= 0.0004 # --> MUESTRAS CON ILUMINACIÓN 


    # Se define el vector asociado a las coordenadas de la máscara de filtrado  --> CALIBRANDO COORDENADAS TRANSFORMADA DE FOURIER
    #coordenadas_MascaraFiltrado = [-0.000928,-0.000665] #coordenada TOMAS EXPERIMENTALES 
    coordenadas_MascaraFiltrado = [-0.0008,-0.000665] #coordenada TOMAS EXPERIMENTALES --> prueba DESPLAZAMIENTO para garantizar 
                                                    #filtro "completo"

    #Creación de la máscara de filtrado para el proceso de reconstrucción
    mascaraReconstruccion_Filtrado = mascaras.funcion_Circulo(radio_pupilaEscaladaPlanoFourier,coordenadas_MascaraFiltrado,
                                                            xx_planoFourier,yy_planoFourier)

    #Se implementa la filtración para obtener el campo filtrado
    campo_Filtrado = mascaraReconstruccion_Filtrado*transformada_FourierReconstruccion

    #Intensidad del campo resultante después del proceso de filtrado
    intensidad_CampoFiltrado = (np.abs(campo_Filtrado))**2

    """ Se vuelve al dominio espacial para determinar el campo complejo asociado al campo filtrado """

    #Calculando la transformada de fourier inversa del campo shifteado
    campo_FiltradoShifteado = np.fft.fftshift(campo_Filtrado)

    #Calculando la transformada de Fourier inversa
    campo_Reconstruccion = np.fft.ifft2(campo_FiltradoShifteado)



    """ Se diseña e implementa una onda plana, inversa al haz de referencia, para despreciar el aporte de la interferencia """

    #Definición de onda plana para eliminar aporte de efectos de interferencia
    #vector_cosenosDirectoresRef =  [-0.04648,-0.0333] #coordenada TOMAS EXPERIMENTALES --> TEST USAF -0.0467,-0.0333
    vector_cosenosDirectoresRef =  [-0.0467,-0.0333] #coordenada TOMAS EXPERIMENTALES --> TEJIDO BUCAL EPITELIAL JOSE

    #Definición de vector de onda asociado al haz de referencia
    vector_ondaRef = [numero_onda_input*vector_cosenosDirectoresRef[0], numero_onda_input*vector_cosenosDirectoresRef[1]]

    # Definición de onda plana INVERSA al haz de referencia  
    onda_PlanaRefInversa = np.exp(-1j*((vector_ondaRef[0]*xx_mascaraReconstruccion)  + (vector_ondaRef[1]*yy_mascaraReconstruccion)))


    # Se multiplica el campo de interés (PLANO MEDICIÓN) con la onda plana generada
    matriz_campoNOContribucionOndaPlana = campo_Reconstruccion*onda_PlanaRefInversa

    #Intensidad del campo resultante de la reconstrucción
    intensidad_matrizCampoNOcontribucionOndaPlana = (np.abs(matriz_campoNOContribucionOndaPlana))**2
    
    ' ################ FIN SECCIÓN DE RECONSTRUCCIÓN ################## '

    return {"Mascara": mascaraReconstruccion,"Intensidad_TransformadaFourier": intensidad_TransformadaFourierRecosntruccionNormalizada,
            "Mascara_ReconstruccionFiltrado":mascaraReconstruccion_Filtrado,"Intensidad_CampoFiltrado":intensidad_CampoFiltrado,
            "Matriz_NOOndaPlana": matriz_campoNOContribucionOndaPlana,"Intensidad_matrizNOOndaPlana":intensidad_matrizCampoNOcontribucionOndaPlana}



""" Llamando función de reconstrucción para obtener resultado de HOLOGRAMA de interés """

objeto_Reconstruido = reconstruccion_Holograma(ruta_imagen)

objeto_Reconstruido_Desfasado = objeto_Reconstruido["Matriz_NOOndaPlana"]*np.exp(1j*((np.pi))/(2)) # TEJIDO EPITELIAL 4.0
#objeto_Reconstruido_Desfasado = objeto_Reconstruido["Matriz_NOOndaPlana"] # TEJIDO EPITELIAL 3.0


""" Llamando función de reconstrucción para obtener resultado de la REFERENCIA --> USAF"""

#referencia_Reconstruida = reconstruccion_Holograma(ruta_imagenReferencia)

' ################ FIN SECCIÓN DE IMPORTACIÓN Y PREPARACIÓN DE HOLOGRAMAS ################## '


' ################ EMPIEZA SECCIÓN DE GRAFICACIÓN ################## '

""" Graficando máscara de transmitancia"""

graph.graficar_intensidad(objeto_Reconstruido["Mascara"],ancho_SensorInput,alto_SensorInput,"Holograma de entrada para reconstrucción")



""" Graficando intensisdad de la transformada de Fourier del holograma"""

graph.graficar_intensidad(np.log(objeto_Reconstruido["Intensidad_TransformadaFourier"]),anchoX_PlanoFourier,altoY_PlanoFourier,
                             "Intensidad normalizada de la Transformada de Fourier del Holograma",1,1)



""" Graficando máscara de filtración para obtener orden +1 de TF"""

graph.graficar_transmitancia(objeto_Reconstruido["Mascara_ReconstruccionFiltrado"],anchoX_PlanoFourier,altoY_PlanoFourier,"Mascara Filtrado")



""" Graficando intensidad del objeto reconstruido """

graph.graficar_intensidad(objeto_Reconstruido["Intensidad_matrizNOOndaPlana"],ancho_SensorInput,alto_SensorInput,
                          "Intensidad asociada al Campo óptico del objeto reconstruido")



""" Graficando información de fase del objeto reconstruido"""

graph.graficar_fase(np.angle(objeto_Reconstruido["Matriz_NOOndaPlana"]),ancho_SensorInput,alto_SensorInput,
                    "Distribución de fase Campo óptico objeto")



""" Graficando información de fase del objeto reconstruido desfasado --> TEJIDO EPITELIAL"""

graph.graficar_fase(np.angle(objeto_Reconstruido_Desfasado),ancho_SensorInput,alto_SensorInput,
                    "Distribución de fase Campo óptico objeto con desfase")


""" Graficando información sobre REFERENCIA --> USAF"""

#graph.graficar_fase(np.angle(referencia_Reconstruida["Matriz_NOOndaPlana"]),ancho_SensorInput,alto_SensorInput,
#                    "Distribución de fase referencia")


' ################ FIN SECCIÓN DE GRAFICACIÓN ################## '


' ################ EMPIEZA SECCIÓN DE CREACIÓN SUBMATRIZ PARA ANALIZAR UNA SECCIÓN DE LA MUESTRA --> TEJIDO EPITELIAL ################## '

# Coordenada espacial X en metros --> EJE CENTRAL DEL RECTANGULO DE SELECCIÓN DE ÁREA DE INTERÉS
#X_m_central = 0.0002  #MUESTRAS PARTE SUPERIOR 4.0
X_m_central = 0.000795  #MUESTRAS PARTE INFERIOR 4.0
#X_m_central = 0.0008  #MUESTRAS 3.0

#Ancho del rectangulo
#ancho_rectangulo = 0.0016 #MUESTRAS PARTE SUPERIOR 4.0
#ancho_rectangulo = 0.00233 #MUESTRAS PARTE INFERIOR 4.0
ancho_rectangulo = 0.0013 #MUESTRAS 3.0

# Definiendo límites de la región cuadrada --> LIMITES EN EJE VERTICAL DE SELECCIÓN DE ÁREA DE INTERÉS
Y_min = 0.0008  # Límite inferior en metros
Y_max = 0.002  # Límite superior en metros


# Convertir coordenadas espaciales a índices de matriz
col_min_idx = int((X_m_central - (ancho_rectangulo/2) + (ancho_SensorInput / 2)) * (resolucion_anchoSensorInput / ancho_SensorInput))  
col_max_idx = int((X_m_central + (ancho_rectangulo/2) + (ancho_SensorInput / 2)) * (resolucion_anchoSensorInput / ancho_SensorInput))  
row_min_idx = int(((alto_SensorInput / 2) - Y_max) * (resolucion_altoSensorInput / alto_SensorInput))
row_max_idx = int(((alto_SensorInput / 2) - Y_min) * (resolucion_altoSensorInput / alto_SensorInput))

#Función para graficar sección seleccionada y para generar submatrizz (también es graficada esta nueva matriz de estudio)
def submatriz(objeto_Reconstruido_Desfasado):


    # Graficar la imagen original con la región marcada
    plt.figure(figsize=(6,5))
    plt.imshow(np.angle(objeto_Reconstruido_Desfasado), 
               extent=[-ancho_SensorInput/2, ancho_SensorInput/2,
                       -alto_SensorInput/2, alto_SensorInput/2], 
               cmap='inferno')

    # Dibujar el rectángulo de la región seleccionada
    rect = patches.Rectangle((X_m_central - (ancho_rectangulo/2), Y_min), ancho_rectangulo, Y_max - Y_min, 
                             linewidth=1.5, edgecolor='red', facecolor='none', linestyle="--")
    plt.gca().add_patch(rect)

    # Configuración de la gráfica original
    plt.title("Región seleccionada en la imagen")
    plt.colorbar(label="Fase del objeto")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.show()

    # Extraer la submatriz dentro del rectángulo
    submatriz_resultado = objeto_Reconstruido_Desfasado[row_min_idx:row_max_idx, col_min_idx:col_max_idx]

    # Graficar la fase de la submatriz extraída
    plt.figure(figsize=(6,5))
    plt.imshow(np.angle(submatriz_resultado), 
                extent=[X_m_central - (ancho_rectangulo/2), X_m_central + (ancho_rectangulo/2), Y_min, Y_max], 
                cmap='inferno')

    # Configuración de la gráfica de la submatriz
    plt.title("Fase de la región recortada")
    plt.colorbar(label="Fase")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.show()

    

    return submatriz_resultado  # Retornar la submatriz extraída

#Asociando la submatriz seleccionada en una variable
submatriz_muestraDesfasada = submatriz(objeto_Reconstruido_Desfasado)

' ################ FIN SECCIÓN DE CREACIÓN SUBMATRIZ PARA ANALIZAR UNA SECCIÓN DE LA MUESTRA --> TEJIDO EPITELIAL ################## '


' ################ EMPIEZA SECCIÓN DE MEDICIÓN RELATIVA DE FASE --> USAF ################## '

# Se calcula la diferencia de fase entre el holograma de interés y la referencia pre establecida
#diferencia_FaseHologramaReferencia = (objeto_Reconstruido["Matriz_NOOndaPlana"])/(referencia_Reconstruida["Matriz_NOOndaPlana"])

# Graficando resultados de la medida de diferencia de Fase
#graph.graficar_fase(np.angle(diferencia_FaseHologramaReferencia),ancho_SensorInput,alto_SensorInput,
#                    "Fase relativa de la muestra") 

' ################ FIN SECCIÓN DE MEDICIÓN RELATIVA DE FASE --> USAF  ################## '


' ################ EMPIEZA SECCIÓN DE CÁLCULO DE LONGITUD DE CAMINO ÓPTICO RELATIVO ################## '

#Calculando la longitud de camino óptico de la muestra --> TEJIDO EPITELIAL
longitud_caminoOpticoMuestra = np.angle(submatriz_muestraDesfasada)/numero_onda_input

#Calculando altura de la muestra NORMALIZADA --> TEJIDO EPITELIAL
longitud_caminoOpticoMuestraNormalizada = longitud_caminoOpticoMuestra/np.max(longitud_caminoOpticoMuestra)

#Graficando el camino óptico de la muestra --> TEJIDO EPITELIAL
plt.figure(figsize=(6,5))
plt.imshow(longitud_caminoOpticoMuestra, 
            extent=[X_m_central - (ancho_rectangulo/2), X_m_central + (ancho_rectangulo/2), Y_min, Y_max], 
            cmap='winter')

# Configuración de la gráfica de la submatriz
plt.title("Longitud de camino óptico de la muestra")
plt.colorbar(label="Longitud camino óptico [m]")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()

#Definiendo el índice de refracción asociado al arreglo --> USAF
#diferencia_indiceRefraccion = 1.52 - 1 #PRIMERO: Indice de refracción del acrilato/porta muestras
                                       #SEGUNDO: Indice de refracción de la REFERENCIA (aire)

#Calculando el camino óptico asociado a la medida relativa de fase --> USAF
#altura_muestra = np.angle(diferencia_FaseHologramaReferencia)/(numero_onda_input*diferencia_indiceRefraccion)

#Calculando altura de la muestra NORMALIZADA --> USAF
#altura_muestraNormalizada = altura_muestra/np.max(altura_muestra)

#Graficando altura de la muestra --> USAF
#graph.graficar_altura(altura_muestra,ancho_SensorInput,alto_SensorInput,"Altura de la muestra")

' ################ FIN SECCIÓN DE CÁLCULO DE LONGITUD DE CAMINO ÓPTICO RELATIVO ################## '


' ################ EMPIEZA SECCIÓN DE CÁLCULO GRÁFICOS 1D ################## '

# Coordenada espacial X en m que queremos analizar --> Recta en vertical de la cual se quiere conocer perfil de altura/camino óptico
X_m = 0.0006 


#GRAFICANDO RECTA EN IMAGEN DE ALTURA --> PARA IDENTIFICAR PERFIL DE GRAFICACIÓN  USAF
# plt.imshow(altura_muestra, extent=[-ancho_SensorInput/2, ancho_SensorInput/2,
#                             -alto_SensorInput/2, alto_SensorInput/2], 
#                             cmap='winter')
# plt.axvline(X_m,color = 'red',linestyle = "--",label = f"X = {X_m} m")
# plt.title("Altura de la muestra")
# plt.colorbar(label="Altura [m]")
# plt.xlabel("X (m)")
# plt.ylabel("Y (m)")
# plt.show()

#GRAFICO PERFIL 1D de la muestra --> TEJIDO EPITELIAL
plt.figure(figsize=(6,5))
plt.imshow(longitud_caminoOpticoMuestra, 
            extent=[X_m_central - (ancho_rectangulo/2), X_m_central + (ancho_rectangulo/2), Y_min, Y_max], 
            cmap='winter')
plt.axvline(X_m,color = 'red',linestyle = "--",label = f"X = {X_m} m")
# Configuración de la gráfica de la submatriz
plt.title("Longitud de camino óptico de la muestra")
plt.colorbar(label="Longitud camino óptico [m]")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()


# Definiendo límites verticales de la recta
Y_min =  0.0009 # Límite inferior Y en metros 
Y_max = 0.0019   # Límite superior Y en metros  


# Convertir la coordenada espacial en índice de matriz
#col_idx = int((X_m + ancho_SensorInput / 2) * (resolucion_anchoSensorInput / ancho_SensorInput)) #USAF
col_idx = int((X_m + (ancho_rectangulo/2)) * ( (col_max_idx - col_min_idx) / ancho_rectangulo)) #TEJIDO EPITELIAL


# Generar las coordenadas Y en metros 
#y_values = np.linspace(-alto_SensorInput/2, alto_SensorInput/2, resolucion_altoSensorInput) #USAF
y_values = np.linspace(Y_min, Y_max, row_max_idx - row_min_idx) #TEJIDO EPITELIAL

# Extraer la columna de la imagen
#column_values = altura_muestra[:, col_idx] #USAF
column_values = longitud_caminoOpticoMuestra[:, col_idx] #TEJIDO EPITELIAL

# Filtrar valores dentro del rango deseado
mask = (y_values >= Y_min) & (y_values <= Y_max)
y_values_recortado = y_values[mask]
column_values_recortado = column_values[mask]

# Graficar la columna en función de Y con límites específicos 
plt.plot(y_values_recortado, column_values_recortado)
plt.xlabel("Posición Y [m]")

# plt.ylabel("Altura [m]") #USAF
# plt.title(f"Perfil de Altura en X = {X_m} m (Limitado a [{Y_min}, {Y_max}])") #USAF


plt.ylabel("Longitud camino óptico [m]") #TEJIDO EPITELIAL
plt.title(f"Perfil de longitud camino óptico en X = {X_m} m (Limitado a [{Y_min}, {Y_max}])") #TEJIDO EPITELIAL

plt.show()

' ################ FIN SECCIÓN DE CÁLCULO GRÁFICOS 1D ################## '


' ################ EMPIEZA SECCIÓN DE CÁLCULO GRÁFICOS 3D ################## '

# Crear figura y ejes 3D
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Crear las coordenadas X e Y --> USAF
# x = np.linspace(-ancho_SensorInput / 2, ancho_SensorInput / 2, resolucion_anchoSensorInput)
# y = np.linspace(-alto_SensorInput / 2, alto_SensorInput / 2, resolucion_altoSensorInput)

# Crear las coordenadas X e Y --> TEJIDO EPITELIAL
x = np.linspace(X_m_central - (ancho_rectangulo/2), X_m_central + (ancho_rectangulo/2), col_max_idx - col_min_idx)
y = np.linspace(Y_min, Y_max, row_max_idx - row_min_idx)

#Creando malla de puntos
X, Y = np.meshgrid(x, y)

# Graficar la superficie en 3D
#ax.plot_surface(X, Y, -altura_muestraNormalizada, cmap='viridis', edgecolor='none') #USAF
ax.plot_surface(X, Y, longitud_caminoOpticoMuestraNormalizada, cmap='viridis', edgecolor='none') #TEJIDO EPITELIAL

# Etiquetas de los ejes
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')

# ax.set_zlabel('Altura Normalizada') #USAF
ax.set_zlabel('Longitud camino óptico normalizado') #TEJIDO EPITELIAL

ax.set_title('Reconstrucción 3D de la muestra')

# Permitir interacción para rotar
plt.show()

' ################ FIN SECCIÓN DE CÁLCULO GRÁFICOS 3D ################## '

