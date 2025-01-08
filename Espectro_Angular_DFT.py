#Impresión para verificar correcto funcionamiento entorno de programación
print("Inicializando entorno de programación Espectro Angular...")

import numpy as np
import matplotlib.pyplot as plt
import Mascaras_Transmitancia as m


def producto_espacio_frecuencia(longitud_onda,z,resolucion,longitud_Arreglo):

    delta_input = longitud_Arreglo/resolucion #Definiendo separación entre número total de muestras

    delta_output = (1)/(resolucion*delta_input) #Definiendo separación entre muestras
                                                            #en plano de medición

    delta_muestreo = [delta_input,delta_output]

    return delta_muestreo

def ecuacion_DFT(resolucion):
    #NOTA --> Se hace ingreso de los términos numero imaginario y pi porque estos
    #se van a definir una única vez en el código (para optimizar el cálculo), asociados a constantes.

    #A continuación se define el término de la exponencial BASE para cálculo de DFT
    termino_Fourier = np.exp((-1j*2*np.pi)/resolucion)

    #A continuación se crea una Matriz llena de 1 (tamaño resolución x resolución)
    matrizDFT = np.ones((resolucion,resolucion),dtype=np.complex64)  

    #Establecimiento primera sumatoria
    for n in range(0,resolucion,1):

        #Estableciendo segunda sumatoria
        for m in range(0,resolucion,1):

            matrizDFT[n,m] = (termino_Fourier**n)**m #Se asigna el valor correspondiente
                                                          #a cada coordenada de la matriz


    return matrizDFT

def calculoDFT1(vector):

    resolucion = len(vector) #Calcular la longitud del vector de entrada

    DFT1_OUT =  ecuacion_DFT(resolucion) @ vector   #Realizando el cálculo de la DFT

    return DFT1_OUT

def calculoDFT2(campo):

    #Calculando las dimensiones de fila y columna del campo de entrada
    fila, columna = campo.shape  
    
    # Calculando la DFT para las filas

    ##  Crea un array con la misma forma y tipo de datos del campo de entrada
    filasDFT = np.zeros_like(campo, dtype=np.complex64)  

    # Ciclo con longitud igual al número de filas del arreglo
    for i in range(0,fila,1):  
        filasDFT[i, :] = calculoDFT1(campo[i, :])  # Calcula las DFT de las filas
    
    # Calculando la DFT para las columnas
    DFT_2D = np.zeros_like(filasDFT, dtype=np.complex64)
    for j in range(0,columna,1):
        DFT_2D[:, j] = calculoDFT1(filasDFT[:, j])    # Calcula las DFT sobre las columnas
    
    return DFT_2D

def funcionShift_DFT(campo):

    #Se define variales a las cuales asignar la forma de la "imagen/función" a shiftear
    M, N = campo.shape
    
    #Primero se divide el campo a evaluar en cuatro cuadrantes
    campoShifted = np.zeros_like(campo)

    # Cuadrante superior izquierdo -> cuadrante inferior derecho
    campoShifted[:M//2, :N//2] = campo[M//2:, N//2:]
    
    # Cuadrante superior derecho -> cuadrante inferior izquierdo
    campoShifted[:M//2, N//2:] = campo[M//2:, :N//2]
    
    # Cuadrante inferior izquierdo -> cuadrante superior derecho
    campoShifted[M//2:, :N//2] = campo[:M//2, N//2:]
    
    # Cuadrante inferior derecho -> cuadrante superior izquierdo
    campoShifted[M//2:, N//2:] = campo[:M//2, :N//2]

    return campoShifted

def ecuacion_DFT_inversa(resolucion):

    #NOTA --> Se hace ingreso de los términos numero imaginario y pi porque estos
    #se van a definir una única vez en el código (para optimizar el cálculo), asociados a constantes.

    #A continuación se define el término de la exponencial BASE para cálculo de DFT
    termino_Fourier = np.exp((1j*2*np.pi)/resolucion)

    #A continuación se crea una Matriz llena de 1 (tamaño resolución x resolución)
    matrizDFT = np.ones((resolucion,resolucion),dtype=np.complex64)  

    #Establecimiento primera sumatoria
    for n in range(0,resolucion,1):

        #Estableciendo segunda sumatoria
        for m in range(0,resolucion,1):

            valorSumatoria = 0 #Variable para cargar resultado de sumatorias

            matrizDFT[n,m] = (termino_Fourier**n)**m #Se asigna el valor correspondiente
                                                          #a cada coordenada de la matriz

            #Se establecen condicones para sumatoria en dominio de variables de salida
            """for xx_input in range(0,resolucion,1):

                for yy_input in range(0,resolucion,1):

                    #A continuación se define el término de la exponencial para el calculo de la DFT
                    termino_Fourier = np.exp((-numero_imaginario*2*numero_pi*(((xx_input*n)+(yy_input*m))/((deltas[0])*(deltas[1]))))/resolucion)

                    #Se actualiza la variable de la suma
                    valorSumatoria += campo*termino_Fourier

            #Declaración de resultado definitivo
            resultado_DFT = valorSumatoria"""

    return matrizDFT
    
def calculoDFT1_inversa(vector):

    resolucion = len(vector) #Calcular la longitud del vector de entrada

    DFT1_OUT =  ecuacion_DFT_inversa(resolucion) @ vector   #Realizando el cálculo de la DFT

    return DFT1_OUT

def calculoDFT2_inversa(campo):

    #Calculando las dimensiones de fila y columna del campo de entrada
    fila, columna = campo.shape  
    
    # Calculando la DFT para las filas

    ##  Crea un array con la misma forma y tipo de datos del campo de entrada
    filasDFT = np.zeros_like(campo, dtype=np.complex64)  

    # Ciclo con longitud igual al número de filas del arreglo
    for i in range(0,fila,1):  
        filasDFT[i, :] = calculoDFT1_inversa(campo[i, :])  # Calcula las DFT de las filas
    
    # Calculando la DFT para las columnas
    DFT_2D = np.zeros_like(filasDFT, dtype=np.complex64)
    for j in range(0,columna,1):
        DFT_2D[:, j] = calculoDFT1_inversa(filasDFT[:, j])    # Calcula las DFT sobre las columnas
    
    return DFT_2D


""" Definición de variables input para creación malla de puntos máscara difractiva"""
resolucion_input = 20
longitud_arreglo_input = 0.0025 #UNIDADES: m

""" Definición de variables input para creación abertura circular"""
radio_input = 0.0005 #UNIDADES: m
centro_input = None 

""" Creación de malla de puntos máscara difractiva"""
xx,yy = m.malla_Puntos(resolucion_input,longitud_arreglo_input)

""" Creación de abertura circular """
mascara_circular = m.funcion_Circulo(radio_input,centro_input,xx,yy)

""" Parámetros de configuración del arreglo difractivo """
longitud_onda_input = 633E-9 #UNIDADES: m
z_input = 0.015 #UNIDADES: m
numero_im = 1j
numero_pi = np.pi
numero_onda = ((2*numero_pi)/longitud_onda_input)

""" Definición de variables input para creación malla de puntos plano medición"""

#Llamando a la función para determinar los deltas de input y output
delta_muestreo = producto_espacio_frecuencia(longitud_onda_input,z_input,resolucion_input,longitud_arreglo_input)

resolucion_medicion = resolucion_input 
longitud_arreglo_medicion = (delta_muestreo[1])*resolucion_medicion #UNIDADES: m

""" Creación de malla de puntos plano de medición"""
xx_espectro,yy_espectro = m.malla_Puntos(resolucion_medicion,longitud_arreglo_medicion)

""" Definición de ecuación Espectro Angular"""

#Definiendo término de espectro angular en máscara difractiva

## Definiendo término correspondiente a "campo entrada"

#### NOTA: (campo entrada) es resultado de la interacción del campo incidente con la transmitancia
#### de la máscara difractiva
campo_incidente = 1 #Se considera por facilidad esta característica en campo incidente 

campo_entrada = campo_incidente*mascara_circular

## Definiendo termino correspondiente a transformada de Fourier del campo de entrada
transformada_campo = calculoDFT2(campo_entrada)
transformada_campo_centrada = funcionShift_DFT(transformada_campo)

## ---- ##
espectro_angular_plano_mascara = ((delta_muestreo[0])**2)*transformada_campo_centrada

#Definiendo término de espectro angular en plano de medición

##Definiendo término de condición de evanescencia

### Definiendo término exponente que determina condición de evanescencia
cond_evanescencia = (np.sqrt(1-((longitud_onda_input**2)*(((xx_espectro)**2)+((yy_espectro)**2)))))

termino_cond_evanescencia = np.exp(numero_im*z_input*numero_onda*cond_evanescencia)

## ---- ##
espectro_angular_plano_medicion = termino_cond_evanescencia*espectro_angular_plano_mascara

#Definiendo término de transformada de Fourier inversa del espectro angular en el plano de medición
transformada_inversa_espectro_angular = calculoDFT2_inversa(espectro_angular_plano_medicion)

#Definiendo resultado de campo difractado en plano de medición
campo_difractado = ((delta_muestreo[1])**2)*transformada_inversa_espectro_angular

""" Calculando irradiancia del campo difractado """
irradiancia_campo_difractado = ((np.abs(campo_difractado))**2)

""" Graficando el resultado de difracción por espectro angular """

#GRAFICANDO MASCARA DIFRACTIVA DE ABERTURA CIRCULAR
plt.imshow(mascara_circular, extent=[-longitud_arreglo_input/2, longitud_arreglo_input/2, -longitud_arreglo_input/2, longitud_arreglo_input/2], cmap='gray')
plt.title("Máscara Circular")
plt.colorbar(label="Transmitancia")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()

#Patrón de difracción de Fresnel
plt.imshow(irradiancia_campo_difractado, extent=[-longitud_arreglo_medicion/2, longitud_arreglo_medicion/2, -longitud_arreglo_medicion/2, longitud_arreglo_medicion/2], cmap='gray', vmax= 1*(np.max(irradiancia_campo_difractado)))  # log1p para mejor visualización
plt.title("Patrón de Difracción")
plt.colorbar(label="Amplitud²")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()