""" Se genera un print de mensaje inicial para verificar correcto funcionamiento del entorno """
print("Inicializando entorno de programación implementación 01 Matrices ABCD transferencia de rayos...")



""" Importando librerias y documentos necesarios """

import Mascaras_Transmitancia as mascaras
import Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import matplotlib.pyplot as plt
import Funciones_importantes as function



""" Definiendo parámetros de máscara difractiva """
resolucion = 3000  # Número de puntos en la malla
longitud_Arreglo = 0.2  # Tamaño físico del área 
radio = 0.02  # Radio del círculo 
centro = None  # El centro será el origen si es None


""" Definición de distancias del arreglo """
distancia_focal01 = 0.1

distancia_focal02 = 0.15

distancia_propagacion_d = 0.05



""" Definiendo parámetros de fuente """

#Definiendo la longitud de onda asociada a la fuente en consideración
longitud_onda_input = 533E-9 #UNIDADES: m

#Calculando el número de onda
numero_onda_input = (2*np.pi)/longitud_onda_input



""" Creando máscara de transmitancia asociada a abertura circular """

# Crear la malla de puntos
xx_mascara, yy_mascara = mascaras.malla_Puntos(resolucion, longitud_Arreglo)

# Crear la máscara circular
mascara = mascaras.funcion_Circulo(radio, centro, xx_mascara, yy_mascara)



""" Graficando máscara de transmitancia asignaada a abertura circular"""

plt.imshow(mascara, extent=[-longitud_Arreglo/2, longitud_Arreglo/2, -longitud_Arreglo/2, longitud_Arreglo/2], cmap='gray')
plt.title("Máscara Circular")
plt.colorbar(label="Transmitancia")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()



""" Se calculan las matrices necesarias para estudiar el arreglo difractivo """

#Primer tramo --> Propagación una distancia focal 01
matriz_propagacion01 = matriz.propagacion_MedioHomogeneo(distancia_focal01)

#Segundo tramo --> Paso a través de una lente delgada 01
matriz_lente01 = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focal01)

#Tercer tramo --> Propagación una distancia focal 01
matriz_propagacion02 = matriz.propagacion_MedioHomogeneo(distancia_focal01)

#Cuarto tramos --> Propagación una distancia genérica "d"
matriz_propagacion03 = matriz.propagacion_MedioHomogeneo(distancia_propagacion_d)

#Quinto tramo --> Paso a través de una lente delgada 02
matriz_lente02 = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focal02)

#Sexto tramo --> Propagación a una distancia focal 02
matriz_propagacion04 = matriz.propagacion_MedioHomogeneo(distancia_focal02)



""" Calculando la matriz del sistema """

#Se crea una lista de matrices para abordar todo el proceso difractivo
lista_matricesABCD = [matriz_propagacion01,matriz_lente01,matriz_propagacion02,matriz_propagacion03,
                      matriz_lente02,matriz_propagacion04]

#Se llama función para calcular la matriz del sistema
matriz_Sistema = matriz.matriz_Sistema(lista_matricesABCD)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL """

camino_optico_central = matriz.camino_Optico(lista_matricesABCD)



""" Calculando la malla de puntos del plano de medición """

#Se llama función para determinar los deltas de muestreo
deltas = function.producto_espacio_frecuencia_TransformadaFresnel(longitud_onda_input,camino_optico_central,
                                                                  resolucion,longitud_Arreglo)

#Se calcula el ancho de la ventana del plano de medición 
ancho_VentanaPlanoMedicion = resolucion*deltas[1]

#Se calcula la malla de puntos asociada al plano de medición
xx_PlanoMedicion, yy_PlanoMedicion = mascaras.malla_Puntos(resolucion, ancho_VentanaPlanoMedicion)



""" Se calcula el resultado del proceso difractivo """

#Se calcula el campo de salida/en plano de medición --> Campo resultante de la difracción 
campo_PlanoMedicion = matriz.matriz_ABCD_Difraccion(camino_optico_central,mascara,matriz_Sistema[0,0],
                                                    matriz_Sistema[0,1],matriz_Sistema[1,1],xx_mascara,
                                                    yy_mascara,xx_PlanoMedicion,yy_PlanoMedicion,numero_onda_input)

#Se calcula la amplitud del campo de salida
amplitud_campoPlanoMedicion = np.abs(campo_PlanoMedicion)

#Se calcula la intensidad del campo de salida
intensidad_campoPlanoMedicion = amplitud_campoPlanoMedicion**2



""" Graficando la intensidad del campo de salida """

plt.imshow(intensidad_campoPlanoMedicion, extent=[-ancho_VentanaPlanoMedicion/2, ancho_VentanaPlanoMedicion/2, -ancho_VentanaPlanoMedicion/2, ancho_VentanaPlanoMedicion/2], 
           cmap='inferno')
plt.title("Intensidad")

plt.colorbar(label="Intensidad")

plt.xlabel("X (m)")
plt.ylabel("Y (m)")

plt.show()