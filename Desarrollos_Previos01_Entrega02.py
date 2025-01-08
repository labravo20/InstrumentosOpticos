""" Se genera un print de mensaje inicial para verificar correcto funcionamiento del entorno """
print("Inicializando entorno de programación implementación 01 Matrices ABCD transferencia de rayos...")



""" Importando librerias y documentos necesarios """

import Mascaras_Transmitancia as mascaras
import Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import matplotlib.pyplot as plt



""" Definiendo parámetros de máscara difractiva """
resolucion = 3000  # Número de puntos en la malla
longitud_Arreglo = 0.5  # Tamaño físico del área 
radio = 0.002  # Radio del círculo 
centro = None  # El centro será el origen si es None


""" Definición de distancias del arreglo """
distancia_focal01 = 10
distancia_focal02 = 10
distancia_propagacion = 10



""" Creando máscara de transmitancia asociada a abertura circular """

# Crear la malla de puntos
xx, yy = mascaras.malla_Puntos(resolucion, longitud_Arreglo)

# Crear la máscara circular
mascara = mascaras.funcion_Circulo(radio, centro, xx, yy)




""" Graficando máscara de transmitancia asignaada a abertura circular"""

plt.imshow(mascara, extent=[-longitud_Arreglo/2, longitud_Arreglo/2, -longitud_Arreglo/2, longitud_Arreglo/2], cmap='gray')
plt.title("Máscara Circular")
plt.colorbar(label="Transmitancia")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()




