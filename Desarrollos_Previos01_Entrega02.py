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
distancia_propagacion_d = 10



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

#CALCULANDO LA MATRIZ DEL SISTEMA (tiene en cuenta contribuciones de cada tramo del proceso)

### Se crea una lista de matrices para abordar todo el proceso difractivo
lista_matricesABCD = [matriz_propagacion01,matriz_lente01,matriz_propagacion02,matriz_propagacion03,
                      matriz_lente02,matriz_propagacion04]