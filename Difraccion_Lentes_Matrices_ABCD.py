""" DESCRIPCIÓN GENERAL DEL PROYECTO:
El presente proyecto tiene el objetivo de desarrollar un código base de funcionamiento para 
el proceso de difracción usando lentes a partir de matrices ABCD --> El arreglo a analizar
consta de un plano del objeto, una lente convergente delgada (y también un diafragma para delimitar
el trabajo en un dominio finito) y un plano de medición asociado 
a la imagen resultante. """



""" Se genera un print de mensaje inicial para verificar correcto funcionamiento del entorno """
print("Inicializando entorno de programación Matrices ABCD transferencia de rayos...")


""" Anotaciones importantes magnificación de imágenes """
#NOTA SOBRE CÁLCULO DISTANCIA IMAGEN Y OBJETO: (1/f) = (1/I) + (1/O)

#NOTA SOBRE MAGNIFICACIÓN (para calcular tamaño aproximado de imagen en función de tamaño objeto):

### M = I/O --> Tamaño de imagen será: (TAMAÑO OBJETO)*M



""" Importando librerias y documentos necesarios """

import LIBRERIA_Mascaras_Transmitancia as mascaras
import Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import matplotlib.pyplot as plt
import LIBRERIA_Funciones_importantes as function
from PIL import Image


""" Definiendo parámetros de máscara difractiva """

resolucion_Input = 3000  # Número de puntos en la malla --> Asociado a comparación con una referencia de cámara
longitud_ArregloInput = 1  #Tamaño físico del área de la ventana
radio = 0.05  # Radio del círculo 
centro = None  # El centro será el origen si es None



""" Definiendo parámetro para el tamaño del diafragma """

radio_diafragmaInput = 0.3 #Se define variable asociada al radio de la abertura circular que representará
                           # el diafragma.



""" Definición de distancias del arreglo """

distancia_focal = 0.07   

distancia_propagacionAribitraria = 0.3 #Se define una distancia de propagación arbitraria

distancia_imagen = 0.09 #Se calcula haciendo uso de la ecuación de distancias para formación de imágenes



""" Definiendo parámetros de fuente """

#Definiendo la longitud de onda asociada a la fuente en consideración
longitud_onda_input = 533E-9 #UNIDADES: m

#Calculando el número de onda
numero_onda_input = (2*np.pi)/longitud_onda_input



""" Creando máscara de transmitancia asociada a abertura circular """

# Crear la malla de puntos
xx_mascara, yy_mascara = mascaras.malla_Puntos(resolucion_Input, longitud_ArregloInput)

# Crear la máscara circular
#mascara = mascaras.funcion_Circulo(radio, centro, xx_mascara, yy_mascara)
#mascara = mascaras.funcion_Rectangulo(radio,radio,centro,xx_mascara,yy_mascara)
#mascara = mascaras.funcion_Corazon(centro,xx_mascara,yy_mascara,radio)

# Cargar la imagen PNG como máscara de transmitancia
ruta_imagen_png = "/home/labravo/Downloads/Ruido_E03.png"  # Especifica la ruta de tu imagen
mascara = function.cargar_imagen_png(ruta_imagen_png, resolucion_Input, longitud_ArregloInput)

""" Se calculan las matrices necesarias para estudiar el PRIMER TRAMO del arreglo difractivo
    OBJETO --> *Propagación* --> LENTE """

#Se calcula la matriz asociada al proceso de propagación desde el plano objeto
#hasta el plano de la lente
matriz_propagacionPrimerTramo = matriz.propagacion_MedioHomogeneo(distancia_propagacionAribitraria)



""" Calculando la matriz del sistema asociada al PRIMER TRAMO """

#NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
#en el arreglo.

#Se define la lista de matrices
lista_matricesPrimerTramoInvertida = [matriz_propagacionPrimerTramo]

#En este caso la matriz del sistema es equivalente a la única matriz presente
matriz_SistemaPrimerTramo = matriz.matriz_Sistema(lista_matricesPrimerTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del PRIMER TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del sistema
camino_opticoCentralPrimerTramo = matriz.camino_Optico(lista_matricesPrimerTramoInvertida)



""" Creando malla de puntos plano lente (en este PRIMER TRAMO se asocia al plano de medición)"""

#Se llama función para determinar los deltas de muestreo del tramo
deltas_tramoObjetoLente = function.producto_espacio_frecuencia_TransformadaFresnel(longitud_onda_input,
                                                                                   matriz_SistemaPrimerTramo[0,1],
                                                                                   resolucion_Input,
                                                                                   longitud_ArregloInput)

#Se calcula el ancho de la ventana del plano de la lente
ancho_VentanaPlanoLente = resolucion_Input*deltas_tramoObjetoLente[1]

#Se calcula la malla de puntos asociada al plano de la lente
xx_PlanoLente, yy_PlanoLente = mascaras.malla_Puntos(resolucion_Input, ancho_VentanaPlanoLente)



""" Se calcula el resultado del proceso difractivo del PRIMER TRAMO"""

#Se calcula el campo de salida/en plano de la lente --> Campo resultante del primer tramo 
campo_PlanoLente = matriz.matriz_ABCD_Difraccion(camino_opticoCentralPrimerTramo,mascara,
                                                 matriz_SistemaPrimerTramo[0,0],
                                                 matriz_SistemaPrimerTramo[0,1],
                                                 matriz_SistemaPrimerTramo[1,1],xx_mascara,yy_mascara,
                                                 xx_PlanoLente,yy_PlanoLente,numero_onda_input,
                                                 deltas_tramoObjetoLente)



""" Se crea diafragma para delimitar dominio del arreglo asociado a lentes FINITAS
NOTA: El diafragma va a construirse a partir de la malla de puntos asociada al plano de la Lente."""

#Creación de máscara circular que representará el diafragma 
diafragma = mascaras.funcion_Circulo(radio_diafragmaInput, centro, xx_PlanoLente, yy_PlanoLente)



""" Se calcula el campo de entrada para el SEGUNDO TRAMO del arreglo 
NOTA: El campo que llega a la lente e interactua con la máscara asociada al diafragma
será el campo de entrada para el segundo tramo."""

#Calculando el campo de entrada al SEGUNDO TRAMO del arreglo
campo_entradaSegundoTramo = campo_PlanoLente*diafragma



""" Se calculan las matrices necesarias para estudiar el SEGUNDO TRAMO del arreglo difractivo
    LENTE --> *Efecto Lente* --> *Propagación* --> MEDICIÓN """

#Se calcula la matriz asociada a la interacción con la lente
matriz_lente = matriz.lente_DelgadaConociendoDistanciaFocal(distancia_focal)

#Se calcula la matriz asociada al proceso de propagación desde el plano de la lente
#hasta el plano de medición
matriz_propagacionSegundoTramo = matriz.propagacion_MedioHomogeneo(distancia_imagen)



""" Calculando la matriz del sistema asociada al SEGUNDO TRAMO 
NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ubicación real
en el arreglo."""

#Creando lista de matrices que describe el arreglo del SEGUNDO TRAMO
lista_matricesSegundoTramoInvertida = [matriz_propagacionSegundoTramo,matriz_lente]

#Se calcula la matriz del sistema
matriz_SistemaSegundoTramo = matriz.matriz_Sistema(lista_matricesSegundoTramoInvertida)



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del SEGUNDO TRAMO """

#Se calcula el camino óptico central a partir de la lista de matrices del Segundo tramo
camino_opticoCentralSegundoTramo = matriz.camino_Optico(lista_matricesSegundoTramoInvertida)



""" Calculando la malla de puntos del plano de medición """

#Se llama función para determinar los deltas de muestreo asociados al SEGUNDO TRAMO
deltas_tramoLenteMedicion = function.producto_espacio_frecuencia_TransformadaFresnel(longitud_onda_input,
                                                                                     matriz_SistemaSegundoTramo[0,1],
                                                                                     resolucion_Input,
                                                                                     ancho_VentanaPlanoLente)

#Se calcula el ancho de la ventana del plano de medición 
ancho_VentanaPlanoMedicion = resolucion_Input*deltas_tramoLenteMedicion[1]

#Se calcula la malla de puntos asociada al plano de medición
xx_PlanoMedicion, yy_PlanoMedicion = mascaras.malla_Puntos(resolucion_Input, ancho_VentanaPlanoMedicion)


""" Se calcula el resultado del proceso difractivo del SEGUNDO TRAMO """

#Se calcula el campo de salida/en plano de medición --> Campo resultante de la difracción 
campo_PlanoMedicion = matriz.matriz_ABCD_Difraccion_Shift(camino_opticoCentralSegundoTramo,
                                                    campo_entradaSegundoTramo,
                                                    matriz_SistemaSegundoTramo[0,0],
                                                    matriz_SistemaSegundoTramo[0,1],
                                                    matriz_SistemaSegundoTramo[1,1],xx_PlanoLente,
                                                    yy_PlanoLente,xx_PlanoMedicion,yy_PlanoMedicion,
                                                    numero_onda_input,deltas_tramoLenteMedicion)

#Se calcula la amplitud del campo de salida
amplitud_campoPlanoMedicion = np.abs(campo_PlanoMedicion)

#Se calcula la intensidad del campo de salida
intensidad_campoPlanoMedicion = amplitud_campoPlanoMedicion**2



""" Graficando máscara de transmitancia asignada a abertura circular"""

plt.imshow(mascara, extent=[-longitud_ArregloInput/2, longitud_ArregloInput/2,
                             -longitud_ArregloInput/2, longitud_ArregloInput/2], 
                             cmap='gray')
plt.title("Máscara")
plt.colorbar(label="Transmitancia")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()



""" Graficando la intensidad del campo de salida """

plt.imshow(intensidad_campoPlanoMedicion, extent=[-ancho_VentanaPlanoMedicion/2, 
                                                  ancho_VentanaPlanoMedicion/2, 
                                                  -ancho_VentanaPlanoMedicion/2, 
                                                  ancho_VentanaPlanoMedicion/2], 
           cmap='gray')
plt.title("Intensidad")
plt.colorbar(label="Intensidad")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()

