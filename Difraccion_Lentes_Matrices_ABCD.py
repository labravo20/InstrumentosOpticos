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

import Mascaras_Transmitancia as mascaras
import Matrices_ABCD_Transferencia_rayos as matriz
import numpy as np
import matplotlib.pyplot as plt
import Funciones_importantes as function



""" Definiendo parámetros de máscara difractiva """

resolucion_Input = 4012  # Número de puntos en la malla --> Asociado a comparación con una referencia de cámara
longitud_ArregloInput = 0.0043   #Tamaño físico del área ---> Resultado de calcular con datos de cámara... 
                        # --> Imagen con tamaño de 0.00775 (longitud arreglo) se usa en condiciones de 
                        # producto espacio frecuencia de transformada de Fresnel
radio = 0.003  # Radio del círculo 
centro = None  # El centro será el origen si es None



""" Definición de distancias del arreglo """
distancia_focal01 = 0.07  #DISTANCIA FOCAL MÁXIMA 0.125 

distancia_focal02 = 0.15

distancia_propagacion_d = 0.315



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
mascara = mascaras.funcion_Rectangulo(radio,radio,centro,xx_mascara,yy_mascara)



""" Se calculan las matrices necesarias para estudiar el PRIMER TRAMO del arreglo difractivo
    OBJETO --> *Propagación* --> LENTE """

#Se calcula la matriz asociada al proceso de propagación desde el plano objeto
#hasta el plano de la lente
matriz_propagacionPrimerTramo = matriz.propagacion_MedioHomogeneo(distancia_focal01)



""" Calculando la matriz del sistema asociada al PRIMER TRAMO """

#NOTA IMPORTANTE: La lista de matrices se debe poner en orden inverso a su ibicación real
#en el arreglo.

#En este caso la matriz del sistema es equivalente a la única matriz presente
matriz_SistemaPrimerTramo = matriz_propagacionPrimerTramo



""" Calculando el camino óptico central --> Asociado a la distancia de propagación TOTAL del PRIMER TRAMO """

#En este caso la lista de la cual se calcula el camino óptico central se conforma por la única matriz
#considerada para este tramo
camino_opticoCentralPrimerTramo = matriz.camino_Optico(matriz_propagacionPrimerTramo)



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



""" Se calcula el resultado del proceso difractivo """

#Se calcula el campo de salida/en plano de la lente --> Campo resultante del primer tramo 
campo_PlanoLente = matriz.matriz_ABCD_Difraccion(camino_opticoCentralPrimerTramo,mascara,
                                                 matriz_SistemaPrimerTramo[0,0],
                                                 matriz_SistemaPrimerTramo[0,1],
                                                 matriz_SistemaPrimerTramo[1,1],xx_mascara,yy_mascara,
                                                 xx_PlanoLente,yy_PlanoLente,numero_onda_input,
                                                 deltas_tramoObjetoLente)


""" Se agrega diafragma para delimitar dominio del arreglo asociado a lentes FINITAS
NOTA: El campo que llega a la lente e interactua con la máscara asociada al diafragma
será el campo de entrada para el segundo tramo. """

#Creación de máscara circular que representará el diafragma 
diafragma = mascaras.funcion_Circulo(radio_diafragmaInput, centro, xx_diafragma, yy_diafragma)





""" Se calcula el resultado del proceso difractivo """

#Se calcula el campo de salida/en plano de medición --> Campo resultante de la difracción 
campo_PlanoMedicion = matriz.matriz_ABCD_Difraccion(camino_optico_central,mascara,matriz_Sistema[0,0],
                                                    matriz_Sistema[0,1],matriz_Sistema[1,1],xx_mascara,
                                                    yy_mascara,xx_PlanoMedicion,yy_PlanoMedicion,numero_onda_input,
                                                    deltas)


#Se calcula la amplitud del campo de salida
amplitud_campoPlanoMedicion = np.abs(campo_PlanoMedicion)

#Se calcula la intensidad del campo de salida
intensidad_campoPlanoMedicion = amplitud_campoPlanoMedicion**2



""" Calculando la malla de puntos del plano de medición """

#Se llama función para determinar los deltas de muestreo
deltas = function.producto_espacio_frecuencia_TransformadaFresnel(longitud_onda_input,matriz_Sistema[0,1],
                                                                  resolucion,longitud_Arreglo)



#Se calcula el ancho de la ventana del plano de medición 
ancho_VentanaPlanoMedicion = resolucion*deltas[1]


#Se calcula la malla de puntos asociada al plano de medición
xx_PlanoMedicion, yy_PlanoMedicion = mascaras.malla_Puntos(resolucion, ancho_VentanaPlanoMedicion)



""" Graficando máscara de transmitancia asignaada a abertura circular"""

plt.imshow(mascara, extent=[-longitud_Arreglo/2, longitud_Arreglo/2, -longitud_Arreglo/2, longitud_Arreglo/2], cmap='gray')
plt.title("Máscara Circular")
plt.colorbar(label="Transmitancia")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()



""" Graficando la intensidad del campo de salida """

plt.imshow(intensidad_campoPlanoMedicion, extent=[-ancho_VentanaPlanoMedicion/2, ancho_VentanaPlanoMedicion/2, -ancho_VentanaPlanoMedicion/2, ancho_VentanaPlanoMedicion/2], 
           cmap='gray')
plt.title("Intensidad")

plt.colorbar(label="Intensidad")

plt.xlabel("X (m)")
plt.ylabel("Y (m)")

plt.show()