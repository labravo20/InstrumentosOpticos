import numpy as np


""" Función producto espacio frecuencia en caso TRANSFORMADA FRESNEL """

def producto_espacio_frecuencia(longitud_onda,z,resolucion,longitud_Arreglo):

    delta_input = longitud_Arreglo/resolucion #Definiendo separación entre número total de muestras

    delta_output = (longitud_onda*z)/(resolucion*delta_input) #Definiendo separación entre muestras
                                                            #en plano de medición

    delta_muestreo = [delta_input,delta_output]

    return delta_muestreo
