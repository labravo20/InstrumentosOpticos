import numpy as np


""" Función producto espacio frecuencia en caso TRANSFORMADA FRESNEL """

def producto_espacio_frecuencia_TransformadaFresnel(longitud_onda,distancia_Propagacion,resolucion,longitud_Arreglo):

    """
    Función para calcular los deltas de muestreo asociadas a las condiciones de Transformada de Fresnel.

    FUNCIÓN RECIBE:
    
    - Longitud de onda de la fuente lumínica utilizada
    - Distancia de propagación
    - Resolución del input 
    - Ancho de la ventana o longitud de arreglo de entrada

    FUNCIÓN RETORNA: Lista con los deltas de muestreo

    """

    #Definiendo separación entre número total de muestras --> En plano de entrada
    delta_input = longitud_Arreglo/resolucion 

    #Definiendo separación entre muestras en plano de medición
    delta_output = (longitud_onda*distancia_Propagacion)/(resolucion*delta_input) 

    #Creando una lista que guarde información de los deltas de muestreo asociados al plano de entrada
    #y de salida del arreglo
    delta_muestreo = [delta_input,delta_output]

    return delta_muestreo


""" Función producto espacio frecuencia en caso TRANSFORMADA FRESNEL para aplicación en SENSORES """

def producto_espacio_frecuencia_TransformadaFresnel_Sensor(resolucion_AnchoSensor, ancho_Sensor, distancia_Propagacion, longitud_Onda, resolucion_AltoSensor = None, alto_Sensor=None):

    ''' 
    Función para calcular los deltas de muestreo asociadas a las condiciones de Transformada de Fresnel
    teniendo en cuenta características específicas de un sensor.

    FUNCIÓN RECIBE:
    
    -Resolucion ancho sensor = Cantidad de muestras/puntos en el ancho del sensor
    -Ancho sensor = Longitud física del ancho del sensor
    -Distancia de propagación entre plano de entrada y plano de salida
    -Longitud de onda de la fuente lumínica utilizada
    -Resolucion alto sensor (OPCIONAL) = Cantidad de muestras/puntos en el alto del sensor
    -Alto sensor (OPCIONAL) = Longitud física del alto del sensor

    FUNCIÓN RETORNA: Lista con los deltas de muestreo del plano de entrada.

    '''


    #Se verifica si se ingresó información sobre el alto del sensor
    if alto_Sensor is None: 

        #En caso de no ser ingresado se establece por default que sea igual al ancho del sensor
        alto_Sensor = ancho_Sensor 

    #Se verifica si se ingresó información sobre la resolucion en el alto
    if resolucion_AltoSensor is None: 
        
        #En caso de no ser ingresado se establece por default que sea igual a la resolucion del ancho del sensor
        resolucion_AltoSensor = resolucion_AnchoSensor
        

    #Se calculan los deltas asociados a las características y en el plano del sensor
    delta_XSensor = ancho_Sensor / resolucion_AnchoSensor
    delta_YSensor = alto_Sensor / resolucion_AltoSensor

    #Usando los deltas del sensor se calculan los deltas asociados al plano de entrada
    delta_EntradaX = (longitud_Onda * distancia_Propagacion) / (resolucion_AnchoSensor * delta_XSensor)
    delta_EntradaY = (longitud_Onda * distancia_Propagacion) / (resolucion_AltoSensor * delta_YSensor)

    #Se crea una lista donde se almacena la información sobre los deltas de muestreo en la entrada
    muestreo_Entrada = [delta_EntradaX, delta_EntradaY]

    return muestreo_Entrada
    