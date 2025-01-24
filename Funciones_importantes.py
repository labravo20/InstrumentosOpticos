import numpy as np
from PIL import Image
import pandas as pd

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
    

""" Función para cargar la imagen PNG y ajustarla al tamaño de la malla """

def cargar_imagen_png(ruta_imagen, resolucion_ancho, resolucion_alto = None):
    """
    Carga una imagen PNG y la ajusta a las dimensiones de la malla.

    Args:
        ruta_imagen (str): Ruta de la imagen PNG.
        resolucion (int): Resolución de la malla (número de puntos).
        longitud_ventana (float): Tamaño físico del área de la ventana.

    Returns:
        numpy.ndarray: Imagen escalada a la resolución de la malla.
    """

    #Inicialmente se verifica si si incluyó información sobre la resolución del alto de la imagen
    if resolucion_alto == None:

        #En caso de NO recibir información se define por caso default una resolución cuadrada
        resolucion_alto = resolucion_ancho

    # Cargar la imagen en escala de grises
    imagen = Image.open(ruta_imagen).convert("L")
    
    # Redimensionar la imagen a la resolución deseada
    imagen = imagen.resize((resolucion_ancho, resolucion_alto), Image.Resampling.LANCZOS)
    
    # Normalizar los valores a un rango de 0 a 1
    imagen_normalizada = np.array(imagen) / 255.0
    
    return imagen_normalizada


""" Función para cargar documento tipo CSV """

def cargar_documento_csv(ruta_documento_csv):

    # Cargar el archivo CSV en un DataFrame
    df = pd.read_csv(ruta_documento_csv)
    
    # Si el CSV tiene datos en columnas o filas que representan la malla de puntos, convertir a un array de NumPy
    # Suponiendo que los datos estén en una forma 2D en el archivo CSV (como una imagen en escala de grises)
    mascara = df.to_numpy()
    
    # Asegúrate de que la máscara esté en la forma correcta (2D)
    dimensiones_documento = mascara.shape  # Verifica las dimensiones

    #impresión para verificación
    #print(dimensiones_documento)

ruta = "/home/labravo/Downloads/MuestraBio_E03.csv"

cargar_documento_csv(ruta)
