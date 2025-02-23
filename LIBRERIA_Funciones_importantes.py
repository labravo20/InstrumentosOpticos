import numpy as np
from PIL import Image
from scipy.ndimage import zoom

""" Función producto espacio frecuencia ESPECTRO ANGULAR """

def producto_espacio_frecuencia(resolucion,longitud_Arreglo):

    delta_input = longitud_Arreglo/resolucion #Definiendo separación entre número total de muestras

    delta_output = (1)/(resolucion*delta_input) #Definiendo separación entre muestras
                                                            #en plano de medición

    delta_muestreo = [delta_input,delta_output]

    return delta_muestreo

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
    muestreo_Entrada = {"deltaPlanoEntrada_X":delta_EntradaX, "deltaPlanoEntrada_Y":delta_EntradaY}

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

def cargar_imagen_png_recortada(ruta_imagen, resolucion_ancho, resolucion_alto=None):
    """
    Carga una imagen PNG de cualquier tamaño, la redimensiona a 3899x3899,
    recorta el área central de 2848x2848 y la ajusta a la resolución deseada.

    Args:
        ruta_imagen (str): Ruta de la imagen PNG.
        resolucion_ancho (int): Resolución final deseada en ancho.
        resolucion_alto (int, opcional): Resolución final deseada en alto. 
                                         Si no se especifica, se usa una imagen cuadrada.

    Returns:
        numpy.ndarray: Imagen procesada (recortada y escalada), normalizada entre 0 y 1.
    """
    if resolucion_alto is None:
        resolucion_alto = resolucion_ancho

    # Cargar la imagen en escala de grises
    imagen = Image.open(ruta_imagen).convert("L")

    # Redimensionar la imagen a 3899x3899
    imagen_redimensionada = imagen.resize((3262, 3262), Image.Resampling.LANCZOS)

    # Definir dimensiones del recorte (2848x2848 centrado)
    ancho_recorte = 2848
    alto_recorte = 2848
    x_inicio = (3262 - ancho_recorte) // 2
    y_inicio = (3262 - alto_recorte) // 2
    x_fin = x_inicio + ancho_recorte
    y_fin = y_inicio + alto_recorte

    # Recortar la imagen
    imagen_recortada = imagen_redimensionada.crop((x_inicio, y_inicio, x_fin, y_fin))

    # Redimensionar la imagen a la resolución deseada
    imagen_final = imagen_recortada.resize((resolucion_ancho, resolucion_alto), Image.Resampling.LANCZOS)

    # Normalizar los valores a un rango de 0 a 1
    imagen_normalizada = np.array(imagen_final) / 255.0

    return imagen_normalizada


""" Función para cargar documento tipo CSV """

def cargar_documento_csv(ruta_csv,resolucion_anchoSensorInput,resolucion_altoSensorInput):

    # Leer el archivo y reemplazar 'i' con 'j'
    with open(ruta_csv, "r") as file:
        contenido = file.read().replace("i", "j")
    
    # Guardar el archivo actualizado
    with open(ruta_csv, "w") as file:
        file.write(contenido)
    
    datos_csv = np.genfromtxt(ruta_csv, delimiter=',', dtype=complex)
    resolucion_x_csv, resolucion_y_csv = datos_csv.shape
    
    # Factores de escala para ajustar a las dimensiones del plano de la máscara
    factor_x = resolucion_anchoSensorInput / resolucion_x_csv
    factor_y = resolucion_altoSensorInput / resolucion_y_csv
    
    #Interpolar el CSV para ajustar su tamaño
    datos_csv_ajustados = zoom(datos_csv, (factor_y, factor_x), order=1)  

    return datos_csv_ajustados


""" Función para cargar documento tipo CSV usando padding --> NO INTERPOLACIÓN"""

def cargar_documento_csv_OPTION02(archivo_csv):

    #Se define las condiciones de muestreo deseadas
    dimensiones_resolucion = (2048, 2048)
    
    #Anexando los datos asociados al archivo CSV
    datos_CSV = np.loadtxt(archivo_csv, delimiter=',', dtype=str)

    print(np.shape(datos_CSV))

    #Para buena interpretación de la información de cambian los términos apropiados para asignar 
    #al término imaginario
    datosCSV_compleja = np.vectorize(lambda x: complex(x.replace(' ', '').replace('i', 'j')))(datos_CSV)

    #Se determina el ancho (COORDENADA X)
    ancho_pad = (dimensiones_resolucion[0] - datosCSV_compleja.shape[0]) // 2

    #Se determina el tamaño que se debe aumentar al arreglo original
    ancho_pad_adicional = (dimensiones_resolucion[0] - datosCSV_compleja.shape[0]) % 2

    #Se determina el alto (COORDENADA Y)
    alto_pad = (dimensiones_resolucion[1] - datosCSV_compleja.shape[1]) // 2

    #Se determina el tamaño que se debe aumentar al arreglo original
    alto_pad_adicional = (dimensiones_resolucion[1] - datosCSV_compleja.shape[1]) % 2

    #Se crea matriz con toda la información de edición del arreglo
    array_padded = np.pad(datosCSV_compleja, ((ancho_pad, ancho_pad + ancho_pad_adicional), 
                                              (alto_pad, alto_pad + alto_pad_adicional)),
                                              mode='constant', constant_values=0)

    return array_padded



""" Función para interferir dos ondas planas """
def interferencia_ondas_malla(X, Y, k1, k2, A1=1, A2=1, fase1=0, fase2=0):
    """
    Calcula el patrón de interferencia entre dos ondas planas sobre una malla dada.

    Args:
        X, Y (numpy.ndarray): Malla de coordenadas espaciales.
        k1 (tuple): Vector de onda (kx1, ky1) de la primera onda.
        k2 (tuple): Vector de onda (kx2, ky2) de la segunda onda.
        A1, A2 (float): Amplitudes de las ondas.
        fase1, fase2 (float): Fases iniciales de las ondas en radianes.

    Returns:
        np.ndarray: Patrón de interferencia (intensidad) sobre la malla.
    """
    # Definir las ondas planas con sus respectivos vectores de onda y fases
    onda1 = A1 * np.exp(1j * (k1[0] * X + k1[1] * Y + fase1))
    onda2 = A2 * np.exp(1j * (k2[0] * X + k2[1] * Y + fase2))

    # Sumar las ondas (interferencia)
    interferencia = np.abs(onda1 + onda2) ** 2  # Intensidad = |E_total|^2

    return interferencia
