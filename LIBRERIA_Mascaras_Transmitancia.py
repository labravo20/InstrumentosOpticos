import numpy as np



''' Definicion de función para crear malla de puntos'''

def malla_Puntos(resolucion_Ancho, ancho_Arreglo, resolucion_Alto = None, alto_Arreglo=None):

    ''' Crea mallas de puntos

    FUNCION RECIBE:

        resolucion_Ancho = Cantidad de muestras/puntos en el ancho de la ventana
        ancho_Arreglo    = Longitud física del ancho de la ventana
        resolucion_Alto (OPCIONAL) = Cantidad de muestras/puntos en el alto de la ventana
        alto_Arreglo  (OPCIONAL)   = Longitud física del alto de la ventana
    
    FUNCIÓN RETORNA:
        xx, yy = malla de puntos bidimensional
        
        '''
    

    #Se verifica si se ingresó el alto del arreglo
    if alto_Arreglo is None:

        #En caso de no ser ingresado se establece por default que sea igual al ancho del arreglo
        alto_Arreglo = ancho_Arreglo

    #Se verifica si se ingresó la resolución asociada al alto del arreglo
    if resolucion_Alto is None:

        #En caso de no ser ingresado se establece por default que sea igual a la resolución del ancho del arreglo 
        resolucion_Alto = resolucion_Ancho
    
    
    #Creando malla de puntos considerando condiciones de ancho y alto del arreglo
    x = np.linspace(-ancho_Arreglo / 2, ancho_Arreglo / 2, resolucion_Ancho) 
    y = np.linspace(-alto_Arreglo / 2, alto_Arreglo / 2, resolucion_Alto) 

    #Creación de malla de puntos bidimensional
    xx, yy = np.meshgrid(x, y)  

    return xx, yy 



''' Definicion de función para máscara circular '''

def funcion_CirculoInvertida(radio, centro, xx, yy): 
    '''
      Crea una máscara invertida con un círculo
    
    FUNCIÓN RECIBE:

        radio  == float 
        centro == type(list) --> [X,Y]
        xx, yy == malla de puntos en la cual se verá el circulo
    
    FUNCIÓN RETORNA:  Máscara circular invertida

    '''

    #Se verifica si se especificó la posición del centro del rectángulo
    if centro is None: 

        #En caso de no ser especificada se ubica el centro en el origen por defecto
        centro = (0, 0)

    #Se calcula la distancia desde el centro de la circunferencia a cada punto
    distancia = (xx - centro[0])**2 + (yy - centro[1])**2 

    #Se ubica las dimensiones del circulo dentro de la malla de puntos predeterminada
    mascara_circular = distancia <= radio**2 

    #Se invierte la máscara para poder generar la máscara circular invertido
    mascara_circularInvertida = ~mascara_circular

    return mascara_circularInvertida 



""" Definición de función para crear máscara invertida con distribución de transmitancia gaussiana """

def funcion_CirculoInvertidoGaussian(radio, centro, xx, yy): 
    '''
      Crea una distribución gaussiana circular
    
    FUNCIÓN RECIBE:

        radio  == float (define la desviación estándar de la gaussiana)
        centro == list o tuple [X, Y]
        xx, yy == malla de puntos en la cual se verá la distribución gaussiana
    
    FUNCIÓN RETORNA:  Máscara circular con distribución gaussiana

    '''
    
    # Se verifica si se especificó la posición del centro
    if centro is None: 
        # En caso de no ser especificada, se ubica el centro en el origen por defecto
        centro = (0, 0)

    # Calculamos la distancia al cuadrado desde el centro (ecuación gaussiana)
    distancia_cuadrada = (xx - centro[0])**2 + (yy - centro[1])**2

    # Calculamos la desviación estándar (sigma) a partir del radio
    sigma = radio / 2  # Relación ajustable entre radio y sigma

    # Generamos la distribución gaussiana
    mascara_gaussiana = np.exp(-distancia_cuadrada / (2 * sigma**2))

    #Invirtiendo el valor de la máscara para que sea una transmitancia invertida
    mascara_gaussianaInvertida = 1 - mascara_gaussiana

    #Ajustando paŕametros de gaussiana para que el valor mínimo sea mayor a cero
    mascara_gaussianaAjustada = (mascara_gaussianaInvertida*0.7) + 0.3

    return mascara_gaussianaAjustada



''' Definicion de función para máscara circular '''

def funcion_Circulo(radio, centro, xx, yy): 
    '''
      Crea una máscara con un círculo
    
    FUNCIÓN RECIBE:

        radio  == float 
        centro == type(list) --> [X,Y]
        xx, yy == malla de puntos en la cual se verá el circulo
    
    FUNCIÓN RETORNA:  Máscara circular

    '''

    #Se verifica si se especificó la posición del centro del rectángulo
    if centro is None: 

        #En caso de no ser especificada se ubica el centro en el origen por defecto
        centro = (0, 0)

    #Se calcula la distancia desde el centro de la circunferencia a cada punto
    distancia = (xx - centro[0])**2 + (yy - centro[1])**2 

    #Se ubica las dimensiones del circulo dentro de la malla de puntos predeterminada
    mascara_circular = distancia <= radio**2 

    return mascara_circular 



""" Definición de función para máscara rectangular """

def funcion_Rectangulo(base, altura, centro, xx, yy): 
    '''
    Crea una máscara con un rectángulo 
    
    FUNCIÓN RECIBE:

        base   == float 
        altura == float 
        centro == type(list) --> [X,Y]
        xx, yy == malla de puntos de salida en la cuál se encontrará el rectángulo
    
    FUNCIÓN RETORNA:  Máscara rectangular
    '''
    
    #Se verifica si se especificó la posición del centro del rectángulo
    if centro is None:

        #En caso de no ser especificada se ubica el centro en el origen por defecto 
        centro = [0, 0] 

    #A continuación se calculan las distancias asociadas a los límites del rectangulo
    x_Min = centro[0] - base / 2
    x_Max = centro[0] + base / 2
    y_Min = centro[1] - altura / 2
    y_Max = centro[1] + altura / 2

    #Se ubica las dimensiones del rectángulo dentro de la malla de puntos predeterminada
    mascara_rectangular = (xx <= x_Max) & (xx >= x_Min) & (yy <= y_Max) & (yy >= y_Min)  

    return mascara_rectangular
 

""" Definición de función para máscara en forma de cruz """

def funcion_Cruz(largo, grosor, centro, xx, yy):
    '''
    Crea una máscara con una cruz
    
    FUNCIÓN RECIBE:
        largo   == float  (longitud de cada brazo de la cruz)
        grosor  == float  (ancho de cada brazo de la cruz)
        centro  == type(list) --> [X, Y] (centro de la cruz)
        xx, yy  == malla de puntos de salida en la cuál se encontrará la cruz
    
    FUNCIÓN RETORNA: Máscara en forma de cruz
    '''
    
    # Se verifica si se especificó la posición del centro de la cruz
    if centro is None:
        centro = [0, 0]  # Se coloca en el origen por defecto
    
    # Definir los límites de los brazos de la cruz
    x_Min_V = centro[0] - grosor / 2  # Vertical
    x_Max_V = centro[0] + grosor / 2
    y_Min_V = centro[1] - largo / 2
    y_Max_V = centro[1] + largo / 2
    
    y_Min_H = centro[1] - grosor / 2  # Horizontal
    y_Max_H = centro[1] + grosor / 2
    x_Min_H = centro[0] - largo / 2
    x_Max_H = centro[0] + largo / 2
    
    # Se ubican las dimensiones de la cruz dentro de la malla de puntos
    mascara_vertical = (xx >= x_Min_V) & (xx <= x_Max_V) & (yy >= y_Min_V) & (yy <= y_Max_V)
    mascara_horizontal = (xx >= x_Min_H) & (xx <= x_Max_H) & (yy >= y_Min_H) & (yy <= y_Max_H)
    
    # La máscara final es la unión de los dos brazos
    mascara_cruz = mascara_vertical | mascara_horizontal
    
    return mascara_cruz


""" Definición de función para máscara corazón """

def funcion_Corazon(centro, xx, yy, escala=1):
    '''
      Crea una máscara con un corazón
    
    FUNCIÓN RECIBE:
        centro == type(list) --> [X, Y]
        xx, yy == malla de puntos en la cual se verá el corazón
        escala == float --> Tamaño del corazón (por defecto 1)
    
    FUNCIÓN RETORNA: Máscara con forma de corazón
    '''
    
    # Si no se especifica el centro, se posiciona en el origen por defecto
    if centro is None:
        centro = (0, 0)
    
    # Reescalamos la malla a las coordenadas ajustadas al centro y la escala
    x = (xx - centro[0]) / escala
    y = (yy - centro[1]) / escala
    
    # Ecuación del corazón: (x^2 + y^2 - 1)^3 - x^2*y^3 <= 0
    mascara_corazon = (x**2 + y**2 - 1)**3 - x**2 * y**3 <= 0
    
    return mascara_corazon


""" Definición de función para máscara de representación Anillos de fase """

def funcion_AnilloFase(radio_interno, radio_externo, fase, xx, yy):
    
    """
    Crea un anillo de fase en el plano definido por las coordenadas xx, yy.
    
    FUNCIÓN RECIBE:
    - radio_interno: Radio interno del anillo.
    - radio_externo: Radio externo del anillo.
    - fase: Retardo de fase introducido por el anillo.
    - xx, yy: Mallas de coordenadas del plano.

    FUNCIÓN RETORNA:
    - anillo_fase: Máscara compleja con el perfil del anillo de fase.
    """
    
    #Se calcula la geometría de un radio dentro de la malla de puntos deseada
    r = np.sqrt(xx**2 + yy**2)

    #Se especifica una geometría de anillo dentro de la malla de puntos
    anillo = (r >= radio_interno) & (r <= radio_externo) 

    #Se aplica la fase dentro de la condición de fase del anillo 
    anillo_fase = np.exp(1j * fase) * anillo  + (1 - anillo) 
    
    return anillo_fase


""" Definición de función para máscara de representación Anillos  """

def funcion_Anillo(radio_interno, radio_externo, xx, yy,transparencia):
    
    """
    Crea un anillo de fase en el plano definido por las coordenadas xx, yy.
    
    FUNCIÓN RECIBE:
    - radio_interno: Radio interno del anillo.
    - radio_externo: Radio externo del anillo.
    - fase: Retardo de fase introducido por el anillo.
    - xx, yy: Mallas de coordenadas del plano.

    FUNCIÓN RETORNA:
    - anillo_fase: Máscara compleja con el perfil del anillo de fase.
    """
    
    #Se calcula la geometría de un radio dentro de la malla de puntos deseada
    r = np.sqrt(xx**2 + yy**2)

    #Se especifica una geometría de anillo dentro de la malla de puntos
    anillo = (r >= radio_interno) & (r <= radio_externo) 

    #Se relaciona la transparencia deseada para el anillo
    anillo = anillo * (1 - transparencia)

    #Se define estructura para anillo invertido
    anillo_invertido = 1-anillo

    return anillo_invertido

