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