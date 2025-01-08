""" Se genera un print de mensaje inicial para verificar correcto funcionamiento del entorno """
print("Inicializando entorno de programación Matrices ABCD transferencia de rayos...")



""" Importando librerias necesarias para el desarrollo del código """
import numpy as np


''' Función para análisis de radio lente '''

def determinacion_Radio(radio):

    """Función retorna el valor al cual se quiere signar el radio de la lente

    FUNCIÓN RECIBE: 

    Caso 01: Valor numérico del radio de interés para la lente
    Caso 02: Un string cualquiera para asignar valor infinito al radio de la lente

    FUNCIÓN RETORNA:

    Caso 01: Valor numérico ingresado para el radio de interés para la lente 
    Caso 02: Valor infinito  
    
    """

    #Inicialmente se verifica si es posible transformar el radio en estudio a un float
    try:

        radio = float(radio)        
        return radio
    
    #Al ingresar un string o valor NO convertibles a float entonces se retorna "infinito" 
    except:

        return np.inf               



''' Matriz de Entrada del sistema '''

#Función para crear una matriz inicial para poder usar como herramienta de cálculo en el desarrollo 
#del código
def matriz_Inicial():

    """ Función para crear una matriz identidad"""
    
    # Se crea y retorna una matriz identidad de tamaño 2x2 --> Este tamaño dada la necesidad
    # de trabajar con matrices ABCD.
    return np.eye(2)       



""" Creación de matriz ABCD para PROPAGACIÓN EN MEDIO HOMOGÉNEO """

def propagacion_MedioHomogeneo(distancia_propagacion):

    '''
    Función para calcular la matriz correspondiente a la propagación en medio homogéneo.
    
    FUNCIÓN RECIBE: Distancia de propagación  (Type: float)
    
    FUNCIÓN RETORNA: Matriz ABCD correspondiente a propagación en medios homogéneos
    '''
    
    # Se crea matriz identidad sobre la cual se calculará la matriz de propagación
    matriz_propagacionMedioHomogeneo = matriz_Inicial()

    # Se asigna el valor de la distancia de propagación en la posicion pre determinada 
    # para una matriz ABCD de propagación en medio homogéneo:
    matriz_propagacionMedioHomogeneo[0,1] = distancia_propagacion # --> Se posiciona el valor de la distancia
                                                                  # de propagación en la primera fila segunda columna.

    ''' MATRIZ RESULTANTE:
    |1   distancia_Propagación|
    |0            1           |
    '''

    return matriz_propagacionMedioHomogeneo



""" Creación de matriz ABCD para REFRACCIÓN"""

def refraccion(n_Incidente, n_Salida):

    '''
    Función para calcular la matriz ABCD para refracción.

    FUNCIÓN RECIBE:

        -índice de refracción medio incidente (Type: float) 
        -índice de refracción medio salida (Type: float)
    
    FUNCIÓN RETORNA:
        
        Matriz ABCD correspondiente a la transferencia de rayos en caso de refracción.
    '''
    
    # Se crea matriz identidad sobre la cual se calculará la matriz para refracción
    matriz_refraccion = matriz_Inicial()
    
    # Se calcula la relación entre los índices de refracción de los medios en consideración
    relacion_indicesRefraccion = n_Incidente/n_Salida
    
    #Se asigna el valor del negativo del poder de convergencia de la lente en la posición
    #pre determinada para una matriz ABCD para lentes delgadas
    matriz_refraccion[1,1] = relacion_indicesRefraccion # --> Se posiciona el valor del poder de
                                                          # la relacion entre los indices de refraccion
                                                          # en la segunda fila segunda columna.

    '''MATRIZ RESULTANTE:
    | 1                    0            | 
    | 0      relacion_indicesRefraccion |
    
    '''
    return matriz_refraccion



""" Creación de matriz ABCD para LENTES DELGADAS """

def lente_Delgada(radio_1, radio_2, n_Incidente, n_Lente, n_Salida):

    '''
    Función para calcular la matriz ABCD de una lente delgada.

    FUNCIÓN RECIBE:

        -radio entrada lente  ((Type: float) para valor establecido o (Type: str) para tendencia a infinito)
        -radio salida lente  ((Type: float) para valor establecido o (Type: str) para tendencia a infinito)
        -índice de refracción medio incidente (Type: float) 
        -índice de refracción del material lente (Type: float) --> Por generalidad se puede considerar 1.5 (vidrio) 
        -índice de refracción medio salida (Type: float)
    
    FUNCIÓN RETORNA:
        
        Matriz ABCD correspondiente a la transferencia de rayos a través de una lente delgada.
    '''
    
    # Se crea matriz identidad sobre la cual se calculará la matriz para lentes delgadas
    matriz_lentesDelgadas = matriz_Inicial()
    
    #Se calcula el término asociado a la ecuación del fabricante de lentes "1/f" siendo f la distancia focal de la lente
    #Este término se asocia al PODER DE CONVERGENCIA de la lente...

    ###Calculando el poder de convergencia de la superficie de entrada:
    poder_convergenciaSuperficieEntrada = ((n_Lente - n_Incidente)/determinacion_Radio(radio_1))

    ###Calculando el poder de convergencia de la superficie de salida:
    poder_convergenciaSuperficieSalida = ((n_Salida - n_Lente)/determinacion_Radio(radio_2))

    ###Calculando el poder de convergencia de la lente:
    poder_convergenciaLente = poder_convergenciaSuperficieEntrada + poder_convergenciaSuperficieSalida
    
    #Se asigna el valor del negativo del poder de convergencia de la lente en la posición
    #pre determinada para una matriz ABCD para lentes delgadas
    matriz_lentesDelgadas[1,0] = -poder_convergenciaLente # --> Se posiciona el valor del poder de
                                                          # convergencia de la lente en la segunda
                                                          # fila primera columna.

    '''MATRIZ RESULTANTE:
    |            1                 0  | 
    |-poder_convergenciaLente      1  |
    
    '''
    return matriz_lentesDelgadas


