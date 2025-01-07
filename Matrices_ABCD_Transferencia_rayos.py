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



""" Creación de matriz ABCD para LENTES DELGADAS """

def lente_Delgada(radio_1, radio_2, n_Incidente, n_Lente, n_Salida, tamaño_Fisico = None):
    '''
    Función para calcular la matriz ABCD de una lente delgada.

    FUNCIÓN RECIBE:

        radio_1     == float si finito, str si infinito
        radio_2     == float si finito, str si infinito
        n_Incidente == float, por default es 1 (aire)
        n_Lente     == float, por default es 1.5 (vidrio)
        n_Salida    == float, por default es 1 (aire)
    
    FUNCIÓN RETORNA:
        
        Matriz ABCD correspondiente
    '''

    ''' Definicion de valores por default para los indices de refraccion, se asume por default que la lente esta hecha de vidrio y que esta 
    inmersa en el aire
    Por tanto:
    n_Incidente = 1 por default
    n_Lente = 1.5 por default
    n_Salida = 1 por default'''
    
    if n_Incidente is None: #establecemos un valor por default para los indices de refraccion
        n_Incidente = 1 #si no se define, entonces el indice de refraccion es el del aire
    
    if n_Salida is None: #Establecemos un valor por default para el indice de refraccion de salida
        n_Salida = 1 #Si no se define, entonces el indice de refraccion de salida es el del aire
    
    if n_Lente is None: #establecemos un valor por default para el indice de refraccion de los lentes en caso de que se deje vacio
        n_Lente = 1.5 #Asumimos por default que las lentes estan hechas de vidrio
    
    matriz = matriz_Inicial()              #Se crea matriz identidad para empezar a trabajar
    
    poder_Lente = ((n_Lente - n_Incidente)/determinacion_Radio(radio_1)) + ((n_Salida - n_Lente)/determinacion_Radio(radio_2)) #se calcula el poder de covergencia de la lente usando la ecuacion del fabricante de lentes
    
    matriz[1,0] = -poder_Lente
    '''
    |1                0  | Poder_lente = 1/f
    |Poder_lente      1  |
    '''
    return matriz


