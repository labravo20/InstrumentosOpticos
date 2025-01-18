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
    
    #Se asigna el valor de la relación entre los índices de refracción en la posición
    #pre determinada para una matriz ABCD para caso de refracción.
    matriz_refraccion[1,1] = relacion_indicesRefraccion # --> Se posiciona el valor de la relación
                                                          # entre los indices de refraccion
                                                          # en la segunda fila segunda columna.

    '''MATRIZ RESULTANTE:
    | 1                    0            | 
    | 0      relacion_indicesRefraccion |
    
    '''
    return matriz_refraccion



""" Creación de matriz ABCD para CURVAS REFRACTIVAS"""

def curva_Refractiva(radio, n_Incidente, n_Salida):

    '''
    Función para calcular la matriz ABCD para curva refractiva.

    FUNCIÓN RECIBE:

        -radio superficie  ((Type: float) para valor establecido o (Type: str) para tendencia a infinito)
        -índice de refracción medio incidente (Type: float) 
        -índice de refracción medio salida (Type: float)
    
    FUNCIÓN RETORNA:
        
        Matriz ABCD correspondiente a la transferencia de rayos en caso curvas refractivas.
    '''
    
    # Se crea matriz identidad sobre la cual se calculará la matriz para curva refractiva
    matriz_curvaRefractiva = matriz_Inicial()
    
    # Se calcula la relación entre los índices de refracción de los medios en consideración
    relacion_indicesRefraccion = n_Incidente/n_Salida

    
    #Se asigna el valor de la relación de los índices de refracción en la posición
    #pre determinada para una matriz ABCD para curvas refracitvas
    matriz_curvaRefractiva[1,1] = relacion_indicesRefraccion # --> Se posiciona el valor de
                                                          # la relacion entre los indices de refraccion
                                                          # en la segunda fila segunda columna.
    
    # Se calcula el término de relación entre índices de refracción considerando el radio de curvatura
    relacion_indicesRefraccionCurvatura = (n_Incidente-n_Salida)/(n_Salida*determinacion_Radio(radio))

    #Se asigna el valor del término de relación entre índices de refracción considerando el radio de
    #curvatura en la posición pre determinada para una matriz ABCD para lentes delgadas
    matriz_curvaRefractiva[1,0] = relacion_indicesRefraccionCurvatura # --> Se posiciona el valor
                                                          # de la relacion entre los indices de 
                                                          # refraccion considerando la curvatura
                                                          # en la segunda fila primera columna.

    '''MATRIZ RESULTANTE:
    |                  1                                     0            | 
    | relacion_indicesRefraccionCurvatura      relacion_indicesRefraccion |
    
    '''
    return matriz_curvaRefractiva



""" Creación de matriz ABCD para CURVAS REFLECTIVAS"""

def curva_Reflectiva(radio):

    '''
    Función para calcular la matriz ABCD para curva reflectiva.

    FUNCIÓN RECIBE:

        -radio superficie  ((Type: float) para valor establecido o (Type: str) para tendencia a infinito)
    
    FUNCIÓN RETORNA:
        
        Matriz ABCD correspondiente a la transferencia de rayos en caso curvas reflectivas.
    '''
    
    # Se crea matriz identidad sobre la cual se calculará la matriz para curvas reflectivas
    matriz_curvaReflectiva = matriz_Inicial()
    
    #Se asigna el valor -1 en la posición
    #pre determinada para una matriz ABCD para curvas reflectivas
    matriz_curvaReflectiva[1,1] = -1 # --> Se posiciona el valor -1
                                     # en la segunda fila segunda columna.
    
    # Se calcula el término relacionado con el aporte del radio de curvatura
    aporte_Radio = 2/determinacion_Radio(radio)

    #Se asigna el valor del término de aporte del radio de curvatura
    #en la posición pre determinada para una matriz ABCD para lentes delgadas
    matriz_curvaReflectiva[1,0] = aporte_Radio # --> Se posiciona el valor
                                                          # del aporte de la curvatura
                                                          # en la segunda fila primera columna.

    '''MATRIZ RESULTANTE:
    |       1            0 | 
    | aporte_Radio      -1 |
    
    '''
    return matriz_curvaReflectiva



""" Creación de matriz ABCD para LENTES DELGADAS sin conocer la distancia focal"""

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



""" Creación de matriz ABCD para LENTES DELGADAS conociendo la distancia focal"""

def lente_DelgadaConociendoDistanciaFocal(distancia_focal):

    '''
    Función para calcular la matriz ABCD de una lente delgada.

    FUNCIÓN RECIBE: distancia focal lente  ((Type: float) 
    
    FUNCIÓN RETORNA:
        
        Matriz ABCD correspondiente a la transferencia de rayos a través de una lente delgada.
    '''
    
    # Se crea matriz identidad sobre la cual se calculará la matriz para lentes delgadas
    matriz_lentesDelgadas = matriz_Inicial()
    
    #Se calcula el término asociado a la ecuación del fabricante de lentes "1/f" siendo f la distancia focal de la lente
    #Este término se asocia al PODER DE CONVERGENCIA de la lente...
    poder_convergenciaLente = 1/distancia_focal
    
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


""" Creación de matriz ABCD para DIFRACCIÓN """
def matriz_ABCD_Difraccion(camino_optico_central, campo_entrada, posicion_A_matriz, posicion_B_matriz, 
                           posicion_D_matriz, xx_entrada, yy_entrada, xx_salida, yy_salida,numero_onda,
                           deltas_muestreo):
    
    '''
    Función para calcular ...

    FUNCIÓN RECIBE: ... 
    
    FUNCIÓN RETORNA:...
    '''
    
    #Calculando término de fase constante:
    fase_constante = np.exp(1j*numero_onda*camino_optico_central)

    #Calculando término de fase parabólica plano medición
    fase_parabolicaPlanoMedicion = np.exp(1j*(numero_onda/(2*posicion_B_matriz))*
                                          posicion_D_matriz*((xx_salida**2)+(yy_salida**2))) 
    
    #Calculando término de fase parabólica plano máscara difractiva
    fase_parabolicaPlanoMascara = np.exp(1j*(numero_onda/(2*posicion_B_matriz))*
                                          posicion_A_matriz*((xx_entrada**2)+(yy_entrada**2))) 
    
    #Calculando función de contribución campo entrada y fase parabólica entrada
    funcion_Entrada = campo_entrada*fase_parabolicaPlanoMascara

    #Calculando la transformada de Fourier de la función asociada al campo de entrada afectada por la
    #fase parabólica de entrada
    #transformada_FourierFuncionEntrada = np.fft.fftshift(np.fft.fft2(funcion_Entrada))
    transformada_FourierFuncionEntrada = (np.fft.fft2(funcion_Entrada))

    #Calculando el campo difractado
    campo_Difractado = deltas_muestreo[0]*deltas_muestreo[1]*fase_constante*fase_parabolicaPlanoMedicion*transformada_FourierFuncionEntrada


    return campo_Difractado


""" Funcion para calcular la matriz del sistema """
def matriz_Sistema(lista_matricesArregloDifractivoInvertida):

    
    # Inicializar el resultado con la primera matriz
    matriz_sistema = lista_matricesArregloDifractivoInvertida[0]
    
    # Multiplicar las matrices en orden
    for matriz in lista_matricesArregloDifractivoInvertida[1:]:
        
        #Se realiza la multiplicación de las matrices para calcular la matriz del sistema 
        matriz_sistema = np.dot(matriz_sistema, matriz)
    
    return matriz_sistema



""" Función para calcular el camino óptico central """
def camino_Optico(lista_matricesArregloDifractivo):

    matriz_inicial = lista_matricesArregloDifractivo[0]

    camino_optico = matriz_inicial[0,1] #Se aigna el valor del término B --> Correspondiente siempre
                                        #a la distancia de propagación...
    
    # Multiplicar las matrices en orden
    for matriz in lista_matricesArregloDifractivo[1:]:

        #Se realiza la multiplicación de las matrices para calcular la matriz del sistema 
        camino_optico = camino_optico + matriz[0,1]


    return camino_optico




