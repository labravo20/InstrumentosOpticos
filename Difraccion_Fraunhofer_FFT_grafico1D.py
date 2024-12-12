import numpy as np
import matplotlib.pyplot as plt


''' Definicion de las difefrentes funciones de transmitancia'''

def malla_Puntos(resolucion, longitud_Arreglo):
    ''' CREACION DE LAS MALLAS DE PUNTOS Y LOS DELTAS PARA LOS PRODUCTOS ESPACIO FRECUENCIA '''
    x = np.linspace(-longitud_Arreglo / 2, longitud_Arreglo / 2, resolucion) #crea las mallas de puntos para el arreglo 
    y = np.linspace(-longitud_Arreglo / 2, longitud_Arreglo / 2, resolucion) 
    xx, yy = np.meshgrid(x, y) #crea una malla de puntos bidimensional 
    return xx, yy #retornamos la malla de puntos


def funcion_Circulo(radio, centro, xx, yy): #definicion de la funcion para hacer circulo transparente
    ''' CREACCION DEL CONJUNTO DE PUNTOS DE LA MASCARA DE DFRACCION '''
    if centro is None: # que pasa si el centro no es definido en la funcion
        centro = (0, 0) #ubica el centro de la circunferencia en el origen
    distancia = (xx - centro[0])**2 + (yy - centro[1])**2 #calculamos la distancia desde el centro de la circunferencia a cada punto
    mascara = distancia <= radio**2 #los puntos de la mascara seran los puntos cuya distancia al centro es menor que el radio
    return mascara #devolvemos los puntos que cumplen la condicion para hacer parte de la mascara



resolucion = 3000  # Número de puntos en la malla
longitud_Arreglo = 0.5  # Tamaño físico del área (10 mm)
radio = 0.002  # Radio del círculo en metros
centro = None  # El centro será el origen si es None


# Crear la malla de puntos
xx, yy = malla_Puntos(resolucion, longitud_Arreglo)

# Crear la máscara circular
mascara = funcion_Circulo(radio, centro, xx, yy)


#Definición de longitud de onda
longitud_onda = 0.000000533 #UNIDADES: m

#Definición de distancia z entre plano de máscara difractiva y plano medición
distancia_entre_planos= 10000000*longitud_onda #UNIDADES: m

#Calculando el término de fase constante
numero_onda = (2*np.pi)/longitud_onda

fase_constante = (np.exp(1j*distancia_entre_planos*numero_onda))/(longitud_onda*distancia_entre_planos*1j)

# Calcular la FFT2
fft2_mascara = np.fft.fft2(mascara)
fft2_mascara_centrada = np.fft.fftshift(fft2_mascara)

# Calcular la magnitud del espectro
amplitud_propagacion_out = fft2_mascara_centrada*fase_constante
amplitud_compleja_out = fft2_mascara_centrada


espectro_magnitud = np.abs(amplitud_compleja_out) # -> Representa el módulo del valor imaginario
                                                  # Para encontrar la irradiancia es necesario modulo²
espectro_modulo_cuadrado = (espectro_magnitud**2)

irradiancia_magnitud = np.abs(amplitud_propagacion_out) # -> Representa el módulo del valor imaginario
                                                  # Para encontrar la irradiancia es necesario modulo²
irradiancia_modulo_cuadrado = (irradiancia_magnitud**2)

# Extraemos el perfil de irradiancia a lo largo del eje x (o eje y)
perfil_irradiancia_x = np.sum(irradiancia_modulo_cuadrado, axis=0)  # Sumamos a lo largo del eje y (corte en x)
perfil_irradiancia_y = np.sum(irradiancia_modulo_cuadrado, axis=1)  # Sumamos a lo largo del eje x (corte en y)

# Gráfica del perfil a lo largo del eje X
plt.plot(xx[0], perfil_irradiancia_x)
plt.title("Perfil de Irradiancia a lo largo del eje X")
plt.xlabel("Posición en X (m)")
plt.ylabel("Irradiancia")
plt.grid(True)
#plt.show()

# Gráfica del perfil a lo largo del eje Y
plt.plot(yy[:, 0], perfil_irradiancia_y)
plt.title("Perfil de Irradiancia a lo largo del eje Y")
plt.xlabel("Posición en Y (m)")
plt.ylabel("Irradiancia")
plt.grid(True)
#plt.show()

