import numpy as np
import matplotlib.pyplot as plt
import Mascaras_Transmitancia as m


resolucion = 3000  # Número de puntos en la malla
longitud_Arreglo = 0.5  # Tamaño físico del área (10 mm)
radio = 0.002  # Radio del círculo en metros
centro = None  # El centro será el origen si es None


# Crear la malla de puntos
xx, yy = m.malla_Puntos(resolucion, longitud_Arreglo)

# Crear la máscara circular
mascara = m.funcion_Circulo(radio, centro, xx, yy)


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

