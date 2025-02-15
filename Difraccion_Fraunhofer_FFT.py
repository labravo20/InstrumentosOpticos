#Mensaje inicial para verificar funcionamiento correcto entorno de programación
print("Inicializando programa simulación Difracción Fraunhofer...")

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

# Visualización
#plt.figure(figsize=(6, 6))

#plt.subplot(2, 2, 1)
plt.imshow(mascara, extent=[-longitud_Arreglo/2, longitud_Arreglo/2, -longitud_Arreglo/2, longitud_Arreglo/2], cmap='gray')
plt.title("Máscara Circular")
plt.colorbar(label="Transmitancia")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()

# Espectro de la FFT2
#plt.subplot(1, 2, 2)
plt.imshow(espectro_magnitud, extent=[-1/(longitud_Arreglo/2), 1/(longitud_Arreglo/2), -1/(longitud_Arreglo/2), 1/(longitud_Arreglo/2)], cmap='inferno', vmax= 0.7*(np.max(espectro_magnitud)))  # log1p para mejor visualización
plt.title("Espectro de Difracción (FFT2)")
plt.colorbar(label="Amplitud del espectro")
plt.xlabel("X (1/m)")
plt.ylabel("Y (1/m)")
plt.show()

# Patrón de difracción de Fraunhofer
plt.imshow(irradiancia_modulo_cuadrado, extent=[-longitud_Arreglo/2, longitud_Arreglo/2, -longitud_Arreglo/2, longitud_Arreglo/2], cmap='gray', vmax= 0.05*(np.max(irradiancia_modulo_cuadrado)))  # log1p para mejor visualización
plt.title("Patrón de Difracción")
plt.colorbar(label="Amplitud²")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()