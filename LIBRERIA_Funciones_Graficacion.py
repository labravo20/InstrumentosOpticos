import matplotlib.pyplot as plt
import numpy as np


""" Función para graficar irradiancias """

def graficar_intensidad(campo,xx_TamañoVentana,yy_TamañoVentana,titulo, valor_min = 1,valor_max = 1):
    
    plt.imshow(campo,extent=[- xx_TamañoVentana/ 2, xx_TamañoVentana/ 2,
                             -yy_TamañoVentana/ 2, yy_TamañoVentana/ 2],
                             cmap="gray",
                             vmin = valor_min*(np.min(campo)),
                             vmax = valor_max*(np.max(campo)))
    
    plt.title(titulo)
    #plt.colorbar(label="Intensidad normalizada")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.show()


""" Función para graficar amplitud """

def graficar_amplitud(campo,xx_TamañoVentana,yy_TamañoVentana,titulo, valor_min = 1,valor_max = 1):
    
    plt.imshow(campo,extent=[- xx_TamañoVentana/ 2, xx_TamañoVentana/ 2,
                             -yy_TamañoVentana/ 2, yy_TamañoVentana/ 2],
                             cmap="gray",
                             vmin = valor_min*(np.min(campo)),
                             vmax = valor_max*(np.max(campo)))
    
    plt.title(titulo)
    plt.colorbar(label="Amplitud")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.show()



""" Función para graficar irradiancias en escala logaritmica """

def graficar_intensidadLog(campo,xx_TamañoVentana,yy_TamañoVentana,titulo, valor_min = 1,valor_max = 1):
    
    plt.imshow(np.log(campo),extent=[- xx_TamañoVentana/ 2, xx_TamañoVentana/ 2,
                             -yy_TamañoVentana/ 2, yy_TamañoVentana/ 2],
                             cmap="gray",
                             vmin = valor_min*(np.min(campo)),
                             vmax = valor_max*(np.max(campo)))
    
    plt.title(titulo)
    plt.colorbar(label="Intensidad")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.show()



""" Función para graficar fase """

def graficar_fase(campo,xx_TamañoVentana,yy_TamañoVentana,titulo):
    
    plt.imshow(campo, extent=[-xx_TamañoVentana/2, xx_TamañoVentana/2,
                              -yy_TamañoVentana/2, yy_TamañoVentana/2], 
                              cmap='inferno')
    plt.title(titulo)
    plt.colorbar(label="Fase")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.show()



""" Función para graficar longitud de camino óptico """

def graficar_longitudCaminoOptico(campo,xx_TamañoVentana,yy_TamañoVentana,titulo):
    
    plt.imshow(campo, extent=[-xx_TamañoVentana/2, xx_TamañoVentana/2,
                              -yy_TamañoVentana/2, yy_TamañoVentana/2], 
                              cmap='winter')
    plt.title(titulo)
    plt.colorbar(label="Longitud de camino óptico [m]")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.show()


""" Función para graficar altura en muestras"""

def graficar_altura(campo,xx_TamañoVentana,yy_TamañoVentana,titulo):
    
    plt.imshow(campo, extent=[-xx_TamañoVentana/2, xx_TamañoVentana/2,
                              -yy_TamañoVentana/2, yy_TamañoVentana/2], 
                              cmap='winter')
    plt.title(titulo)
    plt.colorbar(label="Altura [m]")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.show()


""" Funcion para graficar transmitancias """

def graficar_transmitancia(campo,xx_TamañoVentana,yy_TamañoVentana,titulo):
    
    plt.imshow(campo,extent=[- xx_TamañoVentana/ 2, xx_TamañoVentana/ 2,
                             -yy_TamañoVentana/ 2, yy_TamañoVentana/ 2],
                             cmap="gray")
    
    plt.title(titulo)
    plt.colorbar(label="Transmitancia")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.show()
    