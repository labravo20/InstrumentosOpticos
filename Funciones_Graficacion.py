import matplotlib.pyplot as plt
import numpy as np


""" Función para graficar irradiancias """

def graficar_intensidad(campo,xx_TamañoVentana,yy_TamañoVentana,titulo,titulo_colorBar, titulo_EjeX,titulo_EjeY, valor_min = 1,valor_max = 1):
    
    plt.imshow(campo,extent=[- xx_TamañoVentana/ 2, xx_TamañoVentana/ 2,
                             -yy_TamañoVentana/ 2, yy_TamañoVentana/ 2],
                             cmap="gray",
                             vmin = valor_min*(np.min(campo)),
                             vmax = valor_max*(np.max(campo)))
    
    plt.title(titulo)
    plt.colorbar(label=titulo_colorBar)
    plt.xlabel(titulo_EjeX)
    plt.ylabel(titulo_EjeY)
    plt.show()