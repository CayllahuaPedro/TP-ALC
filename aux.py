import numpy as np
def construye_adyacencia(D,m): 
    # Función que construye la matriz de adyacencia del grafo de museos
    # D matriz de distancias, m cantidad de links por nodo
    # Retorna la matriz de adyacencia como un numpy.
    D = D.copy()
    l = [] # Lista para guardar las filas
    for fila in D: # recorriendo las filas, anexamos vectores lógicos
        l.append(fila<=fila[np.argsort(fila)[m]] ) # En realidad, elegimos todos los nodos que estén a una distancia menor o igual a la del m-esimo más cercano
    A = np.asarray(l).astype(int) # Convertimos a entero
    np.fill_diagonal(A,0) # Borramos diagonal para eliminar autolinks
    return(A)

def calcula_transicion(A):
  n= A.shape[0]
  #calcular A transpuesta
  A_t = np.transpose(A)
  #calcular K
  diagonal = np.sum(A, axis=1)
  K = np.diag(diagonal)
  #calcular K inversa
  K_inv= np.diag(1/diagonal)
  #multiplicacion y obtener C (transicion)
  C = A_t @ K_inv
  return C

#agarra el vector de ranking y devuelve los 3 museos con mayor valor
def calcula_top_3_museos(p): 
    top_3 = np.argsort(p)[::-1][:3]
    return top_3
    



    