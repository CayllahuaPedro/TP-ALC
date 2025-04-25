import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import geopandas as gpd
def construye_adyacencia(D, m):
    n = D.shape[0]
    l = [] # Lista para guardar las filas
    # Ensure m is an integer right at the start
    m_int = int(m) 
    for fila in D: # recorriendo las filas, anexamos vectores lógicos
        num_elementos = len(fila)
        if num_elementos == 0: # Handle empty row explicitly
             l.append(np.zeros(n, dtype=bool)) 
             continue

        # Ensure index is within bounds using the integer m
        # Also handles m_int < 0 by clipping to 0 if necessary
        idx_m_vecino = max(0, min(m_int, num_elementos - 1)) 
        
        # Get the sorted indices (these are integers)
        indices_ordenados = np.argsort(fila)
        
        # Get the index *from the original array* corresponding to the m-th smallest distance
        indice_en_fila_del_m_vecino = indices_ordenados[idx_m_vecino]

        # Now use this index (which MUST be an integer) to get the actual distance value
        # Cast explicitly to int just to be absolutely certain
        distancia_m_vecino = fila[int(indice_en_fila_del_m_vecino)] 
        
        l.append(fila <= distancia_m_vecino) # En realidad, elegimos todos los nodos que estén a una distancia menor o igual a la del m-esimo más cercano
        
    A = np.asarray(l).astype(int) # Convertimos a entero
    np.fill_diagonal(A,0) # Borramos diagonal para eliminar autolinks
    return A

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


def calculaLU(A):
    # matriz es una matriz de NxN
    # Retorna la factorización LU a través de una lista con dos matrices L y U de NxN.
    # Completar! Have fun
    m=A.shape[0]#fila
    n=A.shape[1]
    Ac = A.copy()#U
    #L=np.eye(n)
    if m!=n:
        print('Matriz no cuadrada')
        return
    
    ## desde aqui -- CODIGO A COMPLETAR
    print(n)
    for j in range(n):
        for i in range (j+1,n):
            mult=Ac[i,j]/Ac[j,j] #escalonas/mult es el factor tipo f2-multf1
            Ac[i,j:] = Ac[i,j:]-mult*Ac[j,j:]#resta de finlas// j: dessde j hasta el final
            Ac[i,j]=mult
         #   Ac[i,:]==Ac[i,:]-L[i,j]*Ac[j,:] 
           # cant_op= cant_op+2
          #  return L, Ac, cant_op
                 
    #L = np.tril(Ac,-1) + np.eye(A.shape[0]) #np.eye es la matriz con 1 en diagonal ()
    #U = np.triu(Ac) #CAPTA LA DIAGONAL INFERIOE
         
            
    return Ac

def graph_map_rank (ranking_scores, barrios, museos):
    normalized_ranking_scores = (ranking_scores - ranking_scores.min()) / (ranking_scores.max() - ranking_scores.min())


    # graficar mapa
    fig, ax = plt.subplots(figsize=(10, 10))
    barrios.boundary.plot(color='gray', ax=ax)

    # graficar museos con label de ranking.
    museos.plot(ax=ax, column=normalized_ranking_scores, cmap='viridis', legend=True, markersize=normalized_ranking_scores*100)  # Adjust markersize for better visibility

    # Asociar cada museo con su ranking
    #for x, y, label in zip(museos.geometry.x, museos.geometry.y, normalized_ranking_scores):
    #    ax.annotate(f'{label:.2f}', xy=(x, y), xytext=(3, 3), textcoords="offset points", fontsize=8)

    plt.title("Museos de CABA con Rankings")
    plt.show()

def resolver_LU (L,U,b):
    # Resuelve el sistema Ax=b usando la factorización LU
    # L: matriz triangular inferior
    # U: matriz triangular superior
    # b: vector de términos independientes
    # Retorna el vector solución x
    y = scipy.linalg.solve_triangular(L, b, lower=True)  # Resolvemos Ly = b
    x = scipy.linalg.solve_triangular(U, y)  # Resolvemos Ux = y
    return x

def calcula_matriz_C(A): 
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


def calcula_matriz_C_continua(D): 
    # Función para calcular la matriz de trancisiones C
    # A: Matriz de adyacencia
    # Retorna la matriz C en versión continua
    D = D.copy()
    F = 1/D
    np.fill_diagonal(F,0)
    Kinv = ... # Calcula inversa de la matriz K, que tiene en su diagonal la suma por filas de F 
    C = ... # Calcula C multiplicando Kinv y F
    return C

def calcula_B(C,cantidad_de_visitas):
    # Recibe la matriz T de transiciones, y calcula la matriz B que representa la relación entre el total de visitas y el número inicial de visitantes
    # suponiendo que cada visitante realizó cantidad_de_visitas pasos
    # C: Matirz de transiciones
    # cantidad_de_visitas: Cantidad de pasos en la red dado por los visitantes. Indicado como r en el enunciado
    # Retorna:Una matriz B que vincula la cantidad de visitas w con la cantidad de primeras visitas v
    B = np.eye(C.shape[0])
    for i in range(cantidad_de_visitas-1):
        # Sumamos las matrices de transición para cada cantidad de pasos
        B += np.linalg.matrix_power(C, i + 1)
    return B