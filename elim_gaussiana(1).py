#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Eliminacion Gausianna
"""
import numpy as np

def elim_gaussiana(A):
    cant_op = 0
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
        
    #print(cant_op)
    ## hasta aqui
                 
    L = np.tril(Ac,-1) + np.eye(A.shape[0]) #np.eye es la matriz con 1 en diagonal ()
    U = np.triu(Ac) #CAPTA LA DIAGONAL INFERIOE
         
            
    return L, U, cant_op


    
def main():
    n = 7
    B = np.eye(n) - np.tril(np.ones((n,n)),-1) 
    B[:n,n-1] = 1
    print('Matriz B \n', B)
    
    L,U,cant_oper = elim_gaussiana(B)
    
    print('Matriz L \n', L)
    print('Matriz U \n', U)
    print('Cantidad de operaciones: ', cant_oper)
    print('B=LU? ' , 'Si!' if np.allclose(np.linalg.norm(B - L@U, 1), 0) else 'No!')
    print('Norma infinito de U: ', np.max(np.sum(np.abs(U), axis=1)) )

if __name__ == "__main__":
    main()
    