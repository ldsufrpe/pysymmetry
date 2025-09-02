############## util #######################
__all__ = [ 
 'laplacian2d',
 'laplacian1d',
'advection_diffusion_2d', 
  'to_csr',
  'to_csc',
  'view'
  ] 
from sage.all import *
from sage.libs.gap.element import GapElement
from scipy.sparse import csr_matrix, coo_matrix, csc_matrix
from scipy.linalg import block_diag
from scipy.sparse import spdiags, kronsum, kron, identity


def to_csr(m):

    # Create the appropriate format for the COO format.
    #print(m)
    #m = m.change_ring(RDF)
    term_dict = m.dict()
    data = []#[RDF(1)]*len(term_dict)
    row = []
    col = []
    for k, v in term_dict.items():
        r = k[0]
        c = k[1]
        data.append(v)
        row.append(r)
        col.append(c)
    # Create the COO-matrix
    coo = coo_matrix((data,(row,col)))
    # Let Scipy convert COO to CSR format and return
    return csc_matrix(coo)

def to_csc(m):
    
    # Create the appropriate format for the COO format.    
    
    term_dict = m.dict()
    data = []#[RDF(1)]*len(term_dict)
    row = []
    col = []
    for k, v in term_dict.items():
        r = k[0]
        c = k[1]
        data.append(v)
        row.append(r)
        col.append(c)
    # Create the COO-matrix
    coo = coo_matrix((data,(row,col)))
    # Let Scipy convert COO to CSR format and return
    return csc_matrix(coo)

def view(matrix, latex=True):
    try:
        matrix = matrix.apply_map(lambda x: x.radical_expression())
    except:
        return matrix
    if latex:
        return show(matrix)
    return matrix


def generators1d(n):
    string_reflexao_sigma = ''
    for j in range(1,(n/2+1).floor()):        
            string_reflexao_sigma = string_reflexao_sigma + str((j,n-j+1))            
    return [string_reflexao_sigma]    
 
    


def generators(n):#Nota:temporario

    string_reflexao_sigma = ''
    for j in range(1,(n/2+1).floor()):
        for k in range(0,n):
            string_reflexao_sigma = string_reflexao_sigma + str((k*n+j,k*n+n-j+1))
    
    string_reflexao_miDp = ''
    for i in range(1, n):
        for j in range(i+1, n+1):            
            string_reflexao_miDp = string_reflexao_miDp + str(((i-1)*n+j,(j-1)*n+i))
            
    return [string_reflexao_sigma, string_reflexao_miDp]

import numpy as np
def laplacian2d(n):
    """Return matrix A for 2D Laplace equation using block diagonal
    structure, given the number of unknowns ’n’ in each direction.
    """
    Bdiag = -4 * np.eye(n)
    Bupper = np.diag([1] * (n - 1), 1)
    Blower = np.diag([1] * (n - 1), -1)
    B = Bdiag + Bupper + Blower
    # Creat a list [B,B,B,...,B] with n Bs
    blst = [B] * n
    # Unpack the list of diagonal blocks ’blst’
    # since block_diag expects each block to be passed as separate
    # arguments. This is the same as doing block_diag(B,B,B,...,B)
    A = block_diag(*blst)
    # Upper diagonal array offset by n: we’ve got (n-1) I blocks
    # each containing n ones
    Dupper = np.diag(np.ones(n * (n - 1)), n)
    # Lower diagonal array offset by -n
    Dlower = np.diag(np.ones(n * (n - 1)), -n)
    A += Dupper + Dlower
    return csc_matrix(-A)

from scipy.sparse import spdiags

def laplacian1d(n, h=1):
    x= np.arange(n)
    # h grid width
    #h= #x[width]
    ones=np.ones(n)

    # creating laplace stencil 1/h^2 [1 -2 1] csr format
    data = 1/(h*h)*np.array([ones,-2*ones,ones])
    diags = np.array([-1, 0, 1])
    A = spdiags(data, diags, n, n, format="csc")
    return A



from scipy.sparse import spdiags, kronsum, kron

def advection_diffusion_2d(n, D=1.0, vx=1.0):
    """
    Gera a matriz para a equação de Advecção-Difusão 2D em uma grade n x n.

    - D: Coeficiente de difusão (controla a parte simétrica).
    - vx: Velocidade de advecção na direção x (controla a parte anti-simétrica).    
    """
   
    L = -1 * laplacian2d(n)     
    diagonals_x = [np.ones(n - 1), -np.ones(n - 1)]
    offsets_x = [1, -1]
    Dx_1D = spdiags(diagonals_x, offsets_x, n, n, format="csc")    
    I_y = identity(n, format="csc")
    Ax = kron(I_y, Dx_1D) * 0.5
    A_ad = D * L - vx * Ax
    
    return A_ad