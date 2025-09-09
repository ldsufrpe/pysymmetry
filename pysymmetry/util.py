r"""
Utility functions for sparse linear algebra and simple PDE stencils.

This module provides helpers to convert Sage matrices to SciPy sparse
formats, quick preview utilities, generators for permutation groups in
grid contexts, and standard finite-difference matrices used in PDE
experiments.

Environment versions
- SageMath 10.4
- Python 3.11
- NumPy 1.26
- SciPy 1.14

Notes
- Functions here are numeric and rely on NumPy/SciPy; they are meant to be
  used inside a Sage session.
- See also: pysymmetry.pde for higher-level numerical representation tools.
"""
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
    r"""
    Convert a Sage sparse matrix to a SciPy CSC sparse matrix.

    INPUT:
    - m -- a Sage matrix (typically sparse) with a .dict() method of nonzeros

    OUTPUT:
    - SciPy csc_matrix with the same nonzero pattern and values

    EXAMPLES::

        sage: from sage.all import matrix
        sage: from pysymmetry.util import to_csr
        sage: M = matrix({(0,0): 2, (1,2): -1})
        sage: S = to_csr(M)
        sage: S.get_shape()
        (2, 3)
    """

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
    r"""
    Convert a Sage sparse matrix to a SciPy CSC sparse matrix.

    INPUT:
    - m -- a Sage matrix (typically sparse) with a .dict() method of nonzeros

    OUTPUT:
    - SciPy csc_matrix with the same nonzero pattern and values

    EXAMPLES::

        sage: from sage.all import matrix
        sage: from pysymmetry.util import to_csc
        sage: M = matrix({(0,1): 3, (2,0): -5})
        sage: S = to_csc(M)
        sage: S.get_shape()
        (3, 2)
    """
    
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
    r"""
    Pretty-print a Sage matrix, optionally using LaTeX rendering.

    Attempts to convert symbolic entries to radical expressions for
    readability. If latex=True, shows the matrix via Sage's show().

    INPUT:
    - matrix -- a Sage matrix
    - latex -- bool (default: True); if True, render with show()

    OUTPUT:
    - The displayed object (if latex=True) or the transformed matrix

    EXAMPLES::

        sage: from sage.all import matrix, sqrt
        sage: from pysymmetry.util import view
        sage: M = matrix([[1, sqrt(2)], [0, 1]])
        sage: view(M, latex=False)  # returns a matrix with radical_expression applied
        [1 sqrt(2)]
        [0       1]
    """
    try:
        matrix = matrix.apply_map(lambda x: x.radical_expression())
    except:
        return matrix
    if latex:
        return show(matrix)
    return matrix


def generators1d(n):
    r"""
    Return a list with a single permutation string for 1D reflection symmetry.

    For an integer n, builds the reflection (j, n-j+1) for j=1..floor(n/2),
    concatenated into a single cycle-string suitable for Sage's parser.

    INPUT:
    - n -- positive integer, grid size

    OUTPUT:
    - list -- [sigma_string], where sigma_string encodes the reflection permutation

    EXAMPLES::

        sage: from pysymmetry.util import generators1d
        sage: generators1d(4)
        ['(1,4)(2,3)']
    """
    string_reflexao_sigma = ''
    for j in range(1,(n/2+1).floor()):        
            string_reflexao_sigma = string_reflexao_sigma + str((j,n-j+1))            
    return [string_reflexao_sigma]    
 
    


def generators(n):#Nota:temporario
    r"""
    Build permutation generator strings for an n x n grid (temporary helper).

    Returns two strings:
    - sigma: product of transpositions swapping j with n-j+1 along each row
    - miDp: product of transpositions that reflect across the main diagonal

    INPUT:
    - n -- positive integer, grid side length

    OUTPUT:
    - list -- [sigma_string, miDp_string]

    EXAMPLES::

        sage: from pysymmetry.util import generators
        sage: generators(2)
        ['(1,2)(3,4)', '(2,3)']
    """

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
    r"""
    Return the 2D five-point Laplacian on an n x n grid (CSC sparse).

    Uses a block-diagonal structure with off-diagonal blocks to couple
    vertical neighbors; diagonal blocks contain the standard 1D stencil.

    INPUT:
    - n -- positive integer, grid side length

    OUTPUT:
    - SciPy csc_matrix of shape $(n^2, n^2)$ representing -Δ on the grid

    EXAMPLES::

        sage: from pysymmetry.util import laplacian2d
        sage: A = laplacian2d(3)
        sage: A.get_shape()
        (9, 9)
    """
    Bdiag = -4 * np.eye(n)
    Bupper = np.diag([1] * (n - 1), 1)
    Blower = np.diag([1] * (n - 1), -1)
    B = Bdiag + Bupper + Blower    
    blst = [B] * n
    
    A = block_diag(*blst)    
    Dupper = np.diag(np.ones(n * (n - 1)), n)
    Dlower = np.diag(np.ones(n * (n - 1)), -n)
    A += Dupper + Dlower
    return csc_matrix(-A)

from scipy.sparse import spdiags

def laplacian1d(n, h=1):
    r"""
    Return the 1D three-point Laplacian on a line of length n (CSC sparse).

    The stencil is (1, -2, 1) scaled by 1/h^2.

    INPUT:
    - n -- positive integer, number of grid points
    - h -- grid spacing (default: 1)

    OUTPUT:
    - SciPy csc_matrix of shape (n, n)

    EXAMPLES::

        sage: from pysymmetry.util import laplacian1d
        sage: A = laplacian1d(4)
        sage: A.get_shape()
        (4, 4)
    """
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
    r"""
    Build the 2D advection-diffusion operator on an n x n grid (CSC sparse).

    The operator is A = D*(-Δ) - vx*Dx, where -Δ is the 2D five-point
    Laplacian and Dx is a centered first-difference in the x-direction.

    INPUT:
    - n -- positive integer, grid side length
    - D -- diffusion coefficient (default: 1.0)
    - vx -- advection velocity in x (default: 1.0)

    OUTPUT:
    - SciPy csc_matrix of shape (n^2, n^2)

    EXAMPLES::

        sage: from pysymmetry.util import advection_diffusion_2d
        sage: A = advection_diffusion_2d(3)
        sage: A.get_shape()
        (9, 9)
    """
   
    L = -1 * laplacian2d(n)     
    diagonals_x = [np.ones(n - 1), -np.ones(n - 1)]
    offsets_x = [1, -1]
    Dx_1D = spdiags(diagonals_x, offsets_x, n, n, format="csc")    
    I_y = identity(n, format="csc")
    Ax = kron(I_y, Dx_1D) * 0.5
    A_ad = D * L - vx * Ax
    
    return A_ad

def schrodinger_complex_potential_2d(n, V0=1.0, h=1.0):
    """
    Gera a matriz para a equação de Schrödinger com potencial complexo D4-equivariante.
    O Hamiltoniano é H = -∇² + i*V(x,y), com V(x,y) = V0*(x² + y²).
    Esta matriz é G-equivariante (para D4) e não-hermitiana.

    - n: Tamanho da grade (n x n).
    - V0: Amplitude do potencial.
    - h: Espaçamento da grade.
    """
    L = laplacian2d(n) 

    coords = np.linspace(-(n-1)*h/2, (n-1)*h/2, n)
    x_coords = np.tile(coords, n)
    y_coords = np.repeat(coords, n)
       
    V_values = V0 * (x_coords**2 + y_coords**2)   
    
    V_matrix = diags(1j * V_values, format="csc")   
    H = L + V_matrix
    
    return H