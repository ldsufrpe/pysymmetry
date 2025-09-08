__all__ = [
 'nFiniteGroup',
 'nfinitegroup',  
  'ninner_product',
 'ndegree'] 



from sage.all import *
from sage.interfaces.gap import GapElement
from sage.interfaces.gap import gap


import numpy as np
from scipy.sparse import identity
from scipy.sparse import csr_matrix, coo_matrix, csc_matrix
from collections import OrderedDict

from .parallel import pmap
from .util import to_csr, to_csc
from .pysymmetry import avoid_inverse, check



class gcsr_matrix(csc_matrix):
    r"""
    Sparse matrix wrapper (CSC) with convenience methods for group computations.

    This class inherits from SciPy's csc_matrix and adds helpers
    commonly used in representation-theoretic workflows, such as
    symmetry checks, Sage conversion, and trace/character.

    NOTES:
    - Numeric type: entries are assumed to be numeric (float/int). When converting
      to a Sage matrix with tomatrix(), values are cast to QQ (rational). Adapt
      if you need complex/other exact fields.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def tomatrix(self):
        """
        Convert to a Sage dense/sparse matrix with rational entries.

        OUTPUT:
        - Sage matrix (over QQ) built from the nonzero entries of this sparse matrix.
        """
        dic = self.todok()
        dic_new = dict()
        for key, value in dic.items(): 
            dic_new[key] = QQ(int(value))  # NOTE: adjust for complex if needed
        return matrix(dic_new)

    def issymmetric(self, rtol=1e-10, atol=1e-12):
        """
        Check numerical symmetry A ≈ A^T using NumPy tolerance parameters.

        INPUT:
        - rtol, atol -- relative and absolute tolerances for np.allclose

        OUTPUT:
        - bool -- True if the matrix is numerically symmetric.
        """
        arr = self.toarray()
        return np.allclose(arr, arr.T, rtol=rtol, atol=atol)

    def character(self):
        """
        Return the numerical trace, used as a character value.

        OUTPUT:
        - scalar -- sum of diagonal entries
        """
        return self.diagonal().sum()

class nIsotypicBase(object):
    r"""
    Container for a numerically computed isotypic basis.

    INPUT:
    - base -- list of entries of the form [V0, (degree, multiplicity)] where
      V0 is a sparse matrix whose columns span an isotypic component.

    OUTPUT:
    - An object that stores the list of isotypic components and provides utilities
      to extract blocks of equivariant matrices.

    EXAMPLES::

        sage: # Suppose 'base' was returned by nFiniteGroup.nbase_change_reduction(...)
        sage: # You can list the components and get blocks of an equivariant matrix A
        sage: # base_obj.get_blocks(A)
    """
    def __init__(self, base):
        self._base = list(filter(lambda x: x != None, base))
        self.nisotypic_components = [b[0] for b in self._base]
        self._info = [b[1] for b in self._base]

    def __repr__(self):
        msg = []
        for d, m in self._info:
               msg.append("{} blocks of size  {} x {}\n".format(d, m, m))
            
        
        return ''.join(tuple(msg))

    def __len__(self):
        return len(self._base) 

    
    def get_blocks(self,  matrix_equiv):

        blocks = pmap(lambda b: nget_block(b, matrix_equiv),  self.nisotypic_components)
        return blocks
    
    def list(self):
        return  self.nisotypic_components





def nfactorization(g, words, matrices):    
    r"""
    Factorize a group element g in terms of given generators and map to matrices.

    INPUT:
    - g -- group element (in the permutation group)
    - words -- list of group generators (as GAP/Sage elements)
    - matrices -- list of sparse matrices (same length as words) representing
      the images of the generators.

    OUTPUT:
    - SciPy sparse matrix corresponding to the image of g under the representation
      defined by sending each generator to the corresponding matrix.

    NOTE:
    - Uses GAP through Sage to express g in the free group on the generators, then
      evaluates this word by replacing generators with the provided matrices.
    """
    # g = words[0].parent()(g)
    H = libgap.Group(words)
    ans = H.EpimorphismFromFreeGroup().PreImagesRepresentative(g)
    l1 = str(ans)
    l1 = avoid_inverse(g, words, l1)
    l2 = l1.replace('^', '**')

    def locals2():
        local = dict()
        for i, w_i in enumerate(matrices):
            local['x' + str(i+1)] = w_i
        return local
    return sage_eval(l2, locals=locals2())




class nFiniteGroup(PermutationGroup_generic):
    r"""
    A finite group wrapper optimized for numerical (sparse) computations.

    INPUT:
    - gens -- generators describing the group (e.g., from a Sage permutation group)
    - field -- base field for numeric work (default: RDF)

    This class provides methods to build sparse numerical representations via
    GAP factorization and to compute numerical isotypic bases.
    """

    def __init__(self, gens, field=RDF):
        self.field = field
        super().__init__(gens, gap_group=None, domain=None, canonicalize=True, category=None)

    def nrepresentation(self, gens, image):
        r"""
        Build a dictionary representation g -> matrix using generator images.

        INPUT:
        - gens -- list of group generators (compatible with GAP)
        - image -- list of matrices (SciPy sparse) for the images of gens

        OUTPUT:
        - dict mapping each group element to a gcsr_matrix

        EXAMPLES::

            sage: # Given generators and their matrix images, construct the full rep
            sage: # rep = G.nrepresentation(generators, matrices)
        """
        idty_group = self.identity()
        image = [to_csc(m) for m in image]
        dic = dict()
        _, n = image[0].get_shape()
        elements = self.list()
        for g in elements[1:]:  # except identity
            D = nfactorization(g, gens, image)
            dic[g] = gcsr_matrix(D)
        dic[idty_group] = gcsr_matrix(identity(n, format='csc'))
        return dic

    def nirreducible_representations(self): 
        r"""
        Compute irreducible representations numerically via GAP.

        OUTPUT:
        - (N, irr) where N is the number of irreps and irr is an OrderedDict
          mapping index -> dict representation (g -> gcsr_matrix)

        NOTE:
        - Uses GAP.IrreducibleRepresentations and then converts generator images
          into sparse numerical matrices.
        """
        irr = gap.IrreducibleRepresentations(self)
        gens = gap.GeneratorsOfGroup(self)
        generators = list(gap.List(gens))
        list_irr = OrderedDict()
        for i, s in enumerate(irr):
            image = [matrix(sage_eval(str(gap.Image(s, g)))) for g in generators]
            rep = self.nrepresentation(generators, image)
            list_irr[i] = rep
        N = len(list_irr)
        return N, list_irr

    def nprojection(self, i, j, k, all_irr, right):
        r"""
        Compute the (j,k)-entry of the projection onto the i-th isotypic component.

        INPUT:
        - i -- index of the irreducible component in all_irr
        - j, k -- entry indices for the projected block
        - all_irr -- dictionary of irreducible reps (as from nirreducible_representations)
        - right -- a representation (dict g -> matrix) to be projected

        OUTPUT:
        - sparse matrix corresponding to sum_g left_i(g^{-1})_{j,k} right(g)
        """
        left = all_irr
        s = sum([left[i][g.inverse()][j, k] * right[g] for g in self])
        return s
    
    def nisotypic_component(self, i, left, right):
        r"""
        Compute a basis for the i-th isotypic component of 'right'.

        INPUT:
        - i -- index of irreducible component
        - left -- dictionary of irreducible reps (from nirreducible_representations)
        - right -- target representation (dict g -> matrix)

        OUTPUT:
        - [V0, (degree, multiplicity)] where V0 columns span the isotypic space.
        """
        multiplicity = ninner_product(self, left[i], right)
        if not multiplicity == 0:
            P = self.nprojection(i, 0, 0, left, right)
            P_copy = csr_matrix.copy(P)
            pivots = P.tomatrix().pivots()
            degree = ndegree(self, left[i])
            v0 = P_copy[:, pivots]
            return [v0, (degree, multiplicity)]

    def nbase_change_reduction(self, right):
        r"""
        Build the numerical isotypic basis for a given representation 'right'.

        INPUT:
        - right -- dict representation (g -> matrix)

        OUTPUT:
        - nIsotypicBase object that allows extraction of blocks for equivariant matrices.
        """
        n, l = self.nirreducible_representations()
        r = right
        base = pmap(lambda i: self.nisotypic_component(i, l, r), list(range(n)))
        return nIsotypicBase(base)

nfinitegroup = nFiniteGroup



import numpy as np

def ninner_product(G, left, right):
    r"""
    Inner product of characters using conjugacy class aggregation.

    INPUT:
    - G -- group (nFiniteGroup or Sage permutation group)
    - left -- dictionary representation (g -> sparse matrix) for the left rep
    - right -- dictionary representation (g -> sparse matrix) for the right rep

    OUTPUT:
    - int -- the character inner product <chi_left, chi_right>

    EXAMPLES::

        sage: # For a regular representation 'reg', <reg, reg> equals |G|
        sage: # ninner_product(G, reg, reg)
    """
    s = 0.0  
       
    for conj_class in G.conjugacy_classes():        
        g = conj_class.representative()        
        class_size = len(conj_class)
        char_left_inv = left[g.inverse()].diagonal().sum()
        char_right = right[g].diagonal().sum()        
        s += class_size * char_left_inv * char_right       
    
    return int((1.0 / G.order()) * s)


def ndegree(G, rep):
    r"""
    Degree (matrix size) of a numerical representation.

    INPUT:
    - G -- group
    - rep -- dictionary representation (g -> sparse matrix)

    OUTPUT:
    - int -- the dimension of the representation space
    """
    _, d = rep[G.an_element()].get_shape()       
    return d

# from scipy.linalg import pinv


def nget_block(columm_base, matrix_equiv): 
    r"""
    Compute the block of an equivariant matrix in a given column-space basis.

    Given columns B spanning an isotypic component and an equivariant matrix A,
    returns P * A * B, where P is the pseudoinverse of the dense array of B.

    INPUT:
    - columm_base -- SciPy sparse matrix (columns are basis vectors)
    - matrix_equiv -- SciPy sparse matrix A that commutes with the group action

    OUTPUT:
    - gcsr_matrix -- the block of A restricted to span(B)

    NOTE:
    - This uses a numerical pseudoinverse (NumPy). For exact arithmetic, adapt
      to Sage exact linear algebra if needed.
    """

    columm = columm_base.toarray()
    p = csc_matrix(np.linalg.pinv(columm))#Nota:Aqui da pra configurar pelo Sage pedido o algoritmo numpy, ou chamar outro metodo.
    block = p*matrix_equiv*columm_base
    return gcsr_matrix(block)



