__all__ = [
    'FiniteGroup',
    'group',
    'MapRepresentation',
    'Dg_Linear_Transformation',
    'representation',        
    'get_block'
]

#version 0.2
from sage.all import *
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.groups.perm_gps.permgroup import PermutationGroup_generic
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.libs.gap.libgap import libgap
from sage.interfaces.gap import gap
from sage.libs.gap.element import GapElement
from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_integer_sparse import Matrix_integer_sparse
from sage.misc.sage_eval import sage_eval

from sage.categories.morphism import SetMorphism
from sage.modules.vector_space_morphism import VectorSpaceMorphism
# from numpy import identity
# from numpy import ndarray, array, where, identity


from .util import  view
from .parallel import pmap




################################## FiniteGroup #############################


class FiniteGroup(PermutationGroup_generic):
    def __init__(self, per_group, field=QQbar, matrix=False):

        if isinstance(per_group, PermutationGroup_generic):
            per_group = per_group.gens()

        if matrix:

            check_matrix = all([isinstance(m, (Matrix_integer_dense, Matrix_integer_sparse)) for m in per_group])

            if check_matrix:

                per_group = MatrixGroup(per_group).as_permutation_group().gens()

            else:
                raise TypeError('no matrices')

        super().__init__(per_group,
                         gap_group=None,
                         domain=None,
                         canonicalize=True,
                         category=None)

        self.field = field

    def _hom_(self, n):
        return Hom(self.field**n, self.field**n)    
    
    def _regular_(self, g):
        elements = self.list()
        new_elements = [g*element for element in elements]
        index = [elements.index(element) for element in new_elements]
        base = matrix.identity(self.order())[:, index]
        #base = base  # sage matrix
        H = self._hom_(self.order())        
        return Dg_Linear_Transformation(H, base, g)  # class Representation

    def regular_representation(self):
        r"""
        Return the regular representation defined over the group.

        INPUT:

        - ``self`` -- FiniteGroup ; a  Sage permutation group or a group of the class FiniteGroup.

        OUTPUT: The regular representation defined over the FiniteGroup given by self. 


        EXAMPLES:

        We define the regular representation over the cyclic group of 4 elements ::

            sage: G = FiniteGroup(CyclicPermutationGroup(4)) 
            sage: reg = G.regular_representation()
            sage: [reg(g) for g in G]
            [Linear transformation associated with element g=(), represented by the non-singular matrix:
             [1 0 0 0]
             [0 1 0 0]
             [0 0 1 0]
             [0 0 0 1]
             Representation space: Vector space of dimension 4 over Algebraic Field
             ,
             Linear transformation associated with element g=(1,2,3,4), represented by the non-singular matrix:
             [0 0 0 1]
             [1 0 0 0]
             [0 1 0 0]
             [0 0 1 0]
             Representation space: Vector space of dimension 4 over Algebraic Field
             ,
             Linear transformation associated with element g=(1,3)(2,4), represented by the non-singular matrix:
             [0 0 1 0]
             [0 0 0 1]
             [1 0 0 0]
             [0 1 0 0]
             Representation space: Vector space of dimension 4 over Algebraic Field
             ,
             Linear transformation associated with element g=(1,4,3,2), represented by the non-singular matrix:
             [0 1 0 0]
             [0 0 1 0]
             [0 0 0 1]
             [1 0 0 0]
             Representation space: Vector space of dimension 4 over Algebraic Field
             ]

        We define the group of symmetries of regular hexagon and the regular representation over this group ::
        
            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: G = FiniteGroup(generators)
            sage: reg = G.regular_representation()
            sage: g = G.an_element()
            sage: reg(g)
            Linear transformation associated with element g=(1,3)(4,6), represented by the non-singular matrix:
            [0 0 0 0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 1]
            [0 0 0 0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 0 0 0 1 0]
            [0 1 0 0 0 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0 0 0 0]
            [1 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0 0 0 0]
            [0 0 0 0 0 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0 0 0 0]
            Representation space: Vector space of dimension 12 over Algebraic Field

        We define the regular representation over the symmetric group of 4 simbols ::

            sage: H = SymmetricGroup(4)
            sage: G = FiniteGroup(H);
            sage: reg = G.regular_representation();
            sage: g = G.an_element()
            sage: reg(g)
            Linear transformation associated with element g=(1,3,4), represented by the non-singular matrix:
            24 x 24 dense matrix over Algebraic Field
            Representation space: Vector space of dimension 24 over Algebraic Field

                
        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """                

        
        # usei metodos do orderpy pra deixar os calculos mais rapidos
        #NOTA: precisa verificar se g esta em FiniteGroup
        matrices = [self._regular_(g).matrix() for g in self.gens()]
        M = MatrixGroup(matrices)
        return MapRepresentation(Hom(self, M), self._regular_)

    def irreducible_representations(self, show_table=True): ##Nota: Como documentar saida com o True
        r"""
        Return the number n of irreducible representations of self and the irreducibles representations.


        INPUT:

        - ``self`` -- FiniteGroup ; a  Sage permutation group or a group of the class FiniteGroup.
        - ``show_table`` -- a boolean (default: `True`) ; shows a table of each irreducible representation applied at the generators of self.

        OUTPUT: If show_table=True; return number of irreducible representions and irreducible representions themselves.
                If show_table=False; return the irreducible representions themselves. 


        EXAMPLES:

        We define the symmetric group of 4 simbols and calculate its irreducible representations ::

            sage: G = FiniteGroup(SymmetricGroup(4));
            sage: n, irr = G.irreducible_representations(False) #Irreducible representations
            sage: n
            5
            sage: irr(0)
            Map: 
             From: Permutation FiniteGroup with generators [(1,2), (1,2,3,4)] 
             To: Matrix group over Integer Ring with 2 generators ([1], [1]).
            
            sage: irr(4)
            Map: 
             From: Permutation FiniteGroup with generators [(1,2), (1,2,3,4)] 
             To: Matrix group over Integer Ring with 2 generators (
            [ 0  1  0]  [ 0  0 -1]
            [ 1  0  0]  [ 0  1  0]
            [ 0  0 -1], [ 1  0  0]
            ).

        We define the representation by permutation on the cyclic group calculate its irreducible representations ::

            sage: G = FiniteGroup(CyclicPermutationGroup(6))
            sage: irr = G.irreducible_representations(True) #Irreducible representations
            ||||||SAIDA A SER PENSADA|||||||

        We calculate the irreducible representations of the group of symmetries of a regular tetrahedron ::
        
            sage: G = FiniteGroup(AlternatingGroup(4));
            sage: irr = G.irreducible_representations(True) #Irreducible representations
            sage: for j in range(n):
             ||||||SAIDA A SER PENSADA|||||||
           
        We define the group of symmetries of regular hexagon and calculate its irreducible representations  ::
            
            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: G = FiniteGroup(PermutationGroup(generators))
            sage: n,irr = G.irreducible_representations(False) #Irreducible representations
            sage: for j in range(n):
                    print(irr(j))
            Map: 
             From: Permutation FiniteGroup with generators [(1,2,3,4,5,6), (1,4)(2,3)(5,6)] 
             To: Matrix group over Integer Ring with 2 generators ([1], [1]).
            Map: 
             From: Permutation FiniteGroup with generators [(1,2,3,4,5,6), (1,4)(2,3)(5,6)] 
             To: Matrix group over Integer Ring with 2 generators ([1], [-1]).
            Map: 
             From: Permutation FiniteGroup with generators [(1,2,3,4,5,6), (1,4)(2,3)(5,6)] 
             To: Matrix group over Integer Ring with 2 generators ([-1], [-1]).
            Map: 
             From: Permutation FiniteGroup with generators [(1,2,3,4,5,6), (1,4)(2,3)(5,6)] 
             To: Matrix group over Integer Ring with 2 generators ([-1], [1]).
            Map: 
             From: Permutation FiniteGroup with generators [(1,2,3,4,5,6), (1,4)(2,3)(5,6)] 
             To: Matrix group over Universal Cyclotomic Field with 2 generators (
            [E(3)^2      0]  [0 1]
            [     0   E(3)], [1 0]
            ).
            Map: 
             From: Permutation FiniteGroup with generators [(1,2,3,4,5,6), (1,4)(2,3)(5,6)] 
             To: Matrix group over Universal Cyclotomic Field with 2 generators (
            [-E(3)^2       0]  [ 0 -1]
            [      0   -E(3)], [-1  0]
            ).
               
        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """
        
        
        irr = gap.IrreducibleRepresentations(self)
        gens = gap.GeneratorsOfGroup(self)
        generators = list(gap.List(gens))
        matrices_list = []
        for s in irr:
            matrices = [
                    matrix(sage_eval(str(gap.Image(s, k)))) for k in generators
                    ]
            matrices_list.append(matrices)

        def irr_rep(n):            
            r = representation(generators, matrices_list[n], self.field)
            return r

        N = len(matrices_list)
        if show_table:

            n = self.ngens()
            elements = self.gens()

            matrices_rep = [
                view((irr_rep(k)(x)).matrix(), False) for k in range(N)
                for x in elements
                        ]
            _table = [matrices_rep[i:i + n] for i in range(0, len(matrices_rep), n)]
            header_row = [str(x) for x in elements]
            header_column = [''.join(('$r^', str(k), '$')) for k in range(N)]
            header_column.insert(0, ' ')

            t = table(_table,
                      header_row=header_row,
                      header_column=header_column,
                      align='center')
            show(t)

            return irr_rep

        return N, irr_rep

    def natural_representation(self, field=QQbar):
        r"""
        Return the natural representation of the group.

        The natural representation is defined by the permutation matrices
        associated with each element of the group.

        INPUT:

        - ``field`` -- (default: `QQbar`); The field over which the representation is defined.

        OUTPUT: A representation of the group.

        EXAMPLES:

            sage: G = FiniteGroup(CyclicPermutationGroup(4))
            sage: nat_rep = G.natural_representation()
            sage: g = G.an_element()
            sage: nat_rep(g).matrix()
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            [1 0 0 0]
        """
        gens = self.gens()
        image = [g.matrix() for g in gens]
        return representation(gens, image, field)
    	

    def isotypic_projection(self, right):
        r"""
        Return a list containing the matrices associated to projections operators over the isotypic components(or conglomerates) of right.


        INPUT:

        - ``self`` -- FiniteGroup ; a  Sage permutation group or a group of the class FiniteGroup.
        - ``right`` -- representation ; a representation defined over the same group as self.

        OUTPUT: A list of matrices representing the projections operators over the isotypics subespaces of right.


        EXAMPLES:

        We define the symmetric group of 3 simbols and the isotypic projections associated to the regular representation of this group ::

            sage: G = FiniteGroup(SymmetricGroup(3));
            sage: reg = G.regular_representation();
            sage: G.isotypic_projection(reg)
             [
             [1 1 1 1 1 1]  [ 1  1  1 -1 -1 -1]  [ 2 -1 -1  0  0  0]
             [1 1 1 1 1 1]  [ 1  1  1 -1 -1 -1]  [-1  2 -1  0  0  0]
             [1 1 1 1 1 1]  [ 1  1  1 -1 -1 -1]  [-1 -1  2  0  0  0]
             [1 1 1 1 1 1]  [-1 -1 -1  1  1  1]  [ 0  0  0  2 -1 -1]
             [1 1 1 1 1 1]  [-1 -1 -1  1  1  1]  [ 0  0  0 -1  2 -1]
             [1 1 1 1 1 1], [-1 -1 -1  1  1  1], [ 0  0  0 -1 -1  2]
              ]
        
        We calculate an irreducible representation of G, and the isotypic projection associated ::

            sage: n,irr = G.irreducible_representations(False) #Irreducible representations
            sage: G.isotypic_projection(irr(2))
            [
            [3 0]
            [0 3]
            ]


        We define the representation by permutation on the cyclic group and calculate the isotypic projectors ::

            sage: G = FiniteGroup(CyclicPermutationGroup(4))
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: G.isotypic_projection(rep)
            [
            [1 1 1 1]  [ 1 -1  1 -1]  [   1 -1*I   -1  1*I]  [   1  1*I   -1 -1*I]
            [1 1 1 1]  [-1  1 -1  1]  [ 1*I    1 -1*I   -1]  [-1*I    1  1*I   -1]
            [1 1 1 1]  [ 1 -1  1 -1]  [  -1  1*I    1 -1*I]  [  -1 -1*I    1  1*I]
            [1 1 1 1], [-1  1 -1  1], [-1*I   -1  1*I    1], [ 1*I   -1 -1*I    1]
            ]
        

        We define a representation the group of symmetries of regular hexagon and calculate their isotypic projectors ::
            
            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: rotation_pi_over_3, reflexion_about_x = matrix([[1/2,-sqrt(3)/2],[sqrt(3)/2,1/2]]),matrix([[-1,0],[0,1]])
            sage: matrices = [block_diagonal_matrix(rotation_pi_over_3, rotation_pi_over_3),block_diagonal_matrix(reflexion_about_x, reflexion_about_x)];
            sage: rep = representation(generators, matrices, field=SR)
            sage: G = rep.domain()
            sage: G.isotypic_projection(rep)
            [
            [6 0 0 0]
            [0 6 0 0]
            [0 0 6 0]
            [0 0 0 6]
            ]
                
        TO ADD MORE EXAMPLES
        TO ADD MORE EXAMPLES
        TO ADD MORE EXAMPLES
              
        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """
        n, left = self.irreducible_representations(False)
        proj_lst = []
        for i in range(n):
            if left(i).inner_product(right) != self.field.zero():
                s = sum([conjugate(left(i)(g).character())*right(g).matrix()  for g in self])
                proj_lst.append(s)

        return proj_lst

    def isotypic_base(self, rep, isotypic_components=False):
        r"""
        Return a list with the basis for each isotypic components of rep or the basis change matrix associated to the isotypic decomposition.


        INPUT:

        - ``self`` -- FiniteGroup ; a  Sage permutation group or a group of the class FiniteGroup.
        - ``rep`` -- representation ; a representation defined over the same group as self.
        - ``isotypic_components`` -- a boolean (default: `False`) ; if False returns the change base matrix, if True returns a list with a base for each isotypic component :: 

        OUTPUT: A list with a basis for each isotypic subspace of rep or the change basis matrix to a form in blocks associated to isotypic components of rep.


        EXAMPLES:

        We define the symmetric group of 3 simbols change basis matrix associated to the isotypic decomposition  ::

            sage: G = FiniteGroup(SymmetricGroup(3));
            sage: reg = G.regular_representation();
            sage: G.isotypic_base(reg, isotypic_components=False)
            [ 1| 1| 2 -1  0  0]
            [ 1| 1|-1  2  0  0]
            [ 1| 1|-1 -1  0  0]
            [ 1|-1| 0  0  2 -1]
            [ 1|-1| 0  0 -1  2]
            [ 1|-1| 0  0 -1 -1]

        
        The isotypic basis for each irreducible subspace associated to the regular representation of this group ::

            sage: G.isotypic_base(reg, isotypic_components=True)
            [
            [1]  [ 1]  [ 2 -1  0  0]
            [1]  [ 1]  [-1  2  0  0]
            [1]  [ 1]  [-1 -1  0  0]
            [1]  [-1]  [ 0  0  2 -1]
            [1]  [-1]  [ 0  0 -1  2]
            [1], [-1], [ 0  0 -1 -1]
            ]

        We define the representation by permutation on the cyclic group and calculate its isotypic decompositions ::

            sage: G = FiniteGroup(CyclicPermutationGroup(4))
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: G.isotypic_base(rep, isotypic_components=True)
            [
            [1]  [ 1]  [   1]  [   1]
            [1]  [-1]  [ 1*I]  [-1*I]
            [1]  [ 1]  [  -1]  [  -1]
            [1], [-1], [-1*I], [ 1*I]
            ]

            sage: change_basis = G.isotypic_base(rep, isotypic_components=False); change_basis
            [   1|   1|   1|   1]
            [   1|  -1| 1*I|-1*I]
            [   1|   1|  -1|  -1]
            [   1|  -1|-1*I| 1*I]
            sage: [rep(g).matrix() for g in G]
            [
            [1 0 0 0]  [0 1 0 0]  [0 0 1 0]  [0 0 0 1]
            [0 1 0 0]  [0 0 1 0]  [0 0 0 1]  [1 0 0 0]
            [0 0 1 0]  [0 0 0 1]  [1 0 0 0]  [0 1 0 0]
            [0 0 0 1], [1 0 0 0], [0 1 0 0], [0 0 1 0]
            ]
            sage: [change_basis.inverse()*rep(g).matrix()*change_basis for g in G]
            [
            [1 0 0 0]  [   1    0    0    0]  [ 1  0  0  0]  [   1    0    0    0]
            [0 1 0 0]  [   0   -1    0    0]  [ 0  1  0  0]  [   0   -1    0    0]
            [0 0 1 0]  [   0    0  1*I    0]  [ 0  0 -1  0]  [   0    0 -1*I    0]
            [0 0 0 1], [   0    0    0 -1*I], [ 0  0  0 -1], [   0    0    0  1*I]
            ]
            
        We define a representation the group of symmetries of regular hexagon and calculate their matrix change basis(this representation is irreducible) ::

            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: matrices = [matrix([[1/2,-sqrt(3)/2],[sqrt(3)/2,1/2]]),matrix([[-1,0],[0,1]])];
            sage: rep = representation(generators, matrices)
            sage: G = rep.domain()
            sage: G.isotypic_base(rep, isotypic_components=True)
            [
            [6 0]
            [0 6]
            ]
            
                
        TO ADD MORE EXAMPLES
        TO ADD MORE EXAMPLES
        TO ADD MORE EXAMPLES
              
        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """
        proj_lst = self.isotypic_projection(rep)  #NOTA:devemos usar sempre regular representation(pergunta)??? NAO
        decomposition = [m.matrix_from_columns(m.pivots()) for m in proj_lst]
        if isotypic_components:
            return decomposition
        ncols = len(decomposition)
        return block_matrix(decomposition, ncols=ncols)

    def projection(self, i, j, k, right, left=None):
        r"""
        Return the projection(or transfers operators) associated to irreducible subrepresentations of right.

        INPUT:

        - ``self`` -- FiniteGroup ; a  Sage permutation group or a group of the class FiniteGroup.
        - ``right`` -- representation ; a representation defined over the same group as self.
        - ``i,j,k`` --  integers ; integers numbers representing the projection(or transfer) chosen;
                i- will choose the irreducible representation, and its range is from 0 until the number of irreducibles. ::
                
                j- choose a row of the matrix presentation of the irreducible, its range is the degree of chosen irreducible according i. ::
                
                k- choose a column of the matrix presentation of the irreducible, its range is the degree of chosen irreducible according i. ::
                

        OUTPUT: A matrix associated to the operator that projects over an irreducible subrepresentation of right(if j=k), or is a isomorphism between two equivalent subrepresentations of right.
                Whenever the irreducible, chosen by index i, it is not a subrepresentation of right this matrix will be null.


        EXAMPLES:


        We define the representation by permutation on the cyclic group of four elements and calculate its projectors  ::
        
            sage: G = FiniteGroup(CyclicPermutationGroup(4))
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: G.projection(0,0,0, rep)
            [1 1 1 1]
            [1 1 1 1]
            [1 1 1 1]
            [1 1 1 1]
            sage: n, irr = G.irreducible_representations(False)
            sage: list_of_projectors=[];
            for i in range(n): # Choosing the irreducible
                degree= irr(i).degree() # Calculation the degree of the irreducible
                j=0;# we will use the first row
                for k in range(degree):
                    projector = G.projection(i,j,k, rep);
                    list_of_projectors.append(projector)
            sage: list_of_projectors
            [
            [1 1 1 1]  [ 1 -1  1 -1]  [ 1 -I -1  I]  [ 1  I -1 -I]
            [1 1 1 1]  [-1  1 -1  1]  [ I  1 -I -1]  [-I  1  I -1]
            [1 1 1 1]  [ 1 -1  1 -1]  [-1  I  1 -I]  [-1 -I  1  I]
            [1 1 1 1], [-1  1 -1  1], [-I -1  I  1], [ I -1 -I  1]
            ]


        We define a representation on the group of symmetries of regular hexagon and calculate its projectors ::

            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: matrices = [matrix([[1,0,0],[0,1/2,-sqrt(3)/2],[0,sqrt(3)/2,1/2]]),matrix([[1,0,0],[0,-1,0],[0,0,1]])];
            sage: rep = representation(generators, matrices)
            sage: G = rep.domain()
            sage: G.projection(0,0,0, rep)
            [     12       0       0]
            [      0       0 0.?e-18]
            [      0 0.?e-18       0]
            sage: n, irr = G.irreducible_representations(False)
            sage: list_of_projectors=[];
            for i in range(n): # Choosing the irreducible
                degree= irr(i).degree() # Calculation the degree of the irreducible
                j=0;# we will use the first row
                for k in range(degree):
                    projector = G.projection(i,j,k, rep);
                    list_of_projectors.append(view(projector, latex=False)) # view creates a better visualization
            sage: list_of_projectors
            [
            [12  0  0]  [0 0 0]  [0 0 0]  [0 0 0]  [0 0 0]  [0 0 0]
            [ 0  0  0]  [0 0 0]  [0 0 0]  [0 0 0]  [0 0 0]  [0 0 0]
            [ 0  0  0], [0 0 0], [0 0 0], [0 0 0], [0 0 0], [0 0 0],

            [   0    0    0]  [  0   0   0]
            [   0    3  3*I]  [  0   3 3*I]
            [   0 -3*I    3], [  0 3*I  -3]
            ]


        We define the regular representation over the symmetric group of 3 simbols and we calculate one of its projectors ::

            sage: G = FiniteGroup(SymmetricGroup(3)) 
            sage: reg = G.regular_representation();
            sage: view(G.projection(2,1,1, reg),latex=False) # The function view creates a better visualization
            [                   1  1/2*I*sqrt(3) - 1/2 -1/2*I*sqrt(3) - 1/2                    0                    0                    0]
            [-1/2*I*sqrt(3) - 1/2                    1  1/2*I*sqrt(3) - 1/2                    0                    0                    0]
            [ 1/2*I*sqrt(3) - 1/2 -1/2*I*sqrt(3) - 1/2                    1                    0                    0                    0]
            [                   0                    0                    0                    1 -1/2*I*sqrt(3) - 1/2  1/2*I*sqrt(3) - 1/2]
            [                   0                    0                    0  1/2*I*sqrt(3) - 1/2                    1 -1/2*I*sqrt(3) - 1/2]
            [                   0                    0                    0 -1/2*I*sqrt(3) - 1/2  1/2*I*sqrt(3) - 1/2                    1]

        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """
        
        
        #Nota:   ver como melhorar a posição do left
        #Nota: right pode inclusive ser a irredutivel (pergunta)
        ###Melhorar este nome opcoes: transfer isomorphism 
        if left == None:
            n, left = self.irreducible_representations(False)
        s = sum([matrix(left(i)(g.inverse()))[j][k] * matrix(right(g)) for g in self])  #Nota:
        return s

    def base_to_irreducibles(self, right):
        r"""
        Return the basis change matrix associated to a decomposition of right into irreducibles representations.


        INPUT:

        - ``self`` -- FiniteGroup ; a  Sage permutation group or a group of the class FiniteGroup.
        - ``right`` -- representation ; a representation defined over the same group as self.

        OUTPUT: A change basis matrix that decomposes right in its irreducible components.


        EXAMPLES:


        We define the representation by permutation on the cyclic group of four elements and calculate its decomposition into irreducibles(Note that in this case because multiplicity we get the same result with the isotypic base)  ::
        
            sage: G = FiniteGroup(CyclicPermutationGroup(4))
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: A = G.base_to_irreducibles(rep); A
            [ 1  1  1  1]
            [ 1 -1  I -I]
            [ 1  1 -1 -1]
            [ 1 -1 -I  I]

            sage: [rep(g).matrix() for g in G]
            [
            [1 0 0 0]  [0 1 0 0]  [0 0 1 0]  [0 0 0 1]
            [0 1 0 0]  [0 0 1 0]  [0 0 0 1]  [1 0 0 0]
            [0 0 1 0]  [0 0 0 1]  [1 0 0 0]  [0 1 0 0]
            [0 0 0 1], [1 0 0 0], [0 1 0 0], [0 0 1 0]
            ]

            sage: [A.inverse()*rep(g).matrix()*A for g in G]
            [
            [1 0 0 0]  [ 1  0  0  0]  [ 1  0  0  0]  [ 1  0  0  0]
            [0 1 0 0]  [ 0 -1  0  0]  [ 0  1  0  0]  [ 0 -1  0  0]
            [0 0 1 0]  [ 0  0  I  0]  [ 0  0 -1  0]  [ 0  0 -I  0]
            [0 0 0 1], [ 0  0  0 -I], [ 0  0  0 -1], [ 0  0  0  I]
            ]



        We define a representation on the group of symmetries of regular hexagon and calculate its matrix change basis to decomposing into irreducible blocks ::
            
            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: P = matrix([[5,3,4],[0,3,4],[0,0,4]]);P
            sage: matrices = [P.inverse()*matrix([[1,0,0],[0,1/2,-sqrt(3)/2],[0,sqrt(3)/2,1/2]])*P,P.inverse()*matrix([[1,0,0],[0,-1,0],[0,0,1]])*P];
            sage: rep = representation(generators, matrices, field=QQbar)
            sage: G = rep.domain()
            sage: B = G.base_to_irreducibles(rep); 
            sage: view(B,latex=False) # This function creates a better visualization of the matrix
            [        12       -9/5      -54/5]
            [         0    3*I + 3 -18*I + 18]
            [         0     -9/4*I     27/2*I]
            
            sage: [view(rep(g).matrix(), latex=False) for g in G]
            [
            [1 0 0]  [                 1               9/10 -2/5*sqrt(3) + 6/5]
            [0 1 0]  [                 0  1/2*sqrt(3) - 1/2        4*sqrt(1/3)]
            [0 0 1], [                 0       -3/8*sqrt(3) -1/2*sqrt(3) - 1/2],

            [                 1               9/10  2/5*sqrt(3) + 6/5]
            [                 0 -1/2*sqrt(3) - 1/2       -4*sqrt(1/3)]
            [                 0        3/8*sqrt(3)  1/2*sqrt(3) - 1/2],

            [                 1               3/10 -2/5*sqrt(3) + 2/5]
            [                 0  1/2*sqrt(3) + 1/2        4*sqrt(1/3)]
            [                 0       -3/8*sqrt(3) -1/2*sqrt(3) + 1/2],

            [  1 6/5 8/5]
            [  0  -1   0]
            [  0   0  -1],

            [                 1               3/10  2/5*sqrt(3) + 2/5]
            [                 0 -1/2*sqrt(3) + 1/2       -4*sqrt(1/3)]
            [                 0        3/8*sqrt(3)  1/2*sqrt(3) + 1/2],

            [  1   0   0]
            [  0   1 8/3]
            [  0   0  -1],

            [                 1               9/10 -2/5*sqrt(3) + 6/5]
            [                 0 -1/2*sqrt(3) - 1/2               -4/3]
            [                 0        3/8*sqrt(3)  1/2*sqrt(3) + 1/2],

            [                 1               9/10  2/5*sqrt(3) + 6/5]
            [                 0  1/2*sqrt(3) - 1/2               -4/3]
            [                 0       -3/8*sqrt(3) -1/2*sqrt(3) + 1/2],

            [                 1               3/10 -2/5*sqrt(3) + 2/5]
            [                 0 -1/2*sqrt(3) + 1/2                4/3]
            [                 0        3/8*sqrt(3)  1/2*sqrt(3) - 1/2],

            [   1  6/5  8/5]
            [   0   -1 -8/3]
            [   0    0    1],

            [                 1               3/10  2/5*sqrt(3) + 2/5]
            [                 0  1/2*sqrt(3) + 1/2                4/3]
            [                 0       -3/8*sqrt(3) -1/2*sqrt(3) - 1/2]
            ]

            sage: [view(B.inverse()*rep(g).matrix()*B, latex=False) for g in G]
            [
            [1 0 0]
            [0 1 0]
            [0 0 1],

            [                   1                    0                    0]
            [                   0 -1/2*I*sqrt(3) - 1/2                    0]
            [                   0                    0  1/2*I*sqrt(3) - 1/2],

            [                   1                    0                    0]
            [                   0  1/2*I*sqrt(3) - 1/2                    0]
            [                   0                    0 -1/2*I*sqrt(3) - 1/2],

            [                   1                    0                    0]
            [                   0 -1/2*I*sqrt(3) + 1/2                    0]
            [                   0                    0  1/2*I*sqrt(3) + 1/2],

            [ 1  0  0]
            [ 0 -1  0]
            [ 0  0 -1],

            [                   1                    0                    0]
            [                   0  1/2*I*sqrt(3) + 1/2                    0]
            [                   0                    0 -1/2*I*sqrt(3) + 1/2],

            [  1   0   0]  [                      1                       0                       0]
            [  0   0   6]  [                      0                       0         3*I*sqrt(3) - 3]
            [  0 1/6   0], [                      0 -1/4*I*sqrt(1/3) - 1/12                       0],

            [                     1                      0                      0]  [                      1                       0                       0]
            [                     0                      0       -3*I*sqrt(3) - 3]  [                      0                       0         3*I*sqrt(3) + 3]
            [                     0 1/4*I*sqrt(1/3) - 1/12                      0], [                      0 -1/4*I*sqrt(1/3) + 1/12                       0],

            [   1    0    0]
            [   0    0   -6]
            [   0 -1/6    0],

            [                     1                      0                      0]
            [                     0                      0       -3*I*sqrt(3) + 3]
            [                     0 1/4*I*sqrt(1/3) + 1/12                      0]
            ]




        We define the regular representation over the symmetric group of 3 simbols and we decompose into irreducible components ::

            sage: G = FiniteGroup(SymmetricGroup(3)) 
            sage: reg = G.regular_representation();
            sage: C = G.base_to_irreducibles(reg); 
            sage: view(C, latex=False) # This Function creates a better visualization of the matrix
            [                   1                    1                    1                    0                    0                    3]
            [                   1                    1  1/2*I*sqrt(3) - 1/2                    0                    0 -3/2*I*sqrt(3) - 3/2]
            [                   1                    1 -1/2*I*sqrt(3) - 1/2                    0                    0  3/2*I*sqrt(3) - 3/2]
            [                   1                   -1                    0                    3                    1                    0]
            [                   1                   -1                    0  3/2*I*sqrt(3) - 3/2 -1/2*I*sqrt(3) - 1/2                    0]
            [                   1                   -1                    0 -3/2*I*sqrt(3) - 3/2  1/2*I*sqrt(3) - 1/2                    0]

            sage: g = G.an_element();
            sage: reg(g).matrix()
            [0 0 0 0 0 1]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [0 1 0 0 0 0]
            [0 0 1 0 0 0]
            [1 0 0 0 0 0]

            sage: A = C.inverse()*reg(g).matrix()*C;
            sage: view(A, latex=False)
            [                    1                     0                     0                     0                     0                     0]
            [                    0                    -1                     0                     0                     0                     0]
            [                    0                     0                     0  -3/2*I*sqrt(3) - 3/2                     0                     0]
            [                    0                     0 1/2*I*sqrt(1/3) - 1/6                     0                     0                     0]
            [                    0                     0                     0                     0                     0  -3/2*I*sqrt(3) - 3/2]
            [                    0                     0                     0                     0 1/2*I*sqrt(1/3) - 1/6                     0]


        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """
        
        
        row = 0
        n, left = self.irreducible_representations(False)  #Nota: vou pensar em como implementar o __len__
        base = []
       
        for i in range(n):
            multiplicity = left(i).inner_product(right) 
            if multiplicity != self.field.zero():  #Nota: acho que aqui é um problema em caso de negativa
                P = self.projection(i, row, 0, right, left)  #Nota: conversar com marcelo a melhor forma de deixar row opcional
                pivots = P.pivots()
                for pivot in pivots:
                    v0 = P[:, pivot]
                    degree = left(i).degree()
                    base.append(v0)
                    
                    for k in range(1, degree):
                        v = self.projection(i, row, k, right, left) * v0
                        
                        base.append(v)       
        b = base[0]
        for v in base[1:]:
            b = b.augment(v) 
        
        return b


    def base_equivariant_to_blocks(self, right, row=0):
        r"""
        Return the basis change matrix associated to a symmetry adapted basis to an equivariant operator of right.


        INPUT:

        - ``self`` -- FiniteGroup ; a  Sage permutation group or a group of the class FiniteGroup.
        - ``right`` -- representation ; a representation defined over the same group as self.
        - ``row`` -- an integer (default: 0) ; an integer in the range of the degree of right, this number will choose the set o projectors to be chosen in the contruction of the base :: 

        OUTPUT: A change basis matrix that decomposes the equivariant operator.


        EXAMPLES:


        We define the representation by permutation on the cyclic group of four elements and decompose an equivariant operator under this representation  ::
        
            sage: G = FiniteGroup(CyclicPermutationGroup(4))
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: A = matrix.circulant([1,2,3,4])
            sage: rep.is_equivariant_to(A)
            True
            sage: P = G.base_equivariant_to_blocks(rep); P
            [ 1  1  1  1]
            [ 1 -1  I -I]
            [ 1  1 -1 -1]
            [ 1 -1 -I  I]
            sage: A, P.inverse()*A*P
            (
            [1 2 3 4]  [      10        0        0        0]
            [4 1 2 3]  [       0       -2        0        0]
            [3 4 1 2]  [       0        0 -2*I - 2        0]
            [2 3 4 1], [       0        0        0  2*I - 2]
            )



        We define a representation on the permutation group, and decompose an equivariant operator  ::
            
            sage: G = FiniteGroup(['(2,4)(3,7)(6,8)', '(1,3)(4,6)(7,9)'])
            sage: gens = G.gens()
            sage: matrices = [g.matrix() for g in gens]
            sage: rep = representation(gens, matrices)
            sage: operator = matrix([[ 4, -1, -0, -1, -0, -0, -0, -0, -0],
                                    [-1,  4, -1, -0, -1, -0, -0, -0, -0],
                                    [-0, -1,  4, -0, -0, -1, -0, -0, -0],
                                    [-1, -0, -0,  4, -1, -0, -1, -0, -0],
                                    [-0, -1, -0, -1,  4, -1, -0, -1, -0],
                                    [-0, -0, -1, -0, -1,  4, -0, -0, -1],
                                    [-0, -0, -0, -1, -0, -0,  4, -1, -0],
                                    [-0, -0, -0, -0, -1, -0, -1,  4, -1],
                                    [-0, -0, -0, -0, -0, -1, -0, -1,  4]])
            sage: rep.is_equivariant_to(operator)
            True
            sage: P = G.base_equivariant_to_blocks(rep); P
            [ 2  0  0  2  0  2  0  0  0]
            [ 0  2  0  0  2  0  1  0  4]
            [ 2  0  0 -2  0  0  0  8  0]
            [ 0  2  0  0 -2  0  1  0 -4]
            [ 0  0  8  0  0  0  0  0  0]
            [ 0  2  0  0 -2  0 -1  0  4]
            [ 2  0  0 -2  0  0  0 -8  0]
            [ 0  2  0  0  2  0 -1  0 -4]
            [ 2  0  0  2  0 -2  0  0  0]
            sage: operator, P.inverse()*operator*P
            (
            [ 4 -1  0 -1  0  0  0  0  0]  [ 4 -2  0  0  0  0  0  0  0]
            [-1  4 -1  0 -1  0  0  0  0]  [-2  4 -4  0  0  0  0  0  0]
            [ 0 -1  4  0  0 -1  0  0  0]  [ 0 -1  4  0  0  0  0  0  0]
            [-1  0  0  4 -1  0 -1  0  0]  [ 0  0  0  4  0  0  0  0  0]
            [ 0 -1  0 -1  4 -1  0 -1  0]  [ 0  0  0  0  4  0  0  0  0]
            [ 0  0 -1  0 -1  4  0  0 -1]  [ 0  0  0  0  0  4 -1  0  0]
            [ 0  0  0 -1  0  0  4 -1  0]  [ 0  0  0  0  0 -2  4  0  0]
            [ 0  0  0  0 -1  0 -1  4 -1]  [ 0  0  0  0  0  0  0  4 -1]
            [ 0  0  0  0  0 -1  0 -1  4], [ 0  0  0  0  0  0  0 -2  4]
            )

        TO ADD MORE EXAMPLES
        TO ADD MORE EXAMPLES
        TO ADD MORE EXAMPLES
              

        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """
        ###Mudar Nome, talvez symmetry_adapted_basis
        #Nota: right pode inclusive ser a irredutível???? PODE
        n, left = self.irreducible_representations(False)  #Nota: vou pensar em como implementar o __len__
        base = []

        for i in range(n):

            if left(i).inner_product(right) != self.field.zero():  #Nota: acho que aqui é um problema em caso de negativa
                P = self.projection(i, row, 0, right, left)  #Nota: conversar com marcelo a melhor forma de deixar row opcional. #nota: left???
                pivots = P.pivots()
                geracao1 = []
                for pivot in pivots:
                    v0 = P[:, pivot]
                    geracao1.append(v0)
                base = base + geracao1 # Nota: Traduzir
                degree = left(i).degree()

                for k in range(1, degree):  #Nota: Esse laco roda se o range for vazio? Otimizar
                    geracao2 = []
                    for v0 in geracao1:     # Nota: Traduzir
                        v = self.projection(i, row, k, right, left) * v0
                        geracao2.append(v)    # Nota: Traduzir
                    base = base + geracao2

        b = base[0]
        for v in base[1:]:
            b = b.augment(v)  
        return b 

    #####################################################################################
    def quick_block_prevision(self, right, block_prevision=False):##Nota, consertar row ##Revisar a documentacao do True or False
        r"""
        Return a list with order and multiplicities of blocks to an equivariant operator under right defined on self.

        INPUT:

        - ``self`` -- FiniteGroup ; a  Sage permutation group or a group of the class FiniteGroup.
        - ``right`` -- representation ; a representation defined over an arbitrary group given by self.
        - ``block_prevision`` -- a boolean (default: `False`) ; if set to True then prints a string describing the number and orders of blocks to an equivariant operator under right.

        OUTPUT: A list indicating the degree and multiplicity of the representation defined by right. 


        EXAMPLES:

        We define the regular representation(reg) over the symmetric group of 4 simbols and calculate the structure of an equivariant operator under reg ::

            sage: H = SymmetricGroup(4)
            sage: G = FiniteGroup(H);
            sage: reg = G.regular_representation();
            sage: G.quick_block_prevision(reg)
            [['degree', 'multiplicity'], [1, 1], [1, 1], [2, 2], [3, 3], [3, 3]]

            sage: G.quick_block_prevision(reg, block_prevision=True)
            1 block size 1x1
            1 block size 1x1
            2 block size 2x2
            3 block size 3x3
            3 block size 3x3
            [['degree', 'multiplicity'], [1, 1], [1, 1], [2, 2], [3, 3], [3, 3]]

        We define the representation by permutation on the cyclic group and calculate the structure of an equivariant operator ::

            sage: G = FiniteGroup(CyclicPermutationGroup(4))
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: G.quick_block_prevision(rep, block_prevision=True)
            1 block size 1x1
            1 block size 1x1
            1 block size 1x1
            1 block size 1x1
            [['degree', 'multiplicity'], [1, 1], [1, 1], [1, 1], [1, 1]]

        We calculate the irreducible representations of the group of symmetries of tetrahedron and calculate the structure of an equivariant operator inder its irreducibles ::
        
            sage: G = FiniteGroup(AlternatingGroup(4));
            sage: n,irr = G.irreducible_representations(False) ##Irreducible representations
            sage: for j in range(n):
                    print(G.quick_block_prevision(irr(j),block_prevision=True))
            1 block size 1x1
            [['degree', 'multiplicity'], [1, 1]]
            1 block size 1x1
            [['degree', 'multiplicity'], [1, 1]]
            1 block size 1x1
            [['degree', 'multiplicity'], [1, 1]]
            3 block size 1x1
            [['degree', 'multiplicity'], [3, 1]]

            

        We define two representation rep and rep1 on the group of symmetries of regular hexagon and calculate the structure of an equivariant operator under each one ::
            
            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: matrices = [matrix([[1/2,-sqrt(3)/2],[sqrt(3)/2,1/2]]),matrix([[-1,0],[0,1]])];
            sage: rep = representation(generators, matrices)
            sage: G = rep.domain()
            sage: G.quick_block_prevision(rep,block_prevision=True)
            2 block size 1x1
            [['degree', 'multiplicity'], [2, 1]]

            sage: P = matrix([[2,1],[15,2]]); #A change of basis
            sage: matrices1 = [ block_diagonal_matrix((P.inverse())*A*P, (P.inverse())*A*P) for A in matrices];
            sage: rep1 = representation(generators, matrices1)
            sage: G.quick_block_prevision(rep1,block_prevision=True)
            2 block size 2x2
            [['degree', 'multiplicity'], [2, 2]]

               
        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """        
        n, left = self.irreducible_representations(False) 
        
        info = [['degree', 'multiplicity']]        
        for i in range(n):
        	multiplicity = left(i).inner_product(right)
        	try:
        		multiplicity = multiplicity.radical_expression()
        	except: AttributeError
        	if multiplicity != self.field.zero():
        		degree = left(i).degree()
        		info.append([degree, multiplicity])
        if block_prevision:
        	for k in info[1:]:
        		print(str(k[0]) + ' block size ' + str(k[1]) + 'x' + str(k[1]) )
        return info
           	

      
    def base_change_eigenvalue_reduction_new(self, right, block_prevision=False):
        r"""
        Return part of basis change matrix associated to a symmetry adapted basis to an equivariant operator of right.


        INPUT:

        - ``self`` -- FiniteGroup ; a  Sage permutation group or a group of the class FiniteGroup.
        - ``right`` -- representation ; a representation defined over the same group as self.        
        - ``block_prevision`` -- a boolean (default: `False`) ; if set to True then prints a string describing the number and orders of blocks to an equivariant operator under right.
        
        OUTPUT: A list of basis subespaces, along with degree and multiplicitie(number of equivalent copies).
        The restriction of an equivariant operator to each subespace gives a different block matrix of the equivariant operator relative to the symmetry adapted basis.
        In short, these are parts of the symmetry adapted basis that no generates block repetitions.

        EXAMPLES:


        We define the representation by permutation on the cyclic group of four elements and calculates the subspaces that gives one copy of each block of the equivariant operator under this representation  ::
        
            sage: G = FiniteGroup(CyclicPermutationGroup(4))
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: G.base_change_eigenvalue_reduction_new(rep)
            ([[
            [1]        
            [1]        
            [1]        
            [1], [1, 1]
            ],
              [
            [ 1]        
            [-1]        
            [ 1]        
            [-1], [1, 1]
            ],
              [
            [ 1]        
            [ I]        
            [-1]        
            [-I], [1, 1]
            ],
              [
            [ 1]        
            [-I]        
            [-1]        
            [ I], [1, 1]
            ]],
             [['degree', 'multiplicity'], [1, 1], [1, 1], [1, 1], [1, 1]])

            sage: G.base_change_eigenvalue_reduction_new(rep, block_prevision=True)
            1 block size 1x1
            1 block size 1x1
            1 block size 1x1
            1 block size 1x1
            ([[
            [1]        
            [1]        
            [1]        
            [1], [1, 1]
            ],
              [
            [ 1]        
            [-1]        
            [ 1]        
            [-1], [1, 1]
            ],
              [
            [ 1]        
            [ I]        
            [-1]        
            [-I], [1, 1]
            ],
              [
            [ 1]        
            [-I]        
            [-1]        
            [ I], [1, 1]
            ]],
             [['degree', 'multiplicity'], [1, 1], [1, 1], [1, 1], [1, 1]])



        We define a representation on the permutation group, and calculates parts of the symmetry adapted basis that generates no repetition in the blocks of an equivariant operator  ::

            sage: G = FiniteGroup(['(2,4)(3,7)(6,8)', '(1,3)(4,6)(7,9)'])
            sage: gens = G.gens()
            sage: matrices = [g.matrix() for g in gens]
            sage: rep = representation(gens, matrices)
            sage: G.base_change_eigenvalue_reduction_new(rep)
            ([[
            [2 0 0]        
            [0 2 0]        
            [2 0 0]        
            [0 2 0]        
            [0 0 8]        
            [0 2 0]        
            [2 0 0]        
            [0 2 0]        
            [2 0 0], [1, 3]
            ],
              [
            [ 2]        
            [ 0]        
            [-2]        
            [ 0]        
            [ 0]        
            [ 0]        
            [-2]        
            [ 0]        
            [ 2], [1, 1]
            ],
              [
            [ 0]        
            [ 2]        
            [ 0]        
            [-2]        
            [ 0]        
            [-2]        
            [ 0]        
            [ 2]        
            [ 0], [1, 1]
            ],
              [
            [ 2  0]        
            [ 0  1]        
            [ 0  0]        
            [ 0  1]        
            [ 0  0]        
            [ 0 -1]        
            [ 0  0]        
            [ 0 -1]        
            [-2  0], [2, 2]
            ]],
             [['degree', 'multiplicity'], [1, 3], [1, 1], [1, 1], [2, 2]])
            
        We define the regular representation(reg) over the symmetric group of 3 simbols and calculate the subspaces that gives all the blocks of an equivariant operator without repetition ::

            sage: H = SymmetricGroup(3)
            sage: G = FiniteGroup(H);
            sage: reg = G.regular_representation();
            sage: G.base_change_eigenvalue_reduction_new(reg,block_prevision=True)
            1 block size 1x1
            1 block size 1x1
            2 block size 2x2
            ([[
            [1]        
            [1]        
            [1]        
            [1]        
            [1]        
            [1], [1, 1]
            ],
              [
            [ 1]        
            [ 1]        
            [ 1]        
            [-1]        
            [-1]        
            [-1], [1, 1]
            ],
              [
            [                                         1                                          0]
            [-0.500000000000000? + 0.866025403784439?*I                                          0]
            [-0.500000000000000? - 0.866025403784439?*I                                          0]
            [                                         0                                          1]
            [                                         0 -0.500000000000000? - 0.866025403784439?*I]
            [                                         0 -0.500000000000000? + 0.866025403784439?*I],

            [2, 2]
            ]],
             [['degree', 'multiplicity'], [1, 1], [1, 1], [2, 2]])
           
              

        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """       
        row = 0
        n, left = self.irreducible_representations(False)  #Nota: vou pensar em como implementar o __len__
        base = []
        #base_new = []
        info = [['degree', 'multiplicity']]

        for i in range(n):
            multiplicity = left(i).inner_product(right)
            if multiplicity != self.field.zero():  

                P = self.projection(i, row, 0, right, left)  #Nota: left???
                pivots = P.pivots()
                degree = left(i).degree()                
                v0 = P[:, pivots]
                base.append([v0, [degree, multiplicity]])
                info.append([degree, multiplicity])
        if block_prevision:
            for k in info[1:]:
                print(str(k[0]) + ' block size ' + str(k[1]) + 'x' + str(k[1]) )
        return base, info

group = FiniteGroup

######################################## Morphism ###############################################
   
class IsotypicBase(object):
    """docstring for ClassName"""
    def __init__(self, base):
        #super(IsotypicBase, self).__init__()
        #self.base = list(filter(lambda x: x!=None, base))

        self.isotypic_components = base#[b[0] for b in  self._base]
       

    # def __repr__(self):
    #     msg = []
    #     for d, m in self._info:
    #            msg.append("{} blocks of size  {} x {}\n".format(d, m, m))
            
        
    #     return ''.join(tuple(msg))

    def __len__(self):
        return len(self.base)
    
    def get_blocks(self,  matrix_equiv):

        r"""
    Return the blocks corresponding to each isotypic component using the provided equivariant matrix.

    INPUT:

    - ``matrix_equiv`` -- matrix ; an equivariant operator.

    OUTPUT: A list of blocks, where each block is obtained by applying `get_block` to an isotypic component and the given matrix equivariant.

    EXAMPLES:"""    
        
                
        blocks = pmap(lambda b: get_block(b, matrix_equiv),  self.isotypic_components)
        return blocks
    
    def list(self):
        return  self.isotypic_components

    def matrix(self):
        columns= self.isotypic_components
        c = columns[0]
        for columm in columns[1:]:
            c = c.augment(columm)

        return c




def get_block(columm_base, matrix_equiv): 
    r"""
    Return the matricial representation of matrix_equiv when restricted to the subspace whose basis is columm_base.


    INPUT:

    - ``columm_base`` -- matrix ; a matrix(with full rank) whose columns are part of symmetry adapted basis.
    - ``matrix_equiv`` -- matrix ; a matrix an equivariant operator.

    OUTPUT: A matrix representing the restriction of matrix_equiv to the subespaces given by column_base.

    EXAMPLES:


    We define the representation by permutation on the cyclic group of four elements and calculates the blocks of an equivariant operator under this representation  ::

        sage: G = FiniteGroup(CyclicPermutationGroup(4))
        sage: generators = G.gens()
        sage: matrices = [g.matrix() for g in generators]
        sage: rep = representation(generators, matrices)
        sage: subspaces = G.base_change_eigenvalue_reduction_new(rep) #The Invariant Subspaces
        sage: A = matrix.circulant([1,2,3,4]); #An equivariant Operator 
        sage: rep.is_equivariant_to(A)
        True
        sage: [get_block(k[0],A) for k in subspaces[0]] # The block decomposition of the equivariant operator
        [[10], [-2], [-2 - 2*I], [-2 + 2*I]]


    We define a representation on a permutation group, and calculates the blocks(no repetition) of an equivariant operator  ::

        sage: G = FiniteGroup(['(2,4)(3,7)(6,8)', '(1,3)(4,6)(7,9)'])
        sage: gens = G.gens()
        sage: matrices = [g.matrix() for g in gens]
        sage: rep = representation(gens, matrices)
        sage: operator = matrix([[ 4, -1, -0, -1, -0, -0, -0, -0, -0],
                                [-1,  4, -1, -0, -1, -0, -0, -0, -0],
                                [-0, -1,  4, -0, -0, -1, -0, -0, -0],
                                [-1, -0, -0,  4, -1, -0, -1, -0, -0],
                                [-0, -1, -0, -1,  4, -1, -0, -1, -0],
                                [-0, -0, -1, -0, -1,  4, -0, -0, -1],
                                [-0, -0, -0, -1, -0, -0,  4, -1, -0],
                                [-0, -0, -0, -0, -1, -0, -1,  4, -1],
                                [-0, -0, -0, -0, -0, -1, -0, -1,  4]])  #An equivariant Operator 
        sage: rep.is_equivariant_to(operator)
        True
        sage: subspaces = G.base_change_eigenvalue_reduction_new(rep) #The Invariant Subspaces
        sage: [get_block(k[0],operator) for k in subspaces[0]] # The block decomposition of the equivariant operator
        [
        [ 4 -2  0]                   
        [-2  4 -4]            [ 4 -1]
        [ 0 -1  4], [4], [4], [-2  4]
        ]
        
        #NOTE: It is worth to note, in this case, that the two equal blocks [4], comes from inequivalent subrepresentations.
        In this case, the theory has not way to foresee that.

    We define the regular representation(reg) over the symmetric group of 3 simbols and calculate the blocks of an equivariant operator without repetition ::

        sage: H = SymmetricGroup(3)
        sage: G = FiniteGroup(H);
        sage: reg = G.regular_representation();
        sage: Id = matrix.identity(reg.degree()); # Identity matrix
        sage: reg.is_equivariant_to(Id)
        True
        sage: subspaces = G.base_change_eigenvalue_reduction_new(reg) #The Invariant Subspaces
        sage: [get_block(k[0],Id) for k in subspaces[0]] # The block decomposition of the equivariant operator
        [
        [1], [1],

        [1.000000000000000? + 0.?e-18*I                              0]
        [                             0 1.000000000000000? + 0.?e-18*I]
        ]
      

    REFERENCES:

    For more information, see the following references:

    - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

    - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     

    """
    columm_pinv = columm_base.pseudoinverse()
    block = columm_pinv*matrix_equiv*columm_base
    return block


class MapRepresentation(SetMorphism):


    def __init__(self, H, func):
        self._domain = H.domain()
        self._codomain = H.codomain()
        self.func = func
        

        super().__init__(H, func)

    def _repr_(self):

        return "Map: \n From: {} \n To: {}.".format(self._domain,self._codomain)

    @cached_method
    def inner_product(self, right):
        r"""
        Return the inner product between self and right.

        INPUT:

        - ``self`` -- representation ; a representation defined over an arbitrary group G.
        - ``right`` -- representation ; a representation defined over the same group as self.

        OUTPUT: The inner product between the two representations defined on the group G. 


        EXAMPLES:

        We define the representation by permutation on the cyclic group of 4 elements and calculate its inner products with the trivial representation ::

            sage: G = FiniteGroup(CyclicPermutationGroup(4))
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: n,irr = G.irreducible_representations(False) #Irreducible representations
            sage: rep.inner_product(irr(0))
            1

        We define one representation on the group of symmetries of regular hexagon and calculate the inner product with itself ::

            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: matrices = [matrix([[1/2,-sqrt(3)/2],[sqrt(3)/2,1/2]]),matrix([[-1,0],[0,1]])];
            sage: rep = representation(generators, matrices)
            sage: rep.inner_product(rep)
            1

        We define the regular representation over the symmetric group of 4 simbols and calculate the inner products between the regular representation and the irreducible representations ::

            sage: H = SymmetricGroup(4)
            sage: G = FiniteGroup(H);
            sage: n,irr = G.irreducible_representations(False) #Irreducible representations
            sage: reg = G.regular_representation();
            sage: [reg.inner_product(irr(j)) for j in range(n)] # The inner product between the regular and irreducibles
            [1, 1, 2, 3, 3]
                
        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.        
        """

        
        #NOTA: Verificar se right eh uma representation() (class Representation)
        group = self._domain
        s = 0
        
        # O método conjugacy_classes() retorna uma lista de classes do grupo.
        for conj_class in group.conjugacy_classes():
            # Pega um único representante para a classe.
            g = conj_class.representative()
            
            # O tamanho da classe é o peso do termo na soma.
            class_size = len(conj_class)
            
            term = class_size * self(g.inverse()).character() * right(g).character()
            s += term
            
        if group.field == SR:
            s = s.simplify_full()
            
        return (1 / group.order()) * s
    

    def tensor_product(self, right):
        if not right.is_representation(): # Nota: Em teste, pois pode ser lento em alguns casos.
            msg = '{} is not a representation.'
            raise TypeError(''.join(mgs).format(right))
        left = self
        domain = left._domain
        if not domain==right._domain:
            msg = 'The domains {} and {} is not the same.'#Nota: Melhorar mensagens?
            raise TypeError (''.join(msg).format(domain, right._domain))
        field = left._domain.field
        gens = left._domain.gens()
        image = [(left(g).matrix()).tensor_product(right(g).matrix()) for g in gens]
        # Nota: Vamos pensar melhor como o parametro field é util.
        return representation(gens, image, field)
        

    def direct_sum(self, *reps):
        """
        INPUT: A list of representations.
        OUTPUT: Their direct sum
        """
        reps = list(reps)
        domain = self._domain
        reps.insert(0, self)    
        for k in reps:
            if not k._domain == domain:
                msg = 'The domains {} and {} are not the same.'
                raise TypeError(msg.format(domain, k._domain))            
        field = domain.field
        gens = domain.gens()    
        image = [block_diagonal_matrix([rep(g).matrix() for rep in reps]) for g in gens]
        return representation(gens, image, field)

    def an_element(self):
        r"""
        Return an image of an element of the group under self.

        INPUT:

        - ``self`` -- representation ; a representation defined over an arbitrary group G.

        OUTPUT: A linear transformation corresponding to the image under self of a random element of the group G. 


        EXAMPLES:

        We define the cyclic group of 4 elements and choose an element to show its image under self ::

            sage: G = CyclicPermutationGroup(4) 
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: rep.an_element()
            Linear transformation associated with element g=(1,2,3,4), represented by the non-singular matrix:
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            [1 0 0 0]
            Representation space: Vector space of dimension 4 over Algebraic Field


        We define one representation on the group of symmetries of regular hexagon and show the image of an element of the group ::

            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: matrices = [matrix([[1/2,-sqrt(3)/2],[sqrt(3)/2,1/2]]),matrix([[-1,0],[0,1]])];
            sage: rep = representation(generators, matrices)
            sage: rep.an_element()
            Linear transformation associated with element g=(1,3)(4,6), represented by the non-singular matrix:
            [               -1/2 -0.866025403784439?]
            [-0.866025403784439?                 1/2]
            Representation space: Vector space of dimension 2 over Algebraic Field

        We define the regular representation over the symmetric group of 4 simbols and calculate the image of a random element ::

            sage: H = SymmetricGroup(4)
            sage: G = FiniteGroup(H);
            sage: reg = G.regular_representation();
            sage: reg.an_element()
            Linear transformation associated with element g=(1,3,4), represented by the non-singular matrix:
            24 x 24 dense matrix over Algebraic Field
            Representation space: Vector space of dimension 24 over Algebraic Field


                
        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """                
        g = self._domain.an_element()
        return self.func(g)

    def degree(self):
        r"""
        Return the degree of a representation defined over an arbitrary group G.

        INPUT:

        - ``self`` -- representation; a representation defined over an arbitrary group G.

        OUTPUT: The degree of self, i.e. , the dimension of the subjacent vector space of self. 


        EXAMPLES:

        We define the representation by permutation on the cyclic group of 4 elements ::

            sage: G = CyclicPermutationGroup(4) 
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: rep.degree()
            4

        We define one representation on the group of symmetries of regular hexagon ::

            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: matrices1 = [matrix([[1/2,-sqrt(3)/2],[sqrt(3)/2,1/2]]),matrix([[-1,0],[0,1]])];
            sage: rep = representation(generators, matrices1)
            sage: rep.degree()
            2

        We define the regular representation over the symmetric group of 4 simbols and calculate its degree ::

            sage: H = SymmetricGroup(4)
            sage: G = FiniteGroup(H);
            sage: reg = G.regular_representation();
            sage: reg.degree()
            24

                
        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """        
        d = self.an_element()
        return d.matrix().nrows()
    
    ##Nota: enriquecer documentacao Ex produto tensorial ####Discutir sobre um teste de dominio.
    def is_equivariant_to(self, M):  #Nota: a Testes 
        r"""
        Tests if an operator is equivariant to a representation defined over a an arbitrary group G.

        INPUT:

        - ``M`` -- matrix ; a square matrix associated to an operator over a representation vector space.

        - ``self`` -- representation; a representation defined over an arbitrary group G.

        OUTPUT: True if the operator is equivariant or False if is not equivariant.


        EXAMPLES:

        Define the representation by permutation on the cyclic group of 4 elements ::

            sage: G = CyclicPermutationGroup(4) 
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: v = [1,2,3,4]
            sage: A = matrix.circulant(v); 
            sage: rep.is_equivariant_to(A)
            True
            sage: B = matrix.circulant(SR.var('a b c d'))
            sage: rep.is_equivariant_to(B)
            True
            
        We define a representation over the symmetric group and test if the identity is equivariant ::
        
            sage: G = SymmetricGroup(5)
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices, field=AA)
            sage: rep.is_equivariant_to(identity_matrix(AA, 5))
            True
            sage: C = matrix.hilbert(5).change_ring(AA);
            sage: rep.is_equivariant_to(C);
            False
            
        ENRIQUECER QUANDO DEFINIR PRODUTO TENSORIAL OU SOMA DIRETA.

        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012. 

            """

        for g in self._domain.gens():
            test = self(g).matrix() * M == M * self(g).matrix()
            if not test:
                return False
        return True

    def is_unitary(self):  #Nota: falta testes. Será que basta nos geradores?
        r"""
         Tests if a representation defined over an arbitrary group G is unitary.

         INPUT:

         - ``self`` -- representation; a representation defined over an arbitrary group G.

         OUTPUT: True if the representation is unitary or False if is not unitary.


         EXAMPLES:

         Define the representation by permutation on the cyclic group of 4 elements ::

             sage: G = CyclicPermutationGroup(4) 
             sage: generators = G.gens()
             sage: matrices = [g.matrix() for g in generators]
             sage: rep = representation(generators, matrices)
             sage: rep.is_unitary()
             True
             
         Define two representations on the group of symmetries of regular hexagon ::
         
             sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
             sage: P = matrix([[2,1],[15,2]]); ## A non-unitary base change
             sage: matrices1 = [matrix([[1/2,-sqrt(3)/2],[sqrt(3)/2,1/2]]),matrix([[-1,0],[0,1]])];
             sage: matrices2 = [P.inverse()*A*P for A in matrices1]
             sage: rep1 = representation(generators, matrices1)
             sage: rep2 = representation(generators, matrices2)
             sage: rep1.is_unitary(), rep2.is_unitary()
             (True, False)
 
                
         REFERENCES:

         For more information, see the following references:

         - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

         - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """

        for g in self._domain:

            if self(g).matrix().is_unitary() == False:

                return False
        return True

    def is_irreducible(self):##Nota, Discutir Reducibilidade associada a um corpo.
        r"""
        Tests if a representation defined over an arbitrary group G is irreducible.

        INPUT:

        - ``self`` -- representation; a representation defined over an arbitrary group G.

        OUTPUT: True if the representation is irreducible  or False if is reducible.


        EXAMPLES:

        Define the representation by permutation on the cyclic group of 6 elements ::

            sage: G = CyclicPermutationGroup(6) 
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices)
            sage: rep.is_irreducible()
            False
            
        Define one representation on the group of symmetries of regular hexagon ::
        
            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: matrices = [matrix([[1/2,-sqrt(3)/2],[sqrt(3)/2,1/2]]),matrix([[-1,0],[0,1]])];
            sage: rep = representation(generators, matrices)
            sage: rep.is_irreducible()
            True
            
        We define the regular representation over the symmetric group of 4 simbols and ask if is irreducible ::
        
            sage: H = SymmetricGroup(4)
            sage: G = FiniteGroup(H);
            sage: reg = G.regular_representation();
            sage: reg.is_irreducible()
            False          
            
        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """
        if bool(self.inner_product(self) == self._domain.field.one()):
            return True
        return False


    def is_representation(self):
        """INPUT- A representation.
	    OUTPUT- A boolean True or False acording rep being a representation"""
        domain = self._domain
        for g in domain:
            for h in domain:
                if not self(g*h).matrix()==(self(g).matrix())*(self(h)).matrix(): #Nota: ainda esta lento
                    return False
        return True
        
    
    def is_stable(self, subspace):
        domain = self._domain
        generators = domain.gens()
        space =  self(domain.an_element())._domain
        basis = subspace.basis() 
        if not subspace.is_subspace(space):
            msg = '{} is not a subspace of {}' 
            raise TypeError(''.join(msg).format(subspace, space))
        for g in generators:
            for vector in basis:
                if not self(g).matrix()*vector in subspace:
                    return False
        return True




class Dg_Linear_Transformation(VectorSpaceMorphism):
    
    def __init__(self, H, A, g):

        self._domain = H.domain()
        self._codomain = H.codomain()
        self.g = g
        self.group = g.parent()  # NOTA: talvez nao seja necessario

        super().__init__(H, A)

    def _repr_(self):

        m = self.matrix()
        msg = (
            "Linear transformation associated with element g={}, represented by the non-singular matrix:\n",
            "{!r}\n", "Representation space: {}\n")
        return ''.join(msg).format(self.g, m, self.domain())

    def character(self):
        r"""
        Return the character of a representation defined over an arbitrary group G.

        INPUT:

        - ``self`` -- representation; a representation defined over an arbitrary group G.

        OUTPUT: The character of self. 


        EXAMPLES:

        We define the representation by permutation on the cyclic group of 4 elements and calculate its characters ::

            sage: G = FiniteGroup(CyclicPermutationGroup(4)) 
            sage: generators = G.gens()
            sage: matrices = [g.matrix() for g in generators]
            sage: rep = representation(generators, matrices); rep
            sage: g = G.an_element()
            sage: rep(g).character()
            0
            sage: [rep(g).character() for g in G]
            [4, 0, 0, 0]

        We define one representation on the group of symmetries of regular hexagon and calculate the character of each element in the group ::

            sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
            sage: matrices = [matrix([[1/2,-sqrt(3)/2],[sqrt(3)/2,1/2]]),matrix([[-1,0],[0,1]])];
            sage: rep = representation(generators, matrices)
            sage: G = rep.domain()
            sage: [ [g,rep(g).character()] for g in G]
            [[(), 2],
             [(1,5,3)(2,6,4), -1],
             [(1,3,5)(2,4,6), -1],
             [(1,6,5,4,3,2), 1],
             [(1,4)(2,5)(3,6), -2],
             [(1,2,3,4,5,6), 1],
             [(2,6)(3,5), 0],
             [(1,5)(2,4), 0],
             [(1,3)(4,6), 0],
             [(1,6)(2,5)(3,4), 0],
             [(1,4)(2,3)(5,6), 0],
             [(1,2)(3,6)(4,5), 0]]


        We define the regular representation over the symmetric group of 4 simbols and calculate its characters ::

            sage: H = SymmetricGroup(4)
            sage: G = FiniteGroup(H);
            sage: reg = G.regular_representation();
            sage: [ [g,reg(g).character()] for g in G]
            [[(), 24],
             [(1,4)(2,3), 0],
             [(1,2)(3,4), 0],
             [(1,3)(2,4), 0],
             [(2,4,3), 0],
             [(1,4,2), 0],
             [(1,2,3), 0],
             [(1,3,4), 0],
             [(2,3,4), 0],
             [(1,4,3), 0],
             [(1,2,4), 0],
             [(1,3,2), 0],
             [(3,4), 0],
             [(1,4,2,3), 0],
             [(1,2), 0],
             [(1,3,2,4), 0],
             [(2,4), 0],
             [(1,4,3,2), 0],
             [(1,2,3,4), 0],
             [(1,3), 0],
             [(2,3), 0],
             [(1,4), 0],
             [(1,2,4,3), 0],
             [(1,3,4,2), 0]]


                
        REFERENCES:

        For more information, see the following references:

        - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

        - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
        
        """                
        return self.trace()


def check(s, index, i):
    import copy
    ss = copy.copy(s)
    ss = s +'s'
    sss = ss[index:index+i]
    k = [x for x in sss][-1]    
    try:
        if isinstance(sage_eval(k), sage.rings.integer.Integer):
            return True    
    except: SyntaxError   
    return False
    

def avoid_inverse(g, gens, s):
    import re
    import copy
    def locals():
        local = dict()
        for i, g_i in enumerate(gens):
            local['x'+ str(i+1)] = g_i
        return local
    if len(s)>2:
        ss = copy.copy(s)
        local= locals()        
        keys = local.keys()
        for x_i in keys:
            indexes = [m.start() for m in re.finditer(x_i, ss)]
            g = local[x_i]             
            if isinstance(g, sage.interfaces.gap.GapElement):
                order = gap.Order(g)
            else:                
                order = g.order()            
            for index in indexes:
                i = 5
                test = True
                while test:                    
                    i+=1
                    #print (i)
                    test = check(ss, index, i)
                    if not test:
                        i-=1
                    #print (test)
                s0 = ss[index:index+i]
                #print ('s0',s0)
                if '-' in s0:
                    p = s0[3:]
                    exp_nonnegative = str(order-abs(sage_eval(p)))
                    #print('exp_nonnegative', exp_nonnegative)
                    exp_negative = p
                    #print('exp_negative',exp_negative)
                    s1 = s0.replace(exp_negative, exp_nonnegative)
                    #print('s1',s1)
                    s = s.replace(s0, s1)
    return s

def factorization(g, words, matrices):    
        
    g = words[0].parent()(g)
    H = libgap.Group(words)
    ans = H.EpimorphismFromFreeGroup().PreImagesRepresentative(g)
    l1 = str(ans)
    #l1 = avoid_inverse(g, words, l1)    
    l2 = l1.replace('^', '**')
    
    def locals2():
        local = dict()
        for i, w_i in enumerate(matrices):
            local['x'+ str(i+1)] = w_i        
        return local    
    return sage_eval(l2, locals=locals2())


############################## Reduclible  Representations ################
def representation(generators, matrices, field=QQbar):
    r"""
    Create a representation of an arbitrary group G.

    INPUT:

    - ``generators`` -- list; generators of the group G.

    - ``matrices`` -- list; generators images in the same sequence of generators.

    - ``field`` -- Sage ring (default: `QQbar`); Field of the subjacent representation vector space. Faster computations using the fields `QQbar` or `AA`.

    OUTPUT: Mapping of the group generated by ``generators`` to the group of matrices generated by ``matrices``.


    EXAMPLES:

    This example illustrates the representation by permutation defined on the group of symmetries of the equilateral triangle ::

        sage: G = DihedralGroup(3) 
        sage: generators = G.gens()
        sage: matrices = [g.matrix() for g in generators]
        sage: rep = representation(generators, matrices)
        sage: rep
        Map: 
         From: Permutation FiniteGroup with generators [(1,2,3), (1,3)] 
         To: Matrix group over Algebraic Field with 2 generators (
        [0 1 0]  [0 0 1]
        [0 0 1]  [0 1 0]
        [1 0 0], [1 0 0]
        ).

    Matrix associated to an element g of G ::
    
        sage: g = G.an_element()
        sage: rep(g).matrix()
        [0 1 0]
        [1 0 0]
        [0 0 1]        

    We define the generators of the group of symmetries of an regular hexagon, and the correspondent matrices of symmetry (a rotation through an angle $\pi/3$,
    and a reflexion across the $x$ axis in the cartesian plane) ::

        sage: generators = ["(1,2,3,4,5,6)","(1,4)(2,3)(5,6)"]
        sage: matrices = [matrix([[1/2,-sqrt(3)/2],[sqrt(3)/2,1/2]]),matrix([[-1,0],[0,1]])]
        sage: rep = representation(generators, matrices)
        sage: rep("(1,2,3,4,5,6)")
        Linear transformation associated with element g=(1,2,3,4,5,6), represented by the non-singular matrix:
        [                1/2 -0.866025403784439?]
        [ 0.866025403784439?                 1/2]
        Representation space: Vector space of dimension 2 over Algebraic Field
        
        sage: rep("(1,4)(2,3)(5,6)")
        Linear transformation associated with element g=(1,4)(2,3)(5,6), represented by the non-singular matrix:
        [-1  0]
        [ 0  1]
        Representation space: Vector space of dimension 2 over Algebraic Field

        sage: (rep.domain()).is_isomorphic(DihedralGroup(6))
        True
    
    Next we define the same representation where the field is choosed to be SR ::

        sage: rep = representation(generators, matrices, field=SR)
        sage: rep("(1,2,3,4,5,6)")
        Linear transformation associated with element g=(1,2,3,4,5,6), represented by the non-singular matrix:
        [         1/2 -1/2*sqrt(3)]
        [ 1/2*sqrt(3)          1/2]
        Representation space: Vector space of dimension 2 over Symbolic Ring

    We choose a random element of the group and show its matrix ::

        sage: G = FiniteGroup(generators)
        sage: g = G.an_element()
        sage: rep(g)
        Linear transformation associated with element g=(1,3)(4,6), represented by the non-singular matrix:
        [        -1/2 -1/2*sqrt(3)]
        [-1/2*sqrt(3)          1/2]
        Representation space: Vector space of dimension 2 over Symbolic Ring

    The identity and inverse properties of the representation are illustred ::

        sage: e = G.identity()
        sage: rep(e)
        Linear transformation associated with element g=(), represented by the non-singular matrix:
        [1 0]
        [0 1]
        Representation space: Vector space of dimension 2 over Symbolic Ring
        sage: [rep(g.inverse())==rep(g).inverse() for g in G]
        [True, True, True, True, True, True, True, True, True, True, True, True]

    
    It is an error to choose the group element as a tuple ::

        sage: rep((1,4)(2,3)(5,6))
        TypeError: Traceback (most recent call last):
        ...
        TypeError: 'tuple' object is not callable

        sage: rep((1,2,3,4,5,6))
        TypeError: Traceback (most recent call last):
        ...
        TypeError:  'tuple' object is not callable
    
    The argument can be a string or a list. Be aware that '(1,2,3,4,5,6)' is the usual
    notation of cycle, while [1,2,3,4,5,6] lists the images of the permutation ::

        sage: rep("(1,2,3,4,5,6)")
        Linear transformation associated with element g=(1,2,3,4,5,6), represented by the non-singular matrix:
        [         1/2 -1/2*sqrt(3)]
        [ 1/2*sqrt(3)          1/2]
        Representation space: Vector space of dimension 2 over Symbolic Ring

        sage: rep([1,2,3,4,5,6])
        Linear transformation associated with element g=(), represented by the non-singular matrix:
        [1 0]
        [0 1]
        Representation space: Vector space of dimension 2 over Symbolic Ring

        sage: rep([2,3,4,5,6,1])
        Linear transformation associated with element g=(1,2,3,4,5,6), represented by the non-singular matrix:
        [         1/2 -1/2*sqrt(3)]
        [ 1/2*sqrt(3)          1/2]
        Representation space: Vector space of dimension 2 over Symbolic Ring       

        
    REFERENCES:

    For more information, see the following references:

    - [Ser1977]_Serre, Jean-Pierre. Linear representations of finite groups. Vol. 42. New York: springer, 1977.

    - [Sti2012]_Stiefel, E., and A. Fässler. FiniteGroup theoretical methods and their applications. Springer Science & Business Media, 2012.     
    
    """

    G = FiniteGroup(generators, field)
    identity = G.identity()
    #matrices = [m.change_ring(field) for m in matrices] # NOTA: deixa o codigo lento
    M = MatrixGroup(matrices)
    
    if isinstance(generators[0], str):
        generators = [G(libgap(g)) for g in generators]
        
    n = (matrices[0]).nrows()
    H = Hom(field**n, field**n)

    @cached_function
    def rep(g):
        # Nota: Garantir que g pertence a classe FiniteGroup
        if g == identity:
            return Dg_Linear_Transformation(H, matrix.identity(n), g)
        else:
            A = factorization(g, generators, matrices)

        return Dg_Linear_Transformation(H, A, g)

    return MapRepresentation(Hom(G, M), rep)


