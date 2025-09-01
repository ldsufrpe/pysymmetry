# test.py
import sys
import os
import unittest
import pytest


# Add the parent directory to the path to allow package imports
# This is necessary to run the script directly with 'sage'
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# --- Local Package Import ---
from pysymmetry import *

# --- SageMath Imports (Consolidated) ---
from sage.all import (
    PermutationGroup,
    SymmetricGroup,
    MatrixGroup,
    CyclicPermutationGroup,
    DihedralGroup,
    matrix,
    sqrt,
    QQbar,
    SR,
    ZZ
)
# Direct import to fix the isinstance() TypeError with lazy imports


# --- NumPy Import ---
from numpy import array, array_equal


class TestGroup(unittest.TestCase):
    """
    Consolidated test suite for the pysymmetry package.
    Includes original tests and tests generated from docstrings.
    """

    # --- Original Tests ---

    def test_input(self):
        """Tests the initialization of FiniteGroup with different inputs."""
        # Gap input
        gens = ["(1,2,3,4)", "(1,4)(2,3)", "(8,9)"]
        G = FiniteGroup(gens)
        elements = G.list()
        expected = PermutationGroup(gens).list()
        self.assertEqual(set(elements), set(expected))

        # Matrix input
        mat = [matrix([[0, -1], [1, 0]]), matrix([[0, 1], [1, 0]])]
        G = FiniteGroup(mat, matrix=True)
        elements = G.list()
        expected = MatrixGroup(mat).as_permutation_group().list()
        self.assertEqual(set(elements), set(expected))

    def test_regular(self):
        """Tests the regular representation for a simple group."""
        gens = ["(1,2)", "(3,4)"]
        G = FiniteGroup(gens)
        F = G.field
        elements = G.list()
        g1 = elements[1]
        g2 = elements[2]
        regular = G.regular_representation()
        
        expected1 = matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1],
                            [0, 0, 1, 0]]).change_ring(F)
        expected2 = matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0],
                            [0, 1, 0, 0]]).change_ring(F)

        self.assertEqual(regular(g1).matrix(), expected1)
        self.assertEqual(regular(g2).matrix(), expected2)

    # --- Tests from Docstrings ---

    def test_regular_representation_cyclic(self):
        """Test from the docstring of: regular_representation with CyclicPermutationGroup(4)."""
        G = FiniteGroup(CyclicPermutationGroup(4))
        reg = G.regular_representation()
        
        g1 = G("(1,2,3,4)")
        expected1 = matrix(QQbar, [[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
        self.assertEqual(reg(g1).matrix(), expected1)

        g2 = G("(1,3)(2,4)")
        expected2 = matrix(QQbar, [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
        self.assertEqual(reg(g2).matrix(), expected2)

    def test_irreducible_representations_s4(self):
        """Test from the docstring of: irreducible_representations with SymmetricGroup(4)."""
        G = FiniteGroup(SymmetricGroup(4))
        n, irr = G.irreducible_representations(show_table=False)
        self.assertEqual(n, 5)
        
        irr4 = irr(4)
        g1 = G("(1,2)")
        expected_gen1_img = matrix([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
        self.assertEqual(irr4(g1).matrix(), expected_gen1_img)

    def test_isotypic_projection_s3(self):
        """Test from the docstring of: isotypic_projection with SymmetricGroup(3)."""
        G = FiniteGroup(SymmetricGroup(3))
        reg = G.regular_representation()
        projs = G.isotypic_projection(reg)
        
        # Check proportionality for the first projector
        proj1_calculated = projs[0]
        self.assertNotEqual(proj1_calculated[0,0], 0)
        scalar1 = proj1_calculated[0,0] 
        expected_proj1 = matrix(QQbar, 6, [scalar1]*36)
        self.assertEqual(proj1_calculated, expected_proj1)

        # Check proportionality for the third projector
        proj3_calculated = projs[2]
        self.assertNotEqual(proj3_calculated[0,0], 0)
        scalar3 = proj3_calculated[0,0]
        expected_proj3_base = matrix(QQbar, [
             [ 2, -1, -1,  0,  0,  0], [-1,  2, -1,  0,  0,  0], [-1, -1,  2,  0,  0,  0],
             [ 0,  0,  0,  2, -1, -1], [ 0,  0,  0, -1,  2, -1], [ 0,  0,  0, -1, -1,  2]
        ])
        self.assertEqual(proj3_calculated, (scalar3/2) * expected_proj3_base)

    def test_isotypic_base_s3_regular(self):
        """Test from the docstring of: isotypic_base with S3 and its regular representation."""
        G = FiniteGroup(SymmetricGroup(3))
        reg = G.regular_representation()
        change_basis = G.isotypic_base(reg, isotypic_components=False)
        
        g = G("(1,2,3)")
        Dg = reg(g).matrix()
        A_block = change_basis.inverse() * Dg * change_basis

        self.assertEqual(A_block[0, 2], 0)
        self.assertEqual(A_block[1, 2], 0)
        self.assertEqual(A_block[2, 0], 0)
        self.assertEqual(A_block[2, 1], 0)

    def test_representation_d3(self):
        """
        Test from the docstring of: representation with DihedralGroup(3).
        
        """
        G = DihedralGroup(3)
        generators = G.gens()
        matrices = [g.matrix() for g in generators]
        rep = representation(generators, matrices)

        # Check the type of the codomain by its name, avoiding isinstance() issues
        self.assertIn('MatrixGroup', str(type(rep.codomain())))
        
        g = G.an_element()
        g_img_matrix = rep(g).matrix()
        g_perm_matrix = g.matrix()
        self.assertEqual(g_img_matrix, g_perm_matrix)
        
    def test_representation_hexagon_symbolic(self):
        """Test from the docstring of: representation for the hexagon's symmetry group."""
        generators = ["(1,2,3,4,5,6)", "(1,4)(2,3)(5,6)"]
        matrices = [
            matrix(SR, [[1/2, -sqrt(3)/2], [sqrt(3)/2, 1/2]]),
            matrix(SR, [[-1, 0], [0, 1]])
        ]
        rep = representation(generators, matrices, field=SR)
        
        g1 = rep.domain()("(1,2,3,4,5,6)")
        expected1 = matrix(SR, [[1/2, -1/2*sqrt(3)], [1/2*sqrt(3), 1/2]])
        self.assertEqual(rep(g1).matrix(), expected1)

        g2 = rep.domain()("(1,4)(2,3)(5,6)")
        expected2 = matrix(SR, [[-1, 0], [0, 1]])
        self.assertEqual(rep(g2).matrix(), expected2)
        
        self.assertTrue(rep.domain().is_isomorphic(DihedralGroup(6)))

    def test_base_change_matrix_cyclic4(self):
        """
        Test from the docstring of: base_change_matrix with CyclicPermutationGroup(4).
        This test verifies that the change of basis matrix correctly diagonalizes the representation.
        """
        G = FiniteGroup(CyclicPermutationGroup(4))
        rep = G.natural_representation()
        C = G.base_change_matrix(rep)
        g = G.an_element()
        
        Dg = rep(g).matrix()
        A_block = C.inverse() * Dg * C
        
        # Check if the resulting matrix is diagonal
        self.assertTrue(A_block.is_diagonal())

    def test_base_change_matrix_new_off_filter_optimization_circulant(self):
        """
        Test from the docstring of: base_change_matrix_new_off_filter_optimization
        with a circulant matrix.
        """
        G = FiniteGroup(CyclicPermutationGroup(4))
        rep = G.natural_representation()
        A = matrix.circulant([1,2,3,4])
        self.assertTrue(rep.is_equivariant_to(A))

        P = G.base_change_matrix_new_off_filter_optimization(rep)
        A_block = P.inverse()*A*P
        
        # The expected result is a diagonal matrix
        expected_block = matrix(QQbar, [[10, 0, 0, 0], [0, -2, 0, 0], [0, 0, -2 - 2*sqrt(-1), 0], [0, 0, 0, -2 + 2*sqrt(-1)]])
        
        # We need to sort the diagonal elements for a stable comparison
        self.assertEqual(sorted(A_block.diagonal()), sorted(expected_block.diagonal()))
        
    def test_quick_block_prevision_s4_regular(self):
        """
        Test from the docstring of: quick_block_prevision with SymmetricGroup(4).
        """
        G = FiniteGroup(SymmetricGroup(4))
        reg = G.regular_representation()
        prevision = G.quick_block_prevision(reg)
        
        expected = [['degree', 'multiplicity'], [1, 1], [1, 1], [2, 2], [3, 3], [3, 3]]
        self.assertEqual(prevision, expected)

    def test_inner_product(self):
        """Tests the inner product calculation between representations."""
        
        # Test 1: Regular representation inner product with itself
        G = FiniteGroup(CyclicPermutationGroup(4))
        reg = G.regular_representation()
        self.assertEqual(reg.inner_product(reg), 4)  # Should be equal to group order
        
        # Test 2: Inner product with trivial representation
        n, irr = G.irreducible_representations(show_table=False)
        self.assertEqual(reg.inner_product(irr(0)), 1)  # Regular rep contains trivial rep once
        
        # Test 3: Inner product between different irreducible representations
        G = FiniteGroup(DihedralGroup(3))
        n, irr = G.irreducible_representations(show_table=False)
        # Different irreps should have inner product 0
        self.assertEqual(irr(0).inner_product(irr(1)), 0)
        
        # Test 4: Inner product with symbolic entries
        generators = ["(1,2,3,4,5,6)", "(1,4)(2,3)(5,6)"]
        matrices = [
            matrix(SR, [[1/2, -sqrt(3)/2], [sqrt(3)/2, 1/2]]),
            matrix(SR, [[-1, 0], [0, 1]])
        ]
        rep = representation(generators, matrices, field=SR)
        self.assertEqual(rep.inner_product(rep), 1)  # Irreducible rep has inner product 1 with itself

if __name__ == '__main__':
    unittest.main()
# cd pysymmetry/pysymmetry
# sage test.py -v

# or install pytest
#pytest -v pysymmetry/test.py
#pytest --html=report.html pysymmetry/test.py