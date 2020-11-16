import sys 
sys.path.append('/home/user/')
import unittest
from pysymmetry_01 import *
from sage.all import PermutationGroup, MatrixGroup, matrix, sqrt, ZZ, DihedralGroup, SymmetricGroup
#from parameterized import parameterized
from numpy import array, array_equal
from sage.matrix.matrix_integer_dense import Matrix_integer_dense

# 



class TestGroup(unittest.TestCase):
    
    def test_input(self):

        # Gap input
        gens = ["(1,2,3,4)", "(1,4)(2,3)", "(8,9)"]
        G = Group(gens)
        elements = G.list()
        expected = PermutationGroup(gens).list()
        test_item = set(elements) == set(expected)
        self.assertTrue(test_item)

        # matrix input
        mat = [matrix([[0, -1], [1, 0]]), matrix([[0, 1], [1, 0]])]
        G = Group(mat, matrix=True)
        elements = G.list()
        expected = MatrixGroup(mat).as_permutation_group().list()
        test_item = set(elements) == set(expected)
        self.assertTrue(test_item)
        # PermutationGroup

    def test_regular(self):

        gens = ["(1,2)", "(3,4)"]
        G = Group(gens)
        F = G.field
        elements = G.list()
        g1 = elements[1]
        g2 = elements[2]
        regular = G.regular_representation()
        expected1 = matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1],
                            [0, 0, 1, 0]]).change_ring(F)
        expected2 = matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0],
                            [0, 1, 0, 0]]).change_ring(F)

        self.assertTrue(regular(g1).matrix()==expected1)
        self.assertTrue(regular(g2).matrix()==expected2)
        
    def test_representation(self):
        
        # PermutationGroup
        D = DihedralGroup(3)
        generators = G.gens()
        matrices = [g.matrix() for g in generators]
        rep = representation(gens, matrices)
        
        


if __name__ == '__main__':

    unittest.main()

