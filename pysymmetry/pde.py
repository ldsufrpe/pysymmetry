__all__ = [
 'nGroup',
 'ngroup',  
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



class gcsr_matrix(csc_matrix): # classe array  com atributo index

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        #return csr_matrix(*args, **kwargs)
    def tomatrix(self):#Nota: to rational sparse
    
        dic = self.todok()
        dic_new = dict()
        for key, value in dic.items(): 
            #dic_new[key]=QQ(int(value))
            dic_new[key]=QQ(int(value))# NOta: se for complexo?
            
        return matrix(dic_new)

    def issymmetric(self, rtol=1e-10, atol=1e-12):
        arr = self.toarray()
        return np.allclose(arr, arr.T, rtol=rtol, atol=atol)

    def character(self):
        return self.diagonal().sum()

class nIsotypicBase(object):
    """docstring for ClassName"""
    def __init__(self, base):
        #super(IsotypicBase, self).__init__()
        
        self._base = list(filter(lambda x: x!=None, base))
        self.nisotypic_components = [b[0] for b in  self._base]
        self._info = [b[1] for b in  self._base]

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
        
    #g = words[0].parent()(g)
    H = libgap.FiniteGroup(words)
    ans = H.EpimorphismFromFreeGroup().PreImagesRepresentative(g)
    l1 = str(ans)
    l1 = avoid_inverse(g, words, l1)    
    l2 = l1.replace('^', '**')
    
    def locals2():
        local = dict()
        for i, w_i in enumerate(matrices):
            local['x'+ str(i+1)] = w_i        
        return local    
    return sage_eval(l2, locals=locals2())




class ngroup(PermutationGroup_generic):

    def __init__(self, gens, field=RDF):

       
        self.field = field

        super().__init__(gens, gap_group=None, domain=None, canonicalize=True, category=None)

    def nrepresentation(self, gens, image):
        idty_group = self.identity()
        image = [to_csc(m) for m in image]
        dic = dict()
        _, n  = image[0].get_shape()
        elements = self.list()
        for g in elements[1:]:# Nota exceto o ()     
            D = nfactorization(g, gens, image)
            dic[g] = gcsr_matrix(D) 
        dic[idty_group] = gcsr_matrix(identity(n, format='csc'))
        return dic

    def nirreducible_representations(self): 
          
        irr = gap.IrreducibleRepresentations(self)
        gens = gap.GeneratorsOfGroup(self)
        generators = list(gap.List(gens))
        #generators = [self(g) for g in generators]
        list_irr =  OrderedDict()       
        for i ,s in enumerate(irr):
            image = [matrix(sage_eval(str(gap.Image(s, g)))) for g in generators]#Nota: change_ring verificar se atrapala ou ajuda.
            rep = self.nrepresentation(generators, image)
            list_irr[i] = rep           
        N = len(list_irr)         
        return N, list_irr  
   

    
    def nprojection(self, i, j, k, all_irr, right):
        
        left = all_irr          
        s = sum([left[i][g.inverse()][j,k]*right[g] for g in self])
        
        return s
    
    def nisotypic_component(self, i, left, right):#Nota: implementar a exceção multiplicity.radical_expression()

        multiplicity = ninner_product(self, left[i], right)
            
        if not multiplicity==0:
            
            P = self.nprojection(i, 0, 0, left, right)
            P_copy = csr_matrix.copy(P)               
            pivots = P.tomatrix().pivots()
            degree = ndegree(self, left[i])
            v0 = P_copy[:,pivots]               
            return [v0,(degree, multiplicity)]        


    def nbase_change_reduction(self, right):
        
        n, l = self.nirreducible_representations()
        r = right      
        base = pmap(lambda i: self.nisotypic_component(i, l, r), list(range(n)))
       
        return nIsotypicBase(base)

nGroup = ngroup 



def ninner_product(G, left, right):
    
    l = np.array([left[g.inverse()].diagonal().sum()*right[g].diagonal().sum()  for g in G])
    s = np.sum(l)    
    return int((1./G.order())*s)

def ndegree(G, rep):
    _, d = rep[G.an_element()].get_shape()       
    return d

from scipy.linalg import pinv


def nget_block(columm_base, matrix_equiv): # Nota: traduzir nome das variaveis

    columm = columm_base.toarray()
    p = csc_matrix(np.linalg.pinv(columm))#Nota:Aqui da pra configurar pelo Sage pedido o algoritmo numpy, ou chamar outro metodo.
    block = p*matrix_equiv*columm_base
    return gcsr_matrix(block)



