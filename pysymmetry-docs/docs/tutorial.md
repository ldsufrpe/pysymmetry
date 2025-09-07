# Tutorials

This page provides step-by-step tutorials for solving common problems using `PySymmetry`.

---


## Example 1: A D4-Symmetric System

This tutorial demonstrates how to block-diagonalize any matrix that is equivariant under a group's action. We will use the dihedral group $D_4$, which represents the symmetries of a square.

### Background

Imagine a system with 4 points arranged in a square. Any operator on this system that respects the square's symmetries (rotations and reflections) will have a matrix that commutes with the matrix representation of the $D_4$ group. Such a matrix is called "G-equivariant," and `PySymmetry` can simplify it.

### Implementation with `PySymmetry`

```python
from pysymmetry import FiniteGroup, representation
from sage.all import matrix, QQ

# Step 1: Define the D4 group using its generators
rotation = "(1,2,3,4)"
reflection = "(1,3)"
G = FiniteGroup([rotation, reflection])

# Step 2: Define the Natural Permutation Representation
phi = G.natural_representation()

# Step 3: Define a G-Equivariant Matrix M
a, b, c = 10, 2, 1
M = matrix(QQ, [[a, b, c, b],
                [b, a, b, c],
                [c, b, a, b],
                [b, c, b, a]])

# Step 4: Compute the Symmetry-Adapted Basis
P = G.base_change_matrix(phi)

# Step 5: Perform the Block-Diagonalization
M_block_diagonal = P.inverse() * M * P

print("Original Matrix M:")
show(M)
print("\nSymmetry-Adapted Basis P:")
show(P)
print("\nBlock-Diagonal Matrix P^{-1}MP:")
show(M_block_diagonal)
```

### Results

The output shows that the original matrix `M` is transformed into a block-diagonal matrix containing two 1x1 blocks and one 2x2 block. This decomposition simplifies further analysis, such as finding eigenvalues, which can now be done on the smaller blocks independently.


## Example 1: Molecular Vibrations in Chemistry (GF Method)

This tutorial showcases how `PySymmetry` can be applied to a real-world problem in chemistry: analyzing the vibrational frequencies of a molecule using the GF method.

### Background: The GF Method

The GF method is a classical approach to studying molecular vibrations. The core task is to solve the secular equation $(\mathbf{FG}-\lambda I) = 0$, where:

- **F** is the force-constant matrix (from potential energy).
- **G** is the inverse kinetic energy matrix (from atomic masses and geometry).

When a molecule like water ($H_2O$) has symmetry (in this case, the $C_{2v}$ point group), both **F** and **G** commute with the symmetry operations. This allows us to use `PySymmetry` to find a basis that block-diagonalizes the **FG** matrix, dramatically simplifying the eigenvalue problem.

### Implementation with `PySymmetry`

```python
from pysymmetry import FiniteGroup, representation
from sage.all import matrix, var, show

# Step 1: Define the C2v Group as a Permutation Group
generators = ["(1,2)", "(1)(2)(3)"]
C2v = FiniteGroup(generators)
gens = C2v.gens()

# Step 2: Define the Permutation Representation
C2_matrix = matrix(3, 3, [[0, 1, 0], [1, 0, 0], [0, 0, 1]])
E_matrix = matrix(3, 3, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
rep = representation(gens, [C2_matrix, E_matrix])

# Step 3: Find the Symmetry-Adapted Basis
beta = C2v.base_change_matrix(rep)
print("Symmetry-Adapted Basis (beta):")
show(beta)

# Step 4: Define Symbolic F and G Matrices
f11, f12, f33 = var('f11, f12, f33')
g11, g12, g13, g33 = var('g11, g12, g13, g33')

F_H2O = matrix(3, 3, [[f11, f12, 0], [f12, f11, 0], [0, 0, f33]])
G_H2O = matrix(3, 3, [[g11, g12, g13], [g12, g11, g13], [g13, g13, g33]])

# Step 5: Perform the Block-Diagonalization
FG_original = F_H2O * G_H2O
FG_block_diagonal = beta.inverse() * FG_original * beta

print("\nBlock-Diagonalized FG Matrix (beta^-1 * FG * beta):")
show(FG_block_diagonal)
```

### Results

By transforming to the basis `beta`, the `FG` matrix is converted into a block-diagonal form, separating the problem into a 2x2 block and a 1x1 block. This simplifies the calculation of its eigenvalues, which correspond to the molecule's vibrational frequencies.

---



## Exemplo 3: Exploiting Symmetry to Solve a Physics Problem

This tutorial demonstrates a core use case for `PySymmetry`: simplifying a common eigenvalue problem by exploiting the underlying symmetry of a physical system. We will find the eigenvalues of a 1D Laplacian operator, a task frequently encountered in physics and engineering.

The central idea is that the reflection symmetry of the 1D system allows us to block-diagonalize the Laplacian matrix. This breaks a large problem into smaller, independent ones, making the eigenvalue calculation significantly more efficient.

---

### Step 1: Setup and Problem Definition

First, we import the necessary libraries and define our operator. We will use a helper function from the `util` module to create the Laplacian matrix for a system with 100 points.

```python
# Import necessary components
from pysymmetry import FiniteGroup
from pysymmetry.util import laplacian1d, get_block
from sage.all import *
import numpy as np
```
Define the size of our system and create the Laplacian matrix

```python
n = 100
M = laplacian1d(n)
```


### Step 2: Define the Symmetry Group

A 1D system discretized into `n` points has a reflection symmetry about its center. This symmetry can be described by a group with a single operation: swapping point `j` with point `n - j + 1`. We create this group using its permutation generator.

```python
# Define the generator for the reflection group
def generators1d(n):
    reflection_str = ''
    # Use integer division for compatibility
    for j in range(1, (n // 2 + 1)):
        reflection_str += str((j, n - j + 1))
    return [reflection_str]
```

Create the FiniteGroup object
```python
G = FiniteGroup(generators1d(n))
print(f"Symmetry group created: {G}")
```

### Step 3: Compute the Symmetry-Adapted Basis

Now, we use `PySymmetry` to find the basis that respects the group's structure. This "symmetry-adapted basis" is the key to block-diagonalizing our matrix. We start by computing the natural representation of the group.

```python
# Get the natural (permutation) representation of the group
rep = G.natural_representation()

# Compute the basis that reduces the representation
base_info, _ = G.base_change_eigenvalue_reduction_new(rep)
```

### Step 4: Decompose the Matrix and Find Eigenvalues

With the symmetry-adapted basis, we can project the original large matrix `M` into smaller, independent blocks. The eigenvalues of these small blocks are the same as the eigenvalues of the original matrix, but are much easier to compute.

Use get_block to project M onto the subspaces defined by our new basis:
```python
blocks = [get_block(info[0], M) for info in base_info]
```

Calculate eigenvalues for each small block and then combine them
```python
eigenvalues_from_blocks = sorted(np.concatenate(
    [np.linalg.eigvals(b.toarray()) for b in blocks]
))

print(f"Successfully found {len(eigenvalues_from_blocks)} eigenvalues from the decomposed blocks.")
```

This completes the process. The `eigenvalues_from_blocks` list contains all the eigenvalues of the full `M` matrix. By using symmetry, we avoided direct computation on the large matrix and instead solved the problem on smaller, simpler ones.

