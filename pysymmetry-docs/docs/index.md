# Welcome to PySymmetry

**A Python library for applying group representation theory to problems with symmetry.**

---

### About the Project

`PySymmetry` is a powerful tool built on SageMath designed to simplify complex computational problems by leveraging their underlying symmetries. Many problems in physics, engineering, and mathematics involve systems that are symmetric under certain transformations (like rotations or reflections). These symmetries can be mathematically described using group theory.

The core feature of this library is the **block diagonalization of equivariant operators**. By finding a symmetry-adapted basis, `PySymmetry` can transform a large, complex matrix into a set of smaller, independent block matrices. This decomposition significantly simplifies subsequent calculations, such as finding eigenvalues, and can lead to substantial performance improvements.

### Key Features

* **FiniteGroup and Representation Theory Tools:** Provides an intuitive interface for working with permutation groups and their matrix representations.
* **Automatic Block Decomposition:** Implements algorithms to find the symmetry-adapted basis and automatically block-diagonalize equivariant matrices.
* **Performance:** Accelerates calculations by breaking down large problems into smaller, more manageable ones. The library includes parallel processing capabilities to speed up computations.
* **Built on SageMath:** Leverages the extensive mathematical capabilities of the SageMath ecosystem.

### Quick Example: Eigenvalues of a 1D Laplacian

Let's find the eigenvalues of a 1D Laplacian operator, a common problem in physics. Using `PySymmetry`, we can exploit the reflection symmetry of the system to simplify the calculation.

```python
# 1. Import the library and create the operator
from pysymmetry import FiniteGroup, representation
from pysymmetry.util import laplacian1d, get_block
import numpy as np

n = 100
M = laplacian1d(n)

# 2. Define the symmetry group (in this case, a simple reflection)
def generators1d(n):
    string_reflexao_sigma = ''
    # Note: Integer division for compatibility
    for j in range(1, (n // 2 + 1)):
        string_reflexao_sigma += str((j, n - j + 1))
    return [string_reflexao_sigma]

G = FiniteGroup(generators1d(n))

# 3. Get the natural representation and the symmetry-adapted basis
rep = G.natural_representation()
base_info, _ = G.base_change_eigenvalue_reduction_new(rep)

# 4. Get the blocks from the original matrix
# This is much faster than calculating eigenvalues for the full matrix M
blocks = [get_block(info[0], M) for info in base_info]
eigenvalues_from_blocks = sorted(np.concatenate(
    [np.linalg.eigvals(b.toarray()) for b in blocks]
))
```

The `eigenvalues_from_blocks` will be identical to the eigenvalues of the full `M` matrix, but they are computed more efficiently from the smaller blocks.

### Next Steps

- API Reference: Dive into the details of the available functions and classes.

- Installation: Learn how to install pysymmetry in your own projects.
