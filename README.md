# PySymmetry

PySymmetry is a Sage-based Python package for working with finite groups and their representations. It provides tools for:
- Creating groups from generators.
- Computing regular, natural, and irreducible representations.
- Decomposing representations into isotypic components.
- Calculating projection operators, inner products, and change-of-basis matrices.
- Working with equivariant operators.

## Installation

1. Ensure you have [SageMath](https://www.sagemath.org/) installed.
2. Clone this repository:
   ```
   git clone <repository_url>
   ```
3. Navigate to the project folder:
   ```
   cd /home/leon-denis/JupyterLab/pysymmetry
   ```
4. Start Sage or run scripts using Sage’s environment.

## Usage

Import the package in your Sage worksheet or script:
```python
from pysymmetry import Group, representation, get_block
```

### Example: Regular Representation of a Cyclic Group

```python
# Create a cyclic group of order 4
G = Group(CyclicPermutationGroup(4))
# Compute the regular representation
reg = G.regular_representation()
# Display the matrix of a random element
g = G.an_element()
print(reg(g).matrix())
```

## API Overview

- **Group**: Extends Sage’s permutation group notion to work with representations.
- **representation**: Constructs a representation from a list of generators and their corresponding matrices.
- **MapRepresentation/Dg_Linear_Transformation**: Handle mapping and linear transformations associated with group elements.
- **Projection and Base Change Methods**: Compute isotypic projections, change-of-basis matrices, and block structures for equivariant operators.

For detailed API documentation, refer to inline docstrings in the source code.

## References

- Serre, J.-P. *Linear Representations of Finite Groups*. Springer, 1977.
- Stiefel, E. and Fässler, A. *Group Theoretical Methods and Their Applications*. Springer, 2012.

## License

Distributed under the terms of the MIT License. See `LICENSE` for more information.



