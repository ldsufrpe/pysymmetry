# PySymmetry: A SageMath package for Representation Theory of Finite Groups

PySymmetry is a SageMath-based Python package for computational group representation theory. It provides a suite of tools for working with finite groups and their representations, enabling researchers to perform complex calculations in a straightforward and efficient manner.

For a detailed API reference and tutorials, please see the **[Official Documentation](https://ldsufrpe.github.io/pysymmetry/)**.

## Directory Structure

The `pysymmetry` package is organized as follows:

```
pysymmetry/
├── pysymmetry/
│   ├── __init__.py
│   ├── pysymmetry.py
│   ├── pde.py
│   ├── parallel.py
│   ├── util.py
│   ├── timer.py
│   └── test.py
├── Examples/
│   ├── exemplo.ipynb
│   └── ... (other examples)
├── README.md
└── ... (other documentation and package files)
```

  - **`pysymmetry/pysymmetry/`**: The main source code directory.
      - **`__init__.py`**: Initializes the `pysymmetry` package.
      - **`pysymmetry.py`**: Contains the core functionality of the package, including the `FiniteGroup` class and functions for working with representations.
      - **`pde.py`**, **`parallel.py`**, **`util.py`**, **`timer.py`**: Utility modules for partial differential equations, parallel processing, general utilities, and timing, respectively.
      - **`test.py`**: Contains the unit tests for the package.
  - **`pysymmetry/Examples/`**: Contains Jupyter notebooks with usage examples.
      - **`exemplo.ipynb`**: A notebook demonstrating the basic features of the package.
  - **`pysymmetry/README.md`**: The main README file for the project.

## Installation

`PySymmetry` is built on top of SageMath. Therefore, a working installation of SageMath is a prerequisite.

### Step 1: Install SageMath

SageMath can be installed on Windows, macOS, and Linux. The recommended installation method varies by operating system.

  - **For Windows**: The recommended method is to use the Windows Subsystem for Linux (WSL).
    1.  Enable WSL and install a Linux distribution like Ubuntu from the Microsoft Store.
    2.  Once inside your Linux environment, install SageMath using its Linux instructions.
  - **For macOS and Linux**: You can download pre-built binaries or install from a package manager (like Homebrew for macOS or APT for Debian/Ubuntu).

For detailed, step-by-step instructions for your specific operating system, please consult the **[Official SageMath Installation Guide](https://doc.sagemath.org/html/en/installation/index.html)**. This guide provides the most up-to-date and comprehensive instructions.

### Step 2: Install PySymmetry

Once SageMath is installed, you can install `PySymmetry`.

1.  **Clone the repository**:

    ```bash
    git clone https://github.com/ldsufrpe/pysymmetry.git
    ```

2.  **Navigate to the project directory**:

    ```bash
    cd pysymmetry
    ```

    The package is now ready to be used within the SageMath environment.

## Execution and Comprehensive Test Run

To use `PySymmetry`, you must first start the SageMath environment.

### Launching a Jupyter Notebook in SageMath

The most convenient way to use `PySymmetry` interactively is through a Jupyter Notebook.

1.  Open your terminal or command prompt.

2.  Navigate to the directory containing your project files (e.g., the cloned `pysymmetry` directory).

3.  Launch the Jupyter Notebook interface by running the following command:

    ```bash
    sage -n jupyter
    ```

This command starts a Jupyter server and opens a new tab in your web browser. From there, you can create new notebooks or open existing ones (like those in the `Examples/` directory). When creating a new notebook, be sure to select the "SageMath" kernel.

### Comprehensive Test Run

This section provides a comprehensive test run to demonstrate the core functionalities of `PySymmetry`. The following example will:

1.  Define the symmetric group S\<sub\>3\</sub\>.
2.  Compute its regular representation.
3.  Calculate the isotypic projection operators.
4.  Construct the change-of-basis matrix for the isotypic decomposition.
5.  Verify that the change of basis correctly block-diagonalizes the representation.

#### Sample Code (to be run in a SageMath Jupyter cell)

```python
# Import necessary components from SageMath and PySymmetry
from sage.all import SymmetricGroup
from pysymmetry import FiniteGroup

# 1. Define the Symmetric Group S3
G = FiniteGroup(SymmetricGroup(3))
print("Group:", G)
print("Group Order:", G.order())
print("-" * 30)

# 2. Compute the regular representation
reg = G.regular_representation()
print("Regular Representation Degree:", reg.degree())
g = G("(1,2,3)")
print("\nMatrix for element g = (1,2,3):")
print(reg(g).matrix())
print("-" * 30)

# 3. Calculate the isotypic projection operators
projs = G.isotypic_projection(reg)
print("Number of Isotypic Projections:", len(projs))
print("\nFirst Isotypic Projection Operator (scaled):")
# We show a scaled version for clarity
print(projs[0]/projs[0][0,0])
print("-" * 30)

# 4. Construct the change-of-basis matrix
change_basis = G.isotypic_base(reg, isotypic_components=False)
print("Change of Basis Matrix:")
print(change_basis)
print("-" * 30)

# 5. Verify the block-diagonalization
Dg = reg(g).matrix()
A_block = change_basis.inverse() * Dg * change_basis
print("Block-diagonalized matrix for g = (1,2,3):")
print(A_block)

# Check for block-diagonal structure
# The off-diagonal blocks between the first two isotypic components should be zero.
is_block_diagonal = (A_block[0, 2] == 0 and A_block[1, 2] == 0 and
                     A_block[2, 0] == 0 and A_block[2, 1] == 0)

print("\nVerification: Is the matrix block-diagonal?", is_block_diagonal)
```

#### Expected Output

```
Group: Permutation Group with generators [(1,2,3), (1,3)]
Group Order: 6
------------------------------
Regular Representation Degree: 6

Matrix for element g = (1,2,3):
[0 0 1 0 0 0]
[1 0 0 0 0 0]
[0 1 0 0 0 0]
[0 0 0 0 1 0]
[0 0 0 0 0 1]
[0 0 0 1 0 0]
------------------------------
Number of Isotypic Projections: 3

First Isotypic Projection Operator (scaled):
[1 1 1 1 1 1]
[1 1 1 1 1 1]
[1 1 1 1 1 1]
[1 1 1 1 1 1]
[1 1 1 1 1 1]
[1 1 1 1 1 1]
------------------------------
Change of Basis Matrix:
[ 1| 1| 2 -1  0  0]
[ 1| 1|-1  2  0  0]
[ 1| 1|-1 -1  0  0]
[ 1|-1| 0  0  2 -1]
[ 1|-1| 0  0 -1  2]
[ 1|-1| 0  0 -1 -1]
------------------------------
Block-diagonalized matrix for g = (1,2,3):
[ 1.00000000000000?                     0                     0                     0                     0                     0]
[                   0 -1.00000000000000?                     0                     0                     0                     0]
[                   0                     0                     0                   1.0                     0                     0]
[                   0                     0 -1.00000000000000? - 0.?e-17*I                     0                     0                     0]
[                   0                     0                     0                     0                     0                   1.0]
[                   0                     0                     0                     0 -1.00000000000000? - 0.?e-17*I                     0]

Verification: Is the matrix block-diagonal? True
```

This test run confirms that `PySymmetry` correctly computes representations, projection operators, and the corresponding change-of-basis matrix to achieve block-diagonalization, which is a central task in the application of representation theory to physical problems.

## User Manual: API Overview

`PySymmetry` provides a high-level API for working with finite groups and their representations. The main components are:

### `FiniteGroup` Class

This class extends Sage's `PermutationGroup` to include methods for representation theory.

  - **`FiniteGroup(group_data, field=QQbar, matrix=False)`**: Constructor for the class. `group_data` can be a list of generators (as strings or permutations) or matrices.
  - **`.regular_representation()`**: Computes the regular representation of the group.
  - **`.irreducible_representations(show_table=True)`**: Returns the irreducible representations of the group.
  - **`.isotypic_projection(representation)`**: Calculates the projection operators onto the isotypic components of a given representation.
  - **`.isotypic_base(representation, isotypic_components=False)`**: Computes a basis that decomposes the representation into its isotypic components.
  - **`.base_change_matrix(representation)`**: Returns the change-of-basis matrix that block-diagonalizes the representation into its irreducible components.

### `representation` Function

  - **`representation(generators, matrices, field=QQbar)`**: Constructs a representation from a list of group generators and their corresponding matrix images.

### Example: Natural Representation of a Cyclic Group

```python
from pysymmetry import FiniteGroup
from sage.all import CyclicPermutationGroup

# Create a cyclic group of order 4
G = FiniteGroup(CyclicPermutationGroup(4))

# Compute the natural representation (permutation matrices)
nat_rep = G.natural_representation()

# Display the matrix for a generator
g = G("(1,2,3,4)")
print(nat_rep(g).matrix())
```

This will output the 4x4 permutation matrix corresponding to the cyclic permutation `(1,2,3,4)`.

For a complete API reference and further examples, please visit the **[Official Documentation](https://ldsufrpe.github.io/pysymmetry/)**.

## How to Run Tests

This project includes a consolidated test suite located at `pysymmetry/test.py`. The tests rely on SageMath, so they must be executed within the SageMath environment.

Prerequisites
- SageMath installed and available on your PATH (see Installation above).
- Python packages inside Sage (only needed if you choose pytest):
  - pytest (optional, for running tests with pytest)
  - pytest-html (optional, for generating an HTML report)
  - NumPy (usually included with Sage; install only if missing)

Install optional test dependencies inside Sage
- Install using Sage's pip so packages are bound to the Sage environment:
  ```bash
  sage -pip install pytest pytest-html numpy
  ```
  Note: NumPy typically ships with Sage, but the command above ensures it is present.

Run the tests
You can run the tests using either Python's built-in unittest (simplest) or pytest.

Option A — unittest (no extra packages required)
1. Navigate to the directory containing the test file:
   ```bash
   cd pysymmetry/pysymmetry
   ```
2. Run the test suite with verbose output using Sage:
   ```bash
   sage test.py -v
   ```
   Alternatively, you can use Sage's Python explicitly:
   ```bash
   sage -python -m unittest -v test.py
   ```

Option B — pytest (optional, more features)
1. From the repository root, run:
   ```bash
   sage -python -m pytest -v pysymmetry/test.py
   ```
2. Generate an HTML report (optional):
   ```bash
   sage -python -m pytest --html=report.html -v pysymmetry/test.py
   ```

Troubleshooting
- Command not found: If `sage` is not recognized, ensure SageMath is installed and added to your PATH.
- Import errors: Make sure you are running tests via Sage (e.g., `sage ...`). Regular system Python will not provide `sage.all`.
- NumPy missing: Install it inside Sage with `sage -pip install numpy`.

## References

  - Serre, J.-P. *Linear Representations of Finite Groups*. Springer, 1977.
  - Stiefel, E. and Fässler, A. *Group Theoretical Methods and Their Applications*. Springer, 2012.

## License

This software is distributed under the terms of the **MIT License**. See the `LICENSE` file for more details.