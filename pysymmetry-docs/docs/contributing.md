# Contributing to PySymmetry

We welcome contributions from the community! Whether it's reporting a bug, proposing a new feature, or submitting code, your help is valuable.

## Bug Reports and Feature Requests

Please use the **Issues** tab on our GitHub repository to report bugs or suggest new features. Provide as much detail as possible, including steps to reproduce the bug.

## Development Setup

To contribute code, you'll need to set up a local development environment.

1.  **Fork the repository** on GitHub.
2.  **Clone your fork** locally:
    ```bash
    git clone [https://github.com/YOUR_USERNAME/pysymmetry.git](https://github.com/YOUR_USERNAME/pysymmetry.git)
    cd pysymmetry
    ```
3.  **Set up SageMath**: Ensure you have a working SageMath environment as described in the installation guide.

## Running Tests

Before submitting a contribution, please ensure that all tests pass. The tests for this project are located in `pysymmetry/test.py`.

You can run the test suite using the following command from the project's root directory:

```bash
sage -t pysymmetry/test.py -v
 ```
## Pull Request Process

1.  Create a new branch for your feature or bug fix.
2.  Make your changes and ensure all tests pass.
3.  Add docstrings for any new functions or classes.
4.  Submit a pull request to the `main` branch of the original repository.