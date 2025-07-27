# API Reference

This page provides a detailed reference for the `pysymmetry` library's public API. The content is automatically generated from the docstrings in the source code.

---

## Core Module (`pysymmetry`)

This is the main module containing the core classes and functions for working with group theory and representations.

::: pysymmetry.Group
    options:
      show_root_heading: true
      show_source: false

::: pysymmetry.representation
    options:
      show_root_heading: true

::: pysymmetry.get_block
    options:
      show_root_heading: true

---

## PDE & Numerical Module (`pde`)

This module contains specialized, performance-oriented classes and functions, often used for numerical applications like solving partial differential equations.

::: pysymmetry.pde.nGroup
    options:
      show_root_heading: true
      show_source: false

---

## Utilities Module (`util`)

This module provides helper functions for creating common matrices and other utilities.

::: pysymmetry.util.laplacian1d
    options:
      show_root_heading: true

::: pysymmetry.util.laplacian2d
    options:
      show_root_heading: true

---

## Parallel Processing (`parallel`)

This module contains tools for parallelizing computations to improve performance on multi-core systems.

::: pysymmetry.parallel.pmap
    options:
      show_root_heading: true

::: pysymmetry.parallel.pmap_reduce
    options:
      show_root_heading: true
