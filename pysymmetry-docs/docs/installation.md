# Installation Guide

To use `PySymmetry`, you must have a working installation of **SageMath**. This guide provides detailed steps for setting up the environment.

---

## Step 1: Install SageMath

SageMath is a free, open-source mathematics software system that `PySymmetry` relies on. The installation process varies depending on your operating system.

### Windows

The recommended method for Windows is to use the **Windows Subsystem for Linux (WSL)**.

1.  **Enable WSL**: Open PowerShell as an administrator and run:
    ```bash
    wsl --install
    ```
    This command will enable the necessary features, download the latest Linux kernel, and install a Linux distribution (Ubuntu by default).

2.  **Set up your Linux environment**: Once installed, launch your Linux distribution (e.g., Ubuntu) from the Start Menu and follow the on-screen instructions to create a user account and password.

3.  **Install SageMath**: Inside your Linux terminal, you can install SageMath using its official binaries or from a package manager.

### macOS & Linux

For macOS and Linux, you can use pre-built binaries, which is the easiest method.

1.  **Download SageMath**: Visit the official **[SageMath download page](https://www.sagemath.org/download.html)** and download the binary for your operating system.
2.  **Extract the files**: Unpack the downloaded archive to a location of your choice.
3.  **Run SageMath**: You can run SageMath directly from the extracted folder.

For the most up-to-date and detailed instructions, always refer to the **[Official SageMath Installation Guide](https://doc.sagemath.org/html/en/installation/index.html)**.

---

## Step 2: Install PySymmetry

Once SageMath is installed and working, you can set up `PySymmetry`.

1.  **Clone the Repository**: Open a terminal and run the following command to clone the project from GitHub:
    ```bash
    git clone https://github.com/ldsufrpe/pysymmetry.git
    ```

2.  **Navigate to the Directory**:
    ```bash
    cd pysymmetry
    ```
    The library is now ready to be used within the SageMath environment.

---

## Step 3: Launch and Use

The best way to use `PySymmetry` is through a Jupyter Notebook within the SageMath environment.

1.  **Start Jupyter**: In your terminal, from the `pysymmetry` project directory, run:
    ```bash
    sage -n jupyter
    ```
2.  **Create or Open a Notebook**: This command will open a new tab in your browser. You can create a new notebook (be sure to select the "SageMath" kernel) or open one of the examples from the `Examples/` directory.