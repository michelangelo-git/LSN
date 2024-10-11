# LSN - Numerical Simulation Lab

### Overview
This repository contains exercises assigned for the **Laboratorio di Simulazione Numerica** (Numerical Simulation Lab) course, held by Professor Galli at the Department of Physics, Università degli Studi di Milano.

### Repository Structure

#### Exercises 1, 2, 3, 8, 9, and 10
- Each exercise has its own directory (e.g., `es_1`, `es_2`, etc.).
- Source code is stored under the `src/` subdirectory: `es_n/src`.
- Header files are located in the `include/` subdirectory: `es_n/include`.
- Produced data is stored in `OUTPUT/`: `es_n/OUTPUT`.
- Input files are located in `INPUT/`: `es_n/INPUT`.
- A `Makefile` is provided for compilation:
  - Run `make` in the source directory to compile the code.
  - Run `make run` to execute the code.
  - For **Exercise 9**, there are two configurations:
    - Run `make circle` to use the circular starting configuration.
    - Run `make square` to use the square starting configuration.
    
#### Exercises 4, 6, and 7 (Molecular Dynamics)
- These exercises are grouped under the `NSL_SIMULATOR/` directory.
- To configure each simulation:
  - Modify the input configuration in `./INPUT/input.dat`.
  - Specify the property to measure during the simulation in `./INPUT/properties.dat`.
  - Define the output directory in `./INPUT/odir.dat` (e.g., `../OUTPUT/es_4/solid`).
  - To restart the system, specify the input directory in `./INPUT/idir.dat`.

#### Jupyter Notebooks
- Each exercise has its own Jupyter Notebook, stored in the respective exercise directory (e.g., `es_n.ipynb`).

#### Exercises 11 and 12 (Machine Learning & Neural Networks)
- These exercises consist solely of Jupyter Notebooks, written in Python. No additional source code is provided.

### Missing Exercise
**Exercise 5** is not included in this repository, as it is only assigned to Master’s degree students.
