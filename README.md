# cyclic-gMISpy

A Python package for calculating Genetic Minimal Intervention Sets (gMIS) and Minimal Cut Sets (MCS) in metabolic networks with regulatory cycles.

## Overview

cyclic-gMISpy is a computational tool designed for identifying minimal genetic interventions in metabolic networks that consider regulatory cycles and gene-protein-reaction (GPR) rules. The package extends traditional minimal cut set analysis by incorporating regulatory network dynamics and cyclic dependencies.

## Features

- **Minimal Cut Set (MCS) Calculation**: Find minimal sets of reactions that disrupt metabolic objectives
- **Genetic Minimal Intervention Sets (gMIS)**: Identify minimal gene knockouts considering GPR rules
- **Regulatory Cycles Support**: Handle cyclic dependencies in gene regulatory networks
- **Parallel Computing**: Multi-threaded calculations for improved performance
- **Multiple Solvers**: Support for Gurobi and CPLEX optimization solvers

## Installation

### Prerequisites

- Python ≥ 3.7
- COBRApy
- NumPy
- SciPy
- pandas
- tqdm
- bidict
- bonesis
- mpbn

### Required Optimization Solvers

At least one of the following optimization solvers:
- **Gurobi** (recommended for performance)
- **CPLEX**


### Install Package

```bash
# Clone the repository
git clone https://github.com/your-username/cyclic-gMISpy.git
cd cyclic-gMISpy

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

Or use pip install via piptest
```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ cyclic-gmispy
```



## Quick Start

### Basic MCS Calculation

```python
import cobra
from cyclic_gmis import calculateMCS

# Load your metabolic model
model = cobra.io.read_sbml_model("your_model.xml")
```



### gMIS with Regulatory Cycles

```python
import pandas as pd
from cyclic_gmis import calculateGMIS

# Load metabolic model
model = cobra.io.read_sbml_model("your_model.xml")

# Load regulatory network data
regulatory_df = pd.read_csv("regulatory_network.csv")

from cyclic_gmis import calculateParallelGMIS

# Run parallel gMIS calculation
parallel_results = calculateParallelGMIS(
    cobraModel=model,
    regulatory_dataframe=regulatory_df,
    numWorkers=8,
    maxKOLength=3
)
```

## Core Components

### Main Functions

- `calculateParallelGMIS()`: Calculate genetic minimal intervention sets with regulatory cycles
- `calculateReadjustmentParallel()`: Calculate readjustment cycles efficiently

### Utilities

- `buildGMatrix()`: Build gene-protein-reaction matrix
- `checkGMCS()`: Validate computed genetic minimal cut sets
- `saveSolutions()`: Save results in various formats
- `setSolver()`: Configure optimization solver settings

### Problem Definitions

- Support for different optimization problem formulations
- Flexible constraint handling
- Multiple objective functions

## Configuration Options



## Input Data Formats

### Metabolic Models
- SBML format (via COBRApy)
- JSON format
- MAT files

### Regulatory Networks
- CSV files with regulator-target relationships
- Tab-separated files
- Custom pandas DataFrame format

## Citing

If you use cyclic-gMISpy in your research, please cite:

```
[Pending]
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions and support:
- Open an issue on GitHub
- Contact: [cjrodriguezf@unav.es]

## Acknowledgments

- COBRApy community for metabolic modeling tools
- Bonesis framework for regulatory network analysis
- Optimization solver providers (Gurobi, CPLEX, SCIP)
