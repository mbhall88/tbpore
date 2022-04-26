# TBpore

*Mycobacterium tuberculosis* genomic analysis from Nanopore sequencing data

[![Python CI](https://github.com/mbhall88/tbpore/actions/workflows/ci.yaml/badge.svg)](https://github.com/mbhall88/tbpore/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/mbhall88/tbpore/branch/main/graph/badge.svg)](https://codecov.io/gh/mbhall88/tbpore)
[![PyPI](https://img.shields.io/pypi/v/tbpore)](https://pypi.org/project/tbpore/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/tbpore)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**⚠️WORK IN PROGRESS⚠️**

# Installation

## conda

`conda install tbpore`

## local

### Dependencies
* `git`
* `conda`

### Walkthrough

```
git clone https://github.com/mbhall88/tbpore
cd tbpore
conda env create -f environment.yaml && conda activate tbpore  # install dependencies
just install  # install tbpore
just check  # checks installation is fine
scripts/run_sample_example.sh  # if you want to run tbpore in a sample example
```

