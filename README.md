# LabPod

LabPod is a small Python toolkit for lab data analysis - inspired by Lab work on Imperial's Medical Biosciences BSc.  
Currently includes `protein_quantifier()` for BSA standard curve fitting + protein concentration calculation from your fixed Excel layout - example below.

[![PyPI version](https://img.shields.io/pypi/v/LabPod.svg)](https://pypi.org/project/LabPod/)
[![Downloads](https://static.pepy.tech/badge/LabPod)](https://pepy.tech/project/LabPod)

## Install

```bash
pip install LabPod
```

## How to use
In python file run: 
```python
protein_quantifier("path to you excel file", number of protein sample)
```
Function works for 96-well plates loaded in the manner.
![alt text](image.png)

## Output

The `protein_quantifier()` function generates a BSA standard curve using linear regression and calculates the concentration of unknown protein samples from their absorbance values.

### Example Results

The function returns:
- The fitted regression equation (slope and intercept)
- The coefficient of determination (R²)
- The calculated protein concentration for the unknown sample - assuming you sample is the sam econcetration as your stock.

Example console output:

```text
Slope: 0.00124
Intercept: 0.0321
R²: 0.998

Protein concentration (unknown sample): 1.37 µg/mL
