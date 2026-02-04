# LabPod

LabPod is a small Python toolkit for lab data analysis - inspired by lab work on Imperial's Medical Biosciences BSc.  
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
from labpod import protein_quantifier
protein_quantifier("path to you excel file", number of protein sample)
```
Function works for 96-well plates loaded in the manner.
![alt text](image.png)

## Output

The `protein_quantifier()` function generates a BSA standard curve using linear regression and calculates the concentration of unknown protein samples from their absorbance values.

### Example Results

The function returns:
- The fitted regression equation 
- The coefficient of determination (R²)
- A plotted BSA curve
- The calculated protein concentration for each unknown sample (labeled as in the image above) - assuming you sample is the same econcetration as your stock.

