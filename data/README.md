# Data folder:

This folder contains the data used for the analysis (no data included since it cannot be shared) and the file `create_data.R` which transforms the datasets used into the following format:


1. Early-insulin-treated patients:

| pardm<br>numeric | agerec<br>numeric | hba1c<br>numeric | agedx<br>numeric | sex<br>numeric | bmi<br>numeric | C<br>numeric | A<br>numeric |
|---------|---------|---------|---------|---------|---------|---------|---------|
| 1 - positive<br>0 - negative | \>1 or \<35 | \>3% or \< 15% | \>1 or \<120 | 1 - male<br>2 - female | \>14 or \<70 | 1 - positive<br>0 - negative | 1 - positive<br>0 - negative |


2. Non-early-insulin-treated patients:

| pardm<br>numeric | agerec<br>numeric | hba1c<br>numeric | agedx<br>numeric | sex<br>numeric | bmi<br>numeric | insoroha<br>numeric |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| 1 - positive<br>0 - negative | \>1 or \<35 | \>3% or \<15% | \>1 or \<120 | 1 - male<br>2 - female | \>14 or \<70   | 1 - positive<br>0 - negative |
