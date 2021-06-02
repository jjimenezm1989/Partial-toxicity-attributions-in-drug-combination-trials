
—————————————————————————————————————————————————————————————————————————

TITLE: CANCER PHASE I TRIAL DESIGN USING DRUG COMBINATIONS WHEN A FRACTION OF DOSE LIMITING TOXICITIES IS ATTRIBUTABLE TO ONE OR MORE AGENTS
AUTHORS: JOSE L JIMENEZ, MOURAD TIGHIOUART, MAURO GASPARINI.
E-MAILS: jose.jimenez@polito.it, Mourad.Tighiouart@cshs.org, gasparini@calvino.polito.it
AUTHOR RESPONSIBLE FOR THE CODE: Jose L Jimenez

—————————————————————————————————————————————————————————————————————————

VERSION INFORMATION


R version 3.3.1 (2016-06-21)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.12.6 (Sierra)

locale:
[1] es_ES.UTF-8/es_ES.UTF-8/es_ES.UTF-8/C/es_ES.UTF-8/es_ES.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] rjags_4-6   coda_0.18-1

loaded via a namespace (and not attached):
[1] tools_3.3.1     grid_3.3.1      lattice_0.20-33

—————————————————————————————————————————————————————————————————————————

FILE LIST:

JimenezTighiouartGasparini_BiometricalJournalScript1.R - R CODE
JimenezTighiouartGasparini_BiometricalJournalScript2.R - R CODE
JimenezTighiouartGasparini_BiometricalJournalScript3.R - R CODE
JimenezTighiouartGasparini_BiometricalJournalScript4.R - R CODE
JimenezTighiouartGasparini_BiometricalJournalScript5.R - R CODE
JimenezTighiouartGasparini_BiometricalJournalScript6.R - R CODE
JimenezTighiouartGasparini_BiometricalJournalScript7.R - R CODE
JimenezTighiouartGasparini_BiometricalJournalScript8.R - R CODE


—————————————————————————————————————————————————————————————————————————


R CODE: JimenezTighiouartGasparini_BiometricalJournalScript1.R


- This file simulates 1000 trials and generates the posterior estimates of the model parameters with continuous doses.
- This files saves the posterior estimates and produces the results of Table 2 in the manuscript.
- In the R file there are notes explaining the propose of each function and other notes.


—————————————————————————————————————————————————————————————————————————


R CODE: JimenezTighiouartGasparini_BiometricalJournalScript2.R


- This file is an intermediate step and calculates the posterior probability of DLT at the MTD curve. It is saved in a .csv file.
- This file needs the saved files from JimenezTighiouartGasparini_BiometricalJournalScript1.R.
- In the R file there are notes explaining the propose of each function and other notes.


—————————————————————————————————————————————————————————————————————————


R CODE: JimenezTighiouartGasparini_BiometricalJournalScript3.R


- This file rounds the MTD curve to the nearest discrete dose combination and produces Table 4.
- This file needs the saved files from JimenezTighiouartGasparini_BiometricalJournalScript1.R and JimenezTighiouartGasparini_BiometricalJournalScript2.R.
- In the R file there are notes explaining the propose of each function and other notes.


—————————————————————————————————————————————————————————————————————————


R CODE: JimenezTighiouartGasparini_BiometricalJournalScript4.R

- This program runs the design and computes operating characteristics under model misspecification. It produces Tables 6 and 7.
- It also produces Table S2 and S3 for the supplementary material.
- In the R file there are notes explaining the propose of each function and other notes.


—————————————————————————————————————————————————————————————————————————


R CODE: JimenezTighiouartGasparini_BiometricalJournalScript5.R

- This program generates Figures 2 and 3.
- In the R file there are notes explaining the propose of each function and other notes.


—————————————————————————————————————————————————————————————————————————


R CODE: JimenezTighiouartGasparini_BiometricalJournalScript6.R

- This program generates Figure 4.
- In the R file there are notes explaining the propose of each function and other notes.


—————————————————————————————————————————————————————————————————————————


R CODE: JimenezTighiouartGasparini_BiometricalJournalScript7.R

- This program generates Figure 5.
- In the R file there are notes explaining the propose of each function and other notes.


—————————————————————————————————————————————————————————————————————————


R CODE: JimenezTighiouartGasparini_BiometricalJournalScript8.R

- This program runs the design and computes operating characteristics under attribution misspecification.
- It produces Tables S4 and S5.
- In the R file there are notes explaining the propose of each function and other notes.


—————————————————————————————————————————————————————————————————————————


- Tables 1, 3, 5 and S1 do not require any R code to be produced.

- Since Figure 1 does not show any result, we do not provide R code.

- Since Figures 6, 7 and S1 only present a graphical representation of Tables 3, 5 and S1 respectively, we do not provide any R code. 
We obtained these Figures using the “persp” function following the information available in https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/persp.html





