README file for NIF_L1A
LAS -- 2023-07-24

Set of MATLAB files to: 

(1) Run the 2011 and 2012 NIF over L1A drainages:

2011: run2011_47_6z_20kmB_surfacecrack.m
2012: run2012_47_6z_20kmB_surfacecrack.m

Additional .m files needed for the NIF are in ../NIF_L1A/software (all /software files from Jeff McGuire).

N.B. Calculates and saves MLE_47_all_6z_2_20kmB_surfacecrack.mat and MLE2012_47_range7_gc_20km_surfacecrack.mat. 



(2) Calculate Green function matrices for positions on half-space surface:

run2011_47_6z_20kmB_Gsurface_surfacecrack.m

N.B. Calculates and saves Gsurface500_strain_surfacecrack.mat.



(3) Plotting script for JGR:ES Figure 2: displacement time series for 2012 L1A drainages. 

paperfig2_Stevens2015FigS6_2012FLOW_20230717.m

