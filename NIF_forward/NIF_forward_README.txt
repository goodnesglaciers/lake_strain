README file for NIF_forward
LAS -- 2023-06-30

Set of MATLAB files to: 

(1) multiply Green function matrices for half-space surface by NIF opening and slip from to forward model displacement, strain, and stress at the surface of the elastic half space. 

2011: Nsurface_displ_strain_stress_2011_20230531.m 
2012: Nsurface_displ_strain_stress_2012_20230619.m

Plots a helpful figure of (columns) surface displacement, strain, and stress for (rows) xx, yy, and zz for [opening and slip] along the basal and vertical planes. Saves figure to directory.

(2) plotting script for JGR:ES Figures 2 and S1: principal stresses for winter velocities and L1A lake-drainage deformation.

2011: paperfig2_nif2011_max20kmB_200kPa_20230621.m
2012: paperfigS1_nif2012_max20kmB_200kPa_20230621.m

(3) plotting script for JGR:ES Figure 6: stress decomposition time series for 2011 and 2012 L1A drainages. Calculates and saves n2011_allfour.mat and n2012_allfour.mat.

paperfig6_decomp_elastic_2011_2012_maxHF_20230619.m

(4) plotting script for JGR:ES Figure 7: regress stress decomposition against \sigma_{1}.

paperfig7_stressdecomposition_regression_20230619.m 

(5) plotting script for JGR:ES Figure 8: stress rotation between \sigma_{1,winter} and \sigma_{1}(t_{3}) for 2011 and 2012 L1A drainages.

paperfig8_roseplot360_princstress_2011_2012_20230619.m

(6) .m files needed for the NIF are in /software (all files from Jeff McGuire)
