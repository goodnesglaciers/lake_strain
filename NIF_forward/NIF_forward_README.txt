README file for NIF_forward
LAS -- 2023-06-30

Set of MATLAB files to: 

(1) multiply [Green function matrices for positions on half-space surface] by [NIF opening and slip] to forward model displacement, strain, and stress at the surface of the elastic half space. 

2011: Nsurface_displ_strain_stress_surfC_2011_20230721.m 
2012: Nsurface_displ_strain_stress_surfC_2012_20230724.m

N.B. Requires Gsurface500_strain_surfacecrack.mat, MLE_47_all_6z_2_20kmB_surfacecrack.mat, and MLE2012_47_range7_gc_20km_surfacecrack.mat from ../NIF_L1A. Additional .m files needed for the NIF are in ../NIF_L1A/software (all /software files from Jeff McGuire).

Plots a helpful figure of (columns) surface displacement, strain, and stress for (rows) xx, yy, and zz for [opening and slip] along the basal and vertical planes. Saves figure to directory.



(2) plotting script for JGR:ES Figure 7: stress decomposition time series for 2011 and 2012 L1A drainages. 

paperfig7_decomp_elastic_2011_2012_maxHF_surfC_20230724.m

N.B. Calculates and saves n2011_allfour_surfC.mat and n2012_allfour_surfC.mat.



(3) plotting script for JGR:ES Figures 3 and S1: principal stresses for winter velocities and L1A lake-drainage deformation.

2011: paperfig3_nif2011_max20kmB_200kPa_surfC_20230721.m
2012: paperfigS1_nif2012_max20kmB_200kPa_surfC_20230724.m

N.B. Requires Gsurface500_strain_surfacecrack.mat from ../NIF_L1A.



(4) plotting script for JGR:ES Figure 8: regress stress decomposition against \sigma_{1}.

paperfig8_stressdecomposition_regression_surfC_20230721.m 

N.B. Requires Gsurface500_strain_surfacecrack.mat from ../NIF_L1A.



(5) plotting script for JGR:ES Figure 9: stress rotation between \sigma_{1,winter} and \sigma_{1}(t_{3}) for 2011 and 2012 L1A drainages.

paperfig9_roseplot360_princstress_2011_2012_20230724.m

N.B. Requires Gsurface500_strain_surfacecrack.mat from ../NIF_L1A. Stylistic modifications made in Keynote.



(6) plotting scripts for JGR:ES Movies M1 and M1.

2011: papermovieM1_princstress_dial_allfour_2012_20230724.m
2012: papermovieM2_princstress_dial_allfour_2012_20230724.m

N.B. Requires n2011_allfour_surfC.mat and n2012_allfour_surfC.mat.
