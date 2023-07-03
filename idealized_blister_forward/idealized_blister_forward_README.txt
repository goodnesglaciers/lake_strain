README file for idealized_blister_forward
LAS -- 2023-06-30

Set of MATLAB files to: 

(1) calculate idealised blister geometry given an initial lake volume following Lai et al. (2021). Save blister opening and a slip patch in format needed for okada85. 

within /make_idealized_blisters: run idealized_blister.m 

Plots some helpful figures of blister geometry and a Y/N indicator of where slip occurs.

(2) populate Green function matrices for slip and opening along a 20 x 20 km horizontal plane located at some depth within an elastic halfspace. 

within /forward_okada85: run run_Gmatrix_20km_bed_thickness.m

** This run file stops after Green function matrices are populated. Does not complete inversion. **

(3) multiply Green function matrices from (2) by idealized blister opening and slip from (1) to forward model displacement, strain, and stress at the surface of the elastic half space. N.B. Only considers slip and opening along the bed (no mode-1 opening of a vertical fracture). 

within /forward_okada: Nsurface_displ_strain_stress_*m.m 

Plots a helpful figure of (columns) surface displacement, strain, and stress for (rows) xx, yy, and zz for [opening], [slip], and [opening and slip] defined along the basal plane. Saves files to /forward_okada/Forward_*m_check_plots.

(4) plotting script for JGR:ES Figure 3: idealized blister example. 

within /forward_okada: paperfig3_idealized_example_15GPa_20230630.m 

(5) plotting script for JGR:ES Figure 10: idealized blisters over parameter range. Vary material properties here to create Figures S5 and S6.

within /forward_okada: paperfig10_idealized_15GPa_20230630.m 
within /forward_okada: paperfigS5_idealized_032GPa_20230630.m 
within /forward_okada: paperfigS6_idealized_39GPa_20230630.m 

(6) .m files needed for the NIF are in /software (all files from Jeff McGuire)

N.B. File/script redundancy due to parameter-space ranges for ice thickness and lake volume. Adjust as needed for the geometry and lake volume one is interested in.

Range of ice thicknesses tested: [500:250:3000]; % [ m ]

Range of lake volumes tested: [0.001, 0.005, 0.01, 0.05, 0.1]; % [ km^3 ] 