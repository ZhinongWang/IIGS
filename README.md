This is the supporting information for the manuscript "Simultaneous Estimation of P- and S-wave Velocities by Integrated Inversion of Guided P- and Surface-wave Dispersion Curves".

1. The folder named “Individual inversion” contains the codes for implementing the individual inversion method. 
1.1 “main_individual.m” is the main function;
1.2 “R0.mat” and “R1.mat” are the observed surface-wave dispersion data;
1.3 “process.m” is used to arrange the input observed data;
1.4 “ragleigh_forward.m”, “ragleigh_fzero_brent.m” and “rayleighphase” are the functions used to solve the forward problem of surface-wave dispersion curves;
1.5 “rayleigh_jacobian_vs.m” is the program for calculating the Jacobian matrix of individual inversion method;
1.6 “reduce_delta.m”, “reduce_delta_dfda.m”, “reduce_delta_dffb.m”, “reduce_delta_dfdc.m” and "reduce_delta_dfdh.m" are the programs used to claculate the dispersion function and partial derivatives;
1.7 “bos2vp.m” is the function used to calculate the Vp from Vs and Poisson's ratio;
1.8 “plot_vf.m” and “plotvsthk” are the functions for plot the results.

2. The folder named “Integrated inversion” contains the codes for implementing the proposed integrated inversion method. 
2.1 “main_integrated.m” is the main function;
2.2 “R0.mat”, “R1.mat”, “G1.mat” and “G2.mat” are the observed surface- and guided-P dispersion data;
2.3 “process1.m” and “process2.m” are used to arrange the input observed surface- and guided-P dispersion data;
2.4 “guided_forward_vsvp.m”, “FastDelta2.m” and “guidedphase.m” are the functions used to solve the forward problem of guided-P wave dispersion curves;
2.5 “guided_jacobian_vsvp.m” is the program for calculating the Jacobian matrix of integrated inversion method;
2.6 “reduce_delta.m”, “reduce_delta_dfda.m”, “reduce_delta_dffb.m”, “reduce_delta_dfdc.m” and “reduce_delta_dfdh.m” are the programs used to claculate the disprsion function and partial derivatives;
2.7 "bos2vp.m" is the function used to calculate the Vp from Vs and Poisson's ratio;
2.8 "plot_vf.m" and "plotvsthk" are the functions for plot the results.

3. The folder named “Field data” contains the field seismic data used in this mauscript. 
