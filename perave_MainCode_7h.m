%%%% PBPL FEL simulation code %%%%
%%% Input deck compatible with wafFEL, 1D period average, and GENESIS %%%
clear all
close all
%% Load physical constants
physical_constants        
%% Load the User Determined initial conditions
Perave_User_Input_7h
%% Calculate 1-D FEL parameters
calculate_FEL_parameters
%% Compute the undulator field
compute_undulator_field_v7h
%% Initialize particle distribution
generate_perave_particles_v7h
%% Run the main integration routine
    t0 = tic;            
    perave_core_v7h;        
    disp(['Simulation time = ',num2str(toc(t0)./60),' min'])
%% Post-process output
perave_postprocessor_7h    
if param.saveoutput
save_perave_output 
end