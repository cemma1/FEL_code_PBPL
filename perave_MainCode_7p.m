%%%% 1d FEL simulation code %%%%
% Check for version compatibility
if verLessThan('matlab','9.1')
warning('Including functions in scripts requires MATLAB R2016b or later')
warning('Functions defined inside perave_postprocessor script need to be defined as separate functions')
end
clear all 
close all
%% Load physical constants
physical_constants        
%% Load the User Determined initial conditions
Perave_User_Input_7p
%% Calculate 1-D FEL parameters
calculate_FEL_parameters_p
%% Compute the undulator field
compute_undulator_field_v7p
%% Initialize particle distribution
generate_perave_particles_v7p
%% Run the main integration routine
    t0 = tic;            
    perave_core_v7p;        
    disp(['Simulation time = ',num2str(toc(t0)),' sec'])
%% Post-process output
perave_postprocessor_7p    
if param.saveoutput
save_perave_output 
end