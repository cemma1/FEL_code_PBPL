% Compute_undulator_field

% Works with perave_core_v7
% Defines the resonant phase as a function of z 
% The Perave code will calculate the undulator field K(z) from res_phase(z) at each
% integration step

% This is if you wanna define the undulator K step-by-step using a
% pre-defined resonant phase
%{
startindex=floor(param.z0/param.stepsize);

if param.tapering==0
    res_phase=zeros(1,param.Nsnap);
else
    res_phase(1:startindex)=0;
    res_phase(startindex:param.Nsnap)=param.psir;
end
Kz(1)=param.K;
%plot([1:1:param.Nsnap]*param.stepsize,res_phase)
%xlabel('z [m]');ylabel('\psi_r');enhance_plot
%}

% Import external field
%Kz = importdata('external_magnetic_field.mat');

% Load in a polynomial taper profile

if param.tapering==0
    param.z0=lwig;
end
for i=1:param.Nsnap
    zval=i*param.stepsize;
    if zval>param.z0                
        Kz(i)=param.K*(1-param.ratio*(zval-param.z0)^param.order);
        csiz=Kz(i)^2/(4+2*Kz(i)^2);
        fbess1(i)=besselj(0,csiz)-besselj(1,csiz);                                               
    else
        Kz(i)=param.K;
        fbess1(i)=JJ1;
    end
end
clearvars zval 
input.Kz = Kz;
input.fbess1 = fbess1;