% Compute_undulator_field

% Works with perave_core_v7
% Defines the resonant phase as a function of z 
% The Perave code will calculate the undulator field K(z) from res_phase(z) at each
% integration step

startindex=floor(param.z0/param.stepsize);
res_phase = [];
if param.tapering==0
    res_phase=zeros(1,param.Nsnap);
else
    res_phase(1:startindex)=0;
    res_phase(startindex:param.Nsnap)=param.psir;% For const psir taper        
end
Kz(1)=param.K;
%plot([1:1:param.Nsnap]*param.stepsize,res_phase)
%xlabel('z [m]');ylabel('\psi_r');enhance_plot
