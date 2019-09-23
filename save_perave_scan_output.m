%save_perave_output
    if ~exist('Psir_scan')      
    system('mkdir Psir_scan');
    
    end
dirname=['Psir_scan/Psir_',num2str(psirvalues(psirindex)),'/'];

system(['mkdir ',dirname]);
save([dirname,'/simulation_parameters'],'param')
averagepower=mean(power,2);save([dirname,'/averagepower'],'averagepower')
save([dirname,'/average_bunching'],'bunch')
%save([dirname,'/average_energy'],'meanenergy')
save([dirname,'/undulator_field'],'Kz')

if param.itdp==0   
   save([dirname,'/thetap'],'thetap') 
   save([dirname,'/gammap'],'gammap')
   save([dirname,'/radfield'],'radfield')
end
save([dirname,'/radfield'],'radfield')
save([dirname,'/res_phase'],'res_phase')
if param.itdp
[powerspec,omega]=spectrum_calc(radfield(param.Nsnap,:),param.lambda0,param.zsep);
save([dirname,'/output_spectrum'],'powerspec')
save([dirname,'/output_frequency'],'omega')
end
function [power_spectrum,omega,sigma]=spectrum_calc(field,xlamds,zsep)
% Do the FFT
    ft=fftshift(fft(field));
    power_spectrum=abs(ft).^2;
    
% Calculate the frequency range    
	nslice=size(field,2);
	omegas=2*pi/(xlamds/3e8);
    dt = nslice*zsep*xlamds/3e8;
    df = 2*pi/dt;
    omega=df*(1:length(ft))';
	omega=omega-median(omega)+omegas;
	omega=omega/omegas;
    %Center the spectrum so you have delta omega/omega on x-axis
    omega=omega-1; 
end
