%save_perave_output

%dirname=['Psir_',num2str(psirvalues(psirindex)),'/'];
dirname = 'Simulation_output/';
system(['mkdir ',dirname]);
save([dirname,'/simulation_parameters'],'param')
averagepower=mean(power,2);save([dirname,'/averagepower'],'averagepower')
save([dirname,'/average_bunching'],'bunch')
%save([dirname,'/average_energy'],'meanenergy')
save([dirname,'/undulator_field'],'Kz')
if param.itdp
[powerspec,omega]=spectrum_calc(radfield(end,:),param.lambda0,param.zsep);
save([dirname,'/output_spectrum'],'powerspec')
save([dirname,'/output_frequency'],'omega')

output_field=radfield(end,:);save([dirname,'/output_field'],'output_field')

output_power=power(end,:);save([dirname,'/output_power'],'output_power')
output_time=[1:1:size(power,2)]*param.zsep*param.lambda0;save([dirname,'/output_time'],'output_time')
end
if param.itdp==0   
   save([dirname,'/thetap'],'thetap') 
   save([dirname,'/gammap'],'gammap')
   save([dirname,'/radfield'],'radfield')
end
save([dirname,'/radfield'],'radfield')
