function [power_spectrum,omega,sigma]=g3spectrum2(field,xlamds,zsep)

	nslice=size(field,2);
	%nz=size(field,2);
	omegas=2*pi/(xlamds/3e8);
	df=2*pi/nslice/zsep/(xlamds/3e8);
	
	%ft=zeros(size(field,1),1);
    ft=fftshift(fft(field));

	%omega=df*(1:size(ft,2))';
    omega=df*(1:length(ft))';
	omega=omega-median(omega)+omegas;
	omega=omega/omegas;
	power_spectrum=abs(ft).^2;
	
% 	ssum=sum(power_spectrum);
% 	mu=omega'*power_spectrum./ssum;
% 	for i=1:nz
% 		sigma(1,i)=sqrt( ((omega'-mu(i)).^2) * power_spectrum(:,i) /ssum(i));
%     end
    
omega=omega-1;
%spec=smooth(spec);