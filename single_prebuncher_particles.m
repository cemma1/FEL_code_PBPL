function [thetanew,etanew] = single_prebuncher_particles(thetaold,etaold,param)

ps1(:,1)=thetaold+pi;

ps1(:,2)=(etaold-param.gamma0)./std((etaold-param.gamma0))-mean((etaold-param.gamma0)./std((etaold-param.gamma0)));
%plot(ps1(:,1),ps1(:,2));
%% 4 step transformation
%Abh=param.Abh;% This is if the bucket height value is given at the correct psir
Abh=param.Abh*sqrt(cos(param.psir)+param.psir*sin(param.psir)-pi/2*sin(param.psir));% This to convert from psir=0 to actual bucket height

% Numerical value from N. Sudar optimization (see paper ref)

A2=1.3268590514556047+0.4035665622224646*Abh;

B2=0.2011812564233884+0.0006328650045684311*exp(71.01657740614692/(8.488984985248903+Abh));

%Modulator
ps2(:,1)=ps1(:,1);
ps2(:,2)=ps1(:,2)+A2*sin(ps1(:,1));
%Chicane
ps3(:,1)=ps2(:,1)+B2*ps2(:,2);
ps3(:,2)=ps2(:,2);

thetanew=ps3(:,1)+param.psir-pi/2;
etanew=ps3(:,2)*param.deltagamma+param.gamma0;

%Plot the transformation if you want
%figure(2)
%h1=plot(thetaold,etaold)
%hold on
%h2=plot(thetanew,etanew,'.r')
%xlim([-pi,pi])
%legend([h1,h2],{'No bunching',sprintf(['B=',num2str(abs(mean(exp(1i.*thetanew))))])});
%legend boxoff
%xlabel('\theta');ylabel('\gamma')
