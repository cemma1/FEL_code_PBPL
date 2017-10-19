function [thetanew,etanew] = prebunch_particles(thetaold,etaold,param)

% Apply transformations to phase space (theta,gamma)
ps1(:,1)=thetaold+pi;
ps1(:,2)=(etaold-param.gamma0)./std((etaold-param.gamma0))-mean((etaold-param.gamma0)./std((etaold-param.gamma0)));
%plot(ps1(:,1),ps1(:,2));
%% 4 step transformation
%Abh=param.Abh;% This is if the bucket height value is given at the correct psir
Abh=param.Abh*sqrt(cos(param.psir)+param.psir*sin(param.psir)-pi/2*sin(param.psir));% This to convert from psir=0 to actual bucket height
% Numerical coefficiencts are fit for optimal buncher settings
% Procedure from Sudar et. al. paper http://www.sciencedirect.com/science/article/pii/S0168900217301924

A2=0.614+0.601*Abh-0.003*Abh^2;

A1=A2^2/(0.028*A2^2+3.811*A2-6.309);

B1=(0.086*A2^2+11.974*A2-19.822)/(1.056*A2^2+5.436*A2);

B2=(pi*A2)/(1.965*A2^2+9.486);

% Modulation
ps2(:,1)=ps1(:,1);
ps2(:,2)=ps1(:,2)+A1*sin(ps1(:,1));

% Compression
ps3(:,1)=ps1(:,1)+B1*ps2(:,2);
ps3(:,2)=ps1(:,2)+A1*sin(ps1(:,1));

% Second modulation
ps4(:,1)=ps3(:,1);
ps4(:,2)=ps3(:,2)+A2*sin(ps3(:,1));

% Compression
ps5(:,1)=ps4(:,1)+B2*ps4(:,2);
ps5(:,2)=ps4(:,2);
% Plot the different steps in the transformation if you want
% figure(1)
% subplot(2,2,1)
% plot(ps2(:,1),ps2(:,2),'.')
% legend(sprintf(['B=',num2str(abs(mean(exp(1i.*ps2(:,1)))))]))
% ylabel('\Delta \gamma/\sigma_\gamma')
% subplot(2,2,2)
% plot(ps3(:,1),ps3(:,2),'.')
% legend(sprintf(['B=',num2str(abs(mean(exp(1i.*ps3(:,1)))))]))
% ylabel('\Delta \gamma/\sigma_\gamma')
% subplot(2,2,3)
% plot(ps4(:,1),ps4(:,2),'.')
% legend(sprintf(['B=',num2str(abs(mean(exp(1i.*ps4(:,1)))))]))
% ylabel('\Delta \gamma/\sigma_\gamma')
% subplot(2,2,4)
% plot(ps5(:,1),ps5(:,2),'.')
% legend(sprintf(['B=',num2str(abs(mean(exp(1i.*ps5(:,1)))))]))
% ylabel('\Delta \gamma/\sigma_\gamma')
% set(findobj('type','axes'),'xlim',[0 2*pi],'xgrid','on')

thetanew=ps5(:,1)+param.psir-pi/2;
etanew=ps5(:,2)*param.deltagamma+param.gamma0;

