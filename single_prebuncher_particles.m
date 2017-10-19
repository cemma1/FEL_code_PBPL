function [thetanew,etanew] = single_prebuncher_particles(thetaold,etaold,param)
%% Make a normal distribution
%clear all
%close all
%npart=1e4;
%gamma0=1e4;
%sigmagamma=1;

ps1(:,1)=thetaold+pi;
%ps1(:,2)=(normrnd(gamma0,sigmagamma,paramnpart,1)-gamma0)/sigmagamma;
ps1(:,2)=(etaold-param.gamma0)./std((etaold-param.gamma0))-mean((etaold-param.gamma0)./std((etaold-param.gamma0)));
%plot(ps1(:,1),ps1(:,2));
%% 4 step transformation
%Abh=param.Abh;% This is if the bucket height value is given at the correct psir
Abh=param.Abh*sqrt(cos(param.psir)+param.psir*sin(param.psir)-pi/2*sin(param.psir));% This to convert from psir=0 to actual bucket height
%disp(Abh)

A2=1.3268590514556047+0.4035665622224646*Abh;

B2=0.2011812564233884+0.0006328650045684311*exp(71.01657740614692/(8.488984985248903+Abh));

%Modulator
ps2(:,1)=ps1(:,1);
ps2(:,2)=ps1(:,2)+A2*sin(ps1(:,1));
%Chicane
ps3(:,1)=ps2(:,1)+B2*ps2(:,2);
ps3(:,2)=ps2(:,2);

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

thetanew=ps3(:,1)+param.psir-pi/2;
etanew=ps3(:,2)*param.deltagamma+param.gamma0;
%etanew=ps5(:,2).*std((etaold-param.gamma0))+mean((etaold-param.gamma0)./std((etaold-param.gamma0)));

%Plot stuff if you want
figure(2)
h1=plot(thetaold,etaold)
hold on
h2=plot(thetanew,etanew,'.r')
xlim([-pi,pi])
legend([h1,h2],{'No bunching',sprintf(['B=',num2str(abs(mean(exp(1i.*thetanew))))])});
legend boxoff
xlabel('\theta');ylabel('\gamma')
pause

%% Bunch stuff old code
% Abh=10;
% 
% A2=0.614+0.601*Abh-0.003*Abh^2;
% 
% A1=A2^2/(0.028*A2^2+3.811*A2-6.309);
% 
% B1=(0.086*A2^2+11.974*A2-19.822)/(1.056*A2^2+5.436*A2);
% 
% B2=(pi*A2)/(1.965*A2^2+9.486);
% 
% etanew = etaold+A1*sin(thetaold)+...
%        A2*sin(thetaold+B1*(etaold+A1*sin(thetaold)));
%    
% thetanew = thetaold+(B1+B2)*etaold+(B1*A1+B2*A2)*sin(thetaold)+...
%            B2*A2*sin(thetaold+B1*etaold+B1*A1*sin(thetaold));
%        
%        figure
%        subplot(2,1,1)
%        set(gca,'FontSize',16)
%        plot(thetaold,etaold,'.');
%        xlim([0,2*pi])
%        xlabel('\theta');ylabel('p=\Delta \gamma/\sigma_\gamma')
%        subplot(2,1,2)
%        set(gca,'FontSize',16)
%        plot(thetanew,etanew,'r.');
%        xlim([0,2*pi])
%        xlabel('\theta');ylabel('p=\Delta \gamma/\sigma_\gamma')
%        disp(sprintf(['Bunching Factor = ',num2str(abs(mean(exp(1i.*thetanew))))]))
%        legend(sprintf(['A=',num2str(Abh)]));legend boxoff
%% Look at the double buncher scalings       

% Abh=linspace(5,50);
% 
% A2=0.614+0.601.*Abh-0.003.*Abh.^2;
% 
% A1=A2.^2./(0.028.*A2.^2+3.811.*A2-6.309);
% 
% B1=(0.086.*A2.^2+11.974.*A2-19.822)./(1.056.*A2.^2+5.436.*A2);
% 
% B2=(pi.*A2)./(1.965.*A2.^2+9.486);
% %
% for i=1:length(Abh)
%     thnew=thetaold+(B1(i)+B2(i))*etaold+(B1(i)*A1(i)+B2(i)*A2(i))*sin(thetaold)+...
%            B2(i)*A2(i)*sin(thetaold+B1(i)*etaold+B1(i)*A1(i)*sin(thetaold));
%        bunching(i)=abs(mean(exp(1i.*thnew)));
% end
% figure
% subplot(2,2,1)
% plot(A2,A2./A1);xlabel('A2');ylabel('A2/A1')
% xlim([A2(1),A2(end)])
% enhance_plot
% subplot(2,2,2)
% plot(A2,pi./(A1.*B1));xlabel('A2');ylabel('\pi/A1B1')
% xlim([A2(1),A2(end)])
% enhance_plot
% subplot(2,2,3)
% plot(A2,pi./(A2.*B2));xlabel('A2');ylabel('\pi/A2B2')
% xlim([A2(1),A2(end)])
% enhance_plot
% subplot(2,2,4)
% plot(A2,bunching);xlabel('A2');ylabel('BUNCHING')
% xlim([A2(1),A2(end)])
% enhance_plot
% 
