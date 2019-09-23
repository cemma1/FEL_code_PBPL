%perave_postprocessor
close all
set(groot, 'defaultTextInterpreter','tex')
kw=2*pi/param.lambdau;
hbar=6.582e-16;
%% Field z vs s
figure;imagesc(abs(radfield));set(gca,'YDir','normal')
xlabel('s');ylabel('z');h1 = colorbar
title(h1,'|E|');enhance_plot('FontSize',20)

%% Spectrum as a function of z
zlocations=linspace(param.stepsize,lwig,30);fundpower=[];sidebandpower=[];
zidx=round(zlocations/param.stepsize);
if param.itdp
    dt = param.zsep*param.lambda0/c;
    tposition = [1:size(power,2)]*dt*1e15;
omegamin=-7e-3;omegamax=7e-3; % For sideband filtering 
h=figure(1);
%set(h, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

for n=1:length(zidx)

[powerspec,omega]=spectrum_calc(radfield(zidx(n),:),param.lambda0,param.zsep);
sidebandindex=omega>omegamin & omega<omegamax;
fundspectrum=powerspec(sidebandindex);
fundpower(n)=trapz(fundspectrum)/trapz(powerspec);
figure(1)
subplot(1,2,1)
%plot((omega+1)*hbar*2*pi*c/param.lambda0,powerspec)%Energy spectrum
semilogy(omega,abs(powerspec))
%plot(omega,abs(powerspec))
xlabel('\delta\omega/\omega ','FontSize',16)
    ylabel('P (\omega) [arb. units]','FontSize',16)    
    xlim([-100,100].*rho1D)    
    set(gca,'FontSize',16)
    legend(sprintf(['z / L_u =',num2str(zlocations(n)/lwig)]));
    
subplot(1,2,2)
tailslice = param.Nsnap+1;
if param.currprofile
    yyaxis left
    plot(tposition,power(zidx(n),:))
    ylabel('Output Radiation Power [W]','FontSize',16)
    yyaxis right
    plot(tposition,param.Iprofile(tailslice:end).*1e-3)
    ylabel('Current [kA]','FontSize',16)
else
 plot(tposition,power(zidx(n),:))   
 ylabel('Output Radiation Power [W]','FontSize',16)
end    
xlim([0,tposition(end)])
xlabel('t [fs]','FontSize',16)
set(gca,'FontSize',16)
drawnow

end
end
%% Radiation Power and spectrum at exit
%power3(:,:) = abs(radfield_3rd(:,:)).^2/377*param.A_e;
zpos= [1:param.Nsnap]*param.stepsize;
zoverlg=zpos./Lgain;
% if ~param.itdp
% Lgfit=fit_gainlength(Lgain,zpos,mean(power,2));
% end
figure(2);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
title('Simulation Output')
subplot(2,3,1)
%semilogy(zpos./Lgain,mean(power,2)/param.Ee/param.I)
plot(zpos./Lgain,mean(power,2)/param.Ee/param.I)
hold on
% if ~param.itdp
% Lgfit=gainlength_fit(Lgain,zpos,mean(power,2));
% end
xlim([0,zpos(end)]/Lgain)
%semilogy(ij*param.stepsize,mean(power3,2),'r')
xlabel('z/L_g')
ylabel('P/P_{beam}')
legend(['\eta_{max}=',num2str(max(mean(power,2)/param.Ee/param.I*100),'%.2f'),'%'],'location','NorthWest'); legend boxoff
enhance_plot
if param.itdp
subplot(2,3,2)
plot(tposition,power(end,:))
xlim([0,tposition(end)])
xlabel('t [fs]')
ylabel('Power [W]')
tcoh=sqrt(10*pi)/3/rho1D*param.lambda0/2/pi/3e8;
legend(sprintf(['t_c ~ ',num2str(tcoh*1e15,'%.2f'),' fs']))
enhance_plot
[powerspec,omega]=spectrum_calc(radfield(end,:),param.lambda0,param.zsep);
subplot(2,3,3)
%plot((omega+1)*hbar*2*pi*c/param.lambda0,powerspec,'b');    
semilogy((omega+1)*hbar*2*pi*c/param.lambda0,powerspec,'b');    
xlim([omega(1)+1,omega(end)+1]*hbar*2*pi*c/param.lambda0)    
    xlabel('Photon Energy [eV]')
    ylabel('P (\omega) [arb. units]')
    enhance_plot
    title('Output Spectrum')

% Plot undulator field again
%     subplot(2,3,6)
%     plot([1:1:param.Nsnap]*param.stepsize,Kz)
% xlim([1,param.Nsnap]*param.stepsize)
% xlabel('z [m]')
% ylabel('Undulator K')
% if param.tapering
% legend(sprintf(['a_w(z)=a_0*(1-c*(z-z_0)^d)\n z_0=',num2str(param.z0),'\n dK/K = ',...
%     num2str(param.ratio*(lwig-param.z0)^param.order),'\n d=',num2str(param.order)]),'location','SouthWest')
% ylim([(1-(lwig-param.z0)^param.order*param.ratio),1]*param.K)
% end
% enhance_plot
end
%% Bunching factor and energy loss (takes a while to calculate)

figure(2)
subplot(2,3,4)
plot([1:param.Nsnap-1]*param.stepsize/Lgain,bunch)
xlim([0,param.Nsnap*param.stepsize]/Lgain)
xlabel('z/L_g')
ylabel('Bunching Factor')
enhance_plot
for ij=1:param.Nsnap
meanenergy(ij)=mean(abs(mean(gammap(ij,:,:),3)));    
end
subplot(2,3,5)
plot([1:1:param.Nsnap]*param.stepsize/Lgain,(meanenergy/meanenergy(1)-1))
xlabel('z/L_g')
ylabel('\Delta\gamma/\gamma_0')
xlim([0,param.Nsnap*param.stepsize]/Lgain)
enhance_plot
    
% Sideband Stuff (for seeded FEL)
if param.itdp
figure(22)
semilogy(zlocations,fundpower,'r')
hold on
semilogy(zlocations,1-fundpower,'b')
legend('Fundamental Power (Fractional)', 'Sideband Power (Fractional)','location','NorthWest')
ylim([0,1]);
set(gca,'FontSize',14)
xlabel('z [m]','FontSize',16)
ylabel('P/P_{total}','FontSize',16)    
end

%% Tapering plots 
if param.tapering
meanfield=mean(radfield(:,:),2);
dKdz=zeros(param.Nsnap,1);
dKdz(2:end)=abs(diff(Kz)./param.delz./param.lambdau);
psir=asin(const_resp.*(dKdz./abs(meanfield)));
psirend=asin(const_resp.*(dKdz(end)./abs(radfield(end,:))));
gammarofz=sqrt((param.lambdau/2/param.lambda0).*(1+Kz.^2));
dgammardz=-(param.lambdau/2/param.lambda0).*Kz./gammarofz/2.*dKdz';

options = optimset('Display','off');
for i=1:length(psir)
  if psir(i)>0
    psr=psir(i);

  psi2(i) = pi-psr;
  psi1(i) = fsolve(@(x) cos(x)+x*sin(psr)-cos(psi2(i))-psi2(i)*sin(psr),-pi,options);
 dX=[psi1(i):pi/50:psi2(i)];
 Y=cos(psr)+cos(dX)-(pi-psr-dX)*sin(psr);
 %alpha(i)=sqrt(2)/8*trapz(dX,Y);
 alpha(i)=(1-sin(psr))/(1+sin(psr));%this is an approx. by S.Y. Lee 
  else
      psi2(i)=pi;
      psi1(i)=-pi;
  end
end
figure(2)
    subplot(2,3,6)
       %semilogy([1:1:param.Nsnap]*param.stepsize,sqrt(cos(psir).*mean(abs(radfield),2)*e0.*param.K.*param.ku./(m_e.*param.gamma0.^2))./c) 
       ksynch=param.ku*sqrt(cos(psir).*mean(param.chi2/param.k*abs(radfield),2)*Kz/(1+Kz.^2));
       semilogy(zpos./Lgain,2*pi./ksynch) 
   xlim([1,param.Nsnap]*param.stepsize./Lgain)
   xlabel('z/L_g','FontSize',12)
   ylabel('\lambda_{synch} [m] ','FontSize',12)
   enhance_plot
%% Tapering Plots cont'd

figure(3);
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Tapering Plots', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',30,'Fontname','Times')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
c1 = -0.685;
c2 = -20.462;
c3 = 0.0172;
c4 = -0.128;
chengyingfit = pi/2+c1+c2*atan(c3*exp(c4*2*param.ku*rho1D.*zpos));
subplot(2,3,1)
plot(zoverlg,psir*180/pi)
hold on
plot(zoverlg,chengyingfit.*180/pi,'r--')
%legend('Simulation','CY fit'); legend boxoff
%plot(zoverlg,res_phase*180/pi,'r--')% To compare the computed resonant phase with the real one fed to the code, they should be the same
xlabel('z/L_g');ylabel('\Psi_R [degree]');enhance_plot;xlim([0,zoverlg(end)])
subplot(2,3,2)
plot(zoverlg,abs(meanfield)*1e-12)
xlabel('z/L_g');ylabel('<E_{rad}> [TV/m]');enhance_plot;xlim([0,zoverlg(end)])
subplot(2,3,3)
plot(zoverlg,dgammardz*0.511)
xlabel('z/L_g');ylabel('d\gamma_r/dz [MeV/m]');enhance_plot;xlim([0,zoverlg(end)])
%subplot(2,3,5)
%plot(zoverlg,(psi2-psi1)./2/pi)
%hold on
%plot(zoverlg,alpha,'r')
%ylim([0,1])
%xlim([0,zoverlg(end)])
%enhance_plot
%xlabel('z/L_g')
%legend('f_b','\alpha')
if param.itdp
subplot(2,3,4)
plot(tposition,psirend*180/pi)
xlim([0,tposition(end)])
xlabel('t [fs]');ylabel('\Psi_R [degree]');enhance_plot;
area = sqrt(abs(meanfield).*Kz').*(1-sin(psir))./(1+sin(psir));
subplot(2,3,5)
plot(zoverlg,area./area(2))
xlim([0,zoverlg(end)])
xlabel('z/L_g');ylabel('Bucket Area [arb. units]');enhance_plot; 
else
    area = sqrt(abs(meanfield).*Kz').*(1-sin(psir))./(1+sin(psir));
subplot(2,3,4)
plot(zoverlg,area./area(2))
xlim([0,zoverlg(end)])
xlabel('z/L_g');ylabel('Bucket Area [arb. units]');enhance_plot;    
end
subplot(2,3,6)
plot(zoverlg,Kz)
xlabel('z/L_g');ylabel('Undulator K (rms)');enhance_plot;xlim([0,zoverlg(end)])
end
%% Optical Guiding plots
%{
realpartn = param.chi1(1).*Kz./gammarofz./param.k./abs(meanfield').*cos(res_phase);
imagpartn = param.chi1(1).*Kz./gammarofz/param.k./abs(meanfield').*sin(res_phase);
n = param.chi1(1).*Kz./gammarofz./param.k./abs(meanfield').*exp(1i.*res_phase);

Vsquared = param.k^2*param.sigmax^2.*(n.^2-1);

figure
subplot(2,1,1)
semilogy(zoverlg,abs(n))
hold on
semilogy(zoverlg,realpartn,'r')
hold on
semilogy(zoverlg,imagpartn,'k')
legend('abs(n)-1','Re(n)-1','Im(n)');
subplot(2,1,2)
semilogy(zoverlg,abs(Vsquared))
%}
% Calculate the refractive index as a function of z and s
for n=1:size(radfield,2)% Loop over each slice
    % Calculate the refractive index as a function of z and s
fieldangle(:,n)=unwrap(angle(radfield(:,n)));
Re_n(:,n)= 1+gradient(fieldangle(:,n))/param.k;
Im_n(:,n) = gradient(abs(radfield(:,n)))./abs(radfield(:,n))/param.k; 
%{
meanphase=mean(mean(squeeze(thetap(n,:,:)),2)+fieldangle(n,:),2)+pi/2;
meanenergy=mean(squeeze(gammap(n,:,:)),2); 

% Average energy and phase across arrays - assumes const current, otherwise
% change param.chi1 to match the value for each slice
Re_n(n,:)=1+param.chi1(1)./param.k*Kz(n)./abs(radfield(n,:)).*(cos(meanphase)./meanenergy)';
Im_n(n,:)=param.chi1(1)./param.k*Kz(n)./abs(radfield(n,:)).*(sin(meanphase)./meanenergy)';
%}
end
figure;
subplot(1,2,1);surface(Re_n);shading interp;subplot(1,2,2);surface(Im_n);shading interp
%% eSASE plots
zlocations=linspace(param.stepsize,lwig,30);
zidx=round(zlocations/param.stepsize);
for n=1:length(zidx)
    
    tp=squeeze(thetap(zidx(n),:,:))+pi/2;% You can add the phase of the field if you want
    gp=squeeze(gammap(zidx(n),:,:))+pi/2;% You can add the phase of the field if you want
    
    slice_bunching = abs(mean(exp(1j.*tp),2));
    slice_energy = mean(gp,2);
    
    tposition = [1:size(power(zidx(n),:),2)]*param.zsep*param.lambda0/c;
    tailslice = param.Nsnap+1;    
            
    if param.currprofile
        filename = 'eSASE_movie'
    figure(23) 
    subplot(3,1,1)
        yyaxis left
        plot(tposition/tcoh,power(zidx(n),:)/rho1D/param.I/param.Ee)
        ylabel('P/\rho P_{beam}','FontSize',16)
        yyaxis right
        plot(tposition/tcoh,param.Iprofile(tailslice:end).*1e-3)
        ylabel('Current [kA]','FontSize',16)
        
    subplot(3,1,2)
        yyaxis left
        plot(tposition/tcoh,slice_bunching)% You need to find where B(s) is stored...
        ylabel('Bunching Factor','FontSize',16)
        yyaxis right
        plot(tposition,param.Iprofile(tailslice:end).*1e-3)
        ylabel('Current [kA]','FontSize',16)
    
    subplot(3,1,3)
        yyaxis left
        plot(tposition/tcoh,(slice_energy-param.gamma0)/rho1D/param.gamma0)% You need to find where B(s) is stored...
        ylabel('\Delta \gamma/\rho\gamma_0','FontSize',16)
        yyaxis right
        plot(tposition,param.Iprofile(tailslice:end).*1e-3)
        ylabel('Current [kA]','FontSize',16)
        xlabel('t/t_c','FontSize',16)
        
              frame = getframe(23);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
    end
end
%% Phasespace movie
zlocations=linspace(param.stepsize,lwig,50);
zidx=round(zlocations/param.stepsize);
if param.tapering
bucketheight = sqrt(param.chi2/param.ku.*Kz./gammarofz.^2.*abs(mean(radfield(:,:),2))');
else
    bucketheight = sqrt(param.chi2/param.ku.*Kz'./param.gamma0.^2.*abs(mean(radfield(:,:),2)));
end

if param.phasespacemovie
    filename='particle_movie.gif';
    figure(10);
    
for i=1:length(zidx)    
    x=linspace(-pi,pi,1e4);
    if param.tapering
        %sepa=bucketheight(zindices(i))'.*separatrix(linspace(-pi,pi,1e4),res_phase(zindices(i)))/rho1D+(gammarofz(zindices(i))-param.gamma0)/param.gamma0/rho1D;
        %sepamin=(gammarofz(zindices(i))-param.gamma0)/param.gamma0-bucketheight(zindices(i))'.*separatrix(linspace(-pi,pi,1e4),res_phase(zindices(i)));
        sepa=bucketheight(zidx(i))'.*separatrix(x,res_phase(zidx(i)))+(gammarofz(zidx(i))-param.gamma0)/param.gamma0;
        sepamin=(gammarofz(zidx(i))-param.gamma0)/param.gamma0-bucketheight(zidx(i))'.*separatrix(x,res_phase(zidx(i)));
    else
    sepa=bucketheight(zidx(i))'.*separatrix(x,res_phase(zidx(i)));
    sepamin=-bucketheight(zidx(i))'.*separatrix(x,res_phase(zidx(i)));
    gammarofz=param.gamma0.*ones(length(zidx));
    end
fieldphase(i)=mean(angle(radfield(zidx(i),:)));
tp=squeeze(thetap(zidx(i),:,:))+fieldphase(i)+pi/2;% You can add the phase of the field if you want
gp=squeeze(gammap(zidx(i),:,:));    
tresh=reshape(tp,[1,size(tp,1)*size(tp,2)]);gresh=reshape(gp,[1,size(gp,1)*size(gp,2)]);
%tresh=squeeze(thetap(zindices(i),param.Nsnap/2,:));gresh=squeeze(gammap(zindices(i),param.Nsnap/2,:));
%CALCULATE THE TRAPPING FRACTION
%[psi1,psi2]=bucket_parameters(param.psir);% This is for const res phase
[psi1,psi2]=bucket_parameters(res_phase(zidx(i)));% This is for the general case
indi=(mod(tp,2*pi)-pi)>psi2 & (mod(tp,2*pi)-pi)<psi1;
g2=(gp(indi)./(meanenergy(1))-1);
%g2=(gp(indi)./(param.gamma0)-1);% This is ~ the same as the above line
angoli=(mod(tp(indi),2*pi)-pi);
indi2=g2<(bucketheight(zidx(i))'.*separatrix(angoli,param.psir)+(gammarofz(zidx(i))-param.gamma0)/param.gamma0);
ftrap(i)=numel(g2(indi2))/numel(gp);

ind=x<psi1 & x>psi2;

if param.itdp
plot((mod(tresh,2*pi)-pi)./pi,(gresh./(meanenergy(1))-1)*100,'*k','MarkerSize',1)
else     
    subplot(1,2,1)
    semilogy([1:zidx(i)]*param.stepsize/Lgain,power(1:zidx(i))/param.Ee/param.I,'k')    
    xlim([0,lwig]./Lgain);ylim([min(power),2*max(power)]./param.Ee/param.I)
    %plot([1:zindices(i)]*param.stepsize/Lgain,power(1:zindices(i))/param.Ee/param.I,'k')
    %xlim([0,lwig]./Lgain);ylim([min(power),1.1*max(power)]./param.Ee/param.I)
    xlabel('z/L_g');ylabel('P/P_{beam}')
    set(gca,'YTick',logspace(-4,-1,4))
    enhance_plot('Times',20)
    hold on
    legend(['z/L_u = ',num2str(zlocations(i)/lwig,'%.2f')],'location','SouthEast')
    legend boxoff
    subplot(1,2,2)
    % Without the separatrix
    %plot((mod(tp,2*pi)-pi)./pi,(gp./(meanenergy(1))-1),'.k','MarkerSize',1)        
    plot((mod(tp,2*pi)-pi)./pi,(gp./(meanenergy(1))-1)*100,'.k',x(ind)/pi,sepa(ind)*100,'r',x(ind)/pi,sepamin(ind)*100,'r','LineWidth',2)
    set(gca,'FontSize',20)
    xlabel('\Psi/pi');ylabel('\Delta \gamma/\gamma_0');%enhance_plot('FontSize',16)
    drawnow
end
xlim([-1,1])
set(gca,'FontSize',20)
xlabel('\Psi/pi');ylabel('\Delta \gamma/\gamma_0');%enhance_plot('FontSize',16)
% If you want to save to GIF

%       drawnow
%       frame = getframe(10);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if i == 1;
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%       end

end
figure(2)
subplot(2,3,5)
plot(zlocations./param.lambdau,ftrap)
xlabel('z/\lambda_u');ylabel('f_t [Calculated]');
%legend(sprintf(['f_t [Theory]=',num2str((psi1-psi2)/2/pi)]))
enhance_plot
end
%%
if param.tapering
figure(50)
plot(zoverlg,res_phase*180/pi)
hold on
plot(zoverlg,chengyingfit.*180/pi,'r--');xlabel('z/L_g');ylabel('\Psi_R [degree]');enhance_plot;xlim([0,zoverlg(end)])
end
%% Functions to calculate bucket parameters and fit the gain length

function [ sep ] = separatrix( psi,psir )
% Calculate height of separatrix at phase psi for resonant phase psir
% Detailed explanation goes here
sep = sqrt(abs(cos(psi)+cos(psir)+(psi+psir-pi*sign(psir))*sin(psir)));
end

function [Lgfit] = fit_gainlength(Lgaintheory,z1,pow);

zmin=5*Lgaintheory;zmax=zmin+10*Lgaintheory;
pgain=pow(z1>zmin & z1<zmax);
zgain=z1(z1>zmin & z1<zmax);

Lgfit=(zgain(end)-zgain(1))/(log(pgain(end)/pgain(1)));
disp(sprintf(['Lg theory = ' num2str(Lgaintheory)]));

% figure
% semilogy(linspace(zmin,zmax,length(pgain)),pgain)
% hold on
% semilogy(linspace(zmin,zmax),pgain(1)*exp((linspace(zmin,zmax)-zmin)./Lgfit),'r--')

end

function [ psi1, psi2, bucket_height, capture_fraction, bucket_area, bunching, sinaverage] = bucket_parameters(psir)
% Calculate various bucket parameters as a function of psir
 
psi1 = pi - psir;
pond = @(psi,psir) cos(psi1)+psi1*sin(psir)-(cos(psi)+psi*sin(psir));
psi2 = fsolve(@(psi) pond(psi,psir),[-pi:psi1*0.99]);
psi2 = psi2(1);
bucket_height = sqrt(cos(psir)-(pi/2-psir)*sin(psir));
bucket_area = (1-sin(psir)) / (1+sin(psir));
capture_fraction = (psi1 - psi2)/2/pi;
potential = @(psi,psir) cos(psi)+psi*sin(psir);
bucket_sep = @(psi) sqrt( (-potential(psi1,psir)+potential(psi,psir))/2);
bucket_bun = @(psi) bucket_sep(psi).*exp(1i*psi);
bucket_sin = @(psi) bucket_sep(psi).*sin(psi);
bunching = abs(integral(bucket_bun,psi2,psi1) / integral(bucket_sep,psi2,psi1));
sinaverage = abs(integral(bucket_sin,psi2,psi1) / integral(bucket_sep,psi2,psi1));
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

