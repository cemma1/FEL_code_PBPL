close all
% Plot saved simulation output
IA = 17045;                                                       % Alfven current
%% Import stuff
path=['Simulation_output/']

param=importdata([path,'simulation_parameters.mat']);
if param.itdp
time=importdata([path,'output_time.mat']);
spectrum=importdata([path,'output_spectrum.mat']);
outputpower=importdata([path,'output_power.mat']);
omega=importdata([path,'output_frequency.mat'])-5e-4;
end
%Time independent quantities
avpower=importdata([path,'averagepower.mat']);
%meanenergy=importdata([path,'average_energy.mat']);
avbunching=importdata([path,'average_bunching.mat']);
%% Calculate stuff
rho1D = 1/param.gamma0*(1/8*param.I/IA*param.K.^2./param.sigmax^2/param.ku^2)^(1/3);
Lgain = param.lambdau/(4*sqrt(3)*pi*rho1D);
%% Plot stuff
figure(1)
subplot(2,1,1)
if param.itdp
semilogy(omega,abs(spectrum))
xlabel('\delta\omega/\omega ','FontSize',16)
    ylabel('P (\omega) [arb. units]','FontSize',16)    
    xlim([-20,20].*1e-3)    
    set(gca,'FontSize',16) 
    enhance_plot
    subplot(2,1,2)
plot([1:1:size(outputpower,2)]*param.zsep*param.lambda0*1e15/3e8,outputpower)
xlim([1,size(outputpower,2)]*param.zsep*param.lambda0*1e15/3e8)
xlabel('t [fs]','FontSize',16)
ylabel('Output Radiation Power [W]','FontSize',16)
set(gca,'FontSize',16)
enhance_plot
end
figure(2)
subplot(1,3,1)
plot([1:param.Nsnap]*param.stepsize/Lgain,avpower/param.Ee/param.I)
xlim([0,param.Nsnap*param.stepsize]/Lgain)
hold on
xlabel('z/L_g')
ylabel('P/P_{beam}')
enhance_plot
    subplot(1,3,2)
plot([1:1:param.Nsnap-1]*param.stepsize/Lgain,avbunching)
xlim([0,param.Nsnap*param.stepsize]/Lgain)
xlabel('z/L_g')
ylabel('Bunching Factor')
enhance_plot
% subplot(1,3,3)
% plot([1:1:param.Nsnap-1]*param.stepsize/Lgain,(meanenergy/meanenergy(1)-1)/rho1D)
% xlabel('z/L_g')
% ylabel('\Delta\gamma/\rho\gamma')
% xlim([0,param.Nsnap*param.stepsize]/Lgain)
% enhance_plot
%% Comparison plot
path=['Simulation_output/']

param=importdata([path,'simulation_parameters.mat']);
time=importdata([path,'output_time.mat']);
spectrum=importdata([path,'output_spectrum.mat']);
outputpower=importdata([path,'output_power.mat']);
omega=importdata([path,'output_frequency.mat'])-5e-4;
avpower=importdata([path,'averagepower.mat']);

avbunching=importdata([path,'average_bunching.mat']);

figure(1)
subplot(2,1,1)
hold on
semilogy(omega,abs(spectrum),'r')
xlabel('\delta\omega/\omega ','FontSize',16)
    ylabel('P (\omega) [arb. units]','FontSize',16)            
    set(gca,'FontSize',16) 
    enhance_plot
    subplot(2,1,2)
    hold on
semilogy([1:1:size(outputpower,2)]*param.zsep*param.lambda0*1e15/3e8,outputpower,'r')
xlim([1,size(outputpower,2)]*param.zsep*param.lambda0*1e15/3e8)
xlabel('t [fs]','FontSize',16)
ylabel('Output Radiation Power [W]','FontSize',16)
set(gca,'FontSize',16)
enhance_plot


figure(2)
subplot(1,3,1)
plot([1:param.Nsnap]*param.stepsize/Lgain,avpower/param.Ee/param.I,'r*')
xlim([0,param.Nsnap*param.stepsize]/Lgain)
hold on
xlabel('z/L_g')
ylabel('P/P_{beam}')
enhance_plot
    subplot(1,3,2)
    hold on
plot([1:1:param.Nsnap-1]*param.stepsize/Lgain,avbunching,'r*')
xlim([0,param.Nsnap*param.stepsize]/Lgain)
xlabel('z/L_g')
ylabel('Bunching Factor')
enhance_plot
% subplot(1,3,3)
% hold on
% plot([1:1:param.Nsnap-1]*param.stepsize/Lgain,(meanenergy/meanenergy(1)-1)/rho1D,'r*')
% xlabel('z/L_g')
% ylabel('\Delta\gamma/\rho\gamma')
% xlim([0,param.Nsnap*param.stepsize]/Lgain)
%enhance_plot