function [newphasespace,newevalue]=push_FEL_particles_RK4p(phasespace,evalue,param,kvalue,fbess1,islice)
   
gammar_sq=param.lambdau/(2*param.lambda0)*(1+kvalue^2/2);
khat = kvalue * fbess1;
% RK-4 for the particles

k1theta=param.stepsize*(param.ku*(1-(gammar_sq./phasespace(:,2).^2)));
k1gamma=param.stepsize*(param.chi2*(khat./phasespace(:,2)).*...
    real(evalue*exp(1i*phasespace(:,1))));

k2theta=param.stepsize*(param.ku*(1-(gammar_sq./(phasespace(:,2)+0.5*k1gamma).^2)));
k2gamma=param.stepsize*(param.chi2*(khat./(phasespace(:,2)+0.5*k1gamma)).*...
    real(evalue*exp(1i*(phasespace(:,1)+0.5*k1theta))));

k3theta=param.stepsize*(param.ku*(1-(gammar_sq./(phasespace(:,2)+0.5*k2gamma).^2)));
k3gamma=param.stepsize*(param.chi2*(khat./(phasespace(:,2)+0.5*k2gamma)).*...
    real(evalue*exp(1i*(phasespace(:,1)+0.5*k2theta))));

k4theta=param.stepsize*(param.ku*(1-(gammar_sq./(phasespace(:,2)+k3gamma).^2)));
k4gamma=param.stepsize*(param.chi2*(khat./(phasespace(:,2)+k3gamma)).*...
    real(evalue*exp(1i*(phasespace(:,1)+k3theta))));

newphasespace(:,1)=phasespace(:,1)+1/6*(k1theta+2*k2theta+2*k3theta+k4theta);
newphasespace(:,2)=phasespace(:,2)+1/6*(k1gamma+2*k2gamma+2*k3gamma+k4gamma);

% Add LSC energy loss along bunch - this adds zero if LSC is disabled
newphasespace(:,2)=newphasespace(:,2)+param.LSCEloss(islice); 

% Predictor-Corrector method for the field

f1=param.chi1(islice)*khat*...
     mean(exp(-1i*phasespace(:,1))./phasespace(:,2));
f1star=param.chi1(islice)*khat*...
     mean(exp(-1i*newphasespace(:,1))./newphasespace(:,2));

newevalue=evalue-param.stepsize/2.*(f1+f1star);
 
