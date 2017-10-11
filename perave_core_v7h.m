% This version is for a KMR style 1-D TESSA where you choose the undulator
% K iteratively by specifying the resonant phase
%% Solve system of equations

% define initial conditions
%options = odeset('RelTol', param.accuracy,'OutputFcn',@odeplot,'OutputSel',[1]); %,'Stats','on');   % set solver options
%options = odeset('RelTol', param.accuracy); %,'Stats','on'); 

total_simtime = 0;
%hl = 0;
% Constant for the resonant phase calculation   
const_resp=1/param.chi2*(param.lambdau/2/param.lambda0);

if(param.itdp)
for ij = 1:param.Nsnap-1  % takes Nsnap snapshots along length of undulator
    tstart = tic;
    
firstslice=ij;
    
    for islice = firstslice:param.nslices       
     gammaf = squeeze(gammap(ij,islice,:));
     thetaf = squeeze(thetap(ij,islice,:));
     E_q0 = radfield(ij,islice);
     
% RK4th order integration     
     phasespaceold = [thetaf,gammaf];
     evaluesold = E_q0;

     [phasespacenew,evaluesnew]=push_FEL_particles_RK4(phasespaceold,evaluesold,param,Kz(ij));       

     thetap(ij+1,islice,:) = phasespacenew(:,1);
     gammap(ij+1,islice,:) = phasespacenew(:,2);
     radfield(ij+1,islice) = evaluesnew;
    end
% Slippage of the radiation field          
            %Fundamental
            B = radfield(ij+1,:).';
            B = circshift(B,1);
            radfield(ij+1,:) = B;
            radfield(ij+1,1) = param.E0;
            
            % Compute bunching 
            bunch(ij)=mean(abs(mean(exp(1j.*thetap(ij,:,:)),3)));
            
            % Compute undulator field at next step (constant res phase)
            Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*abs(mean(radfield(ij,param.Nsnap:param.nslices),2)).*sin(res_phase(ij));
            
            telapsed = toc(tstart);
    total_simtime = total_simtime+telapsed;
    
    formatSpec = '%.3f sec from z = %.3f to z = %.3f, total length = %.3f \n';
    fprintf(formatSpec, telapsed, param.stepsize*(ij-1), param.stepsize*ij, param.Nsnap*param.stepsize);
end
    else
        for ij = 1:param.Nsnap-1  % takes Nsnap snapshots along length of undulator
     gammaf = squeeze(gammap(ij,1,:));
     thetaf = squeeze(thetap(ij,1,:));
     E_q0 = radfield(ij,1);   
     
     % RK4th order integration     
     phasespaceold = [thetaf,gammaf];
     evaluesold = E_q0;
     [phasespacenew,evaluesnew]=push_FEL_particles_RK4(phasespaceold,evaluesold,param,Kz(ij));       
     %[phasespacenew,evaluesnew,evaluesnew3rdharm]=push_FEL_particles_3rdharm(phasespaceold,Kz(ij),fbess1(ij),fbess3(ij),evaluesold,0,param);       
     thetap(ij+1,1,:) = phasespacenew(:,1);
     gammap(ij+1,1,:) = phasespacenew(:,2);
     radfield(ij+1,1) = evaluesnew;          
     % Compute bunching and mean energy
    bunch(ij)=mean(abs(mean(exp(1j.*thetap(ij,:,:)),3))); 
    
     % Compute undulator field at next step (constant res phase)     
            Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*mean(abs(radfield(ij,:)),2).*sin(res_phase(ij));           
        end
    end                    

% Remove slices within one total slippage length
 if(param.itdp)
 radfield(:,1:param.Nsnap)=[];
 gammap(:,1:param.Nsnap,:)=[];
 thetap(:,1:param.Nsnap,:)=[];
 end
% Calculate radiation power 
 power(:,:) = abs(radfield(:,:)).^2/377*param.A_e;