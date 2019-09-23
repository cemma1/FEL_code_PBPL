% This version is for a KMR style 1-D taper where you compute the undulator
% K step-by-step by specifying a constant resonant phase a priori
total_simtime = 0;
% Time dependent case
if(param.itdp)
for ij = 1:param.Nsnap-1  % takes Nsnap snapshots along length of undulator
    tstart = tic;
    firstslice=ij;
    
    for islice = firstslice:param.nslices       
            % RK4th order integration     
     [phasespacenew,radfield(ij+1,islice)]=...
         push_FEL_particles_RK4([squeeze(thetap(ij,islice,:)),squeeze(gammap(ij,islice,:))],radfield(ij,islice),param,Kz(ij));       
     thetap(ij+1,islice,:) = phasespacenew(:,1);
     gammap(ij+1,islice,:) = phasespacenew(:,2);     
    end
            % Slippage of the radiation field          
            radfield(ij+1,:) = circshift(radfield(ij+1,:).',1);
            radfield(ij+1,1:(ij-1)) = 0; % Set field which slips into the window to zero
            
            % Compute bunching 
            bunch(ij)=mean(abs(mean(exp(1j.*thetap(ij,:,:)),3)));
            
            % Compute undulator field at next step (constant res phase)
            Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*abs(mean(radfield(ij,param.Nsnap:param.nslices),2)).*sin(res_phase(ij));
            
            % Print time taken                
    formatSpec = '%.3f sec from z = %.3f to z = %.3f, total length = %.3f \n';
    fprintf(formatSpec, toc(tstart), param.stepsize*(ij-1), param.stepsize*ij, param.Nsnap*param.stepsize);
end
    else
        % Time independent case
        for ij = 1:param.Nsnap-1  % takes Nsnap snapshots along length of undulator
    % RK4th order integration                     
     [phasespacenew,radfield(ij+1,1)]=...
         push_FEL_particles_RK4([squeeze(thetap(ij,1,:)),squeeze(gammap(ij,1,:))],radfield(ij,1),param,Kz(ij));       
     thetap(ij+1,1,:) = phasespacenew(:,1);
     gammap(ij+1,1,:) = phasespacenew(:,2);
              
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
