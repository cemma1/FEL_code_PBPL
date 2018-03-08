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
            if kmrtaper
            Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*abs(mean(radfield(ij,param.Nsnap:param.nslices),2)).*sin(res_phase(ij));
            end
            % Compute res phase at next step to preserve bucket area
            if constareataper
                if ij>100 % Choose the start point of the constant area taper (in units of integration steps)
                psvals = linspace(eps,pi/2-eps,100);                
                for nn=1:length(psvals)
                    Kzguess(nn)=Kz(ij)-param.stepsize/const_resp*mean(abs(radfield(ij,param.Nsnap:param.nslices)),2).*sin(psvals(nn));                                          
                    alpha2 = (1-sin(psvals(nn)))/(1+sin(psvals(nn)));                    
                    Area2(nn) = alpha2.*sqrt((mean(abs(radfield(ij+1,:)).*Kzguess(nn))));                                        
                    areaconst(nn) = (abs((Area2(nn)-1*Area1)))/Area1;% For constant changes to bucket area                     
                    %areaconst(nn)=(abs((Area2(nn)-(1-(ij-z0)/(10*param.Nsnap))*Area1)))/Area1; % For variable changes to bucket area z0 is the start point of const area taper
                end

                [mina,ind]=min(abs(areaconst));
                res_phase(ij+1)=psvals(ind);% Sets the psir value to that which keeps the area constant                
                Kz(ij+1)=Kzguess(ind);% Sets the K value to that which keeps the area constant
                    %{ 
                figure% Some diagnostics for the constant-area code in case it breaks for some reason                
                subplot(1,2,1)
                plot(psvals,abs(areaconst));xlabel('\psi_r');ylabel('\delta A/A');                
                subplot(1,2,2)
                plot(psvals,Kzguess);
                disp(res_phase(ij))
                disp(psvals(ind))
                disp(Kz(ij))
                disp(Kzguess(ind))
                disp(Area2(ind))
                disp(Area1)
                ac(ij) = areaconst(ind);
                               pause
                %} 
                                
                else
                Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*mean(abs(radfield(ij,param.Nsnap:param.nslices)),2).*sin(res_phase(ij));  
                alpha1 = (1-sin(res_phase(ij)))/(1+sin(res_phase(ij)));
                Area1 = alpha1*sqrt((mean(abs(radfield(ij,:)).*Kz(ij))));
                end
            end
            %Linear taper if you want
            %{
        if lineartaper % Here add the increasing psir taper
            if ij>34
                res_phase(ij+1)=res_phase(ij)+psirgradient*param.stepsize;
                Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*mean(abs(radfield(ij,param.Nsnap:param.nslices)),2).*sin(res_phase(ij));  
            else
                Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*abs(mean(radfield(ij,param.Nsnap:param.nslices),2)).*sin(res_phase(ij));
            end
        end
        %}
        
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
     %       Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*mean(abs(radfield(ij,:)),2).*sin(res_phase(ij));  
     
     % Compute res phase at next step to preserve bucket area  
            if constareataper
                if ij>50 % Choose the start point of the constant area taper (in units of integration steps)
                psvals = linspace(eps,pi/2-eps,100);                
                for nn=1:length(psvals)
                    Kzguess(nn)=Kz(ij)-param.stepsize/const_resp*mean(abs(radfield(ij,:)),2).*sin(psvals(nn));                                          
                    alpha2 = (1-sin(psvals(nn)))/(1+sin(psvals(nn)));                    
                    Area2(nn) = alpha2.*sqrt((mean(abs(radfield(ij+1,:)).*Kzguess(nn))));                                        
                    areaconst(nn) = (abs((Area2(nn)-1*Area1)))/Area1;% For constant changes to bucket area                     
                    %areaconst(nn)=(abs((Area2(nn)-(1-(ij-z0)/(10*param.Nsnap))*Area1)))/Area1; % For variable changes to bucket area z0 is the start point of const area taper
                end

                [mina,ind]=min(abs(areaconst));
                res_phase(ij+1)=psvals(ind);% Sets the psir value to that which keeps the area constant                
                Kz(ij+1)=Kzguess(ind);% Sets the K value to that which keeps the area constant
                    %{ 
                figure% Some diagnostics for the constant-area code in case it breaks for some reason                
                subplot(1,2,1)
                plot(psvals,abs(areaconst));xlabel('\psi_r');ylabel('\delta A/A');                
                subplot(1,2,2)
                plot(psvals,Kzguess);
                disp(res_phase(ij))
                disp(psvals(ind))
                disp(Kz(ij))
                disp(Kzguess(ind))
                disp(Area2(ind))
                disp(Area1)
                ac(ij) = areaconst(ind);
                               pause
                %} 
                                
                else
                Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*mean(abs(radfield(ij,:)),2).*sin(res_phase(ij));  
                alpha1 = (1-sin(res_phase(ij)))/(1+sin(res_phase(ij)));
                Area1 = alpha1*sqrt((mean(abs(radfield(ij,:)).*Kz(ij))));
                end
            end
            if lineartaper % Here add the increasing psir taper
                res_phase(ij+1)=res_phase(ij)+psirgradient*param.stepsize;
                Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*mean(abs(radfield(ij,:)),2).*sin(res_phase(ij));  
            end
            
            if kmrtaper
                Kz(ij+1)=Kz(ij)-param.stepsize/const_resp*mean(abs(radfield(ij,:)),2).*sin(res_phase(ij));  
            end
            
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
