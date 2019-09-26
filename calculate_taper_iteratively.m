function [Kznew,psirnew] = calculate_taper_iteratively(Kz,param,zstep,const_resp,radfield,res_phase)

firststep = round(param.z0/param.stepsize);

% Time independent case
if ~(param.itdp)
    switch param.tapering    
        case 0
            Kznew = Kz(zstep);
            psirnew = 0.0;
            
        case 1 % KMR-style with constant resonant phase
            Kznew = Kz(zstep)-param.stepsize/const_resp*mean(abs(radfield(zstep,:)),2).*sin(res_phase(zstep));
            psirnew = res_phase(zstep+1);

        case 2 % Linearly increasing resonant phase     
            if zstep>firststep
            psirnew = res_phase(zstep)+param.psirgradient*param.stepsize;
            Kznew = Kz(zstep)-param.stepsize/const_resp*mean(abs(radfield(zstep,:)),2).*sin(res_phase(zstep)); 
            else
            psirnew = res_phase(zstep+1);
            Kznew = Kz(zstep);
            end

        case 3 % Constant bucket area                                  
                   if zstep>firststep % Choose the start point of the constant area taper (in units of integration steps)
                    alpha1 = (1-sin(res_phase(firststep)))/(1+sin(res_phase(firststep)));
                    Area1 = alpha1*sqrt((mean(abs(radfield(firststep,:)).*Kz(firststep))));   
                    psvals = linspace(eps,pi/2-eps,100);                

                    for nn=1:length(psvals)
                        Kzguess(nn)=Kz(zstep)-param.stepsize/const_resp*mean(abs(radfield(zstep,:)),2).*sin(psvals(nn));                                          
                        alpha2 = (1-sin(psvals(nn)))/(1+sin(psvals(nn)));                    
                        Area2(nn) = alpha2.*sqrt((mean(abs(radfield(zstep+1,:)).*Kzguess(nn))));                                        
                        areaconst(nn) = (abs((Area2(nn)-1*Area1)))/Area1;% For constant changes to bucket area                     
                        %areaconst(nn)=(abs((Area2(nn)-(1-(zstep-z0)/(10*param.Nsnap))*Area1)))/Area1; % For variable changes to bucket area z0 is the start point of const area taper
                    end
                    [minareadiff,ind]=min(abs(areaconst));
                    psirnew=psvals(ind);% Sets the psir value to that which keeps the area constant                
                    Kznew=Kzguess(ind);% Sets the K value to that which keeps the area constant

                    else
                    Kznew=Kz(zstep)-param.stepsize/const_resp*mean(abs(radfield(zstep,:)),2).*sin(res_phase(zstep));  
                    alpha1 = (1-sin(res_phase(zstep)))/(1+sin(res_phase(zstep)));
                    Area1 = alpha1*sqrt((mean(abs(radfield(zstep,:)).*Kz(zstep))));
                    psirnew = res_phase(zstep+1);
                   end
                   
        case 4 % Polynomial tapering
            if zstep > firststep
                Kznew = param.K*(1-param.ctaper*(zstep-firststep).^param.dtaper*param.stepsize^param.dtaper);
                psirnew = res_phase(zstep+1);
            else
                Kznew = Kz(zstep);
                psirnew = res_phase(zstep+1);
            end
    end
    
else
% Time dependent case    
    switch param.tapering
        
        case 0 % No taper
            Kznew = Kz(zstep);
            psirnew = 0.0;
            
        case 1 % KMR
            Kznew = Kz(zstep)-param.stepsize/const_resp*mean(abs(radfield(zstep,param.Nsnap:param.nslices)),2).*sin(res_phase(zstep));
            psirnew = res_phase(zstep+1);            
            
        case 2 % Linear Taper
            if zstep>firststep
            psirnew = res_phase(zstep)+param.psirgradient*param.stepsize;
            Kznew=Kz(zstep)-param.stepsize/const_resp*abs(mean(radfield(zstep,param.Nsnap:param.nslices),2)).*sin(res_phase(zstep));
            else                
            psirnew = res_phase(zstep+1);     
            Kznew = Kz(zstep);
            end
            
                            
        case 3 % Constant Bucket Area           
                if zstep>firststep % Choose the start point of the constant area taper (in units of integration steps)
                alpha1 = (1-sin(res_phase(firststep)))/(1+sin(res_phase(firststep)));
                Area1 = alpha1*sqrt((mean(abs(radfield(firststep,:)).*Kz(firststep))));
                psvals = linspace(eps,pi/2-eps,100);                
                
                for nn=1:length(psvals)
                    Kzguess(nn)=Kz(zstep)-param.stepsize/const_resp*mean(abs(radfield(zstep,param.Nsnap:param.nslices)),2).*sin(psvals(nn));                                          
                    alpha2 = (1-sin(psvals(nn)))/(1+sin(psvals(nn)));                    
                    Area2(nn) = alpha2.*sqrt((mean(abs(radfield(zstep+1,:)).*Kzguess(nn))));                                        
                    areaconst(nn) = (abs((Area2(nn)-1*Area1)))/Area1;% For constant changes to bucket area                     
                    %areaconst(nn)=(abs((Area2(nn)-(1-(zstep-z0)/(10*param.Nsnap))*Area1)))/Area1; % For variable changes to bucket area z0 is the start point of const area taper
                end

                [minareadiff,ind]=min(abs(areaconst));
                psirnew=psvals(ind);% Sets the psir value to that which keeps the area constant                
                Kznew=Kzguess(ind);% Sets the K value to that which keeps the area constant
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
                Kznew=Kz(zstep)-param.stepsize/const_resp*mean(abs(radfield(zstep,param.Nsnap:param.nslices)),2).*sin(res_phase(zstep));  
                alpha1 = (1-sin(res_phase(zstep)))/(1+sin(res_phase(zstep)));
                Area1 = alpha1*sqrt((mean(abs(radfield(zstep,:)).*Kz(zstep))));
                psirnew = res_phase(zstep+1);
                end
        
        case 4 % Polynomial tapering
            if zstep > firststep
                Kznew = param.K*(1-param.ctaper*(zstep-firststep).^param.dtaper*param.stepsize^param.dtaper);
                psirnew = res_phase(zstep+1);
            else
                Kznew = Kz(zstep);
                psirnew = res_phase(zstep+1);               
            end
    end
end