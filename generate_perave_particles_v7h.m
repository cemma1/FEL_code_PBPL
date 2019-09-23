%% Initialize current profile if set in user input
if param.currprofile
    dt = param.zsep*param.lambda0/c;
    tvector = [1:param.nslices].*dt-round((param.nslices+param.Nsnap)/2.5).*dt;
    param.Iprofile = param.I.*exp(-tvector.^2/2/param.sigmat^2);    
else
    param.Iprofile = param.I.*ones(1,param.nslices);
end
param.chi1=mu0*c/2.*param.Iprofile./param.A_e; % Simplifying constant for current density
%% initialize phase space (Quiet - start problem )
tic
disp('Loading particles ...');
nbins = 64;
mpart = param.Np/nbins;
n_electron = param.I*param.lambda0*param.zsep/e0/c;
p1 = zeros(param.Np,1);    

%radfield=param.E0*ones(param.Nsnap,param.nslices);
thetap = zeros(param.Nsnap,param.nslices,param.Np);
gammap=zeros(param.Nsnap,param.nslices,param.Np);

for islice = 1:param.nslices
X0 = hammersley(2,param.Np);
gammap(1,islice,:) = param.gamma0+param.deltagamma*X0(1,:);

auxtheta1 = hammersley(1,mpart)'*2*pi/nbins-pi;

for jbin = 1:nbins
    for ipart = 1:mpart
        thetap(1,islice,ipart+(jbin-1)*mpart)=auxtheta1(ipart)+2*(jbin-1)*pi/nbins;
    end
end

if(param.shotnoise)
    an = 2*sqrt(-log(rand(1))/n_electron);    
    phin = rand(1)*2*pi;
    for ipart = 1:param.Np;
    thetap(1,islice,ipart) = thetap(1,islice,ipart)-an*sin(thetap(1,islice,ipart)+phin);
    end    
end

if (param.prebunching)
   % Double buncher scheme a la N. Sudar http://www.sciencedirect.com/science/article/pii/S0168900217301924  
    [thetap(1,islice,:),gammap(1,islice,:)]=prebunch_particles(squeeze(thetap(1,islice,:)),squeeze(gammap(1,islice,:)),param);    
   % Single buncher scheme
    %[thetap(1,islice,:),gammap(1,islice,:)]=single_prebuncher_particles(squeeze(thetap(1,islice,:)),squeeze(gammap(1,islice,:)),param);    
end
bunching(islice) = (sum(exp(1i.*thetap(1,islice,:))/param.Np));
end

disp('Particles loaded successfully');
toc
