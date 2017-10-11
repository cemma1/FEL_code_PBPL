%% initialize phase space (Quiet - start problem )
tic
disp('Loading particles ...');
Np = param.Np;
nbins = 64;
mpart = Np/nbins;
n_electron = param.I*param.lambda0*param.zsep/e0/c;
p1 = zeros(Np,1);    

radfield=ones(param.Nsnap,param.nslices)*param.E0;
%radfield_3rd=ones(param.Nsnap,param.nslices)*param.E03rdharm;
thetap = zeros(param.Nsnap,param.nslices,Np);
gammap=zeros(param.Nsnap,param.nslices,Np);

for islice = 1:param.nslices
X0 = hammersley(2,Np);
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
    for ipart = 1:Np;
    thetap(1,islice,ipart) = thetap(1,islice,ipart)-an*sin(thetap(1,islice,ipart)+phin);
    end    
end

if (param.prebunching)
   % thetap(1,islice,:) = thetap(1,islice,:)-2.*param.bunch*sin(thetap(1,islice,:)+param.bunchphase);
   % Double buncher scheme a la Nick   
    [thetap(1,islice,:),gammap(1,islice,:)]=prebunch_particles(squeeze(thetap(1,islice,:)),squeeze(gammap(1,islice,:)),param);    
   % Single buncher scheme
    %[thetap(1,islice,:),gammap(1,islice,:)]=single_prebuncher_particles(squeeze(thetap(1,islice,:)),squeeze(gammap(1,islice,:)),param);    
end
bunching(islice) = (sum(exp(1i.*thetap(1,islice,:))/Np));
end

disp('Particles loaded successfully');
toc