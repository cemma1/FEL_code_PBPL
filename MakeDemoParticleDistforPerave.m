% Generates a mock Gaussian phase space distribution, slices it up and stores it as input for Perave
charge = param.I*sqrt(2*pi*param.sigmat^2);% Must be specified for the current profile calculation
Ntotal = 1e6;% Decide how many particles to fill your distribution with
tGauss = randn(Ntotal,1)*param.sigmat;
gammaGauss = param.gamma0 + randn(Ntotal,1)*param.deltagamma;
% Cut particles at +/- 5 sigma
tmaxGauss = 5*param.sigmat;
idx = abs(tGauss)<tmaxGauss;
% PS stores your phase space co-ordinates
ps(:,5) = tGauss(idx);
ps(:,6) = gammaGauss(idx);
figure;histogram2(ps(:,5),ps(:,6),'DisplayStyle','Tile');xlabel('t [as]');ylabel('\gamma');
%% Load mock particle dist into Perave theta and gamma arrays
NPART = param.Np;
XLAMDS = param.lambda0;
NSLICE = round(3e8*(range(ps(:,5)))/XLAMDS);
zsep = 1;

thetap = zeros(param.Nsnap,param.nslices,param.Np);
gammap=zeros(param.Nsnap,param.nslices,param.Np);

% Load the particles slice by slice
slice_ps_cut = zeros(NPART,6);
t0 = -1*tmaxGauss;
for islice = 1:NSLICE
   clearvars slice_ps slice_ps_cut
   tmin = t0+(islice-1)*XLAMDS/3e8;
   tmax = t0+islice*XLAMDS/3e8;
   
   idslice = tmin<ps(:,5) & ps(:,5)<tmax;
   npart_per_slice(islice) = sum(idslice);
   
   jump = floor(npart_per_slice(islice)/NPART);
   slice_ps(1:sum(idslice),1:6)= ps(idslice,:);

   % Cut the particles to NPART instead of however many there are in the slice   
   if jump<1       
       tminslice = tmin*3e8*2*pi/XLAMDS-2*pi*(islice-1);
       dummycoords = repmat([0.0,0.0,0.0,0.0,tmin*3e8*2*pi/XLAMDS-2*pi*(islice-1),mean(slice_ps(:,6))],NPART-npart_per_slice(islice),1);                
       dummycoords(:,5) = tmin+ (tmax-tmin)*rand(size(dummycoords,1),1);% Randomly assign a t value ( to be converted to theta later in the code)
       slice_ps_cut(1:NPART,1:6) = cat(1,slice_ps(:,1:6),dummycoords);% If there's less more than NPART just fill with dummy coordinates                              
   else  
        idxp = randi(size(slice_ps,1),NPART,1);
        slice_ps_cut(1:NPART,1:6) = slice_ps(idxp,1:6);
   end
                      
   % Store average slice quantities for beamfile
   zpos(islice) = 3e8*tmin;
   gamma0(islice) = mean(slice_ps(:,6));
   delgam(islice) = std(slice_ps(:,6));
         
   % Check if any quantities are NaN set them to previous slice value
   if any(isnan([gamma0(islice),delgam(islice)]))
       gamma0(islice) = mean(ps(:,6)); delgam(islice) = 0.0;
   end
    
    % Load LPS into Perave theta and gamma 
    thetap(1,islice,:) = slice_ps_cut(:,5)*3e8*2*pi/XLAMDS-2*pi*(islice-1);%theta, Note: Don't forget to go from t to theta
    gammap(1,islice,:) = slice_ps_cut(:,6);%gamma
    bunching(islice) = abs(sum(exp(1i.*thetap(1,islice,:))/NPART));           
        
end
% Calculate the current profile
curpeak = npart_per_slice*charge/trapz(linspace(0,tmax,NSLICE),npart_per_slice);
param.Iprofile = curpeak;

if length(param.Iprofile)~=param.nslices
   error(['Current profile slices not equal to nslices in the simulation.'...
       'Set them equal by changing param.nslices in the input file']) 
end
% Plot current profile and bunching
figure;yyaxis left;plot(bunching);yyaxis right;plot(curpeak)

