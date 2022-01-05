% Perave_code_user input
%%%%User entered parameters%%%%%
%% Undulator parameters
param.lambdau = 3.0e-2;                                       % undulator period
param.K = 3.5; %e0*Bfield*/me/c/ku;                           % Peak undulator parameter
param.ku = 2.*pi./param.lambdau;                              % undulator wavenumber
lwig=param.lambdau*30;                                        % Undulator length m    
% Tapering options
param.tapering = 0;                                           % tapering (0 no tapering ; 1 decelation)    
param.z0 = param.lambdau;
param.psir = pi/4;
% Polynomial Taper options
param.order = 1;
param.ratio = 0.03/(lwig-param.z0)^param.order;
%% Simulation control options
param.phasespacemovie=0;
param.itdp = 1;
param.saveoutput=0;
% Set simulation length and # of snapshots
param.delz=1;
param.stepsize = param.lambdau*param.delz;
param.Nsnap = round(lwig/param.stepsize);                     % Num of snapshots to take over the length of the undulator
param.shotnoise = 1;
param.zsep = param.delz;                                                              
if(~param.itdp)
    param.nslices = 1;
    param.shotnoise =0;                                       % Note if you want to model time independent start-up from noise set P0 = pnoise
else
    param.nslices = round(6*param.Nsnap);                     % Num of slices: Note you want more than 1 slippage length Nsnap
end
param.addLSC = 1;                                             % LSC on/off
%% radiation parameters
param.lambda0 = 2.5e-9;                                       % Seed wavelength
param.k = 2*pi/param.lambda0;                                 % wavenumber in free space
P0 = 1e5; param.P0=P0;                                        % Seed power (W) 
zr = 5;                                                       % Rayleigh length of seed
param.waist = sqrt(zr*param.lambda0/pi);
A_mode = pi*param.waist^2/2;
param.E0 = sqrt(P0/c/eps0/A_mode);                            % Assume circular polarization  
%% Electron beam parameters
param.gamma0 = sqrt(param.k/2/param.ku*(1+param.K^2/2));      % relativistic gamma factor
param.Np = 8192;                                              % # of macroparticles (500-1000 well) 
param.Ee = param.gamma0*me*c^2/e0;                            % Total e-beam energy (eV)
energyspread = 10;                                            % Absolute energy spread MeV
param.deltagammarel = energyspread/param.gamma0/0.511;        % Relative energy spread dgamma/gamma
param.deltagamma = param.gamma0*param.deltagammarel;
param.prebunching = 0;                                        % set to 1 to start from a pre-bunched beam. 
if param.prebunching
    % Parameter definition in C. Emma et al PRAB 20, 110701 (2017)
    param.Abh = 20;                                           % Initial bunching modulation amplitude
    param.bunchphase = param.psir;                            % Initial bunching phase
    P0required=(param.Abh*param.deltagamma)^4*(0.511e6*param.ku/param.K/2)^2*A_mode/377; 
    param.P0=P0required;                                      % Set the seed power to the power required for pre-bunching
    param.E0 = sqrt(param.P0/c/eps0/A_mode);                  % Assume circular polarization
end
betax=1;                                                      % Beta function
emitx=0.4e-6;                                                 % Transverse emittance
param.I = 8e5;                                                % beam current 
param.sigmax = sqrt(betax*emitx/param.gamma0);                % beam radius
param.equivEnSpread = emitx^2*param.lambdau/(4*param.lambda0*param.sigmax^2*param.gamma0);  % Equivalent energy spread from emittance Maier in Synchrotron Light sources and FELs, units of gamma
param.A_e = 2*pi*param.sigmax^2;                              % beam cross section 
bunchlength=param.nslices*param.zsep*param.lambda0/c;
% To include a non-uniform current distribution 
param.currprofile = 1;                                        % 0 = uniform current profile 1 = gaussian 2 = trapezoid
param.sigmat = 20e-18;                                        % beam sigma [s], valid only for Gaussian Current param.currprofile =1
param.currgradient = param.I;                                 % Slope of trapezoidal current distribution
%% Simplifying constants
param.chi2 = e0/me/c^2/2; % These are the planar ones
% Constant for the resonant phase calculation   
const_resp=1/param.chi2*(param.lambdau/2/param.lambda0);
% This is not correct for the planar case it has to be multipled by 
% 2*JJ(z), see eq 5 in my section of Halavanau J. Synchrotron Rad. (2019). 26, 635–646

% These are the corresponding constatns for a helical undulator
%param.chi2 = e0/me/c^2; 
%param.chi1=mu0*c/2*param.I/param.A_e; 
%const_resp=1/param.chi2*(param.lambdau/2/param.lambda0);



