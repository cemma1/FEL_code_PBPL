% Perave_code_user input
%%%%User entered parameters%%%%%
%% Undulator parameters
param.lambdau = 3.0e-2;                                     % undulator period
param.K = 2.475; %e0*Bfield*/me/c/ku;                       % RMS undulator parameter
param.ku = 2.*pi./param.lambdau;                            % undulator wavenumber
lwig=param.lambdau*1e3;                                     % Undulator length m    
% Tapering options
param.tapering = 0;                                         % tapering (0 no tapering ; 1 decelation)    
param.z0 = param.lambdau*10;
param.psir = pi/4;
%% Simulation control options
param.phasespacemovie=1;
param.itdp = 0;
param.saveoutput=1;
% Set simulation length and # of snapshots
param.delz=1;
param.stepsize = param.lambdau*param.delz;
param.Nsnap = round(lwig/param.stepsize);                    % number of snapshots to take over the length of the undulator
param.shotnoise = 1;                                % Only matters if param.itdp = 1;
param.zsep = param.delz;                                                              
if(~param.itdp)
    param.nslices = 1;
    param.shotnoise =0;                            % Note if you want to model time independent start-up from noise set P0 = pnoise
else
    param.nslices = round(5*param.Nsnap);          % Num of slices: Note you want more than 1 slippage length Nsnap
end
%% radiation parameters
param.lambda0 = 1.5*1e-10;                                    % Seed wavelength (m)
param.k = 2*pi/param.lambda0;                                 % wavenumber in free space
P0 = 1e4; param.P0=P0;                                        % Seed power (W) 
zr = 5;                                                       % Rayleigh length of seed (m)
param.waist = sqrt(zr*param.lambda0/pi);
A_mode = pi*param.waist^2/2;
param.E0 = sqrt(2*P0/c/eps0/A_mode/2);                        % Assume circular polarization  
%% Electron beam parameters
param.gamma0 = sqrt(param.k/2/param.ku*(1+param.K^2));        % relativistic gamma factor
param.Np = 512;                                               % # of macroparticles (500-1000 well) 
param.Ee = param.gamma0*me*c^2/e0;                            % Total e-beam energy (eV)
energyspread = 1.5;                                           % Absolute energy spread MeV
param.deltagammarel = energyspread/param.gamma0/0.511;        % Relative energy spread dgamma/gamma
param.deltagamma = param.gamma0*param.deltagammarel;
param.prebunching = 0;                                        % set to 1 to start from a pre-bunched beam. 
if param.prebunching
    % Parameter definition in C. Emma et al PRAB 20, 110701 (2017)
    param.Abh = 20;                                               % Initial bunching modulation amplitude
    param.bunchphase = param.psir;                                % Initial bunching phase
    P0required=(param.Abh*param.deltagamma)^4*(0.511e6*param.ku/param.K/2)^2*A_mode/377; 
    param.P0=P0required;                                          % Set the seed power to the power required for pre-bunching
    param.E0 = sqrt(param.P0/c/eps0/A_mode);                      % Assume circular polarization
end
betax=10;                                                     % Beta function (m)
emitx=0.4e-6;                                                 % Transverse emittance
param.I = 4000;                                               % beam current 
param.sigmax = sqrt(betax*emitx/param.gamma0);                % beam radius
param.A_e = 2*pi*param.sigmax^2;                              % beam cross section 
bunchlength=param.nslices*param.zsep*param.lambda0/c;
%% Simplifying constants
param.chi2 = e0/me/c^2;
param.chi1=mu0*c/2*param.I/param.A_e;
% Constant for the resonant phase calculation   
const_resp=1/param.chi2*(param.lambdau/2/param.lambda0);


