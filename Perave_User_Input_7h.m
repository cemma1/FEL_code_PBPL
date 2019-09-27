
% Perave_code_user input
%%%%User entered parameters%%%%%
%% Undulator parameters
param.lambdau = 3.9e-2;                                     % undulator period
param.K = 5.5/sqrt(2); %e0*Bfield*/me/c/ku;                 % RMS undulator parameter
param.ku = 2.*pi./param.lambdau;                            % undulator wavenumber
lwig=param.lambdau*100;                                     % Undulator length m    
% Tapering options
param.tapering = 0;                                         % tapering (0 no tapering ; 1 KMR; 2 linear; 3 constant area; 4 polynomial)    
param.z0 = param.lambdau;
param.psir = 3*pi/8;
param.psirgradient = 27*1/lwig*pi/180;                      % 0 by default, Only used if param.tapering = 2 (linear psir taper)
param.dtaper = 1;
deltaKoverK = 0.1;
param.ctaper = deltaKoverK/((lwig-param.z0).^param.dtaper);
%param.psir = psirvalues(psirindex);% For scanning
%% Simulation control options
param.phasespacemovie=0;
param.itdp = 1;
param.saveoutput=0;
% Set simulation length and # of snapshots
param.delz=1;   
param.stepsize = param.lambdau*param.delz;
param.Nsnap = round(lwig/param.stepsize);                    % number of snapshots to take over the length of the undulator
param.shotnoise = 1;
param.zsep = param.delz;                                                              
if(~param.itdp)
    param.nslices = 1;
    param.shotnoise =1;                                        % Note if you want to model time independent start-up from noise set P0 = pnoise
else
    param.nslices = round(3.5*param.Nsnap);                    % Note you want more than 1 slippage length (Nsnap)
end
%% radiation parameters
param.lambda0 = 8e-10;                                        % Seed wavelength
param.k = 2*pi/param.lambda0;                                 % wavenumber in free space
P0 = 2.7e6; param.P0=P0;                                      % Seed power (W) 
zr = 5;                                                       % Rayleigh length of seed
param.waist = sqrt(zr*param.lambda0/pi);
A_mode = pi*param.waist^2/2;
param.E0 = sqrt(2*P0/c/eps0/A_mode/2);                        % Assume circular polarization
% To include a non-uniform seed field distribution 
param.fieldprofile = 0;                                       % 0 = uniform current profile 1 = gaussian
param.sigmatfield = 250e-18;                                  % beam sigma [s], only for param.currprofile =1
%% Electron beam parameters
param.gamma0 = sqrt(param.k/2/param.ku*(1+param.K^2));        % relativistic gamma factor
param.chirpslope = 0e-3/param.nslices;                       % dimensionless chirp parameter
param.Np = 512;                                               % # of macroparticles (500-1000 well) 
param.Ee = param.gamma0*me*c^2/e0;                            % Total e-beam energy (eV)
energyspread = 8;                                             % Absolute energy spread MeV
param.deltagammarel = energyspread/param.gamma0/0.511;        % Relative energy spread dgamma/gamma
param.deltagamma = param.gamma0*param.deltagammarel;
param.prebunching = 0;                                        % set to 1 to start from a pre-bunched beam. 
if param.prebunching
param.Abh = 20;                                               % Initial bunching modulation amplitude
param.bunchphase = param.psir;                                % Initial bunching phase
P0required=(param.Abh*param.deltagamma)^4*(0.511e6*param.ku/param.K/2)^2*A_mode/377; 
param.P0=P0required;                                          % Set the seed power to the power required for pre-bunching
param.E0 = sqrt(param.P0/c/eps0/A_mode);                      % Assume circular polarization
end
betax=10;                                                     % Beta function
emitx=50e-9;                                                  % Transverse emittance
param.I = 80e3;                                               % Peak current 
param.sigmax = sqrt(betax*emitx/param.gamma0);                % beam radius
param.A_e = 2*pi*param.sigmax^2;                              % beam cross section 
bunchlength=param.nslices*param.zsep*param.lambda0/c;
% To include a non-uniform current distribution 
param.currprofile = 0;                                        % 0 = uniform current profile 1 = gaussian 2 = trapezoid
param.sigmat = (110/2.355)*1e-18;                             % beam sigma [s], only for gaussian (param.currprofile =1)
param.currgradient = param.I;                                 % Slope of trapezoidal current distribution
%% Simplifying constants
param.chi2 = e0/me/c^2;
% Constant for the resonant phase calculation   
const_resp=1/param.chi2*(param.lambdau/2/param.lambda0);