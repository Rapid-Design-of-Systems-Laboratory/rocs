function [h,p,rho,T] = standardAtmosphere(hG,system)
%
% This function calculates pressure, density, and temperature at the
% appropriate geometric altitude.
%
% Inputs:
%   hG - geometric altitude, m
%   system - system of units (1=metric, 2=English)
%
% Outputs:
%   h - potential altitude, m
%   p - pressure, Pa
%   rho - density, kg/m^3
%   T - temperature, K
%          
%
% Developed by Michael J. Grant on 08/30/12
%

%% Constants
T0 = 288.16; % Surface temperature, K
p0 = 1.01325e5; % Surface pressure, Pa
rho0 = 1.2250; % Surface density, kg/m^3
re = 6.356766e6; % Earth radius, m
h0 = 0; % Surface altitude, m
g0 = 9.806; % Surface acceleration, m/s^2
R = 287; % Specific gas constant, J/kg-K
hGLayers = [11000 25000 47000 53000 79000]; % Layer boundaries, m
layerGradient = [-6.5e-3 inf 3e-3 inf -4.5e-3]; % Temperature change gradient in layer, K/m
ft2m = 0.3048; % feet to meter conversion
m2ft = 1/ft2m; % meter to feet conversion
K2R = 1.8; % Kelvin to Rankine conversion
densM2densEng = 0.00194032033; % kg/m^3 conversion to slug/ft^3
Pa2psf = 0.0208854342; % Pascal to psf conversion

%% Default to metric if no inputs provided
if nargin < 2
  system = 1;
end

%% Convert to metric if in English
if system == 2
  hG = hG*ft2m;
end

%% Determine layer altitude resides in
layer = find(hG > hGLayers,1,'last') + 1; % Add layer for final altitude calculation
if isempty(layer) % Account for first layer
  layer = 1;
end

% Initialize first layer values which are surface values
T1 = T0;
p1 = p0;
rho1 = rho0;
h1 = h0;

% Calculate up from surface to layer altitude resides in
for ctrLayer = 1 : 1 : layer
  
  % Determine geometric altitude of calculation
  if ctrLayer ~= layer % Top of layer
    hGLayer = hGLayers(ctrLayer);
  else
    hGLayer = hG;
  end
  
  % Determine geopotential alitude of calculation
  hLayer = re/(re+hGLayer)*hGLayer;
  
  % If layer is gradient, perform gradient calculation. Otherwise, perform
  % isothermal calculation.
  if layerGradient(ctrLayer) ~= inf % Gradient
    
    % Gradient value
    a = layerGradient(ctrLayer);
    
    % Values at top of gradient layer
    T = T1 + a*(hLayer-h1);
    p = p1*(T/T1)^(-g0/(a*R));
    rho = rho1*(T/T1)^(-(g0/(a*R)+1));
    
    % Assign values as base of next layer
    T1 = T;
    p1 = p;
    rho1 = rho;
    h1 = hLayer;
    
  else % Isothermal
    
    % Values at top of isothermal layer
    T = T1;
    p = p1*exp(-g0/(R*T)*(hLayer-h1));
    rho = rho1*p/p1;
    
    % Assign values as base of next layer
    T1 = T;
    p1 = p;
    rho1 = rho;
    h1 = hLayer;
    
  end
  
end

% With loop finished, calculation performed to requested altitude. Have all
% outputs
h = hLayer;

%% Output

% Convert back to English if requested
if system == 2
  h = h*m2ft;
  p = p*Pa2psf;
  rho = rho*densM2densEng;
  T = T*K2R;
end

return

