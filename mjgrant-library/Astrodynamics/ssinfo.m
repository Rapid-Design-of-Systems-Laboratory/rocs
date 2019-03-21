function [axper,r,mu,a,orbper,e,inc] = ssinfo(body)
%SSINFO Obtain information for a body in the solar system
%  [AXPER,R,MU,A,ORBPER,E,INC] = SSINFO(BODY) returns properties
%  cooresponding to BODY.
%
%  Input:
%    BODY = body name ('sun', 'moon', 'mercury', 'venus', 'earth',
%      'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'pluto')
%      [string]
%
%  Output:
%    AXPER = axial rotational period [rev/day]
%    R = equatorial radius [km]
%    MU = gravitational parameter (G*m) [km^3/s^2]
%    A = semi-major axis of orbit [km]
%    ORBPER = orbital period [sec]
%    E = eccentricity of orbit [nondimensional]
%    INC = inclination of orbit to ecliptic [deg]
%

%   Michael J. Grant / May 28, 2006


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OBTAIN PROPERTIES OF CELESTIAL BODY %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch body
  case 'sun'
    axper = 0.0394011;
    r = 695990;
    mu = 132712200000;
    a = NaN;
    orbper = NaN;
    e = NaN;
    inc = NaN;
  case 'moon'
    axper = 0.036601;
    r = 1738.2;
    mu = 4902.7949;
    a = 384388.174; % With respect to Earth
    orbper = 2360591;
    e = 0.0549; % With respect to Earth
    inc = 5.228388;
  case 'mercury'
    axper = 0.0170513958333;
    r = 2439.7;
    mu = 22032.09;
    a = 57909175.407279;
    orbper = 7600558.6632;
    e = 0.20563069;
    inc = 7.00497999999999;
  case 'venus'
    axper = 0.0041148030556; % Westward
    r = 6051.9;
    mu = 324858.592079;
    a = 108208925.006886;
    orbper = 19414180.7168;
    e = 0.00677323;
    inc = 3.39471;
  case 'earth'
    axper = 1.0027378777778;
    r = (6356.75231424 + 6378.137) / 2; % Average
    mu = 398600.4418;
    a = 149597886.455765;
    orbper = 31558229.5413;
    e = 0.01671022;
    inc = 0.00004999999999998;
  case 'mars'
    axper = 0.9746999527778;
    r = (3375 + 3397) / 2; % Average
    mu = 42828.371901284;
    a = 227936636.175727;
    orbper = 59353406.0271;
    e = 0.09341233;
    inc = 1.85060999999999;
  case 'jupiter'
    axper = 2.4181555555556;
    r = (66854 + 71492) / 2; % Average
    mu = 126686535;
    a = 778412023.132788;
    orbper = 374574964.4705;
    e = 0.04839266;
    inc = 1.3053;
  case 'saturn'
    axper = 2.2522052844444;
    r = (54364 + 60268) / 2; % Average
    mu = 37931284;
    a = 1426725405.91221;
    orbper = 929469271.2322;
    e = 0.05415059999;
    inc = 2.48445999;
  case 'uranus'
    axper = 1.3921113688889;
    r = (24973 + 25559) / 2; % Average
    mu = 5803200;
    a = 2870972206.53583;
    orbper = 2653187400.862;
    e = 0.04716771;
    inc = 0.76986;
  case 'neptune'
    axper = 1.3020833333333;
    r = (24800 + 25269) / 2; % Average
    mu = 6835107;
    a = 4498252889.71578;
    orbper = 5203436337.63;
    e = 0.00858587;
    inc = 1.76917;
  case 'pluto'
    axper = 0.1565666666667;
    r = 1162;
    mu = 1009.076;
    a = 5906376244.79917;
    orbper = 7828996647.593;
    e = 0.24880766;
    inc = 17.1417499999999;
end

return

