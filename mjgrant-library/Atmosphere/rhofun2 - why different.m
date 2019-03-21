function [rho]=rhofun2(h)
%  function [rho]=rhofun2(h)
%  Standard Atmosphere Computations in English Units
%  for the 1976 standard atmosphere up to 230,000 ft.
%  Author: Ilan Kroo (kroo@leland.stanford.edu)  31 Dec 95.
%  Converted to MATLAB by D. Andrisani, 2 Nov 99.
%  Vector input h is geometric altitude in feet.
%  Vector output rho is the same size as vector input h.
%
%  Output   Units
%  rho      slug/ft^3     (density)
%
%  Because h is a vector, the arithmetic operations
%  used below (e.g., .^)are the array operations.
%  The find command is used to determine which region
%  in the atmosphere the altitude data is in, and then
%  the correct equation can be used for each atmospheric region.


RHOSL = 0.00237689;   % slug/ft^3
saSigma =NaN*zeros(size(h)) ; % Initialize saSigma to the
% correct array size and fill the array with elements 
% that are not numbers. NaN is the indicator that MATLAB
% uses for 'not a number'. In this application we fill this
% array initially with NaN so that the output
% array rho will contain NaN corresponding to the 
% altitudes which are out of range (h<232940). MATLAB will not
% allow computations to be done on NaN array elements.

   j1=find(h<232940 & h>=167323);
   if  isempty(j1) 
   else
      saSigma(j1) = ( 0.79899-h(j1)/606330).^11.20114;
   end

   j2=find(h<167323 & h>=154199);
   if   isempty(j2)
   else saSigma(j2) = 0.00116533*exp((h(j2)-154200)/-25992);
   end

   j3=find(h<154199 & h>=104987);
   if isempty(j3)
   else   saSigma(j3) = (0.857003+h(j3)/190115).^-13.20114;
   end

   j4=find(h<104987 & h>=65617);
   if isempty(j4)
   else  saSigma(j4) = (0.978261+h(j4)/659515).^-35.16319;
   end

   j5=find(h<65617 & h>=36089);
%  This is the stratosphere.
   if isempty(j5)
   else   saSigma(j5) = 0.297076 *exp((36089-h(j5))/20806);
   end

   j6=find(h<36089);
%  This is the troposphere.
   if isempty(j6)
   else   saSigma(j6) = (1.0-h(j6)/145442).^4.255876;
   end

rho = RHOSL * saSigma;   % slug/ft^3
