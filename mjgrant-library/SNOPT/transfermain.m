%  function [x,F,INFO] = transfermain()
%function [x,F,INFO] = snoptmain3()
% Defines the NLP problem and calls the mex interface for snopt.
% Some, but not all first derivatives are provided.

transfermain.spc = which('transfermain.spc');

snprint  ( 'transfermain.out' );
snsummary( 'transfermain.sum' );
snspec   (  transfermain.spc  );
snseti   ( 'Verify level     ', 3);
snseti   ( 'Derivative option', 0);

% REMEMBER to define iGfun and jGvar
A = [];
iAfun = [];
jAvar = [];

[x,F,INFO] = snopt(x,xlow,xupp,Flow,Fupp,'transferfun', ...
		   A, iAfun, jAvar, iGfun, jGvar);

snsummary off;
snprint   off; % Closes the file and empties the print buffer
