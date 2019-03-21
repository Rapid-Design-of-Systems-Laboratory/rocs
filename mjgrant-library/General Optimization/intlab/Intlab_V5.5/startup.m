%STARTUP      Dummy routine to call startintlab 
%
%Insert the call to "startintlab" into your personal startup-file.
%

% written  09/10/00     S.M. Rump  
% modified 05/19/08     S.M. Rump  catch error from BLAS_VERSION
%

  try
    startintlab
  catch
    disp('========================================================================')
    disp('========================================================================')
    disp('*** This should not happen. It might be that you used INTLAB         ***')
    disp('***   before and set the environment variable BLAS_VERSION.          ***')
    disp('*** In previous versions of INTLAB this was necessary because        ***')
    disp('***   Intel Math Kernel Library (IMKL) did not work with my          ***')
    disp('***   switching of the rounding mode.                                ***')
    disp('*** Meanwhile this is changed, IMKL should work properly. So please  ***')
    disp('***   delete the environment variable BLAS_VERSION and try again.    ***')
    disp('*** Sorry for inconvenience.                                         ***')
    disp('========================================================================')
    disp('========================================================================')
  end
