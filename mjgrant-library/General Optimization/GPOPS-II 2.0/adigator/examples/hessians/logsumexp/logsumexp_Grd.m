% function [Grd,Fun] = logsumexp_Grd(x)
% 
% Gradient wrapper file generated by ADiGator
% �2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function [Grd,Fun] = logsumexp_Grd(x)
gator_x.f = x;
gator_x.dx = ones(32,1);
y = logsumexp_ADiGatorGrd(gator_x);
Grd = reshape(y.dx,[1 32]);
Fun = y.f;
end