function [solution] = formSolver(expression,value,variable,assumptions,in,oc)
	% This function creates an expression that represents the solution to 'expression = value' for variable.
	% The supplied assumptions are used to try and construct an analytic solution. If no analytic solution
	% can be found, a numerical root solving (fzero) expression is constructed.
	
	% initialGuess = 0; % make this smarter in the future
	
	% Determine if an analytic solution can be found
    expression
    variable
    %%% for h adjoining
    if strcmp(char(variable),'Atrig1')
        solution = {'atan(((Amax*lamGAM*rho0*v*exp(-h/H)*sin(alfa)*(((T4*((y/2 - 1/2)*Mc^2 + 1))/(T0*(((y/2 - 1/2)*v^2)/(R*T0*y) + 1)))^(1/2) - 1))/(2*mass) + (Amax*lamV*rho0*v^2*exp(-h/H)*cos(alfa)*(((T4*((y/2 - 1/2)*Mc^2 + 1))/(T0*(((y/2 - 1/2)*v^2)/(R*T0*y) + 1)))^(1/2) - 1))/(2*mass) + (Amax*T0*cp*lamMASS*rho0*v*exp(-h/H)*(((y/2 - 1/2)*v^2)/(R*T0*y) - (T4*((y/2 - 1/2)*Mc^2 + 1))/T0 + 1))/(2*hpr))/(epsilon2*lamH))'
                    'atan(((Amax*lamGAM*rho0*v*exp(-h/H)*sin(alfa)*(((T4*((y/2 - 1/2)*Mc^2 + 1))/(T0*(((y/2 - 1/2)*v^2)/(R*T0*y) + 1)))^(1/2) - 1))/(2*mass) + (Amax*lamV*rho0*v^2*exp(-h/H)*cos(alfa)*(((T4*((y/2 - 1/2)*Mc^2 + 1))/(T0*(((y/2 - 1/2)*v^2)/(R*T0*y) + 1)))^(1/2) - 1))/(2*mass) + (Amax*T0*cp*lamMASS*rho0*v*exp(-h/H)*(((y/2 - 1/2)*v^2)/(R*T0*y) - (T4*((y/2 - 1/2)*Mc^2 + 1))/T0 + 1))/(2*hpr))/(epsilon2*lamH)) + pi'
                    };
    %%% for m adjoining
    elseif strcmp(char(variable),'Atrig2')
        solution = {'atan(((Amax*lamGAM*rho0*v*exp(-h/H)*sin(alfa)*(((T4*((y/2 - 1/2)*Mc^2 + 1))/(T0*(((y/2 - 1/2)*v^2)/(R*T0*y) + 1)))^(1/2) - 1))/(2*mass) + (Amax*lamV*rho0*v^2*exp(-h/H)*cos(alfa)*(((T4*((y/2 - 1/2)*Mc^2 + 1))/(T0*(((y/2 - 1/2)*v^2)/(R*T0*y) + 1)))^(1/2) - 1))/(2*mass) + (Amax*T0*cp*lamMASS*rho0*v*exp(-h/H)*(((y/2 - 1/2)*v^2)/(R*T0*y) - (T4*((y/2 - 1/2)*Mc^2 + 1))/T0 + 1))/(2*hpr))/(epsilon2*lamMASS))'
                    'atan(((Amax*lamGAM*rho0*v*exp(-h/H)*sin(alfa)*(((T4*((y/2 - 1/2)*Mc^2 + 1))/(T0*(((y/2 - 1/2)*v^2)/(R*T0*y) + 1)))^(1/2) - 1))/(2*mass) + (Amax*lamV*rho0*v^2*exp(-h/H)*cos(alfa)*(((T4*((y/2 - 1/2)*Mc^2 + 1))/(T0*(((y/2 - 1/2)*v^2)/(R*T0*y) + 1)))^(1/2) - 1))/(2*mass) + (Amax*T0*cp*lamMASS*rho0*v*exp(-h/H)*(((y/2 - 1/2)*v^2)/(R*T0*y) - (T4*((y/2 - 1/2)*Mc^2 + 1))/T0 + 1))/(2*hpr))/(epsilon2*lamMASS)) + pi'
                    };
                
    elseif strcmp(char(variable),'Atrig3')
        solution = {'atan(((Amax*lamV*rho0*v^2*exp(-h/H)*cos(alfa)*(((((atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q1)) - atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q2)))*(((y/2 - 1/2)*v^2)/(R*T0*y) - (T4*((y/2 - 1/2)*Mc^2 + 1))/T0 + 1))/pi + (v^2*(y/2 - 1/2))/(R*T0*y) + 1)/(((y/2 - 1/2)*v^2)/(R*T0*y) + 1))^(1/2) - 1))/(2*mass) -(lamMASS*Amax*T0*cp*rho0*v*exp(-h/H)*(atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q1)) - atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q2)))*(((y/2 - 1/2)*v^2)/(R*T0*y) - (T4*((y/2 - 1/2)*Mc^2 + 1))/T0 + 1))/(2*pi*hpr) + (Amax*lamGAM*rho0*v*exp(-h/H)*sin(alfa)*(((((atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q1)) - atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q2)))*(((y/2 - 1/2)*v^2)/(R*T0*y) - (T4*((y/2 - 1/2)*Mc^2 + 1))/T0 + 1))/pi + (v^2*(y/2 - 1/2))/(R*T0*y) + 1)/(((y/2 - 1/2)*v^2)/(R*T0*y) + 1))^(1/2) - 1))/(2*mass))/(lamMASS*epsilon2))'
                    'atan(((Amax*lamV*rho0*v^2*exp(-h/H)*cos(alfa)*(((((atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q1)) - atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q2)))*(((y/2 - 1/2)*v^2)/(R*T0*y) - (T4*((y/2 - 1/2)*Mc^2 + 1))/T0 + 1))/pi + (v^2*(y/2 - 1/2))/(R*T0*y) + 1)/(((y/2 - 1/2)*v^2)/(R*T0*y) + 1))^(1/2) - 1))/(2*mass) -(lamMASS*Amax*T0*cp*rho0*v*exp(-h/H)*(atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q1)) - atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q2)))*(((y/2 - 1/2)*v^2)/(R*T0*y) - (T4*((y/2 - 1/2)*Mc^2 + 1))/T0 + 1))/(2*pi*hpr) + (Amax*lamGAM*rho0*v*exp(-h/H)*sin(alfa)*(((((atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q1)) - atan(qscale*(- 0.5*rho0*exp(-h/H)*v^2 + q2)))*(((y/2 - 1/2)*v^2)/(R*T0*y) - (T4*((y/2 - 1/2)*Mc^2 + 1))/T0 + 1))/pi + (v^2*(y/2 - 1/2))/(R*T0*y) + 1)/(((y/2 - 1/2)*v^2)/(R*T0*y) + 1))^(1/2) - 1))/(2*mass))/(lamMASS*epsilon2)) + pi'
                    };

      elseif strcmp(char(variable),'Tx')
        solution = {'((1/(1+y/(1+x)))*(1/5.2)*(((50.+z)^2.5*sec(gam)*(-1.0*lamV*v*cos(alfamax*sin(alfatrig))+sin(alfamax*sin(alfatrig))*(-1.0*lamGAM*cos(bankmax*sin(banktrig))-1.0*lamPSII*sec(gam)*sin(bankmax*sin(banktrig)))))/mass))^(1/4.2)', 
                    '300',
                    '3420'};
                
    elseif strcmp(char(variable),'Ttrig')
        solution = {'pi/2',
                    '-pi/2',
                   'asin(((((v*(z + 50)^2.5)*(-(1560*lamV*cos(alfamax*sin(alfatrig)))/mass-(1560*lamGAM*sin(alfamax*sin(alfatrig))*cos(bankmax*sin(banktrig)))/(mass*v)-(1560*lamPSII*sin(alfamax*sin(alfatrig))*sin(bankmax*sin(banktrig)))/(mass*v*cos(gam))))/(8112.0*cos(gam)*(5500/(x+1)+4500/(y+1))))^(1/4.2)-1860)/1560)'
                     };
                
                % Txy mast hai
    elseif strcmp(char(variable),'Ttrig1')
        solution = {'pi/2',
                    '-pi/2',
                   'asin(((((v*(z + 50)^2.5)*(-(1560*lamV*cos(alfamax*sin(alfatrig)))/mass-(1560*lamGAM*sin(alfamax*sin(alfatrig))*cos(bankmax*sin(banktrig)))/(mass*v)-(1560*lamPSII*sin(alfamax*sin(alfatrig))*sin(bankmax*sin(banktrig)))/(mass*v*cos(gam))))/(8112.0*cos(gam)*(5.5-x/1000)))^(1/4.2)-1860)/1560)'
                     };               
    % Dont touch below control as it is working
    elseif strcmp(char(variable),'Ttrignew')
         solution = {'pi/2',
                    '-pi/2',
                   'asin(((((v*(z + 50)^2.5)*(-(1560*lamV*cos(alfamax*sin(alfatrig)))/mass-(1560*lamGAM*sin(alfamax*sin(alfatrig))*cos(bankmax*sin(banktrignew)))/(mass*v)-(1560*lamPSII*sin(alfamax*sin(alfatrig))*sin(bankmax*sin(banktrignew)))/(mass*v*cos(gam))))/(8112.0*cos(gam)))^(1/4.2)-1860)/1560)'
                     };
     % Dont touch below control as it is working
     elseif strcmp(char(variable),'Ttrigpopnew')
         solution = {'pi/2',
                    '-pi/2',
                   'asin(((((v*(z + 50)^2.5)*(-(1560*lamV*cos(alfamax*sin(alfatrig)))/mass-(1560*lamGAM*sin(alfamax*sin(alfatrig))*cos(bankmax*sin(banktrig)))/(mass*v)-(1560*lamPSII*sin(alfamax*sin(alfatrig))*sin(bankmax*sin(banktrig)))/(mass*v*cos(gam))))/(8112.0*cos(gam)*(1+y/(x+1))))^(1/4.2)-1860)/1560)'
                     };

                 % New for cda
     elseif strcmp(char(variable),'Ttrigpopnew_x')
         solution = {'pi/2',
                    '-pi/2',
                    'asin((((-((1560*lamV*cos(alfamax*sin(alfatrignew)))/mass + (1560*lamPSII*sin(alfamax*sin(alfatrignew))*sin(bankmax*sin(banktrignew)))/(mass*v*cos(gam)))*(v*(y + 1)*(z + 50)^2.5))^(1/4.2))/(37315200.0*cos(gam))-1860)/1560)'
                     };   
%
    elseif strcmp(char(variable),'banktrig')
          solution = {'-pi/2',
                      'pi/2',
                      '-asin(3.1416/(2*bankmax))',
                      'asin(3.1416/(2*bankmax))'
                      };

                  
                 % cda
%      elseif strcmp(char(variable),'Ttrigpopnew_x')
%          solution = {'pi/2',
%                     '-pi/2',
%                    'asin(((((v*(z + 50)^2.5)*(-(1560*lamV*cos(alfamax*sin(alfatrig)))/mass-(1560*lamGAM*sin(alfamax*sin(alfatrig))*cos(bankmax*sin(banktrig)))/(mass*v)-(1560*lamPSII*sin(alfamax*sin(alfatrig))*sin(bankmax*sin(banktrig)))/(mass*v*cos(gam))))/(8112.0*cos(gam)*(4600/(y+1))))^(1/4.2)-1860)/1560)'
%                      };   

%      elseif strcmp(char(variable),'alfatrig')
%          solution = {'pi/2',
%                     '-pi/2',
%                    'asin(atan((lamPSII*sin(bankmax*sin(banktrig)))/(v*lamV*cos(gam))+(lamGAM*cos(bankmax*sin(banktrig)))/(v*lamV))/alfamax)',
%                    'asin(atan((lamPSII*sin(bankmax*sin(banktrig)))/(v*lamV*cos(gam))+(lamGAM*cos(bankmax*sin(banktrig)))/(v*lamV) + pi)/alfamax)'
%                      };
%                  
%      elseif strcmp(char(variable),'banktrig')
%          solution = {'pi/2',
%                     '-pi/2',
%                     'asin(atan(lamPSII/(cos(gam)*lamGAM))/bankmax)',
%                     'asin((atan(lamPSII/(cos(gam)*lamGAM) + pi))/bankmax)'
%                      };
                                                                
     elseif strcmp(char(variable),'Ttrigpopnew2')
         solution = {'pi/2',
                    '-pi/2',
                   'asin(((((v*(z + 50)^2.5)*(-(1560*lamV*cos(alfamax*sin(alfatrig)))/mass-(1560*lamGAM*sin(alfamax*sin(alfatrig))*cos(bankmax*sin(banktrig)))/(mass*v)-(1560*lamPSII*sin(alfamax*sin(alfatrig))*sin(bankmax*sin(banktrig)))/(mass*v*cos(gam))))/(8112.0*cos(gam)*(5+y/(x+1))))^(1/4.2)-1860)/1560)'
                     };
%     if strcmp(char(variable),'Atrigcomb')
%         solution = {'atan(((((1-w)*Tmax*Amax/2)*(w*rho0*exp(-h/H)*v*Amax/2*(v/gc*(sqrt((T4/T0*(1+(y-1)/2*Mc^2))/(1 + (y-1)/2*(v/(sqrt(y*R*T0)))^2))-1))))*(lamV*cos(alfa)/m + lamGAM*sin(alfa)/(m*v) - lamMASS/(g0*((1-w)*IspR + w*(v/gc*(sqrt((T4/T0*(1+(y-1)/2*Mc^2))/(1 + (y-1)/2*(v/(sqrt(y*R*T0)))^2))-1))*hpr/(g0*cp*T0*((T4/T0*(1+(y-1)/2*Mc^2))-(1 + (y-1)/2*(v/(sqrt(y*R*T0)))^2)))))))/(lamV*epsilon2))',
%                     'atan(((((1-w)*Tmax*Amax/2)*(w*rho0*exp(-h/H)*v*Amax/2*(v/gc*(sqrt((T4/T0*(1+(y-1)/2*Mc^2))/(1 + (y-1)/2*(v/(sqrt(y*R*T0)))^2))-1))))*(lamV*cos(alfa)/m + lamGAM*sin(alfa)/(m*v) - lamMASS/(g0*((1-w)*IspR + w*(v/gc*(sqrt((T4/T0*(1+(y-1)/2*Mc^2))/(1 + (y-1)/2*(v/(sqrt(y*R*T0)))^2))-1))*hpr/(g0*cp*T0*((T4/T0*(1+(y-1)/2*Mc^2))-(1 + (y-1)/2*(v/(sqrt(y*R*T0)))^2)))))))/(lamV*epsilon2))+ pi'
%                     }; 
                
    elseif strcmp(char(variable),'alfadotTtrig')
        solution = {'atan((lamALFA*alfarate)/(lamV*epsilon*cos(alfa)))',
                    'atan((lamALFA*alfarate)/(lamV*epsilon*cos(alfa))) + pi'
                    };
                
    elseif strcmp(char(variable),'alfadotTtrigh')
        solution = {'atan((lamALFATRIG*alfarate)/(lamH*alfamax*epsilon*cos(alfatrig)))',
                    'atan((lamALFATRIG*alfarate)/(lamH*alfamax*epsilon*cos(alfatrig))) + pi'
                    };
                
%     elseif strcmp(char(variable),'ThrottleTtrig')
%         solution = {'atan((Tmax*lamGAM*sin(alfa)/(2*mass*v) + lamV*(Tmax*cos(alfa)/(2*mass)) - Tmax*lamMASS/(2*Isp*g0))/(lamV*epsilon))',
%                     'atan((Tmax*lamGAM*sin(alfa)/(2*mass*v) + lamV*(Tmax*cos(alfa)/(2*mass)) - Tmax*lamMASS/(2*Isp*g0))/(lamV*epsilon)) + pi'
%                     };
                
    elseif strcmp(char(variable),'alfadotTtrighub')
        solution = {'atan((lamALFA*alfarate)/(lamH*alfamax*epsilon*cos(alfa)))',
                    'atan((lamALFA*alfarate)/(lamH*alfamax*epsilon*cos(alfa))) + pi'
                    };
                
    elseif strcmp(char(variable),'ThrottleTtrig')
        solution = {'atan((Tmax*lamGAM*sin(alfamax*sin(alfatrig))/(2*mass*v) + Tmax*lamV*cos(alfamax*sin(alfatrig))/(2*mass) - Tmax*lamMASS/(2*Isp*g0))/(epsilon*lamH))',
                    'atan((Tmax*lamGAM*sin(alfamax*sin(alfatrig))/(2*mass*v) + Tmax*lamV*cos(alfamax*sin(alfatrig))/(2*mass) - Tmax*lamMASS/(2*Isp*g0))/(epsilon*lamH)) + pi'
                    };
                
    elseif strcmp(char(variable),'ATRIG')
        solution = {'atan((Amax*T0*cp*lamMASS*rho0*v*(1 - T4*(Mc^2*(y/2 - 1/2) + 1)/T0 + v^2*(y/2 - 1/2)/(R*T0*y))*exp(-h/H)/(2*hpr) + Amax*lamGAM*rho0*v*(sqrt(T4*(Mc^2*(y/2 - 1/2) + 1)/(T0*(1 + v^2*(y/2 - 1/2)/(R*T0*y)))) - 1)*exp(-h/H)*sin(alfamax*sin(alfatrig))/(2*gc*mass) + Amax*lamV*rho0*v^2*(sqrt(T4*(Mc^2*(y/2 - 1/2) + 1)/(T0*(1 + v^2*(y/2 - 1/2)/(R*T0*y)))) - 1)*exp(-h/H)*cos(alfamax*sin(alfatrig))/(2*gc*mass))/(epsilon*lamH))',
                    'atan((Amax*T0*cp*lamMASS*rho0*v*(1 - T4*(Mc^2*(y/2 - 1/2) + 1)/T0 + v^2*(y/2 - 1/2)/(R*T0*y))*exp(-h/H)/(2*hpr) + Amax*lamGAM*rho0*v*(sqrt(T4*(Mc^2*(y/2 - 1/2) + 1)/(T0*(1 + v^2*(y/2 - 1/2)/(R*T0*y)))) - 1)*exp(-h/H)*sin(alfamax*sin(alfatrig))/(2*gc*mass) + Amax*lamV*rho0*v^2*(sqrt(T4*(Mc^2*(y/2 - 1/2) + 1)/(T0*(1 + v^2*(y/2 - 1/2)/(R*T0*y)))) - 1)*exp(-h/H)*cos(alfamax*sin(alfatrig))/(2*gc*mass))/(epsilon*lamH)) + pi'
                    };  
                
    elseif strcmp(char(variable),'ATRIGqb')
        solution = {'atan( (-Amax*T0*cp*lamMASS*rho0*v*(atan(qscale*(q1 - 0.5*rho0*v^2*exp(-h/H))) - atan(qscale*(q2 - 0.5*rho0*v^2*exp(-h/H))))*(1 - T4*(Mc^2*(y/2 - 1/2) + 1)/T0 + v^2*(y/2 - 1/2)/(R*T0*y))*exp(-h/H)/(2*pi*hpr) + Amax*lamGAM*rho0*v*(sqrt(((atan(qscale*(q1 - 0.5*rho0*v^2*exp(-h/H))) - atan(qscale*(q2 - 0.5*rho0*v^2*exp(-h/H))))*(1 - T4*(Mc^2*(y/2 - 1/2) + 1)/T0 + v^2*(y/2 - 1/2)/(R*T0*y))/pi + 1 + v^2*(y/2 - 1/2)/(R*T0*y))/(1 + v^2*(y/2 - 1/2)/(R*T0*y))) - 1)*exp(-h/H)*sin(alfamax*sin(alfatrig))/(2*gc*mass) + Amax*lamV*rho0*v^2*(sqrt(((atan(qscale*(q1 - 0.5*rho0*v^2*exp(-h/H))) - atan(qscale*(q2 - 0.5*rho0*v^2*exp(-h/H))))*(1 - T4*(Mc^2*(y/2 - 1/2) + 1)/T0 + v^2*(y/2 - 1/2)/(R*T0*y))/pi + 1 + v^2*(y/2 - 1/2)/(R*T0*y))/(1 + v^2*(y/2 - 1/2)/(R*T0*y))) - 1)*exp(-h/H)*cos(alfamax*sin(alfatrig))/(2*gc*mass))/(epsilon*lamH) )',
                    'atan( (-Amax*T0*cp*lamMASS*rho0*v*(atan(qscale*(q1 - 0.5*rho0*v^2*exp(-h/H))) - atan(qscale*(q2 - 0.5*rho0*v^2*exp(-h/H))))*(1 - T4*(Mc^2*(y/2 - 1/2) + 1)/T0 + v^2*(y/2 - 1/2)/(R*T0*y))*exp(-h/H)/(2*pi*hpr) + Amax*lamGAM*rho0*v*(sqrt(((atan(qscale*(q1 - 0.5*rho0*v^2*exp(-h/H))) - atan(qscale*(q2 - 0.5*rho0*v^2*exp(-h/H))))*(1 - T4*(Mc^2*(y/2 - 1/2) + 1)/T0 + v^2*(y/2 - 1/2)/(R*T0*y))/pi + 1 + v^2*(y/2 - 1/2)/(R*T0*y))/(1 + v^2*(y/2 - 1/2)/(R*T0*y))) - 1)*exp(-h/H)*sin(alfamax*sin(alfatrig))/(2*gc*mass) + Amax*lamV*rho0*v^2*(sqrt(((atan(qscale*(q1 - 0.5*rho0*v^2*exp(-h/H))) - atan(qscale*(q2 - 0.5*rho0*v^2*exp(-h/H))))*(1 - T4*(Mc^2*(y/2 - 1/2) + 1)/T0 + v^2*(y/2 - 1/2)/(R*T0*y))/pi + 1 + v^2*(y/2 - 1/2)/(R*T0*y))/(1 + v^2*(y/2 - 1/2)/(R*T0*y))) - 1)*exp(-h/H)*cos(alfamax*sin(alfatrig))/(2*gc*mass))/(epsilon*lamH) ) + pi'
                    };
    
   
      elseif strcmp(char(variable),'alfat')
         solution = {'pi/2',
                    '-pi/2',
                   'asin((0.5*(Cl1*lamGAM/(Cd2*lamV))^2-(((-(Cl1^2 + 2*Cd2*Cd0) + sqrt((Cl1^2 + 2*Cd2*Cd0)^2 - 4*(Cd2^2)*(Cd0^2 - (2*mass*gEarth*gMax*exp(h/H)/(rho0*Aref*v^2))^2)))/(2*(Cd2^2)))+((-(Cl1^2 + 2*Cd2*Cd0) - sqrt((Cl1^2 + 2*Cd2*Cd0)^2 - 4*(Cd2^2)*(Cd0^2 - (2*mass*gEarth*gMax*exp(h/H)/(rho0*Aref*v^2))^2)))/(2*(Cd2^2)))))/(((-(Cl1^2 + 2*Cd2*Cd0) + sqrt((Cl1^2 + 2*Cd2*Cd0)^2 - 4*(Cd2^2)*(Cd0^2 - (2*mass*gEarth*gMax*exp(h/H)/(rho0*Aref*v^2))^2)))/(2*(Cd2^2)))-((-(Cl1^2 + 2*Cd2*Cd0) - sqrt((Cl1^2 + 2*Cd2*Cd0)^2 - 4*(Cd2^2)*(Cd0^2 - (2*mass*gEarth*gMax*exp(h/H)/(rho0*Aref*v^2))^2)))/(2*(Cd2^2)))))'
                     };
      elseif strcmp(char(variable),'CLt')
         solution = {'pi/2',
                    '-pi/2',
                    'asin((((lamGAM*cos(pi*(1+sin(bankt))/4) + lamPSII*sin(pi*(1+sin(bankt))/4)/cos(gam))/(1.86*v*lamV))^(1/0.86)-0.5*(((0.110717 + 0.834519*(h/hs - 1) + 1.213679*(h/hs - 1)^2 - 1.060833*(h/hs - 1)^3)*(b^2*h^2/v^2)+(-0.672677 + 2.73417*(h/hs - 1) - 0.864369*(h/hs - 1)^2 - 12.1*(h/hs - 1)^3)*(b*h/v - (b^2*h^2/v^2))+(0.812241 + 2.337815*(h/hs - 1) + 10.31628*(h/hs - 1)^2 + 22.97486*(h/hs - 1)^3)*(1 - b*h/v - (b*h/v - (b^2*h^2/v^2)))+(-3.151267 - 13.62131*(h/hs - 1) - 40.4855*(h/hs - 1)^2 - 57.83333*(h/hs - 1)^3)*(v/(b*h) - 2 + b*h/v - (1 - b*h/v - (b*h/v - (b^2*h^2/v^2))))+(2.368095 + 19.0734*(h/hs - 1) + 69.86905*(h/hs - 1)^2 + 127.777778*(h/hs - 1)^3)*(v^2/(b*h)^2 - 3*v/(b*h) + 3 - b*h/v -(v/(b*h) - 2 + b*h/v - (1 - b*h/v - (b*h/v - (b^2*h^2/v^2)))))+ delCLH)+CLLB))/(0.5*(((0.110717 + 0.834519*(h/hs - 1) + 1.213679*(h/hs - 1)^2 - 1.060833*(h/hs - 1)^3)*(b^2*h^2/v^2)+(-0.672677 + 2.73417*(h/hs - 1) - 0.864369*(h/hs - 1)^2 - 12.1*(h/hs - 1)^3)*(b*h/v - (b^2*h^2/v^2))+(0.812241 + 2.337815*(h/hs - 1) + 10.31628*(h/hs - 1)^2 + 22.97486*(h/hs - 1)^3)*(1 - b*h/v - (b*h/v - (b^2*h^2/v^2)))+(-3.151267 - 13.62131*(h/hs - 1) - 40.4855*(h/hs - 1)^2 - 57.83333*(h/hs - 1)^3)*(v/(b*h) - 2 + b*h/v - (1 - b*h/v - (b*h/v - (b^2*h^2/v^2))))+(2.368095 + 19.0734*(h/hs - 1) + 69.86905*(h/hs - 1)^2 + 127.777778*(h/hs - 1)^3)*(v^2/(b*h)^2 - 3*v/(b*h) + 3 - b*h/v -(v/(b*h) - 2 + b*h/v - (1 - b*h/v - (b*h/v - (b^2*h^2/v^2)))))+ delCLH)-CLLB)))'
                   }; 
elseif strcmp(char(variable),'CLTrigNew')  
             solution = {'pi/2',
                    '-pi/2',
                    'asin((((lamGAM*cos(pi*(1+sin(bankt))/4) + lamPSII*sin(pi*(1+sin(bankt))/4)/cos(gam))/(1.86*v*lamV))^(1/0.86)+2.5-(((0.110717 + 0.834519*(h/hs - 1) + 1.213679*(h/hs - 1)^2 - 1.060833*(h/hs - 1)^3)*(b^2*h^2/v^2)+(-0.672677 + 2.73417*(h/hs - 1) - 0.864369*(h/hs - 1)^2 - 12.1*(h/hs - 1)^3)*(b*h/v - (b^2*h^2/v^2))+(0.812241 + 2.337815*(h/hs - 1) + 10.31628*(h/hs - 1)^2 + 22.97486*(h/hs - 1)^3)*(1 - b*h/v - (b*h/v - (b^2*h^2/v^2)))+(-3.151267 - 13.62131*(h/hs - 1) - 40.4855*(h/hs - 1)^2 - 57.83333*(h/hs - 1)^3)*(v/(b*h) - 2 + b*h/v - (1 - b*h/v - (b*h/v - (b^2*h^2/v^2))))+(2.368095 + 19.0734*(h/hs - 1) + 69.86905*(h/hs - 1)^2 + 127.777778*(h/hs - 1)^3)*(v^2/(b*h)^2 - 3*v/(b*h) + 3 - b*h/v -(v/(b*h) - 2 + b*h/v - (1 - b*h/v - (b*h/v - (b^2*h^2/v^2)))))+ delCLH)))/2.5)'
                     };

        elseif strcmp(char(variable),'CLg')
         CLG = '(nmax*g/(SA*exp(-h/H)*v^2))';   
         solution = {'pi/2',
                    '-pi/2',
                    ['asin(2/(',CLG,')*(((lamGAM*cos(pi*(1+sin(bankt))/4) + lamPSII*sin(pi*(1+sin(bankt))/4)/cos(gam))/(1.86*v*lamV))^(1/0.86))-1)']
                   };
        elseif strcmp(char(variable),'bankt')
          solution = {'-pi/2',
                      'pi/2',
                      'asin(4*atan(lamPSII/(lamGAM*cos(gam)))/pi-1)',
                      'asin(4*(atan(lamPSII/(lamGAM*cos(gam)))+pi)/pi-1)'
                      };
                  
      elseif strcmp(char(variable),'CLtrig')
         solution = {'pi/2',
                    '-pi/2',
                    'asin((((lamGAM*cos(pi*(1+sin(bankt))/4) + lamPSII*sin(pi*(1+sin(bankt))/4)/cos(gam))/(1.86*v*lamV))^(1/0.86)-0.2)/0.1)'
                   }; 
               
      elseif strcmp(char(variable),'CLw')
         solution = {'pi/2',
                    '-pi/2',
                    'asin((((lamGAM*cos(pi*(1+sin(bankt))/4) + lamPSII*sin(pi*(1+sin(bankt))/4)/cos(gam))/(1.86*v*lamV))^(1/0.86)-0.5*((1-w)*CLUB+w*(((0.110717 + 0.834519*(h/hs - 1) + 1.213679*(h/hs - 1)^2 - 1.060833*(h/hs - 1)^3)*(b^2*h^2/v^2)+(-0.672677 + 2.73417*(h/hs - 1) - 0.864369*(h/hs - 1)^2 - 12.1*(h/hs - 1)^3)*(b*h/v - (b^2*h^2/v^2))+(0.812241 + 2.337815*(h/hs - 1) + 10.31628*(h/hs - 1)^2 + 22.97486*(h/hs - 1)^3)*(1 - b*h/v - (b*h/v - (b^2*h^2/v^2)))+(-3.151267 - 13.62131*(h/hs - 1) - 40.4855*(h/hs - 1)^2 - 57.83333*(h/hs - 1)^3)*(v/(b*h) - 2 + b*h/v - (1 - b*h/v - (b*h/v - (b^2*h^2/v^2))))+(2.368095 + 19.0734*(h/hs - 1) + 69.86905*(h/hs - 1)^2 + 127.777778*(h/hs - 1)^3)*(v^2/(b*h)^2 - 3*v/(b*h) + 3 - b*h/v -(v/(b*h) - 2 + b*h/v - (1 - b*h/v - (b*h/v - (b^2*h^2/v^2)))))+ delCLH))+CLLB))/(0.5*((1-w)*CLUB+w*(((0.110717 + 0.834519*(h/hs - 1) + 1.213679*(h/hs - 1)^2 - 1.060833*(h/hs - 1)^3)*(b^2*h^2/v^2)+(-0.672677 + 2.73417*(h/hs - 1) - 0.864369*(h/hs - 1)^2 - 12.1*(h/hs - 1)^3)*(b*h/v - (b^2*h^2/v^2))+(0.812241 + 2.337815*(h/hs - 1) + 10.31628*(h/hs - 1)^2 + 22.97486*(h/hs - 1)^3)*(1 - b*h/v - (b*h/v - (b^2*h^2/v^2)))+(-3.151267 - 13.62131*(h/hs - 1) - 40.4855*(h/hs - 1)^2 - 57.83333*(h/hs - 1)^3)*(v/(b*h) - 2 + b*h/v - (1 - b*h/v - (b*h/v - (b^2*h^2/v^2))))+(2.368095 + 19.0734*(h/hs - 1) + 69.86905*(h/hs - 1)^2 + 127.777778*(h/hs - 1)^3)*(v^2/(b*h)^2 - 3*v/(b*h) + 3 - b*h/v -(v/(b*h) - 2 + b*h/v - (1 - b*h/v - (b*h/v - (b^2*h^2/v^2)))))+ delCLH))-CLLB)))'
                   }; 


      elseif strcmp(char(variable),'utrig')
         solution = {'pi/2',
                    '-pi/2',
                    'asin(1-3*lamY/8)'
                   };   
               
elseif strcmp(char(variable),'weccontrol')
         solution = {'atan((-gam*(x2 + lamX2/m))/epsilon)',
                    'atan((-gam*(x2 + lamX2/m))/epsilon) + pi'
                   };

elseif strcmp(char(variable),'distcontrol')
         solution = {'atan(((lamX32*(0.8*x31 - (0.8*vol*x32)/(x32*(vol - 1) + 1)))/areb - (lamX17*(0.8*x17 - 0.8*x16 + 0.8*((vol*x17)/(x17*(vol - 1) + 1) - (vol*x18)/(x18*(vol - 1) + 1))))/atray + (lamX2*(0.8*(x1 - x2) - 0.8*((vol*x2)/(x2*(vol - 1) + 1) - (vol*x3)/(x3*(vol - 1) + 1))))/atray + (lamX3*(0.8*(x2 - x3) - 0.8*((vol*x3)/(x3*(vol - 1) + 1) - (vol*x4)/(x4*(vol - 1) + 1))))/atray + (lamX4*(0.8*(x3 - x4) - 0.8*((vol*x4)/(x4*(vol - 1) + 1) - (vol*x5)/(x5*(vol - 1) + 1))))/atray + (lamX5*(0.8*(x4 - x5) - 0.8*((vol*x5)/(x5*(vol - 1) + 1) - (vol*x6)/(x6*(vol - 1) + 1))))/atray + (lamX6*(0.8*(x5 - x6) - 0.8*((vol*x6)/(x6*(vol - 1) + 1) - (vol*x7)/(x7*(vol - 1) + 1))))/atray + (lamX7*(0.8*(x6 - x7) - 0.8*((vol*x7)/(x7*(vol - 1) + 1) - (vol*x8)/(x8*(vol - 1) + 1))))/atray + (lamX8*(0.8*(x7 - x8) - 0.8*((vol*x8)/(x8*(vol - 1) + 1) - (vol*x9)/(x9*(vol - 1) + 1))))/atray + (lamX9*(0.8*(x8 - x9) - 0.8*((vol*x9)/(x9*(vol - 1) + 1) - (vol*x10)/(x10*(vol - 1) + 1))))/atray + (lamX10*(0.8*(x9 - x10) - 0.8*((vol*x10)/(x10*(vol - 1) + 1) - (vol*x11)/(x11*(vol - 1) + 1))))/atray + (lamX11*(0.8*(x10 - x11) - 0.8*((vol*x11)/(x11*(vol - 1) + 1) - (vol*x12)/(x12*(vol - 1) + 1))))/atray + (lamX12*(0.8*(x11 - x12) - 0.8*((vol*x12)/(x12*(vol - 1) + 1) - (vol*x13)/(x13*(vol - 1) + 1))))/atray + (lamX13*(0.8*(x12 - x13) - 0.8*((vol*x13)/(x13*(vol - 1) + 1) - (vol*x14)/(x14*(vol - 1) + 1))))/atray + (lamX14*(0.8*(x13 - x14) - 0.8*((vol*x14)/(x14*(vol - 1) + 1) - (vol*x15)/(x15*(vol - 1) + 1))))/atray + (lamX15*(0.8*(x14 - x15) - 0.8*((vol*x15)/(x15*(vol - 1) + 1) - (vol*x16)/(x16*(vol - 1) + 1))))/atray + (lamX16*(0.8*(x15 - x16) - 0.8*((vol*x16)/(x16*(vol - 1) + 1) - (vol*x17)/(x17*(vol - 1) + 1))))/atray + (lamX18*(0.8*(x17 - x18) - 0.8*((vol*x18)/(x18*(vol - 1) + 1) - (vol*x19)/(x19*(vol - 1) + 1))))/atray + (lamX19*(0.8*(x18 - x19) - 0.8*((vol*x19)/(x19*(vol - 1) + 1) - (vol*x20)/(x20*(vol - 1) + 1))))/atray + (lamX20*(0.8*(x19 - x20) - 0.8*((vol*x20)/(x20*(vol - 1) + 1) - (vol*x21)/(x21*(vol - 1) + 1))))/atray + (lamX21*(0.8*(x20 - x21) - 0.8*((vol*x21)/(x21*(vol - 1) + 1) - (vol*x22)/(x22*(vol - 1) + 1))))/atray + (lamX22*(0.8*(x21 - x22) - 0.8*((vol*x22)/(x22*(vol - 1) + 1) - (vol*x23)/(x23*(vol - 1) + 1))))/atray + (lamX23*(0.8*(x22 - x23) - 0.8*((vol*x23)/(x23*(vol - 1) + 1) - (vol*x24)/(x24*(vol - 1) + 1))))/atray + (lamX24*(0.8*(x23 - x24) - 0.8*((vol*x24)/(x24*(vol - 1) + 1) - (vol*x25)/(x25*(vol - 1) + 1))))/atray + (lamX25*(0.8*(x24 - x25) - 0.8*((vol*x25)/(x25*(vol - 1) + 1) - (vol*x26)/(x26*(vol - 1) + 1))))/atray + (lamX26*(0.8*(x25 - x26) - 0.8*((vol*x26)/(x26*(vol - 1) + 1) - (vol*x27)/(x27*(vol - 1) + 1))))/atray + (lamX27*(0.8*(x26 - x27) - 0.8*((vol*x27)/(x27*(vol - 1) + 1) - (vol*x28)/(x28*(vol - 1) + 1))))/atray + (lamX28*(0.8*(x27 - x28) - 0.8*((vol*x28)/(x28*(vol - 1) + 1) - (vol*x29)/(x29*(vol - 1) + 1))))/atray + (lamX29*(0.8*(x28 - x29) - 0.8*((vol*x29)/(x29*(vol - 1) + 1) - (vol*x30)/(x30*(vol - 1) + 1))))/atray + (lamX30*(0.8*(x29 - x30) - 0.8*((vol*x30)/(x30*(vol - 1) + 1) - (vol*x31)/(x31*(vol - 1) + 1))))/atray + (lamX31*(0.8*(x30 - x31) - 0.8*((vol*x31)/(x31*(vol - 1) + 1) - (vol*x32)/(x32*(vol - 1) + 1))))/atray - (0.8*lamX1*(x1 - (vol*x2)/(x2*(vol - 1) + 1)))/acond)/epsilon)',
                    'atan(((lamX32*(0.8*x31 - (0.8*vol*x32)/(x32*(vol - 1) + 1)))/areb - (lamX17*(0.8*x17 - 0.8*x16 + 0.8*((vol*x17)/(x17*(vol - 1) + 1) - (vol*x18)/(x18*(vol - 1) + 1))))/atray + (lamX2*(0.8*(x1 - x2) - 0.8*((vol*x2)/(x2*(vol - 1) + 1) - (vol*x3)/(x3*(vol - 1) + 1))))/atray + (lamX3*(0.8*(x2 - x3) - 0.8*((vol*x3)/(x3*(vol - 1) + 1) - (vol*x4)/(x4*(vol - 1) + 1))))/atray + (lamX4*(0.8*(x3 - x4) - 0.8*((vol*x4)/(x4*(vol - 1) + 1) - (vol*x5)/(x5*(vol - 1) + 1))))/atray + (lamX5*(0.8*(x4 - x5) - 0.8*((vol*x5)/(x5*(vol - 1) + 1) - (vol*x6)/(x6*(vol - 1) + 1))))/atray + (lamX6*(0.8*(x5 - x6) - 0.8*((vol*x6)/(x6*(vol - 1) + 1) - (vol*x7)/(x7*(vol - 1) + 1))))/atray + (lamX7*(0.8*(x6 - x7) - 0.8*((vol*x7)/(x7*(vol - 1) + 1) - (vol*x8)/(x8*(vol - 1) + 1))))/atray + (lamX8*(0.8*(x7 - x8) - 0.8*((vol*x8)/(x8*(vol - 1) + 1) - (vol*x9)/(x9*(vol - 1) + 1))))/atray + (lamX9*(0.8*(x8 - x9) - 0.8*((vol*x9)/(x9*(vol - 1) + 1) - (vol*x10)/(x10*(vol - 1) + 1))))/atray + (lamX10*(0.8*(x9 - x10) - 0.8*((vol*x10)/(x10*(vol - 1) + 1) - (vol*x11)/(x11*(vol - 1) + 1))))/atray + (lamX11*(0.8*(x10 - x11) - 0.8*((vol*x11)/(x11*(vol - 1) + 1) - (vol*x12)/(x12*(vol - 1) + 1))))/atray + (lamX12*(0.8*(x11 - x12) - 0.8*((vol*x12)/(x12*(vol - 1) + 1) - (vol*x13)/(x13*(vol - 1) + 1))))/atray + (lamX13*(0.8*(x12 - x13) - 0.8*((vol*x13)/(x13*(vol - 1) + 1) - (vol*x14)/(x14*(vol - 1) + 1))))/atray + (lamX14*(0.8*(x13 - x14) - 0.8*((vol*x14)/(x14*(vol - 1) + 1) - (vol*x15)/(x15*(vol - 1) + 1))))/atray + (lamX15*(0.8*(x14 - x15) - 0.8*((vol*x15)/(x15*(vol - 1) + 1) - (vol*x16)/(x16*(vol - 1) + 1))))/atray + (lamX16*(0.8*(x15 - x16) - 0.8*((vol*x16)/(x16*(vol - 1) + 1) - (vol*x17)/(x17*(vol - 1) + 1))))/atray + (lamX18*(0.8*(x17 - x18) - 0.8*((vol*x18)/(x18*(vol - 1) + 1) - (vol*x19)/(x19*(vol - 1) + 1))))/atray + (lamX19*(0.8*(x18 - x19) - 0.8*((vol*x19)/(x19*(vol - 1) + 1) - (vol*x20)/(x20*(vol - 1) + 1))))/atray + (lamX20*(0.8*(x19 - x20) - 0.8*((vol*x20)/(x20*(vol - 1) + 1) - (vol*x21)/(x21*(vol - 1) + 1))))/atray + (lamX21*(0.8*(x20 - x21) - 0.8*((vol*x21)/(x21*(vol - 1) + 1) - (vol*x22)/(x22*(vol - 1) + 1))))/atray + (lamX22*(0.8*(x21 - x22) - 0.8*((vol*x22)/(x22*(vol - 1) + 1) - (vol*x23)/(x23*(vol - 1) + 1))))/atray + (lamX23*(0.8*(x22 - x23) - 0.8*((vol*x23)/(x23*(vol - 1) + 1) - (vol*x24)/(x24*(vol - 1) + 1))))/atray + (lamX24*(0.8*(x23 - x24) - 0.8*((vol*x24)/(x24*(vol - 1) + 1) - (vol*x25)/(x25*(vol - 1) + 1))))/atray + (lamX25*(0.8*(x24 - x25) - 0.8*((vol*x25)/(x25*(vol - 1) + 1) - (vol*x26)/(x26*(vol - 1) + 1))))/atray + (lamX26*(0.8*(x25 - x26) - 0.8*((vol*x26)/(x26*(vol - 1) + 1) - (vol*x27)/(x27*(vol - 1) + 1))))/atray + (lamX27*(0.8*(x26 - x27) - 0.8*((vol*x27)/(x27*(vol - 1) + 1) - (vol*x28)/(x28*(vol - 1) + 1))))/atray + (lamX28*(0.8*(x27 - x28) - 0.8*((vol*x28)/(x28*(vol - 1) + 1) - (vol*x29)/(x29*(vol - 1) + 1))))/atray + (lamX29*(0.8*(x28 - x29) - 0.8*((vol*x29)/(x29*(vol - 1) + 1) - (vol*x30)/(x30*(vol - 1) + 1))))/atray + (lamX30*(0.8*(x29 - x30) - 0.8*((vol*x30)/(x30*(vol - 1) + 1) - (vol*x31)/(x31*(vol - 1) + 1))))/atray + (lamX31*(0.8*(x30 - x31) - 0.8*((vol*x31)/(x31*(vol - 1) + 1) - (vol*x32)/(x32*(vol - 1) + 1))))/atray - (0.8*lamX1*(x1 - (vol*x2)/(x2*(vol - 1) + 1)))/acond)/epsilon) + pi'
                   };
               
 
elseif strcmp(char(variable),'u1')
         solution = {'atan(0.5*(lamX1-lamX2*x2/x1+lamX3*(c1-x3)/x1-lamX4*x4/x1-lamX5*x5/x1)/epsilon1)',
                     'atan(0.5*(lamX1-lamX2*x2/x1+lamX3*(c1-x3)/x1-lamX4*x4/x1-lamX5*x5/x1)/epsilon1) + pi'
                    };
              
elseif strcmp(char(variable),'u2')
         solution = {'atan(0.5*(lamX1-lamX2*x2/x1-lamX3*x3/x1-lamX4*x4/x1+lamX5*(c3-x5)/x1)/epsilon2)',
                     'atan(0.5*(lamX1-lamX2*x2/x1-lamX3*x3/x1-lamX4*x4/x1+lamX5*(c3-x5)/x1)/epsilon2) + pi'
                    };
%       elseif strcmp(char(variable),'alfagtrig')   
%          solution = {'pi/2',
%                     '-pi/2',
%                     'asin((lamGAM*Cl/(2*Cd2*v*lamV) - (0.5*w*((sqrt((-(Cl^2 + 2*Cd2*Cd0)-sqrt((Cl^2 + 2*Cd2*Cd0)^2-4*(Cd2^2)*(Cd0^2 - (2*mass*g*gmax/((rho0*exp(-h/H))*v^2*Aref))^2)))/(2*(Cd2^2))))-(-sqrt((-(Cl^2 + 2*Cd2*Cd0)+sqrt((Cl^2 + 2*Cd2*Cd0)^2-4*(Cd2^2)*(Cd0^2 - (2*mass*g*gmax/((rho0*exp(-h/H))*v^2*Aref))^2)))/(2*(Cd2^2)))))))/((0.5*w*((sqrt((-(Cl^2 + 2*Cd2*Cd0)-sqrt((Cl^2 + 2*Cd2*Cd0)^2-4*(Cd2^2)*(Cd0^2 - (2*mass*g*gmax/((rho0*exp(-h/H))*v^2*Aref))^2)))/(2*(Cd2^2))))+(-sqrt((-(Cl^2 + 2*Cd2*Cd0)+sqrt((Cl^2 + 2*Cd2*Cd0)^2-4*(Cd2^2)*(Cd0^2 - (2*mass*g*gmax/((rho0*exp(-h/H))*v^2*Aref))^2)))/(2*(Cd2^2)))))+20*pi/180*(1-w))))'
%                     };

%       elseif strcmp(char(variable),'alfagtrig')   
%          solution = {'pi/2',
%                     '-pi/2',
%                     'asin((lamGAM*Cl/(2*Cd2*v*lamV))/real((sqrt((-(Cl1^2 + 2*Cd2*Cd0)+(sqrt((Cl1^2 + 2*Cd2*Cd0)^2-4*(Cd2^2)*(Cd0^2 - (2*mass*g*gmax/((rho0*exp(-h/H))*v^2*Aref))^2)))))/(2*(Cd2^2))))))'
%                     };
                
    else
        solution = mathematica_solve(expression,value,variable,assumptions);
        solution
       % keyboard;
    end
	% Check first solution (in case multiple solutions exist for Mathematica inverse function of unknown, user-supplied function)
	if isempty(strfind(char(solution(1)),'InverseFunction')) && isempty(strfind(char(solution(1)),'Solve')) % Analytic solution found, finish
		
    if length(char(solution(1))) > 1
      % solution not singular
      return;
    else % solution singular
      
    end
		
	else % Analytic solution not found, form numerical root solving expression
			
		% solution{1} = ['fzero(@(',char(variable),') ',char(expression),'-',char(value),',',num2str(initialGuess),')'];
		
		% Determine range for interval search. Start by searching for a function that uses variable as input.
		indFunc = 0;
		indInput = 0;
		
		endSearch = false;
		
    if isfield(in,'customFunc')

      while true

        indFunc = indFunc + 1;
        if ~isempty(strfind(char(expression),in.customFunc(indFunc).name{1}))
          % Find index of input variable
          while true

            indInput = indInput + 1;
            if ~isempty(strcmp(in.customFunc(indFunc).input{indInput}{1},char(variable)))

              % Variable found, get range for interval search
              interval = in.customFunc(indFunc).input{indInput}{3};
              break;

            else % Do nothing, loop

            end

          end

          break;

        else % Do nothing, loop

        end

      end

    else
      
      interval = [-48*pi/180 48*pi/180]; %%% NEED TO GENERALIZE
      
    end
		
		% Write first part of solver line
		solution{1} = ['fzero(''',char(variable),'SolutionFunction'',',num2str(mean(interval)),',',in.fzeroOptions];
		% solution{1} = ['fzero(@',char(variable),'SolutionFunction,[',num2str(interval),'],',in.fzeroOptions];
		% solution{1} = ['fzero(@',char(variable),'SolutionFunction,0*pi/180,',in.fzeroOptions];
		% solution{1} = ['zeroByLooping(@',char(variable),'SolutionFunction,[',num2str(interval),']'];
		
		% Write remaining custom function to be passed into the function
		extraInput = [];
		% for ctrInput = 2 : 1 : length(customFunc(indFunc).input)
		% 	extraInput = [extraInput,',',char(customFunc(indFunc).input(ctrInput){1})];
		% end
		
		% Write general variables to pass to root solving function
% 		extraInput = [extraInput,',xAndLambda,const,constraint,bank']; %%%%%%% FIX HARDCODED CONTROL!!!!!!!!!!!!!
    extraInput = [extraInput,',xAndLambda,const,constraint']; %%%%%%% FIX HARDCODED CONTROL!!!!!!!!!!!!!
		
		% Write full solution expression
		solution{1} = [solution{1},extraInput,')'];
		
		% Write function that fzero will call
		fid = fopen([in.autocodeDirTraj,'/',char(variable),'SolutionFunction.m'],'w');
		
		% Header
		fprintf(fid,['function [zeroVal] = ',char(variable),'SolutionFunction(',char(variable),extraInput,')\n']);
		% fprintf(fid,'coder.extrinsic(''keyboard'');\n'); %%%%%%%%%%%%%% DELETE LATER!!!!!!!!!!!!!!!!!!!!!!
		% Assign input values
		if isfield(in,'const')
			writeConstants(fid,in.const,false);
			fprintf(fid,'\n');
		end
	
		if isfield(in,'constraintVal')
			writeConstraints(fid,in.constraintVal,false);
		end
		fprintf(fid,'\n');
	
		% Use real values since these will be fed into the user custom functions. Some Matlab algorithms require real numbers
		% (e.g., interpolation routines).
	  fprintf(fid,'\n');
	  fprintf(fid,'%% States\n');
	  for ctrState = 1 : 1 : oc.num.states
	    fprintf(fid,[char(oc.state.var(ctrState)),' = real(xAndLambda(',int2str(ctrState),',1));\n']);
	  end
	  fprintf(fid,'\n');

	  fprintf(fid,'%% Costates\n');
	  for ctrCostate = 1 : 1 : oc.num.costates
	    fprintf(fid,[char(oc.costate.var(ctrCostate)),' = real(xAndLambda(',int2str(ctrCostate+oc.num.states),',1));\n']);
	  end
	  fprintf(fid,'\n');
	
		% Write derived quantities
		if isfield(in,'quantity')
			writeQuantities(fid,in.quantity);
			fprintf(fid,'\n');
		end
		
		% Write root-solving expression
		fprintf(fid,['zeroVal = real(',char(expression),'-',char(value),');\n']); % real expression required by fzero
		
		% fprintf(fid,'if zeroVal == inf\n');%abs(real(zeroVal)) > 1e8 || abs(imag(zeroVal)) > 0\n');
		% fprintf(fid,'  xAndLambda\n');
		% fprintf(fid,'  alfa\n');
		% fprintf(fid,'  bank\n');
		% fprintf(fid,'  zeroVal = zeroVal\n');
		% 
		% % fprintf(fid,'  alfaSet = linspace(-10*pi/180,10*pi/180,100);\n');
		% % fprintf(fid,'  zeroValSet = NaN(size(alfaSet));\n');
		% % fprintf(fid,'  for ctr = 1 : 1 : length(alfaSet)\n');
		% % fprintf(fid,'    alfa = alfaSet(ctr);\n');
		% % fprintf(fid,['    zeroValSet(ctr) = ',char(expression),'-',char(value),';\n']);
		% % fprintf(fid,'  end\n');
		% % fprintf(fid,'plot(alfaSet,zeroValSet);\n');
		% % 
		% % fprintf(fid,'  keyboard\n');
		% fprintf(fid,'end\n');
		
		% fprintf(fid,'if abs(real(zeroVal)) < 1e-3\n');
		% fprintf(fid,'  alfa\n');
		% fprintf(fid,'  bank\n');
		% fprintf(fid,'  zeroVal = zeroVal\n');
		% fprintf(fid,'  keyboard\n');
		% fprintf(fid,'end\n');
		
		fprintf(fid,'\n');
		fprintf(fid,'return\n');
		fprintf(fid,'\n');

		% Close file
		fclose(fid);
		
	end
	
return

