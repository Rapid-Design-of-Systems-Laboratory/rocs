function [derivative] = formDerivative(expression,variables)
	%
	% This function determines if an analytic derivative can be performed on the expression with respect to variables.
	% An analytic derivative enables multiple solution branches to be propagated forward.
	% If an analytic derivative is not possible, then a derivative is provided that uses complex step.
	%
	
% 	epsilon = 1e-30; % step size
	epsilon = 1e-6; % step size
	
	derivative = sym(NaN(1,length(variables)));
	
	for ctrVariable = 1 : 1 : length(variables)
	
		% Determine if an analytic derivative can be performed
		derivative(ctrVariable) = jacobian(expression,variables(ctrVariable))
	
		if isempty(strfind(char(derivative(ctrVariable)),'D(')) && ...
        isempty(strfind(char(derivative(ctrVariable)),'[1]')) && ...
        isempty(strfind(char(derivative(ctrVariable)),'diff')) % Analytic derivative found, finish. Will usually find [1] in Mathematica numerical D[] result.
		
			continue;
		
		else % Analytic derivative not found, form numerical derivative
			
% 			% Complex step derivative does not work when performing table interpolation, etc.
% 			derivative(ctrVariable) = imag(subs(expression,variables(ctrVariable),variables(ctrVariable)+epsilon*1i)/epsilon);
			
			% Perform forward finite difference
			derivative(ctrVariable) = simplify((subs(expression,variables(ctrVariable),variables(ctrVariable)+epsilon)-expression)/epsilon);
			
			% % Perform central finite difference
% 			derivative(ctrVariable) = simplify((subs(expression,variables(ctrVariable),variables(ctrVariable)+epsilon)-subs(expression,variables(ctrVariable),variables(ctrVariable)-epsilon))/(2*epsilon));
			
		end
	
	end

return

