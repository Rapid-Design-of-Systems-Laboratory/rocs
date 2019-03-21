function [oc2] = singularArcCalc(in,oc)
	if sum(in.oc.sequence.singularArc) > 0 || sum(in.oc.switch.enforced) > 0
	    disp('I am inside singular thing')
	    oc.singular.H = subs(oc.H.unconstrained,in.oc.u,sym(0));
	    oc.singular.dHudt(1,1) = diff(oc.singular.H,sym('Cl'));
    
	    % Take derivatives until control found
	    index = 1;
	    while true
        
	        if isempty(strfind(char(oc.singular.dHudt(index,1)),char(in.oc.u)))
	            % Control not found yet, take another derivative
	            oc.singular.dHudt(index+1,1) = jacobian(oc.singular.dHudt(index,1),[oc.x; oc.lambda])* ...
	                [oc.xDot; oc.lambdaDot.unconstrained];
	            index = index + 1;
	        else
	            % Control found, exit
	            break;
	        end
        
	    end
	    % This part is kind of hard coded, make it general (remove bank)
	    %   bankExpression = simplify(solve(oc.singular.dHudt(end,1) == 0,in.oc.u));
	    oc.bankExpression = solve(oc.singular.dHudt(end,1) == 0,in.oc.u);
	    oc.bankExpressionChoose = oc.bankExpression(1);
	    %   bankCosArg = simplify(cos(bankExpressionChoose));
	    oc.bankCosArg = cos(oc.bankExpressionChoose);
	end
	oc2 = oc;
end