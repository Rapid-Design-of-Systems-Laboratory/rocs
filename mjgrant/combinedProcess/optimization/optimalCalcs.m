function [oc,in_] = optimalCalcs(in)
%
% This function calculates the various necessary conditions associated with a
% given optimal control problem. At the moment, the optimal control problem
% inputs are provided within this function. The optimal control problem should
% be input to the function in the future.
%
% input : in and out [structures]
%         in is obtained from inputsRunTrajectoryProcess.m
%         out is obtained from initialCalculations.m
% output : oc [structure]
% Developed by : Dr. M.J. Grant, Kshitij Mall and Thomas Antony
% Last modified : 21 Feb, 2014


% Initialize in case don't want to run all of optimalCalcs. Need to reorganize
% required calculations even if symbolic math doesn't change (so user doesn't
% have to rederive necessary conditions of optimality)
oc.num.lagrangeMultipliers.interiorPoint = [];

%%%%%%%%%%%%%%%%%%%%%%%%
%% State Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Number of states
oc.num.states = length(in.oc.state(:,1));
oc.num.origStates = oc.num.states;

oc.independentVariable = sym((in.oc.independentVariable(1)));

% Assign symbolic variables for state variables and associated expressions.
% Note, reserved function expressions will cause an error.

oc.state.var = sym(in.oc.state(:,1)); % States

oc.state.rate = sym(in.oc.stateRate(:,1)); % Process equations

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Costate Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

oc.num.costates = oc.num.states;

% Initialize costate symbolic variables
oc.costate.var = sym(nan(oc.num.costates,1));
for ctrCostate = 1 : 1 : oc.num.costates
    oc.costate.var(ctrCostate,1) = sym(['lam',upper(char(oc.state.var(ctrCostate,1)))]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Control Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize control symbolic variables
oc.control.var = sym(in.oc.control(:,1));
oc.num.controls = length(oc.control.var);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Path Cost Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign path cost symbolic variables
oc.cost.path = sym(in.oc.cost.path{1,1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Path Constraint Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine number of path constraints
if isfield(in.oc.constraint,'path')
    oc.num.constraints.path = length(in.oc.constraint.path(:,1));
else
    oc.num.constraints.path = 0;
end

% Scale constraints such that S/val - 1 <= 0 is used rather than S - val <= 0 and vice-versa.
for ctrPathConstraint = 1 : 1 : oc.num.constraints.path
    if strcmp(in.oc.constraint.path{ctrPathConstraint,2},'>')
        % oc.constraint.path{ctrPathConstraint,1} = sym(in.oc.constraint.path{ctrPathConstraint,3})/sym(in.oc.constraint.path{ctrPathConstraint,1}) - sym(1);
        oc.constraint.path{ctrPathConstraint,1} = sym(in.oc.constraint.path{ctrPathConstraint,3}) - sym(in.oc.constraint.path{ctrPathConstraint,1});
    elseif strcmp(in.oc.constraint.path{ctrPathConstraint,2},'<')
        % oc.constraint.path{ctrPathConstraint,1} = sym(1) - sym(in.oc.constraint.path{ctrPathConstraint,3})/sym(in.oc.constraint.path{ctrPathConstraint,1});
        oc.constraint.path{ctrPathConstraint,1} = sym(in.oc.constraint.path{ctrPathConstraint,1}) - sym(in.oc.constraint.path{ctrPathConstraint,3});
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interior Point Constraint Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine number of interior point constraints
if isfield(in.oc.constraint,'interiorPoint')
    oc.num.constraints.interiorPoint = length(in.oc.constraint.interiorPoint(:,1));
else
    oc.num.constraints.interiorPoint = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Point Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial cost
oc.cost.augmented.initial = sym(in.oc.cost.initial{1,1});

% Determine number of initial constraints
if isfield(in.oc.constraint,'initial')
    oc.num.constraints.initial = length(in.oc.constraint.initial(:,1));
else
    oc.num.constraints.initial = 0;
end

% Initial constraint symbolic variable assignment
oc.constraint.initial = sym(in.oc.constraint.initial(:,1));

% Compute augmented initial cost. augmented initial cost = initial cost + Lagrange multiplier*initial constraint
oc.num.lagrangeMultiplier.initial = oc.num.constraints.initial;

for ctrInitialConstraint = 1 : 1 : oc.num.lagrangeMultiplier.initial
    oc.cost.augmented.initial = oc.cost.augmented.initial + ...
        sym(['lmInitial',int2str(ctrInitialConstraint)])*oc.constraint.initial(ctrInitialConstraint);
end

% Compute initial costate conditions
oc.costate.initial = formDerivative(oc.cost.augmented.initial,oc.state.var).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Terminal Point Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Terminal cost
oc.cost.augmented.terminal = sym(in.oc.cost.terminal{1});

% Determine number of terminal constraints
if isfield(in.oc.constraint,'terminal')
    oc.num.constraints.terminal = length(in.oc.constraint.terminal(:,1));
    % Terminal constraint
    oc.constraint.terminal = sym(in.oc.constraint.terminal(:,1));
    % Compute augmented terminal cost. augmented terminal cost = terminal cost + Lagrange multiplier*terminal constraint
    
else
    oc.num.constraints.terminal = 0;
end
oc.num.lagrangeMultiplier.terminal = oc.num.constraints.terminal;
for ctrTerminalConstraint = 1 : 1 : oc.num.lagrangeMultiplier.terminal
    oc.cost.augmented.terminal = oc.cost.augmented.terminal + ...
        sym(['lmTerminal',int2str(ctrTerminalConstraint)])*oc.constraint.terminal(ctrTerminalConstraint);
end

% Compute final lambda conditions
oc.costate.terminal = formDerivative(oc.cost.augmented.terminal,oc.state.var).';

oc.num.parameters = oc.num.constraints.initial + oc.num.constraints.terminal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stop Execution if Not Wanting to Run optimalCalcs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The lines below can be computationally intensive)

if ~in.run.optimalCalcs
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unconstrained Arc Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hamiltonian
oc.hamiltonian.unconstrained.expression = oc.cost.path + oc.costate.var.'*oc.state.rate;

fid_ctrl_c = fopen([in.autocodeDirTraj,'/computeControl.h'],'w');
% Control law has to be computed only if the files are being written
if in.oc.writeEquations
    
    % Optimal control law
    singularArc = false*ones(1,oc.num.controls);
    oc.hamiltonian.unconstrained.controlPartial = formDerivative(oc.hamiltonian.unconstrained.expression,oc.control.var);
    % keyboard
 
    for ctrControl = 1 : 1 : oc.num.controls
%        oc.control.unconstrained.expression{ctrControl} = ...
%        		solve([char(oc.hamiltonian.unconstrained.controlPartial(ctrControl)),'=0'],char(oc.control.var(ctrControl)));
        % NEED TO AUTOMATE TRANSITION TO COUPLED CONTROL EQUATIONS!!!!!!!!!!!!!!!!!!!!!
      oc.control.unconstrained.expression{ctrControl} = formSolver(oc.hamiltonian.unconstrained.controlPartial(ctrControl), ...
        sym('0'),oc.control.var(ctrControl),in.oc.assumptions.control,in,oc);
%         oc.control.unconstrained.expression{ctrControl} = formSolver(oc.hamiltonian.unconstrained.expression, ...
%             sym('0'),oc.control.var(ctrControl),in.oc.assumptions.control,in,oc);
        % keyboard
        % Add protection when answer is a numeric value (only one solution returned by
        % Mathematica)
        try
            if ctrControl == 1 && isnumeric(eval(oc.control.unconstrained.expression{ctrControl}{1}))
                oc.control.unconstrained.expression{ctrControl}{1,1} = '-pi/2'; % Changed from 0
                oc.control.unconstrained.expression{ctrControl}{2,1} = 'pi/2'; % Changed from pi
            end
        catch
        end
        
        % Check for bang/bang, singular arc
        if isempty(oc.control.unconstrained.expression{ctrControl}{1})
            singularArc(ctrControl) = true;
        end
        
        for ctrCoefficientSet = 1 : 1 : length(oc.control.unconstrained.expression{ctrControl})
            unconstrainedCoefficients{ctrControl}{ctrCoefficientSet} = [];
        end
    end
    
    % Reconfigure problem with bang/bang, singular arc capability via smoothing
    % (Silva and Trelat)
    singularArcEquation = cell(1,oc.num.controls);
    if sum(singularArc) > 0
        
        extraControlIndex = [];
        
        % Determine EOMs without control
        for ctrState = 1 : 1 : oc.num.states
            
            controlFound = false;
            
            for ctrControl = 1 : 1 : oc.num.controls
                if ~isempty(strfind(char(oc.state.rate(ctrState)),char(oc.control.var(ctrControl))))
                    %         controlFound = true;
                    if singularArc(ctrControl)
                        singularArcEquation{ctrControl} = [singularArcEquation{ctrControl}, ctrState];
                        controlFound = true;
                    end
                end
            end
            
            if ~controlFound % need to add extra control
                extraControlIndex = [extraControlIndex, ctrState];
            end
            
        end
        
        % Change EOMs
        controlTerm = sym('0');
        for ctrExtraControl = 1 : 1 : length(extraControlIndex)
            
            indState = extraControlIndex(ctrExtraControl);
            oc.state.rate(indState) = oc.state.rate(indState) + ...
                sym(['epsilon*uEpsilon',int2str(indState)]);
            
            % Used in modified control
            controlTerm = controlTerm + sym('epsilon')^2*oc.costate.var(indState,1)^2;
            
        end
        
        % Change control equation for singular control - NOTE: negative sign due to Trelat's paper using
        % Pontryagin's MAXIMUM Principle
        for ctrControl = 1 : 1 : oc.num.controls
            if singularArc(ctrControl)
                oc.control.unconstrained.expression{ctrControl} = ...
                    {char(-oc.hamiltonian.unconstrained.controlPartial(ctrControl)/ ...
                    sqrt(oc.hamiltonian.unconstrained.controlPartial(ctrControl)^2 + ...
                    controlTerm))};
                singularControlIndex = ctrControl;
                unconstrainedCoefficients{ctrControl}{ctrCoefficientSet} = [];
            end
        end
        
        % Add control expressions for extra controls
        for ctrExtraControl = 1 : 1 : length(extraControlIndex)
            
            indState = extraControlIndex(ctrExtraControl);
            
            % Change controls - NOTE: negative sign due to Trelat's paper using
            % Pontryagin's MAXIMUM Principle
            oc.control.unconstrained.expression{oc.num.controls+ctrExtraControl} = ...
                {char(sym('-epsilon')*oc.costate.var(indState,1)/ ...
                sqrt(oc.hamiltonian.unconstrained.controlPartial(singularControlIndex)^2 + ...
                controlTerm))};
            
            oc.control.var(oc.num.controls+ctrExtraControl) = sym(['uEpsilon',int2str(indState)]);
            
            unconstrainedCoefficients{oc.num.controls+ctrExtraControl}{ctrCoefficientSet} = [];
            
        end
        
        oc.num.controls = length(oc.control.var);
        
        % Change Hamiltonian to include extra controls
        oc.hamiltonian.unconstrained.expression = oc.cost.path + oc.costate.var.'*oc.state.rate;
        
    end
    
    % Determine control write order to function. Select expression with fewest control variables first.
    controlsInExpression = zeros(1,oc.num.controls);
    for ctrExpression = 1 : 1 : oc.num.controls
        for ctrControl = 1 : 1 : oc.num.controls
            
            if ~isempty(strfind(char(oc.control.unconstrained.expression{ctrExpression}(1)),char(oc.control.var(ctrControl))))
                controlsInExpression(ctrExpression) = controlsInExpression(ctrExpression) + 1;
            end
            
        end
    end
    
    % Sort controls starting with those with fewest appearances in control equations
    [Y,oc.control.unconstrained.writeOrder] = sort(controlsInExpression);
    
    % Costate derivatives
    oc.costate.rate.unconstrained = formDerivative(-oc.hamiltonian.unconstrained.expression,oc.state.var).';
    
    % Write unconstrained control function
    writeControl(in,oc,oc.control.unconstrained.expression,unconstrainedCoefficients,'computeControlUnconstrained', ...
        oc.hamiltonian.unconstrained.expression,oc.control.unconstrained.writeOrder);
    if in.rootSolving == 1
        fprintf(fid_ctrl_c,'#include "computeControlUnconstrained_.h"\n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constrained Arc Calculations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize costate jump partials
oc.num.lagrangeMultipliers.interiorPoint = [];
oc.units.interiorPoint = {};
oc.constraint.interiorPoint = {};

% Loop through each constraint and perform appropriate calculations upfront
for ctrPathConstraint = 1 : 1 : oc.num.constraints.path
    
    % Initialize constraint values
    index = 1;
    controlFound = false;
    constraintExpression = oc.constraint.path{ctrPathConstraint,1};
    
    % Take derivatives until control found for selected constraint
    while true
        % Check if no controls appear in constraint
        for ctrControl = 1 : 1 : oc.num.controls
            
            if ~isempty(strfind(char(constraintExpression),char(oc.control.var(ctrControl))))
                controlFound = true;
                controlInEquation = oc.control.var(ctrControl);
                controlIndex = ctrControl;
            end
            
        end
      
        % If control found, break while loop. Else take time derivative of constraint to find control.
        if controlFound    
            break;
            
        else
            % Save interior point constraint
            oc.constraint.interiorPoint.expression{ctrPathConstraint,1}(index,1) = constraintExpression;
            
            % Compute partial of interior point constraint with respect to the state
            oc.constraint.interiorPoint.statePartial{ctrPathConstraint}(1:length(oc.state.var),index) = ...
                formDerivative(oc.constraint.interiorPoint.expression{ctrPathConstraint,1}(index,1),oc.state.var).';
            
            % Compute partial of interior point constraint with respect to the independent variable
            oc.constraint.interiorPoint.independentVariablePartial{ctrPathConstraint}(index,1) = ...
                formDerivative(oc.constraint.interiorPoint.expression{ctrPathConstraint,1}(index,1),sym(in.oc.independentVariable{1}));
            
            % Save units
            oc.units.interiorPoint{ctrPathConstraint,1}(index,1) = ...
                sym(in.oc.constraint.path{ctrPathConstraint,4})/sym(in.oc.independentVariable{2})^(index-1);
            
            % Control not found yet, take another derivative.
            index = index + 1;
            constraintExpression = formDerivative(constraintExpression,oc.state.var)*oc.state.rate;
            
        end
        % constraintExpression
        if in.rootSolving == 1
            fprintf(fid_ctrl_c,['#include "computeControlConstraint',num2str(ctrPathConstraint),'_.h"\n']);
        end
    end
    %% Change made in March 2015: If no interior point constraint then initialize its value and take paritals wrt time and state
    if isempty(oc.constraint.interiorPoint)
        oc.constraint.interiorPoint.expression{ctrPathConstraint,1}(index,1) = '0';
        % Compute partial of interior point constraint with respect to the state
        oc.constraint.interiorPoint.statePartial{ctrPathConstraint}(1:length(oc.state.var),index) = ...
            formDerivative(oc.constraint.interiorPoint.expression{ctrPathConstraint,1}(index,1),oc.state.var).';
        
        % Compute partial of interior point constraint with respect to the independent variable
        oc.constraint.interiorPoint.independentVariablePartial{ctrPathConstraint}(index,1) = ...
            formDerivative(oc.constraint.interiorPoint.expression{ctrPathConstraint,1}(index,1),sym(in.oc.independentVariable{1}));
        
        % Save units
        oc.units.interiorPoint{ctrPathConstraint,1}(index,1) = ...
            sym(in.oc.constraint.path{ctrPathConstraint,4})/sym(in.oc.independentVariable{2})^(index-1);
        
        
    end
    
    % Determine number of Lagrange multipliers
    oc.num.lagrangeMultipliers.interiorPoint(ctrPathConstraint) = length(oc.constraint.interiorPoint.expression{ctrPathConstraint,1});
    
    if in.oc.writeEquations
        % Write external function for interior point state partial
        funcName = ['interiorPoint',int2str(ctrPathConstraint),'statePartial'];
        fid = fopen([in.autocodeDirTraj,'/',funcName,'.m'],'w');
        fprintf(fid,['function dNdx = ',funcName,'(']);
        for ctrControl = 1 : 1 : oc.num.controls
            fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
        end
        for ctrState = 1 : 1 : oc.num.states
            fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
            if ctrState ~= oc.num.states
                fprintf(fid,',');
            end
        end
        fprintf(fid,',const,constraint)\n');
        fprintf(fid,'\n');
        
        writeConstants(fid,in.const,false);
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        if isfield(in,'constraintVal')
            writeConstraints(fid,in.constraintVal,false);
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'\n');
        [numRow,numCol] = size(oc.constraint.interiorPoint.statePartial{ctrPathConstraint});
        fprintf(fid,'dNdx = [');
        for row = 1 : 1 : numRow
            for col = 1 : 1 : numCol
                fprintf(fid,'%s',char(oc.constraint.interiorPoint.statePartial{ctrPathConstraint}(row,col)));
                if col ~= numCol
                    fprintf(fid,', ');
                elseif row ~= numRow
                    fprintf(fid,';\n  ');
                end
            end
        end
        fprintf(fid,'];\n\n');
        fprintf(fid,'return\n\n');
        fclose(fid);
        
        % Write external function for interior point independent variable partial
        funcName = ['interiorPoint',int2str(ctrPathConstraint),'independentVariablePartial'];
        fid = fopen([in.autocodeDirTraj,'/',funcName,'.m'],'w');
        fprintf(fid,['function dNdt = ',funcName,'(']);
        for ctrControl = 1 : 1 : oc.num.controls
            fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
        end
        for ctrState = 1 : 1 : oc.num.states
            fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
            if ctrState ~= oc.num.states
                fprintf(fid,',');
            end
        end
        fprintf(fid,',const,constraint)\n');
        fprintf(fid,'\n');
        
        writeConstants(fid,in.const,false);
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        if isfield(in,'constraintVal')
            writeConstraints(fid,in.constraintVal,false);
        end
        fprintf(fid,'\n');
        [numRow,numCol] = size(oc.constraint.interiorPoint.independentVariablePartial{ctrPathConstraint});
        fprintf(fid,'dNdt = [');
        for row = 1 : 1 : numRow
            for col = 1 : 1 : numCol
                fprintf(fid,'%s',char(oc.constraint.interiorPoint.independentVariablePartial{ctrPathConstraint}(row,col)));
                if col ~= numCol
                    fprintf(fid,', ');
                elseif row ~= numRow
                    fprintf(fid,';\n  ');
                end
            end
        end
        fprintf(fid,'];\n\n');
        fprintf(fid,'return\n\n');
        fclose(fid);
        
        % Determine which control appears in constraint equation
        constraintWithControl = constraintExpression;
        % keyboard
        % Collect coefficients of control for polynomial expressions (solution to higher order powers often result in large expressions)
        polyExpression = collect(constraintWithControl,controlInEquation);
        [controlCoefficients,controlExpressions] = coeffs(polyExpression,controlInEquation);
        %% Changes made in March 2015: constants and trigonometric function check made
        numSinTerms = length(coeffs(polyExpression,sin(controlInEquation)))-1; % Take constant term out, there will always be one constant term
        numCosTerms = length(coeffs(polyExpression,cos(controlInEquation)))-1; % Take constant term out, there will always be one constant term
        numTrigTerms = numSinTerms + numCosTerms;
        % keyboard
        % Rewrite control polynomial using coefficients
        %% Made a change here! (if length(controlExpressions)>1 originally).
        %% Changed it to >= instead.
        if (length(controlExpressions) >= 1 && numTrigTerms == 0)  % polynomial most likely found

            constraintCoefficientsCollect{controlIndex} = controlCoefficients(:);
            
            compactControlPoly = 0;
            for ctrCoeff = 1 : 1 : length(controlExpressions)
                compactControlPoly = compactControlPoly + sym(['coeff',int2str(ctrCoeff)])*controlExpressions(ctrCoeff);
            end
            
        else % result is not a polynomial equation. assume trig (cos and sin) used.
            
            % Determine coefficients of terms with cos(control)
            polyExpressionCos = collect(constraintWithControl,cos(controlInEquation));
            [controlCoefficientsCos,controlExpressionsCos] = coeffs(polyExpressionCos,cos(controlInEquation));
            
            % Write out resulting expression. Skip cos expression with 1 to prevent double-counting the constant term that has neither sin nor cos. Use sin expression with 1.
            compactControlPoly = 0;
            index = 0;
            for ctrCoeff = 1 : 1 : length(controlExpressionsCos)
                if controlExpressionsCos(ctrCoeff) ~= 1
                    index = index + 1;
                    compactControlPoly = compactControlPoly + sym(['coeff',int2str(index)])*controlExpressionsCos(ctrCoeff);
                    constraintCoefficientsCollect{controlIndex}(index) = controlCoefficientsCos(ctrCoeff);
                else
                    % Determine coefficients of terms with sin(control)
                    polyExpressionSin = collect(constraintWithControl,sin(controlInEquation));
                    [controlCoefficientsSin,controlExpressionsSin] = coeffs(controlCoefficientsCos(ctrCoeff),sin(controlInEquation));
                end
            end
            
            for ctrCoeff = 1 : 1 : length(controlExpressionsSin)
                index = index + 1;
                compactControlPoly = compactControlPoly + sym(['coeff',int2str(index)])*controlExpressionsSin(ctrCoeff);
                constraintCoefficientsCollect{controlIndex}(index) = controlCoefficientsSin(ctrCoeff);
            end
            
        end
        
        % Solve for control solution(s)
        [constraintControl{controlIndex}] = formSolverCombined(compactControlPoly,sym('0'),controlInEquation, ...
            in.oc.assumptions.control,in,oc);
        % Determine constraint coefficients
        for ctrCoefficientSet = 1 : 1 : length(constraintControl{controlIndex})
            constraintCoefficients{controlIndex}{ctrCoefficientSet} = constraintCoefficientsCollect{controlIndex};
        end
        
        % Add remaining control solution(s) from unconstrained arcs and save corresponding coefficients
        for ctrControl = 1 : 1 : oc.num.controls
            if ctrControl ~= controlIndex % new control
                constraintControl{ctrControl} = oc.control.unconstrained.expression{ctrControl};
                for ctrCoefficientSet = 1 : 1 : length(constraintControl{ctrControl})
                    constraintCoefficients{ctrControl}{ctrCoefficientSet} = unconstrainedCoefficients{ctrControl}{ctrCoefficientSet};
                end
            end
        end
        
        % Save constrained control set
        oc.control.constraint.path{ctrPathConstraint}.expression = constraintControl;
        
        % Determine the write order for constraint controls
        controlsInExpression = zeros(oc.num.controls,1);
        
        for ctrExpression = 1 : 1 : oc.num.controls
            for ctrControl = 1 : 1 : oc.num.controls
                
                if ~isempty(strfind(char(oc.control.constraint.path{ctrPathConstraint}.expression{ctrExpression}(1)), ...
                        char(oc.control.var(ctrControl))))
                    controlsInExpression(ctrExpression) = controlsInExpression(ctrExpression) + 1;
                end
                
            end
        end
        
        % Sort controls starting with those with fewest appearances in control equations
        [Y,oc.control.constraint.path{ctrPathConstraint}.writeOrder] = sort(controlsInExpression);
        
        % Write control options. Note that Hamiltonian along constraint is the same as that along an unconstrained arc.
        coefficients = constraintCoefficients;
        
        writeControl(in,oc,oc.control.constraint.path{ctrPathConstraint}.expression,coefficients,['computeControlConstraint',int2str(ctrPathConstraint)], ...
            oc.hamiltonian.unconstrained.expression,oc.control.constraint.path{ctrPathConstraint}.writeOrder);
        
        % Compute Hamiltonian along constraint
        muConstraint = sym(['mu',int2str(ctrPathConstraint)]);
        oc.hamiltonian.constraint.path{ctrPathConstraint}.expression = oc.cost.path + oc.costate.var.'*oc.state.rate + ...
            muConstraint*constraintExpression;
        
        % Compute partial of Hamiltonian with respect to control
        for ctrControl = 1 : 1 : oc.num.controls
            oc.hamiltonian.constraint.path{ctrPathConstraint}.controlPartial(ctrControl,1) = ...
                diff(oc.hamiltonian.constraint.path{ctrPathConstraint}.expression,oc.control.var(ctrControl));
        end
        
        % Compute Lagrange multiplier derivatives
        oc.costate.rate.constraint.path{ctrPathConstraint} = ...
            -formDerivative(oc.hamiltonian.constraint.path{ctrPathConstraint}.expression,oc.state.var).';
        
        % Calculate constraint Lagrange multiplier
        muConstraint = sym(['mu',int2str(ctrPathConstraint)]);
        oc.lagrangeMultiplier.constraint.path{ctrPathConstraint} = ...
            simplify(solve([char(oc.hamiltonian.constraint.path{ctrPathConstraint}.controlPartial(controlIndex,1)),'=0'], ...
            char(muConstraint)));
    end
end
fclose(fid_ctrl_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Interior Point Cost    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oc.cost.interiorPoint.constraintIndex = in.oc.cost.interiorPoint{3};
if in.oc.cost.interiorPoint{3} > 0
    oc.cost.interiorPoint.expression(1,1) = sym(in.oc.cost.interiorPoint{1});
    
    % Compute partial of interior point cost with respect to the state
    oc.cost.interiorPoint.statePartial(1:length(oc.state.var),1) = ...
        formDerivative(oc.cost.interiorPoint.expression(1,1),oc.state.var).';
    
    % Compute partial of interior point cost with respect to the independent variable
    oc.cost.interiorPoint.independentVariablePartial(1,1) = ...
        formDerivative(oc.cost.interiorPoint.expression(1,1),sym(in.oc.independentVariable{1}));
    
    if in.oc.writeEquations
        % Write out the partial derivatives
        funcName = ['interiorPointCostStatePartial'];
        fid = fopen([in.autocodeDirTraj,'/',funcName,'.m'],'w');
        fprintf(fid,['function dIdx = ',funcName,'(']);
        for ctrControl = 1 : 1 : oc.num.controls
            fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
        end
        for ctrState = 1 : 1 : oc.num.states
            fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
            if ctrState ~= oc.num.states
                fprintf(fid,',');
            end
        end
        fprintf(fid,',const,constraint)\n');
        fprintf(fid,'\n');
        
        writeConstants(fid,in.const,false);
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        if isfield(in,'constraintVal')
            writeConstraints(fid,in.constraintVal,false);
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'\n');
        [numRow,numCol] = size(oc.cost.interiorPoint.statePartial);
        fprintf(fid,'dIdx = [');
        for row = 1 : 1 : numRow
            for col = 1 : 1 : numCol
                fprintf(fid,'%s',char(oc.cost.interiorPoint.statePartial(row,col)));
                if col ~= numCol
                    fprintf(fid,', ');
                elseif row ~= numRow
                    fprintf(fid,';\n  ');
                end
            end
        end
        fprintf(fid,'];\n\n');
        fprintf(fid,'return\n\n');
        fclose(fid);
        
        % Write external function for interior point independent variable partial
        funcName = ['interiorPointCostIndependentVariablePartial'];
        fid = fopen([in.autocodeDirTraj,'/',funcName,'.m'],'w');
        fprintf(fid,['function dIdt = ',funcName,'(']);
        for ctrControl = 1 : 1 : oc.num.controls
            fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
        end
        for ctrState = 1 : 1 : oc.num.states
            fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
            if ctrState ~= oc.num.states
                fprintf(fid,',');
            end
        end
        fprintf(fid,',const,constraint)\n');
        fprintf(fid,'\n');
        
        writeConstants(fid,in.const,false);
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        if isfield(in,'constraintVal')
            writeConstraints(fid,in.constraintVal,false);
        end
        fprintf(fid,'\n');
        [numRow,numCol] = size(oc.cost.interiorPoint.independentVariablePartial);
        fprintf(fid,'dIdt = [');
        for row = 1 : 1 : numRow
            for col = 1 : 1 : numCol
                fprintf(fid,'%s',char(oc.cost.interiorPoint.independentVariablePartial(row,col)));
                if col ~= numCol
                    fprintf(fid,', ');
                elseif row ~= numRow
                    fprintf(fid,';\n  ');
                end
            end
        end
        fprintf(fid,'];\n\n');
        fprintf(fid,'return\n\n');
        fclose(fid);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interior Point Constraint %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize costate jump partials
indexInteriorPoint = oc.num.constraints.path; % Continue counting interior point constraints after last path constraint

% Loop through each constraint and perform appropriate calculations upfront
for ctrInteriorPoint = 1 : 1 : oc.num.constraints.interiorPoint
    
    % Initialize
    indexInteriorPoint = indexInteriorPoint + 1; % index to count from path constraint interior points
    oc.num.lagrangeMultipliers.interiorPoint(indexInteriorPoint) = 1; % only one Lagrange multiplier
    % oc.constraint.interiorPoint.expression{indexInteriorPoint,1}(1,1) = sym(in.oc.constraint.interiorPoint{ctrInteriorPoint,3})/sym(in.oc.constraint.interiorPoint{ctrInteriorPoint,1}) - sym(1);
    oc.constraint.interiorPoint.expression{indexInteriorPoint,1}(1,1) = sym(in.oc.constraint.interiorPoint{ctrInteriorPoint,3}) - sym(in.oc.constraint.interiorPoint{ctrInteriorPoint,1});
    
    % Compute partial of interior point constraint with respect to the state
    oc.constraint.interiorPoint.statePartial{indexInteriorPoint}(1:length(oc.state.var),1) = ...
        formDerivative(oc.constraint.interiorPoint.expression{indexInteriorPoint,1}(1,1),oc.state.var).';
    
    % Compute partial of interior point constraint with respect to the independent variable
    oc.constraint.interiorPoint.independentVariablePartial{indexInteriorPoint}(1,1) = ...
        formDerivative(oc.constraint.interiorPoint.expression{indexInteriorPoint,1}(1,1),sym(in.oc.independentVariable{1}));
    
    if in.oc.writeEquations
        
        % Write external function for interior point state partial
        funcName = ['interiorPoint',int2str(indexInteriorPoint),'statePartial'];
        fid = fopen([in.autocodeDirTraj,'/',funcName,'.m'],'w');
        fprintf(fid,['function dNdx = ',funcName,'(']);
        for ctrControl = 1 : 1 : oc.num.controls
            fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
        end
        for ctrState = 1 : 1 : oc.num.states
            fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
            if ctrState ~= oc.num.states
                fprintf(fid,',');
            end
        end
        fprintf(fid,',const,constraint)\n');
        fprintf(fid,'\n');
        
        writeConstants(fid,in.const,false);
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        if isfield(in,'constraintVal')
            writeConstraints(fid,in.constraintVal,false);
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'\n');
        [numRow,numCol] = size(oc.constraint.interiorPoint.statePartial{indexInteriorPoint});
        fprintf(fid,'dNdx = [');
        for row = 1 : 1 : numRow
            for col = 1 : 1 : numCol
                fprintf(fid,'%s',char(oc.constraint.interiorPoint.statePartial{indexInteriorPoint}(row,col)));
                if col ~= numCol
                    fprintf(fid,', ');
                elseif row ~= numRow
                    fprintf(fid,';\n  ');
                end
            end
        end
        fprintf(fid,'];\n\n');
        fprintf(fid,'return\n\n');
        fclose(fid);
        
        % Write external function for interior point independent variable partial
        funcName = ['interiorPoint',int2str(indexInteriorPoint),'independentVariablePartial'];
        fid = fopen([in.autocodeDirTraj,'/',funcName,'.m'],'w');
        fprintf(fid,['function dNdt = ',funcName,'(']);
        for ctrControl = 1 : 1 : oc.num.controls
            fprintf(fid,'%s,',char(oc.control.var(ctrControl)));
        end
        for ctrState = 1 : 1 : oc.num.states
            fprintf(fid,'%s',char(oc.state.var(ctrState,1)));
            if ctrState ~= oc.num.states
                fprintf(fid,',');
            end
        end
        fprintf(fid,',const,constraint)\n');
        fprintf(fid,'\n');
        
        writeConstants(fid,in.const,false);
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        if isfield(in,'constraintVal')
            writeConstraints(fid,in.constraintVal,false);
        end
        fprintf(fid,'\n');
        [numRow,numCol] = size(oc.constraint.interiorPoint.independentVariablePartial{indexInteriorPoint});
        fprintf(fid,'dNdt = [');
        for row = 1 : 1 : numRow
            for col = 1 : 1 : numCol
                fprintf(fid,'%s',char(oc.constraint.interiorPoint.independentVariablePartial{indexInteriorPoint}(row,col)));
                if col ~= numCol
                    fprintf(fid,', ');
                elseif row ~= numRow
                    fprintf(fid,';\n  ');
                end
            end
        end
        fprintf(fid,'];\n\n');
        fprintf(fid,'return\n\n');
        fclose(fid);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write Files for Indirect Solver %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary Conditions Function (bc.m) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if in.oc.writeEquations
    fid = fopen([in.autocodeDirTraj,'/bc.m'],'w');
    writeBC(fid,in,oc);
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Derivative Function (derivFunc.m) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Write function for single region case (to support BVP4C).
    multipointBVP = false;
    fid = fopen([in.autocodeDirTraj,'/derivFunc.m'],'w');
    writeDerivFunc(fid,in,oc,multipointBVP);
    fclose(fid);
    
    % Write function for multi-region case (to support BVP4C).
    multipointBVP = true;
    fid = fopen([in.autocodeDirTraj,'/derivFuncRegion.m'],'w');
    writeDerivFunc(fid,in,oc,multipointBVP);
    fclose(fid);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jacobian of Derivative Function (derivFuncJac.m) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if in.writeJac
    
    % Used by BVP4C. Takes a while to generate this function, but execution speed
    % is dramatically increased.
    jacY = [oc.x oc.lambda];
    jacP = [];
    fid = fopen([in.autocodeDirTraj,'/derivFuncJac.m'],'w');
    writeDerivFuncJac(fid,in,oc);
    fclose(fid);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jacobian of Boundary Conditions Function (bcJac.m) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Future work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% State Transition Tensor Function (derivFuncSTT.m) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This currently only works for a single unconstrained arc.

% Create symbolic variables
orderSTT = in.cont.orderSTT;

if orderSTT > 0
    
    % Create augmented vector of state and costate equations
    oc.fAug = oc.xDot;
    for ctr = 1 : 1 : oc.numStates
        oc.fAug = [oc.fAug; sym(oc.lambdaDot.unconstrained{ctr})];
    end
    
    % Augmented state vector
    oc.xAug = [oc.x; oc.lambda];
    oc.numAugStates = length(oc.xAug);
    
    % Construct initial values for STT regardless if function is made.
    for ctr = 1 : 1 : orderSTT
        % Determine size of Jacobian values are being assigned to
        sizeJac = oc.numAugStates*ones(1,ctr+1);
        
        % Construct empty cell array to store subscript indices
        sub = cell(1,length(sizeJac));
        
        % Compute Jacobian
        for ind = 1 : 1 : oc.numAugStates^(ctr+1)
            % Construct cell array of subscript elements in Jacobian
            [sub{:}] = ind2sub(sizeJac,ind);
            % Construct initial values for STT. Identity for first order STT and zero
            % for rest.
            if (ctr == 1) && (sum(diff([sub{:}])) == 0)
                val = 1;
            else
                val = 0;
            end
            X0{ctr}(ind,1) = val;
        end
        
    end
    
    % Initial values for state transition tensors
    oc.X0 = X0;
    
    % Write state transition tensors to function
    if in.oc.writeSTTfunc
        
        % Open function that will be used to compute STTs symbolically
        fidSymSTT = fopen([in.autocodeDirTraj,'/',in.robust.funcSymSTT,'.m'],'w');
        fid = fopen([in.autocodeDirTraj,'/derivFuncSTT.m'],'w');
        writeSTT(fid,fidSymSTT,in,oc);
        
        fclose(fid);
        fclose(fidSymSTT);
        
        % Create codegen file to speed up calculations
        curDir = pwd;
        cd([curDir,'/autocode/optimization']);
        codegen derivFuncSTT
        cd(curDir);
        
    end
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Write Custom Fixed-Step Solver %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This is useful in the presence of state rate discontinuities such as those
% % arising from switches in control.
%
% maxPts = 10000;
% fid = fopen([in.autocodeDirTraj,'/ode4.m'],'w');
% fprintf(fid,'function [xOut] = ode4(dt,X0,p,const,constraint,scale,VEH,x0,xf) %%#ok<INUSD,INUSL>\n');
% writeHeader(fid,in,out.VEH);
% fprintf(fid,'assert(isa(dt, ''double''));\n');
% fprintf(fid,'assert(all(size(dt)== [1 1]));\n');
% fprintf(fid,'assert(isa(X0, ''double''));\n');
% fprintf(fid,['assert(all(size(X0)== [',int2str(2*oc.numStates),' 1]));\n']);
% fprintf(fid,'assert(isa(p, ''double''));\n');
% fprintf(fid,['assert(all(size(p)== [',int2str(oc.numP + oc.numNu0 + oc.numNuF),' 1]));\n']);
% fprintf(fid,'assert(isa(x0, ''double''));\n');
% fprintf(fid,['assert(all(size(x0)== [',int2str(oc.numStates),' 1]));\n']);
% fprintf(fid,'assert(isa(xf, ''double''));\n');
% fprintf(fid,['assert(all(size(xf)== [',int2str(oc.numStates),' 1]));\n']);
% fprintf(fid,['xOut = nan(',int2str(1/in.dt+1),',',num2str(2*oc.numStates),');\n']);
% fprintf(fid,'\n');
% fprintf(fid,'xCur = X0;\n');
% fprintf(fid,'xOut(1,:) = xCur'';\n');
% fprintf(fid,'t = 0;\n');
% fprintf(fid,['for ctr = 2 : 1 : ',int2str(1/in.dt+1),'\n']);
% fprintf(fid,'  F1 = derivFunc1(t,xCur,p,const,constraint,scale,VEH,x0,xf);\n');
% fprintf(fid,'  F2 = derivFunc1(t+0.5*dt,xCur+0.5*dt*F1,p,const,constraint,scale,VEH,x0,xf);\n');
% fprintf(fid,'  F3 = derivFunc1(t+0.5*dt,xCur+0.5*dt*F2,p,const,constraint,scale,VEH,x0,xf);\n');
% fprintf(fid,'  F4 = derivFunc1(t+dt,xCur+dt*F3,p,const,constraint,scale,VEH,x0,xf);\n');
% fprintf(fid,'  xCur = xCur + (dt/6)*(F1 + 2*F2 + 2*F3 + F4);\n');
% fprintf(fid,'  xOut(ctr,:) = xCur'';\n');
% fprintf(fid,'  t = t + dt;\n');
% fprintf(fid,'end\n');
% fprintf(fid,'\n');
% fprintf(fid,'return\n');
% fprintf(fid,'\n');
%
% % Create codegen file to speed up calculations
% curDir = pwd;
% cd([curDir,'/autocode/optimization']);
% codegen ode4;
% [SUCCESS,MESSAGE,MESSAGEID] = copyfile('./codegen/mex/ode4/ode4_mex.mexmaci64','.');
% cd(curDir);

return

