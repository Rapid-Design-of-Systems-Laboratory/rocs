function [inMod, setCONT] = runContinuation(in,out)
% This function builds upon the solution obtained by runIndirect.m using
% continuation method.
% input : in and out [structures]
%         in is obtained from inputsRunTrajectoryProcess.m
%         out is obtained from initialCalculations.m,
%         optimalCalcs.m and runIndirect.m
% output : out.setCONT [structure]
% Developed by : Dr. M.J. Grant
% Modified by : Kshitij Mall

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

oc = out.oc;

% Determine number of states and other parameters
numStates = out.oc.num.states; % total number of states
numControl = out.oc.num.controls;
% numNuTotal = out.oc.num.constraints.initial + out.oc.num.constraints.terminal;

if ~isfield(in,'skipUnconvergedSolutions')
    in.skipUnconvergedSolutions = false;
end

in.vars.numContinuations = length(in.CONT);

if isfield(in.oc,'initialGuessFunc') || strcmp(in.oc.guess.mode,'file')==0
    % else
    % Initialize continuation index
    
    % Initialize arc type sequence. 0 implies unconstrained, nonzero values determine index from constraint cell.
    in.vars.arcTypeSequence = 0; % initially unconstrained
    
    if oc.num.constraints.path+oc.num.constraints.interiorPoint > 0
        in.vars.interiorPointConstraintSequence = NaN(oc.num.constraints.path+oc.num.constraints.interiorPoint,1); % initially no interior point constraints
        in.vars.interiorPointConstraintSequence(1,1) = 0;
        in.vars.interiorPointNumLagrangeMultipliers = 0; % initially no interior point Lagrange multipliers
        in.vars.interiorPointNumLagrangeMultipliersCumSum = cumsum(in.vars.interiorPointNumLagrangeMultipliers);
    else
        in.vars.interiorPointConstraintSequence = [];
        in.vars.interiorPointNumLagrangeMultipliers = [];
        in.vars.interiorPointNumLagrangeMultipliersCumSum = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% Preproccess Scaling %%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % [out] = initializeScaling(in,out,numStates);
    [in] = initializeScaling(in,out);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run Continuation Set %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Determine number of continuation sets
    % Continuation cycle will run for both unconstrained and constrained arcs
    % For further information please refer to inputsRunTrajectoryProcess.m bottom-most part
    out.setCONT = struct([]); % initialize out.setCONT structure
    
    if length(in.cont.method) == 1
        in.cont.method = ones(in.vars.numContinuations,1);
    end
end

in.vars.numContinuations = length(in.CONT);


% Initialize temporary variables used to increment continuation parameters
% inTemp = out.IG.in;
inTemp = in;
inMod = in;
% outTemp = out;

sol = out.IG.sol;
% piiIndex = [];

% Optimal control file index
% fileIndex = 1;
contIndex = 1;

% Create handles to mex files
if in.useMex
    derivFunc = str2func('derivFunc_mex');
    bcFunc = str2func('bc_mex');
    derivFuncRegion = str2func('derivFuncRegion_mex');
    derivFuncJac_mex = str2func('derivFuncJac_mex');
else
    derivFunc = str2func('derivFunc');
    bcFunc = str2func('bc');
    derivFuncRegion = str2func('derivFuncRegion');
    derivFuncJac = str2func('derivFuncJac');
end

derivFuncRegionGPU = str2func('derivFuncRegion');

% Initialize new active constraint flag
newActiveConstraint = false;


for setCtr = in.cont.startIndex : 1 : in.vars.numContinuations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize Continuation Run %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Output to screen
    fprintf('\nExecuting Continuation Set #%i:\n',setCtr);
    
    % Determine number of continuation runs
    if in.cont.method(setCtr) == 0
        in.cont.method(setCtr) = 1; % makes fixed step default method
    end
    
    if in.cont.method(setCtr) == 1
        numRuns = in.CONT{setCtr}.numCases; % user defined, doesn't change
    elseif in.cont.method(setCtr) == 2
        numRuns = 1; % modified during continuation until new BC reached
    end
    
    % Initialize continuation variables
    CONT = struct([]);
    
    % Determine how to split trajectory with new constraint
    if newActiveConstraint
        
        % Find indices for arc insertion
        %%%%% MADE CHANGE HERE FOR SURFACE CONSTRAINT
        Iconstraint1 = find(sol.x <= timeNormalizedConstraintVal,1,'last')-1;
        if timeNormalizedConstraintVal == sol.x(end) % protect for constraint at end of unconstrained arc
            Iconstraint1 = Iconstraint1 - 2; % results in constraint and following unconstrained arc to only have two points
        end
        
        % Determine number of original arcs
        numOrigArcs = length(in.vars.arcTypeSequence);
        timeIndex = 2*in.oc.num.origStates;
        
        % Determine arc that new constraint will be inserted into. First arc is considered the zeroth arc.
        arcNum = floor(sol.x(Iconstraint1));
        
        % Find additional indices for arc insertion
        if strcmp(constraintType,'path')
            Iconstraint2 = Iconstraint1 + 1;
        elseif strcmp(constraintType,'interior')
            Iconstraint2 = Iconstraint1; % only one point for interior point constraint
        end
        IoldArcBeforeConstraint = find(sol.x == arcNum,1,'last'); % beginning of segment before constraint
        IoldArcAfterConstraint = find(sol.x == arcNum+1,1,'first'); % end of segment after constraint
        
        % Save times across each arc
        
        if ~in.convertParametersToStates
            independentVariableArcBeforeConstraint = (sol.x(Iconstraint1)-arcNum)*sol.parameters(arcNum+1);
            if strcmp(constraintType,'path')
                independentVariableArcDuringConstraint = (sol.x(Iconstraint2)-sol.x(Iconstraint1))*sol.parameters(arcNum+1);
            end
            independentVariableArcAfterConstraint = (sol.x(IoldArcAfterConstraint)-sol.x(Iconstraint2))*sol.parameters(arcNum+1);
        else
            independentVariableArcBeforeConstraint = (sol.x(Iconstraint1)-arcNum)*sol.y(timeIndex+arcNum+1,1);
            if strcmp(constraintType,'path')
                independentVariableArcDuringConstraint = (sol.x(Iconstraint2)-sol.x(Iconstraint1))*sol.y(timeIndex+arcNum+1,1);
            end
            independentVariableArcAfterConstraint = (sol.x(IoldArcAfterConstraint)-sol.x(Iconstraint2))*sol.y(timeIndex+arcNum+1,1);
        end
        
        % Break trajectory arc into an appropriate number of segments
        
        if strcmp(constraintType,'interior')
            % Save delta time between values at location of interior point constraint
            dtSave = sol.x(Iconstraint2+1) - sol.x(Iconstraint1);
        end
        
        % First segment (old arc)
        sol.x(IoldArcBeforeConstraint:Iconstraint1) = (sol.x(IoldArcBeforeConstraint:Iconstraint1)-sol.x(IoldArcBeforeConstraint)) / ...
            (sol.x(Iconstraint1)-sol.x(IoldArcBeforeConstraint)) + arcNum;
        
        
        if strcmp(constraintType,'path')
            
            % Third segment (old arc)
            sol.x(Iconstraint2:IoldArcAfterConstraint) = (sol.x(Iconstraint2:IoldArcAfterConstraint)-sol.x(Iconstraint2)) / ...
                (sol.x(IoldArcAfterConstraint)-sol.x(Iconstraint2)) + arcNum+2;
            
            % Rest of arcs after third segment (if exist)
            if IoldArcAfterConstraint ~= length(sol.x)
                sol.x(IoldArcAfterConstraint+1:end) = sol.x(IoldArcAfterConstraint+1:end) + 2;
            end
            
            % Second segment (new constrained arc). Consists of only two points that must be repeated. Insert arc into time segment and solution vector.
            sol.x = [sol.x(1:Iconstraint1) sol.x(Iconstraint1:Iconstraint2) sol.x(Iconstraint2:end)];
            sol.y = [sol.y(:,1:Iconstraint1) sol.y(:,Iconstraint1:Iconstraint2) sol.y(:,Iconstraint2:end)];
            
        elseif strcmp(constraintType,'interior')
            
            % Second segment (old arc)
            sol.x(Iconstraint2+1:IoldArcAfterConstraint) = (sol.x(Iconstraint2+1:IoldArcAfterConstraint) - sol.x(Iconstraint2+1) + dtSave) / ...
                (sol.x(IoldArcAfterConstraint) - sol.x(Iconstraint2+1) + dtSave) + arcNum+1;
            
            % Rest of arcs after second segment (if exist)
            if IoldArcAfterConstraint ~= length(sol.x)
                sol.x(IoldArcAfterConstraint+1:end) = sol.x(IoldArcAfterConstraint+1:end) + 1;
            end
            
            % Split original arc into two for both the time segment and solution vector.
            sol.x = [sol.x(1:Iconstraint1) sol.x(Iconstraint2:end)];
            sol.y = [sol.y(:,1:Iconstraint1) sol.y(:,Iconstraint2:end)];
            
        end
        % keyboard
        % Alter parameters vector to include extra arc times and costate jump conditions. Uses original time values.
        parametersOrig = sol.parameters;
        if ~in.convertParametersToStates
            independentVariableSetOrig = parametersOrig(1:numOrigArcs);
            lagrangeMultiplierInitialSetOrig = parametersOrig(numOrigArcs+1:numOrigArcs+oc.num.lagrangeMultiplier.initial);
            lagrangeMultiplierFinalSetOrig = parametersOrig(numOrigArcs+oc.num.lagrangeMultiplier.initial+1: ...
                numOrigArcs+oc.num.lagrangeMultiplier.initial+oc.num.lagrangeMultiplier.terminal);
            lagrangeMultiplierInteriorPointSetOrig = parametersOrig(numOrigArcs+oc.num.lagrangeMultiplier.initial+ ...
                oc.num.lagrangeMultiplier.terminal+1:end);
        else
            independentVariableSetOrig = sol.y(timeIndex+1:timeIndex+numOrigArcs,1);
            lagrangeMultiplierInitialSetOrig = parametersOrig(1:oc.num.lagrangeMultiplier.initial);
            lagrangeMultiplierFinalSetOrig = parametersOrig(oc.num.lagrangeMultiplier.initial+1: ...
                oc.num.lagrangeMultiplier.initial+oc.num.lagrangeMultiplier.terminal);
            lagrangeMultiplierInteriorPointSetOrig = parametersOrig(oc.num.lagrangeMultiplier.initial+ ...
                oc.num.lagrangeMultiplier.terminal+1:end);
        end
        
        % Construct time parameter sequence
        if ~in.convertParametersToStates
            if strcmp(constraintType,'path')
                sol.parameters = [independentVariableSetOrig(1:arcNum);independentVariableArcBeforeConstraint; ...
                    independentVariableArcDuringConstraint;independentVariableArcAfterConstraint;independentVariableSetOrig(arcNum+2:end)];
            elseif strcmp(constraintType,'interior')
                sol.parameters = [independentVariableSetOrig(1:arcNum);independentVariableArcBeforeConstraint; ...
                    independentVariableArcAfterConstraint;independentVariableSetOrig(arcNum+2:end)];
            end
        else
            sol.parameters = [];
            if strcmp(constraintType,'path')
                indepVarValues = [independentVariableSetOrig(1:arcNum);independentVariableArcBeforeConstraint; ...
                    independentVariableArcDuringConstraint;independentVariableArcAfterConstraint;independentVariableSetOrig(arcNum+2:end)];
            else
                indepVarValues = [independentVariableSetOrig(1:arcNum);independentVariableArcBeforeConstraint; ...
                    independentVariableArcAfterConstraint;independentVariableSetOrig(arcNum+2:end)];
            end
            % independentVariableSetOrig
            % indepVarValues
            timeSteps = size(sol.y,2);
            sol.y(timeIndex+1:timeIndex+length(indepVarValues),:) = repmat(indepVarValues,1,timeSteps);
        end
        
        % Construct Lagrange multiplier parameter sequence at initial and terminal points
        sol.parameters = [sol.parameters; lagrangeMultiplierInitialSetOrig; lagrangeMultiplierFinalSetOrig];
        
        % Insert interior point Lagrange multipliers that occur before inserted arc
        if arcNum == 0
            lagrangeMultipliersBeforeArc = [];
        else
            lagrangeMultipliersBeforeArc = lagrangeMultiplierInteriorPointSetOrig(1:in.vars.interiorPointNumLagrangeMultipliersCumSum(arcNum));
        end
        sol.parameters = [sol.parameters; lagrangeMultipliersBeforeArc];
        
        % Add interior point Lagrange multipliers guessed for inserted arc
        if strcmp(constraintType,'path')
            sol.parameters = [sol.parameters; zeros(sum(oc.num.lagrangeMultipliers.interiorPoint(in.CONT{setCtr}.constraint.path)),1)];
            % sol.parameters = [sol.parameters; zeros(sum(oc.num.lagrangeMultipliers.interiorPoint(in.CONT{setCtr}.constraint.path))+2,1)];
        elseif strcmp(constraintType,'interior')
            sol.parameters = [sol.parameters; zeros(sum(oc.num.lagrangeMultipliers.interiorPoint(oc.num.constraints.path+in.CONT{setCtr}.constraint.interiorPoint)),1)];
        end
        
        % Add remaining interior point Lagrange multipliers that occur after inserted arc
        if arcNum == 0
            lagrangeMultipliersAfterArc = lagrangeMultiplierInteriorPointSetOrig(1:in.vars.interiorPointNumLagrangeMultipliersCumSum(end));
        else
            lagrangeMultipliersAfterArc = lagrangeMultiplierInteriorPointSetOrig(in.vars.interiorPointNumLagrangeMultipliersCumSum(arcNum)+1: ...
                in.vars.interiorPointNumLagrangeMultipliersCumSum(end));
        end
        sol.parameters = [sol.parameters; lagrangeMultipliersAfterArc];
        
        % Determine index of active constraint and track order of constrained and unconstrained arcs. 0 implies unconstrained, which these interior point constraints are inserted into.
        interiorPointConstraintInsert = zeros(oc.num.constraints.path+oc.num.constraints.interiorPoint,1);
        interiorPointConstraintInsertZeros = zeros(oc.num.constraints.path+oc.num.constraints.interiorPoint,1);
        
        if strcmp(constraintType,'path')
            
            in.vars.arcTypeSequence = [in.vars.arcTypeSequence(1:arcNum), 0, in.CONT{setCtr}.constraint.path, 0, in.vars.arcTypeSequence(arcNum+2:end)];
            
            interiorPointConstraintInsert(1) = in.CONT{setCtr}.constraint.path;
            in.vars.interiorPointConstraintSequence = [in.vars.interiorPointConstraintSequence(:,1:arcNum), interiorPointConstraintInsertZeros, ...
                interiorPointConstraintInsert, interiorPointConstraintInsertZeros, in.vars.interiorPointConstraintSequence(:,arcNum+2:end)];
            in.vars.interiorPointNumLagrangeMultipliers = [in.vars.interiorPointNumLagrangeMultipliers(1:arcNum), 0, ...
                oc.num.lagrangeMultipliers.interiorPoint(in.CONT{setCtr}.constraint.path), 0, ...
                in.vars.interiorPointNumLagrangeMultipliers(arcNum+2:end)];
            
        elseif strcmp(constraintType,'interior')
            
            in.vars.arcTypeSequence = [in.vars.arcTypeSequence(1:arcNum), 0, 0, in.vars.arcTypeSequence(arcNum+2:end)];
            
            interiorPointConstraintInsert(1:length(in.CONT{setCtr}.constraint.interiorPoint)) = in.CONT{setCtr}.constraint.interiorPoint'+oc.num.constraints.path;
            in.vars.interiorPointConstraintSequence = [in.vars.interiorPointConstraintSequence(:,1:arcNum), interiorPointConstraintInsertZeros, ...
                interiorPointConstraintInsert, in.vars.interiorPointConstraintSequence(:,arcNum+2:end)];
            in.vars.interiorPointNumLagrangeMultipliers = [in.vars.interiorPointNumLagrangeMultipliers(1:arcNum), 0, ...
                sum(oc.num.lagrangeMultipliers.interiorPoint(in.CONT{setCtr}.constraint.interiorPoint+oc.num.constraints.path)), ...
                in.vars.interiorPointNumLagrangeMultipliers(arcNum+2:end)];
            
        end
        
        in.vars.interiorPointNumLagrangeMultipliersCumSum = cumsum(in.vars.interiorPointNumLagrangeMultipliers);
        
        % Update constraint value in input
        for ctrConstraint = 1 : 1 : numConstraints
            
            if strcmp(constraintType,'path')
                inTemp.constraintVal.(in.oc.constraint.path{in.CONT{setCtr}.constraint.path(ctrConstraint)}){1} = ...
                    constraintVal(ctrConstraint);
            elseif strcmp(constraintType,'interior')
                inTemp.constraintVal.(in.oc.constraint.interiorPoint{in.CONT{setCtr}.constraint.interiorPoint(ctrConstraint)}){1} = ...
                    constraintVal(ctrConstraint);
            end
            
        end
        
        
        %     in.vars.interiorPointConstraintSequence
        %
        %     lagrangeMultiplierInteriorPointSetOrig
        
        %     keyboard
        
    end
    
    
    
    
    
    % Chris' mods
    out.numArcs = length(in.vars.arcTypeSequence);
    
    if in.cont.method(setCtr) == 1
        % Initial state
        x0Save = sol.y(1:numStates,1);
        
        % Final state
        xfSave = sol.y(1:numStates,end);
        
        dx = struct;
        
        if isfield(in.CONT{setCtr},'constraint')
            if isfield(in.CONT{setCtr}.constraint,'initial')
                for ctrState = 1 : 1 : numStates
                    if isfield(in.CONT{setCtr}.constraint.initial,inTemp.oc.state{ctrState,1})
                        if numRuns == 1
                            dx.initial.(inTemp.oc.state{ctrState,1}) = 0;
                        else
                            if length(in.CONT{setCtr}.constraint.initial.(inTemp.oc.state{ctrState,1})) == 1
                                dx.initial.(inTemp.oc.state{ctrState,1}) = linspace(0,in.CONT{setCtr}.constraint.initial.(inTemp.oc.state{ctrState,1}) - x0Save(ctrState,1),numRuns);
                            else
                                dx.initial.(inTemp.oc.state{ctrState,1}) = in.CONT{setCtr}.constraint.initial.(inTemp.oc.state{ctrState,1});
                            end
                        end
                    end
                end
            end
        end
        
        if isfield(in.CONT{setCtr},'constraint')
            if isfield(in.CONT{setCtr}.constraint,'terminal')
                for ctrState = 1 : 1 : numStates
                    if isfield(in.CONT{setCtr}.constraint.terminal,inTemp.oc.state{ctrState,1})
                        if numRuns == 1
                            dx.terminal.(inTemp.oc.state{ctrState,1}) = 0;
                        else
                            if length(in.CONT{setCtr}.constraint.terminal.(inTemp.oc.state{ctrState,1})) == 1
                                dx.terminal.(inTemp.oc.state{ctrState,1}) = linspace(0,in.CONT{setCtr}.constraint.terminal.(inTemp.oc.state{ctrState,1}) - xfSave(ctrState,1),numRuns);
                            else
                                dx.terminal.(inTemp.oc.state{ctrState,1}) = in.CONT{setCtr}.constraint.terminal.(inTemp.oc.state{ctrState,1});
                            end
                        end
                    end
                end
            end
        end
        
        
        if isfield(in.CONT{setCtr},'constraint')
            if isfield(in.CONT{setCtr}.constraint,'path')
                for ctrConstraint = 1 : 1 : numConstraints
                    if isfield(in.CONT{setCtr}.constraint.path,inTemp.oc.constraint.path{ctrConstraint,1})
                        if numRuns == 1
                            dx.path.(inTemp.oc.constraint.path{ctrConstraint,1}) = 0;
                        else
                            if length(in.CONT{setCtr}.constraint.path.(inTemp.oc.constraint.path{ctrConstraint,1})) == 1
                                dx.path.(inTemp.oc.constraint.path{ctrConstraint,1}) = linspace(0,in.CONT{setCtr}.constraint.path.(inTemp.oc.constraint.path{ctrConstraint,1}) - constraintVal,numRuns);
                            else
                                dx.path.(inTemp.oc.constraint.path{ctrConstraint,1}) = in.CONT{setCtr}.constraint.path.(inTemp.oc.constraint.path{ctrConstraint,1});
                            end
                        end
                    end
                end
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run Continuation Set %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    continuationParamSave = struct();
    
    ctrContinuationChange = 1; % continuation counter
    
    %for ctrContinuationChange = 1 : 1 : numRuns
    while ctrContinuationChange <= numRuns
        %         if in.cont.method(setCtr) == 1
        %             if ctrContinuationChange == numRuns
        %                 in.gpuSolve.TimeSteps = 256*2;
        %                 in.gpuSolve.ArcSegments = 512;
        %                 in.gpuSolve.Tol = 1e-6;
        %             	  in.bvpOptions.AbsTol = 1e-6;
        %                 in.bvpOptions.RelTol = 1e-6;
        %             else
        %                 in.gpuSolve.TimeSteps = 256;
        %                 in.gpuSolve.Tol = 1e-4;
        %                 in.gpuSolve.TimeSteps = 256;
        %                   in.bvpOptions.AbsTol = 1e-4;
        %                   in.bvpOptions.RelTol = 1e-4;
        %             end
        %         end
        
        % Initialize temporary structures used to modify inputs during continuation
        % 		inTemp = in;
        inTemp.oc.arcTypeSequnce = in.vars.arcTypeSequence;
        outTemp = out;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Save Beginning of Continuation Data %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save continuation parameters from which a delta is computed in current
        % continuation sequence.
        
        % Check if first run of current continuation set
        if ctrContinuationChange == 1
            
            % Initial state
            x0Save = sol.y(1:numStates,1);
            
            % Final state
            xfSave = sol.y(1:numStates,end);
            
            % Continuation variables
            for continuationVar = in.continuationVarSet
                if isfield(in,continuationVar)
                    names = fieldnames(inTemp.(continuationVar{:}));
                else
                    names = {};
                end
                for fieldCtr = 1 : 1 : length(names)
                    continuationParamSave.(names{fieldCtr}) = inTemp.(continuationVar{:}).(names{fieldCtr}){1};
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Change Continuation Parameters %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if (in.cont.method(setCtr) == 1) % Manually changing continuation parameters
            
            x0 = x0Save;
            xf = xfSave;
            
            % Assign new initial and terminal conditions
            if isfield(in.CONT{setCtr},'constraint')
                if isfield(in.CONT{setCtr}.constraint,'initial')
                    for ctrState = 1 : 1 : numStates
                        if isfield(in.CONT{setCtr}.constraint.initial,inTemp.oc.state{ctrState,1})
                            x0(ctrState,1) = x0Save(ctrState,1) + dx.initial.(inTemp.oc.state{ctrState,1})(ctrContinuationChange);
                        end
                    end
                end
            end
            
            if isfield(in.CONT{setCtr},'constraint')
                if isfield(in.CONT{setCtr}.constraint,'terminal')
                    for ctrState = 1 : 1 : numStates
                        if isfield(in.CONT{setCtr}.constraint.terminal,inTemp.oc.state{ctrState,1})
                            xf(ctrState,1) = xfSave(ctrState,1) + dx.terminal.(inTemp.oc.state{ctrState,1})(ctrContinuationChange);
                        end
                    end
                end
            end
            
            if isfield(in.CONT{setCtr},'constraint')
                if isfield(in.CONT{setCtr}.constraint,'path')
                    for ctrConstraint = 1 : 1 : numConstraints
                        if isfield(in.CONT{setCtr}.constraint.path,inTemp.oc.constraint.path{ctrConstraint,1})
                            inTemp.constraintVal.(inTemp.oc.constraint.path{ctrConstraint,1}){1} = constraintVal0(ctrConstraint) + dx.path.(inTemp.oc.constraint.path{ctrConstraint,1})(ctrContinuationChange);
                        end
                    end
                end
            end
            
            %[xf xfSave]
            inMod.CONT{setCtr}.dx = dx;
            
            % Combine general code below from Grant with above code from Dapkus
            
            % Input deck variables
            for continuationVar = in.continuationVarSet
                
                % Determine fieldnames of structure
                if isfield(in,continuationVar)
                    names = fieldnames(inTemp.(continuationVar{:}));
                else
                    names = {};
                end
                
                % Update values
                for fieldCtr = 1 : 1 : length(names)
                    
                    % Only need to increment variable if exists in continuation structure
                    if isfield(in.CONT{setCtr},continuationVar{:}) && ...
                            isfield(in.CONT{setCtr}.(continuationVar{:}),names{fieldCtr})
                        inTemp.(continuationVar{:}).(names{fieldCtr}){1} = continuationParamSave.(names{fieldCtr}) + ...
                            in.CONT{setCtr}.(continuationVar{:}).(names{fieldCtr})(ctrContinuationChange);
                    else
                        % Keep same value
                        inTemp.(continuationVar{:}).(names{fieldCtr}){1} = continuationParamSave.(names{fieldCtr});
                    end
                    
                end
                
            end
            
        elseif (in.cont.method(setCtr) == 2) % Auto changing continuation parameters
            
            %       numRuns;
            
            x0 = x0Save;
            xf = xfSave;
            
            if ctrContinuationChange <= 3
                numDataPoints = 2;
            else
                numDataPoints = 3;
            end
            
            if isfield(in.CONT{setCtr}.constraint,'initial') == 1
                for ctrState = 1 : 1 : numStates
                    if isfield(in.CONT{setCtr}.constraint.initial,inTemp.oc.state{ctrState,1}) == 1
                        if ctrContinuationChange < 3
                            dx(1:2) = linspace(0,in.CONT{setCtr}.step1,2);
                        else
                            asdf
                        end
                        x0(ctrState,1) = x0Save(ctrState,1) + dx(ctrContinuationChange);
                    end
                end
            elseif isfield(in.CONT{setCtr}.constraint,'terminal') == 1
                for ctrState = 1 : 1 : numStates
                    if isfield(in.CONT{setCtr}.constraint.terminal,inTemp.oc.state{ctrState,1}) == 1
                        in.CONT{setCtr}.newValue = in.CONT{setCtr}.constraint.terminal.(inTemp.oc.state{ctrState,1});
                        contMarker = ctrState;
                        if ctrContinuationChange < 3
                            dx.terminal.(inTemp.oc.state{ctrState,1})(1:2) = linspace(0,in.CONT{setCtr}.step1,2);
                        else
                            stepPredict = curvePrediction(numDataPoints,inMod,out,CONT,setCtr,0);
                            if abs(stepPredict) > abs(dxRemaining(ctrContinuationChange-1))
                                stepPredict = dxRemaining(ctrContinuationChange-1);
                            end
                            %               [stepPredict dxRemaining(ctrContinuationChange-1)];
                            dx.terminal.(inTemp.oc.state{ctrState,1})(ctrContinuationChange) = dx.terminal.(inTemp.oc.state{ctrState,1})(ctrContinuationChange-1) + stepPredict;
                        end
                        xf(ctrState,1) = xfSave(ctrState,1) + dx.terminal.(inTemp.oc.state{ctrState,1})(ctrContinuationChange);
                    end
                end
            elseif isfield(in.CONT{setCtr}.constraint,'path') == 1
                for ctrConstraint = 1 : 1 : numConstraints
                    if isfield(in.CONT{setCtr}.constraint.path,inTemp.oc.constraint.path{ctrConstraint,1}) == 1
                        in.CONT{setCtr}.newValue = in.CONT{setCtr}.constraint.path.(inTemp.oc.constraint.path{ctrConstraint,1});
                        contMarker = inTemp.oc.constraint.path{ctrConstraint,1};
                        if ctrContinuationChange < 3
                            dx.path.(inTemp.oc.constraint.path{ctrConstraint,1})(1:2) = linspace(0,in.CONT{setCtr}.step1,2);
                        else
                            stepPredict = curvePrediction(numDataPoints,inMod,out,CONT,setCtr,0);
                            if abs(stepPredict) > abs(dxRemaining(ctrContinuationChange-1))
                                stepPredict = dxRemaining(ctrContinuationChange-1);
                            end
                            %               [stepPredict dxRemaining(ctrContinuationChange-1)];
                            dx.path.(inTemp.oc.constraint.path{ctrConstraint,1})(ctrContinuationChange) = dx.path.(inTemp.oc.constraint.path{ctrConstraint,1})(ctrContinuationChange-1) + stepPredict;
                        end
                        inTemp.constraintVal.(inTemp.oc.constraint.path{ctrConstraint,1}){1} = constraintVal0(ctrConstraint) + dx.path.(inTemp.oc.constraint.path{ctrConstraint,1})(ctrContinuationChange);
                    end
                end
                
            end
            
            %             dx;
            %             [xf xfSave];
            inMod.CONT{setCtr}.dx = dx;
            
            % 			% Input deck variables
            % 			for continuationVar = in.continuationVarSet
            %
            % 				% Determine fieldnames of structure
            % 				if isfield(in,continuationVar)
            % 					names = fieldnames(inTemp.(continuationVar{:}));
            % 				else
            % 					names = {};
            % 				endcontinuationParamSave
            %
            % 				% Update values
            % 				for fieldCtr = 1 : 1 : length(names)
            %
            % 					% Only need to increment variable if exists in continuation structure
            % 					if isfield(in.CONT{setCtr},continuationVar{:}) && ...
            % 						isfield(in.CONT{setCtr}.(continuationVar{:}),names{fieldCtr})
            % 						inTemp.(continuationVar{:}).(names{fieldCtr}){1} = continuationParamSave.(names{fieldCtr}) + ...
            % 						in.CONT{setCtr}.(continuationVar{:}).(names{fieldCtr})(ctrContinuationChange);
            %
            % 					else
            % 						% Keep same value
            % 						inTemp.(continuationVar{:}).(names{fieldCtr}){1} = continuationParamSave.(names{fieldCtr});
            %
            % 					end
            % 				end
            % 			end
            
        end
        % if setCtr == 2
        % 	x0
        % 	xf
        % 	keyboard;
        % end
        %%%%%%%%%%%%%
        %% Scaling %%
        %%%%%%%%%%%%%
        
        if in.autoScale
            
            % Compute scaling values
            % [scaleVal] = computeScalingValues(in,out,numStates,sol);
            [scaleVal] = computeScalingValues(in,out,sol);
            
            % Obtain fieldnames of constants
            if isfield(in,'const')
                names = fieldnames(in.const);
            else
                names = {};
            end
            
            % Determine appropriate scaling
            for ctrConst = 1 : 1 : length(names)
                scaleFn = in.vars.scaleFn.const{ctrConst};
                [scale] = scaleFn(scaleVal{:,2});
                
                inTemp.const.(char(names{ctrConst})){1} = inTemp.const.(char(names{ctrConst})){1}/scale;
                CONT(ctrContinuationChange).scale.const.(char(names{ctrConst})) = scale;
            end
            
            % Obtain fieldnames of constraints
            if isfield(in,'constraintVal')
                names = fieldnames(in.constraintVal);
            else
                names = {};
            end
            
            % Determine appropriate scaling for constraints
            for ctrConstraint = 1 : 1 : length(names)
                scaleFn = in.vars.scaleFn.constraintVal{ctrConstraint};
                [scale] = real(scaleFn(scaleVal{:,2}));
                
                inTemp.constraintVal.(char(names{ctrConstraint})){1} = inTemp.constraintVal.(char(names{ctrConstraint})){1}/scale;
                CONT(ctrContinuationChange).scale.constraintVal.(char(names{ctrConstraint})) = scale;
            end
            
            % Scale states, x0, and xf
            for ctrState = 1 : 1 : numStates
                scaleFn = in.vars.scaleFn.state{ctrState};
                [scale] = scaleFn(scaleVal{:,2});
                
                sol.y(ctrState,:) = sol.y(ctrState,:)/scale;
                x0(ctrState) = x0(ctrState)/scale;
                xf(ctrState) = xf(ctrState)/scale;
                CONT(ctrContinuationChange).scale.state(ctrState) = scale;
            end
            
            % Scale costates (units = COST/STATE)
            for ctrCostate = 1 : 1 : numStates
                scaleFn = in.vars.scaleFn.costate{ctrCostate};
                scale = scaleFn(scaleVal{:,2});
                sol.y(numStates+ctrCostate,:) = sol.y(numStates+ctrCostate,:)/scale;
                CONT(ctrContinuationChange).scale.costate(ctrCostate) = scale;
            end
            
            % Scale parameters
            
            % Determine number of time parameters
            if ~in.convertParametersToStates
                numTimeParam = length(in.vars.arcTypeSequence);
            else
                numTimeParam = 0;
            end
            
            % First set of parameters consists of arc independent variables
            
            % 			units = in.oc.independentVariable{2};
            if ~in.convertParametersToStates
                for ctrParam = 1 : 1 : numTimeParam
                    scaleFn = in.vars.scaleFn.parameters.independentVariable;
                    scale = scaleFn(scaleVal{:,2});
                    sol.parameters(ctrParam) = sol.parameters(ctrParam)/scale;
                end
                CONT(ctrContinuationChange).scale.parameters.independentVariable = scale;
            else
                for ctrParam = 1 : 1 : length(in.vars.arcTypeSequence)
                    scaleFn = in.vars.scaleFn.parameters.independentVariable;
                    scale = scaleFn(scaleVal{:,2});
                    sol.y(2*numStates+ctrParam,:) = sol.y(2*numStates+ctrParam,:)/scale;
                end
                CONT(ctrContinuationChange).scale.parameters.independentVariable = scale;
            end
            
            % Scale nu parameters (units = COST/CONSTRAINT)
            for ctrParam = 1 : 1 : oc.num.constraints.initial
                scaleFn = in.vars.scaleFn.lagrangeMultipliers.initial{ctrParam};
                scale = scaleFn(scaleVal{:,2});
                sol.parameters(numTimeParam+ctrParam) = sol.parameters(numTimeParam+ctrParam)/scale;
                CONT(ctrContinuationChange).scale.lagrangeMultipliers.initial(ctrParam) = scale;
            end
            
            for ctrParam = 1 : 1 : oc.num.constraints.terminal
                scaleFn = in.vars.scaleFn.lagrangeMultipliers.terminal{ctrParam};
                scale = scaleFn(scaleVal{:,2});
                sol.parameters(numTimeParam+oc.num.constraints.initial+ctrParam) = ...
                    sol.parameters(numTimeParam+oc.num.constraints.initial+ctrParam)/scale;
                CONT(ctrContinuationChange).scale.lagrangeMultipliers.terminal(ctrParam) = scale;
            end
            
            % Scale pii parameters from interior point constraints
            interiorPointIndex = numTimeParam+oc.num.constraints.initial+oc.num.constraints.terminal;
            index = 0;
            
            [m,n] = size(in.vars.interiorPointConstraintSequence);
            
            for ctrInteriorPointConstraintCol = 1 : 1 : n
                for ctrInteriorPointConstraintRow = 1 : 1 : m
                    interiorPointType = in.vars.interiorPointConstraintSequence(ctrInteriorPointConstraintRow,ctrInteriorPointConstraintCol);
                    if interiorPointType == 0 || isnan(interiorPointType)
                        continue;
                    end
                    for ctrInteriorPoint = 1 : 1 : oc.num.lagrangeMultipliers.interiorPoint(interiorPointType)
                        scaleFn = in.vars.scaleFn.lagrangeMultipliers.interiorPoint{interiorPointType}{ctrInteriorPoint};
                        scale = scaleFn(scaleVal{:,2});
                        interiorPointIndex = interiorPointIndex + 1;
                        sol.parameters(interiorPointIndex) = sol.parameters(interiorPointIndex)/scale;
                        index = index + 1;
                        CONT(ctrContinuationChange).scale.lagrangeMultipliers.interiorPoint(index) = scale;
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%
        %% Solve BVP %%
        %%%%%%%%%%%%%%%
        
        % Use prior solution as initial guess
        solInit = sol;
        
        % Set timer
        tStart = tic;
        
        % Set options
        options_bvp = in.bvpOptions;
        
        % Save values from constants
        if isfield(in,'const')
            names = fieldnames(in.const);
        else
            names = {};
        end
        
        for ctr1 = 1 : 1 : length(names)
            const.(char(names{ctr1})) = inTemp.const.(char(names{ctr1})){1};
        end
        
        % Save values from constraints
        if isfield(in,'constraintVal')
            names = fieldnames(in.constraintVal);
        else
            names = {};
        end
        
        constraintVal = struct; % intitialize empty structure
        for ctr1 = 1 : 1 : length(names)
            constraintVal.(char(names{ctr1})) = inTemp.constraintVal.(char(names{ctr1})){1};
        end
        
        % Determine mex function used
        strInd = int2str(contIndex);
        
        % Solve BVP based on solver type
        if in.skipUnconvergedSolutions % skipping activated
            try
                errorFlag = false;
                if (in.rootSolving == 0) % using BVP4C
                    % Determine if use and include Jacobian file
                    if in.useJac
                        options_bvp = bvpset(options_bvp,'FJacobian',derivFuncJac_mex(strInd));
                    end
                    % Solve BVP
                    if length(in.vars.arcTypeSequence) == 1
                        sol = bvp4c(derivFunc,bcFunc,solInit,options_bvp,const,constraintVal,in.vars.arcTypeSequence,in.vars.interiorPointConstraintSequence,oc.num.lagrangeMultipliers.interiorPoint,x0,xf);
                        % sol = bvp4c_rdsl(derivFunc,bcFunc,solInit,options_bvp,const,constraintVal,in.vars.arcTypeSequence,in.vars.interiorPointConstraintSequence,oc.num.lagrangeMultipliers.interiorPoint,x0,xf);
                    else
                        sol = bvp4c(derivFuncRegion,bcFunc,solInit,options_bvp,const,constraintVal,in.vars.arcTypeSequence,in.vars.interiorPointConstraintSequence,oc.num.lagrangeMultipliers.interiorPoint,x0,xf);
                        % sol = bvp4c_rdsl(derivFuncRegion,bcFunc,solInit,options_bvp,const,constraintVal,in.vars.arcTypeSequence,in.vars.interiorPointConstraintSequence,oc.num.lagrangeMultipliers.interiorPoint,x0,xf);
                    end
                    sol.stats; % this line is used to catch situations where a warning is provided by bvp4c that the node limit is reached
                elseif (in.rootSolving == 1) % using GPU
                    % Solve BVP
                    [sol, MNPhi] = bvpgpu(derivFuncRegionGPU,bcFunc,solInit,in,out,const,constraintVal,in.vars.arcTypeSequence,in.vars.interiorPointConstraintSequence,oc.num.lagrangeMultipliers.interiorPoint,real(x0),real(xf));
                    sol.MNPhi = MNPhi;
                else
                    error('Incorrect value for root solver variable, in.rootSolving\n');
                end
            catch
                % skip over solution
                errorFlag = true;
                sol = solInit;
            end
        else % skipping disabled
            errorFlag = false;
            if (in.rootSolving == 0) % using BVP4C
                % Determine if use and include Jacobian file
                if in.useJac
                    options_bvp = bvpset(options_bvp,'FJacobian',derivFuncJac);
                end
                % Solve BVP
                if length(in.vars.arcTypeSequence) == 1
                    sol = bvp4c(derivFunc,bcFunc,solInit,options_bvp,const,constraintVal,in.vars.arcTypeSequence,in.vars.interiorPointConstraintSequence,oc.num.lagrangeMultipliers.interiorPoint,x0,xf);
                % sol = bvp4c_rdsl(derivFunc,bcFunc,solInit,options_bvp,const,constraintVal,in.vars.arcTypeSequence,in.vars.interiorPointConstraintSequence,oc.num.lagrangeMultipliers.interiorPoint,x0,xf);
                else
                    sol = bvp4c(derivFuncRegion,bcFunc,solInit,options_bvp,const,constraintVal,in.vars.arcTypeSequence,in.vars.interiorPointConstraintSequence,oc.num.lagrangeMultipliers.interiorPoint,x0,xf);
                    % sol = bvp4c_rdsl(derivFuncRegion,bcFunc,solInit,options_bvp,const,constraintVal,in.vars.arcTypeSequence,in.vars.interiorPointConstraintSequence,oc.num.lagrangeMultipliers.interiorPoint,x0,xf);
                end
                % try
                % 	sol.stats
                % catch
                % 	keyboard
                % end
            elseif (in.rootSolving == 1) % using GPU
                % Solve BVP
                [sol, MNPhi] = bvpgpu(derivFuncRegionGPU,bcFunc,solInit,in,out,const,constraintVal,in.vars.arcTypeSequence,in.vars.interiorPointConstraintSequence,oc.num.lagrangeMultipliers.interiorPoint,real(x0),real(xf));
                sol.MNPhi = MNPhi;
            else
                error('Incorrect value for root solver variable, in.rootSolving\n');
            end
        end
        
        % Assign error flag
        sol.errorFlag = errorFlag;
        
        % Find execution time
        runTime = toc(tStart);
        
        %%%%%%%%%%%%%%%%%%%%%
        %% Unscale Problem %%
        %%%%%%%%%%%%%%%%%%%%%
        
        if in.autoScale
            
            % Obtain fieldnames of constants
            if isfield(in,'const')
                names = fieldnames(in.const);
            else
                names = {};
            end
            
            % Determine appropriate scaling
            for ctrConst = 1 : 1 : length(names)
                % 				units = in.const.(char(names{ctrConst})){2};
                scale = CONT(ctrContinuationChange).scale.const.(char(names{ctrConst}));
                inTemp.const.(char(names{ctrConst})){1} = inTemp.const.(char(names{ctrConst})){1}*scale;
            end
            
            % Obtain fieldnames of constraints
            if isfield(in,'constraintVal')
                names = fieldnames(in.constraintVal);
            else
                names = {};
            end
            
            % Determine appropriate scaling
            for ctrConstraint = 1 : 1 : length(names)
                % 				units = in.constraintVal.(char(names{ctrConstraint})){2};
                scale = CONT(ctrContinuationChange).scale.constraintVal.(char(names{ctrConstraint}));
                inTemp.constraintVal.(char(names{ctrConstraint})){1} = inTemp.constraintVal.(char(names{ctrConstraint})){1}*scale;
            end
            
            % Scale states
            for ctrState = 1 : 1 : numStates
                % 				units = in.oc.state{ctrState,2};
                scale = CONT(ctrContinuationChange).scale.state(ctrState);
                sol.y(ctrState,:) = sol.y(ctrState,:)*scale;
                x0(ctrState) = x0(ctrState)*scale;
                xf(ctrState) = xf(ctrState)*scale;
            end
            
            % Scale costates (units = COST/STATE)
            for ctrCostate = 1 : 1 : numStates
                % 				units = ['(',in.oc.cost.path{2},')/(',in.oc.state{ctrCostate,2},')'];
                scale = CONT(ctrContinuationChange).scale.costate(ctrCostate);
                sol.y(numStates+ctrCostate,:) = sol.y(numStates+ctrCostate,:)*scale;
            end
            
            % First set of parameters consists of arc times
            if ~in.convertParametersToStates
                % 				units = in.oc.independentVariable{2};
                for ctrParam = 1 : 1 : numTimeParam
                    scale = CONT(ctrContinuationChange).scale.parameters.independentVariable;
                    sol.parameters(ctrParam) = sol.parameters(ctrParam)*scale;
                end
            else
                % 				units = in.oc.independentVariable{2};
                for ctrParam = 1 : 1 : length(in.vars.arcTypeSequence)
                    scale = CONT(ctrContinuationChange).scale.parameters.independentVariable;
                    sol.y(2*numStates+ctrParam,:) = sol.y(2*numStates+ctrParam,:)*scale;
                end
            end
            
            % Scale nu parameters (units = COST/CONSTRAINT)
            for ctrParam = 1 : 1 : oc.num.constraints.initial
                scale = CONT(ctrContinuationChange).scale.lagrangeMultipliers.initial(ctrParam);
                sol.parameters(numTimeParam+ctrParam) = sol.parameters(numTimeParam+ctrParam)*scale;
            end
            
            for ctrParam = 1 : 1 : oc.num.constraints.terminal
                scale = CONT(ctrContinuationChange).scale.lagrangeMultipliers.terminal(ctrParam);
                sol.parameters(numTimeParam+oc.num.constraints.initial+ctrParam) = ...
                    sol.parameters(numTimeParam+oc.num.constraints.initial+ctrParam)*scale;
            end
            
            % Scale pii parameters
            numPii = length(sol.parameters(numTimeParam+oc.num.constraints.initial+oc.num.constraints.terminal+1:end));
            for ctrPii = 1 : 1 : numPii
                
                scale = CONT(ctrContinuationChange).scale.lagrangeMultipliers.interiorPoint(ctrPii);
                sol.parameters(numTimeParam+oc.num.constraints.initial+oc.num.constraints.terminal+ctrPii) = ...
                    sol.parameters(numTimeParam+oc.num.constraints.initial+oc.num.constraints.terminal+ctrPii)*scale;
                
            end
        end
        
        %%%%%%%%%%%%
        %% Output %%
        %%%%%%%%%%%%
        
        % Save values from constants
        if isfield(in,'const')
            names = fieldnames(in.const);
        else
            names = {};
        end
        
        for ctrConst = 1 : 1 : length(names)
            const.(char(names{ctrConst})) = inTemp.const.(char(names{ctrConst})){1};
        end
        
        % Save values from constraints
        if isfield(in,'constraintVal')
            names = fieldnames(in.constraintVal);
        else
            names = {};
        end
        
        constraintVal = struct; % intitialize empty structure
        for ctrConstraint = 1 : 1 : length(names)
            constraintVal.(char(names{ctrConstraint})) = inTemp.constraintVal.(char(names{ctrConstraint})){1};
            %             sol.(char(names{ctrConstraint})) = constraintVal.(char(names{ctrConstraint}));
        end
        
        % Compute control history and unconstrained sufficiency condition
        % Compute control history
        sol.control = NaN(numControl,length(sol.x));
        % sol.d2Hdu2 = NaN(numControl,length(sol.x));
        sol.d2Hdu2 = NaN(1,length(sol.x));
        controlOutput = cell(1,1+numControl);
        
        % Find repeating elements in sol.x to find arc-junctions
        repVals = [0 find(diff(sol.x) == 0)];
        
        for arcCtr = 1:length(repVals)
            if in.vars.arcTypeSequence(arcCtr) == 0
                controlFnStr = 'computeControlUnconstrained';
            else
                controlFnStr = ['computeControlConstraint',num2str(in.vars.arcTypeSequence(arcCtr))];
            end
            if in.useMex
                controlFnStr = [controlFnStr,'_mex'];
            end
            controlFn = str2func(controlFnStr);
            startIndex = repVals(arcCtr)+1;
            if arcCtr == length(repVals)
                endIndex = length(sol.x);
            else
                endIndex = repVals(arcCtr+1);
            end
            for timeCtr = startIndex : 1 : endIndex
                
                % Scale parameters, necessary for evaluation of bang/bang and singular
                % arc logic. Need to generalize.
                
                
                
                
                % Obtain fieldnames of constants
                if isfield(in,'const')
                    names = fieldnames(in.const);
                else
                    names = {};
                end
                
                % Determine appropriate scaling
                constControl = struct;
                for ctrConst = 1 : 1 : length(names)
                    constControl.(names{ctrConst}) = const.(names{ctrConst})/ ...
                        CONT(ctrContinuationChange).scale.const.(char(names{ctrConst}));
                end
                
                % Obtain fieldnames of constraints
                if isfield(in,'constraintVal')
                    names = fieldnames(in.constraintVal);
                else
                    names = {};
                end
                
                % Determine appropriate scaling for constraints
                constraintValControl = struct;
                for ctrConstraint = 1 : 1 : length(names)
                    
                    constraintValControl.(char(names{ctrConstraint})) = constraintVal.(char(names{ctrConstraint}))/ ...
                        CONT(ctrContinuationChange).scale.constraintVal.(char(names{ctrConstraint}));
                end
                
                % Scale states, x0, and xf
                solControl = [];
                for ctrState = 1 : 1 : numStates
                    
                    solControl.y(ctrState,:) = sol.y(ctrState,:)/ ...
                        CONT(ctrContinuationChange).scale.state(ctrState);
                    
                end
                
                % Scale costates (units = COST/STATE)
                for ctrCostate = 1 : 1 : numStates
                    
                    solControl.y(numStates+ctrCostate,:) = sol.y(numStates+ctrCostate,:)/ ...
                        CONT(ctrContinuationChange).scale.costate(ctrCostate);
                    
                end
                
                if ~in.convertParametersToStates
                    [controlOutput{:}] = controlFn(solControl.y(:,timeCtr),constControl,constraintValControl,length(in.vars.arcTypeSequence));
                else
                    % [controlOutput{:}] = controlFn([solControl.y(:,timeCtr);zeros(in.maxNumArcs-length(in.vars.arcTypeSequence),1)],constControl,constraintValControl,length(in.vars.arcTypeSequence));
                    [controlOutput{:}] = controlFn([solControl.y(:,timeCtr);zeros(in.maxNumArcs,1)],constControl,constraintValControl,length(in.vars.arcTypeSequence));
                end
                for ctr1 = 1 : 1 : numControl
                    sol.control(ctr1,timeCtr) = controlOutput{ctr1};
                end
                % sol.d2Hdu2(:,timeCtr) = controlOutput{numControl+1};
            end
        end
        % Continuation data
        %     in = inTemp;
        CONT(ctrContinuationChange).in = inTemp;
        % clear outTemp.setCONT;
        % CONT(ctrContinuationChange).out = outTemp;
        CONT(ctrContinuationChange).sol = sol;
        
        % Adjust solution to NaNs if error. Reduces saved file size and will help
        % with plotting logic.
        if errorFlag
            CONT(ctrContinuationChange).sol.x = NaN;
            CONT(ctrContinuationChange).sol.y = NaN(2*numStates,1);
            CONT(ctrContinuationChange).sol.parameters = NaN(size(sol.parameters));
        end
        
        % Run time
        CONT(ctrContinuationChange).runTime = runTime;
        
        % Output to screen
        if ~sol.errorFlag
            fprintf('Continuation solution %d/%d converged in %g seconds.\n\n',ctrContinuationChange,numRuns,runTime);
        elseif sol.errorFlag
            fprintf('Solver failed for solution %d/%d. Time spent is %g seconds.\n\n',ctrContinuationChange,numRuns,runTime);
        end
        % save -v7.3 data/results.mat in out;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Automated Continuation BC Check %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if (in.cont.method(setCtr) == 2)
            if isfield(in.CONT{setCtr}.constraint,'initial')
                currentBC = x0(contMarker,1);
            elseif isfield(in.CONT{setCtr}.constraint,'terminal')
                currentBC = xf(contMarker,1);
            elseif isfield(in.CONT{setCtr}.constraint,'path')
                currentBC = constraintVal.(contMarker);
            end
            
            dxRemaining(ctrContinuationChange) = in.CONT{setCtr}.newValue - currentBC;
            
            if abs(dxRemaining(ctrContinuationChange)/dxRemaining(1)) > in.gpuSolve.Tol
                numRuns = numRuns + 1;
            end
            
            remainingPercentage = ((1-dxRemaining(ctrContinuationChange)/dxRemaining(1))*100);
            remainingPercentageString = num2str(remainingPercentage);
            
            stop_case = 1000;  % safety/test value to limit # of steps in automation
            
            if numRuns > stop_case
                numRuns = stop_case;
            end
            
        end
        
        % Output to screen
        if in.cont.method(setCtr) == 1
            % fprintf('Continuation solution %d/%d converged in %g seconds\n\n', ...
            % ctrContinuationChange,numRuns,runTime);
        elseif in.cont.method(setCtr) == 2
            fprintf('Continuation solution is %s percent converged in %g seconds\n\n', ...
                remainingPercentageString,runTime);
        end
        % step counter
        
        ctrContinuationChange = ctrContinuationChange + 1;
        
    end % current continuation set
    
    clearvars dx
    
    if in.cont.method(setCtr) == 2
        fprintf('Continuation process completed in %d steps\n\n', ...
            ctrContinuationChange);
    end
    
    % Output of full continuation process
    out.setCONT(setCtr).CONT = CONT;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute Maximum Value of Constraint at End of Continuation Process %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    newActiveConstraint = false;
    if (setCtr < in.cont.indexContBeforeTrades) && (in.cont.method(setCtr) == 1) && setCtr ~=in.vars.numContinuations && ...
            isfield(in.CONT{setCtr+1},'constraint') && ...
            (isfield(in.CONT{setCtr+1}.constraint,'path') || isfield(in.CONT{setCtr+1}.constraint,'interiorPoint'))
        
        % Determine type of constraint
        if isfield(in.CONT{setCtr+1}.constraint,'path')
            numConstraints = length(in.CONT{setCtr+1}.constraint.path);
            direction = cell(numConstraints,1);
            constraintExpression = cell(numConstraints,1);
            for ctrConstraint = 1 : 1 : numConstraints
                constraintType = 'path';
                direction{ctrConstraint} = inTemp.oc.constraint.path{in.CONT{setCtr+1}.constraint.path(ctrConstraint),2};
                constraintExpression{ctrConstraint} = inTemp.oc.constraint.path{in.CONT{setCtr+1}.constraint.path(ctrConstraint),3};
            end
        elseif isfield(in.CONT{setCtr+1}.constraint,'interiorPoint')
            numConstraints = length(in.CONT{setCtr+1}.constraint.interiorPoint);
            direction = cell(numConstraints,1);
            constraintExpression = cell(numConstraints,1);
            
            for ctrConstraint = 1 : 1 : numConstraints
                constraintType = 'interior';
                direction{ctrConstraint} = inTemp.oc.constraint.interiorPoint{in.CONT{setCtr+1}.constraint.interiorPoint(ctrConstraint),2};
                constraintExpression{ctrConstraint} = inTemp.oc.constraint.interiorPoint{in.CONT{setCtr+1}.constraint.interiorPoint(ctrConstraint),3};
            end
        end
        
        % Compute all constants
        names = fieldnames(inTemp.const);
        for ctrNames = 1 : 1 : length(names)
            eval([char(names{ctrNames}),'= inTemp.const.',char(names{ctrNames}),'{1};']);
        end
        
        % Determine maximum/minimum constraint value. If multiple constraints, must be in same direction.
        maximize = NaN(numConstraints,1);
        for ctrConstraint = 1 : 1 : numConstraints
            if strcmp(direction{ctrConstraint},'>')
                maximize(ctrConstraint) = true;
                constraintValueTotal = -inf;
            elseif strcmp(direction{ctrConstraint},'<')
                maximize(ctrConstraint) = false;
                constraintValueTotal = inf;
            end
        end
        
        for ctrConstraint = 1 : 1 : numConstraints-1
            if maximize(ctrConstraint) ~= maximize(ctrConstraint+1)
                error('When specifying multiple constraints simultaneously, direction (maximize or minimize) must be the same.');
            end
        end
        
        timeNormalizedConstraintVal = 0;
        
        for trajCtr = 1 : 1 : length(CONT(end).sol.y(1,:))
            
            % Compute states
            for ctrStates = 1 : 1 : length(inTemp.oc.state)
                eval([inTemp.oc.state{ctrStates,1},' = CONT(end).sol.y(ctrStates,trajCtr);']);
            end
            
            for ctrControls = 1 : 1 : size(inTemp.oc.control,1)
                eval([inTemp.oc.control{ctrControls,1},' = CONT(end).sol.control(ctrControls,trajCtr);'])
            end
            
            % Compute constraint value for each constraint and aggregate
            val = NaN(numConstraints,1);
            for ctrConstraint = 1 : 1 : numConstraints
                val(ctrConstraint) = eval(constraintExpression{ctrConstraint});
            end
            valTotal = sum(val);
            
            if maximize
                if valTotal > constraintValueTotal
                    constraintValueTotal = valTotal;
                    constraintVal = val;
                    timeNormalizedConstraintVal = sol.x(trajCtr);
                end
            elseif ~maximize
                if valTotal < constraintValueTotal
                    constraintValueTotal = valTotal;
                    constraintVal = val;
                    timeNormalizedConstraintVal = sol.x(trajCtr);
                end
            end
            
        end

        newActiveConstraint = true;
        if strcmp(constraintType,'path')
            
            % Due to discrete approximation of trajectory, sometimes qdotMaxTraj is
            % directly at a point of tangency. Need to artifically reduce this to get a
            % trajectory point on the arc.
            if maximize
                constraintVal
                constraintVal = constraintVal - 0.00001*sign(constraintVal)*constraintVal % 0.01% reduction
                
            elseif ~maximize
                constraintVal
                %keyboard
                constraintVal = constraintVal + 0.00001*sign(constraintVal)*constraintVal % 0.01% addition
            end
            constraintVal0 = constraintVal;
            
        end
        
    end
    
end

setCONT = out.setCONT;
return

