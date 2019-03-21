function [min_step_predict] = curvePrediction(numData,in,out,CONT,setCtr,postCalc)

warning off all

numSTT = 1; % highest order of STT information to fit

if numData < 4
    doEnhancePolyfit = true; % run extra fit type, only for 2 or 3 data points
else
    doEnhancePolyfit = false;
end

if postCalc == 1
    close all; clc;
    format long
    postData = true; % load data file
    doPlots = true; % display fits
    predict_range = 0.6; % determines how far post data to predict (1 unit = data set range)
else
    postData = false; % use input data
    doPlots = false; % don't generate plots
    predict_range = 1.0;
end

min_error_scale = 1;
max_error_scale = 30;

%%%%%%%%%%%% 
%% Inputs %% 
%%%%%%%%%%%% 
    
ctrSet = setCtr; % looks at current CONT set only

% Load Data
if postData == true
    data = load('data/results.mat'); % import data file
else
    data.in = in; % uses function inputs
    data.out = out;
    data.out.setCONT(ctrSet).CONT = CONT;
end

% Clear variables 
clear scale constraint const traj VEH CONT ID setCONT 

% Checks current system variables
numTimesteps = data.in.gpuSolve.TimeSteps; % number of timesteps per arc
arcSteps = numTimesteps + 1; % adds the 'zero' value to arc size
numArcs = data.out.setCONT(ctrSet).CONT(1).out.numArcs; % number of arcs
numStates = data.in.oc.num.states; % number of state
numParameters = data.in.oc.num.parameters; % number of parameters (related to fixed initial & terminal states only)
allStates = data.in.oc.state; % list of all states
if numArcs > 1
    allConstraintVariables = data.in.oc.constraintVariable; % list of all constraint variables, doesn't exist if no path constraints are applied
    allVariables = [allStates; allConstraintVariables]; % variable labels
else
    allVariables = allStates;
end
fixedInitialStates = data.in.oc.constraint.initial; % list of fixed initial states
fixedTerminalStates = data.in.oc.constraint.terminal; % list of fixed final states
markerInitialFixTF = ismember(allStates,fixedInitialStates); % compare all and fixed states for T/F values
markerTerminalFixTF = ismember(allStates,fixedTerminalStates); % compare all and fixed states for T/F values
fixedInitialStatesCount = data.in.oc.num.constraints.initial; % number of fixed initial states
fixedTerminalStatesCount = data.in.oc.num.constraints.terminal; % number of fixed terminal states
numPathConstraints = data.in.oc.num.constraints.path; % number of path constraints
%fixedInitialCostatesCount = max(size(markerInitialFixTF(:,1))) - fixedInitialStatesCount; % number of fixed initial costates

numVariables = numStates + numPathConstraints; % adds in path constraints variables, like qdot, if present
numLagrangeMultipliers = numParameters + numPathConstraints; % !!!! Replace with number of PIs, currently assumes 1 per path constraint

% Pull information from step vector name
if isfield(data.in.CONT{ctrSet}.dx,'initial')
    contVariableLocation = 'initial';
elseif isfield(data.in.CONT{ctrSet}.dx,'path')
    contVariableLocation = 'path';
elseif isfield(data.in.CONT{ctrSet}.dx,'interior')
    contVariableLocation = 'interior';
elseif isfield(data.in.CONT{ctrSet}.dx,'terminal')
    contVariableLocation = 'terminal';
end

% Pull variable name and index, index only applies for state variable stepping
for i = 1 : 1 : numVariables
    if isfield(data.in.CONT{ctrSet}.dx.(contVariableLocation),allVariables{i})
        contVariable = allVariables{i};
        contVariableIndex = i;
    end
end

% Determines if the variable being stepped is a state
if ismember(contVariable,allStates)
    contVariableType = 'state';
else
    contVariableType = 'other';
end

% Generate Fitting functions
quadfunname = ['STT1_Quad',num2str(numData),'.m']; % store function file name
quadfun = ['STT1_Quad',num2str(numData)]; % creates function handle
exist(quadfunname); % check to see if curvefit function file exists
if exist(quadfunname) == 0
    STT1_quad(numData);
end
exist(quadfunname);

if doEnhancePolyfit == true
    epolyfunname = ['Enhance_Polyfit',num2str(numData),num2str(numSTT),'.m'];
    epolyfun = ['Enhance_Polyfit',num2str(numData),num2str(numSTT)];
    exist(epolyfunname);
    if exist(epolyfunname) == 0
        enhance_polyfit(numData,numSTT);
    end
    exist(epolyfunname);
end

%%%%%%%%%%%%%%%%% 
%% Assign Data %% 
%%%%%%%%%%%%%%%%% 

ctrCont_max = length(data.out.setCONT(ctrSet).CONT); % Determine number of continuation steps performed

j=1; % counter for # of solutions
for ctrCont = (ctrCont_max - numData + 1) : 1 : ctrCont_max 
          
    i=1; % counter for # of trends
    
    % Assign states and costates 
    for ctrArc = 1 : 1 : numArcs
        setCONT(ctrSet).CONT(j).x0_unscaled((i):(numStates*2+numArcs+i-1)) = data.out.setCONT(ctrSet).CONT(ctrCont).sol.y(:,1+arcSteps*(ctrArc-1));  
        setCONT(ctrSet).CONT(j).xf_unscaled((i):(numStates*2+numArcs+i-1)) = data.out.setCONT(ctrSet).CONT(ctrCont).sol.y(:,arcSteps*ctrArc);  
        scales(j,(i):(numStates+i-1)) = data.out.setCONT(ctrSet).CONT(ctrCont).scale.state(:);
        scales(j,(numStates+i):(numStates*2+i-1)) = data.out.setCONT(ctrSet).CONT(ctrCont).scale.costate(:);
        scales(j,(numStates*2+i):(numStates*2+numArcs+i-1)) = data.out.setCONT(ctrSet).CONT(ctrCont).scale.parameters.independentVariable;
        i=i+2*numStates+numArcs;
    end
    
    % Grab paramaters value and scaling factor
    setCONT(ctrSet).CONT(j).x0_unscaled(i:i+numLagrangeMultipliers-1) = data.out.setCONT(ctrSet).CONT(ctrCont).sol.parameters;
    setCONT(ctrSet).CONT(j).xf_unscaled(i:i+numLagrangeMultipliers-1) = data.out.setCONT(ctrSet).CONT(ctrCont).sol.parameters;
    scales(j,i:i+fixedInitialStatesCount-1) = data.out.setCONT(ctrSet).CONT(ctrCont).scale.lagrangeMultipliers.initial(:);
    scales(j,i+fixedInitialStatesCount:i+numParameters-1) = data.out.setCONT(ctrSet).CONT(ctrCont).scale.lagrangeMultipliers.terminal(:);
    if numArcs > 1
        scales(j,i+numParameters:i+numLagrangeMultipliers-1) = data.out.setCONT(ctrSet).CONT(ctrCont).scale.lagrangeMultipliers.interiorPoint(:);
    end
    i=i+numLagrangeMultipliers;
    
    % Grab path constraint value
    if strcmp(contVariableType,'other') == 1
        constraint_x0_unscaled(j) = data.out.setCONT(ctrSet).CONT(ctrCont).sol.(contVariable);
    end
    
    % M+N*Phi data, generated from scaled states, costates, and parameters
    MNPhi(:,:,j) = data.out.setCONT(ctrSet).CONT(ctrCont).sol.MNPhi;
    
    % Assign variable being stepped
    contValues(j) = data.in.CONT{ctrSet}.dx.(contVariableLocation).(contVariable)(ctrCont);
    if strcmp(contVariableType,'other') == 1
        contValuesScales(j) = data.out.setCONT(ctrSet).CONT(ctrCont).scale.constraintVal.(contVariable);
    end

    j=j+1;
          
end

% Translates problem to start at 0 and scales stepping values
contValuesOffset = contValues(1);
contValues = contValues - contValuesOffset;
if strcmp(contVariableType,'other') == 1
    contValues = contValues.*contValuesScales;
else
    contValues = contValues./scales(:,contVariableIndex)';
end


%%%%%%%%%%%%%%%%%%%%%%% 
%% Data Calculations %% 
%%%%%%%%%%%%%%%%%%%%%%%

% Puts all data into scaled values
for i = 1 : numData 
    A(i,:) = setCONT(ctrSet).CONT(i).x0_unscaled;
    B(i,:) = setCONT(ctrSet).CONT(i).xf_unscaled;
    M(i,:) = setCONT(ctrSet).CONT(i).x0_unscaled./scales(i,:);
end

if strcmp(contVariableType,'state') == 1
    predictionScaleV = ((A(:,contVariableIndex)+B(:,contVariableIndex))/2)./scales(:,contVariableIndex);
elseif strcmp(contVariableType,'other') == 1
    predictionScaleV = constraint_x0_unscaled./contValuesScales;
end

predictionScaleV;
predictionScale = abs(sum(predictionScaleV)/numData);

if strcmp(contVariableType,'other') == 1
    predictionUnscale = abs(sum(constraint_x0_unscaled)/numData);
end

% Search method to determine which BC is changed by stepping
initialBCcount=0;
terminalBCcount=0;

for i = 1 : numStates
   initialBCcount = initialBCcount + markerInitialFixTF(i,1);
   terminalBCcount = terminalBCcount + markerTerminalFixTF(i,1);
   if strcmp(allStates(i,1),contVariable) == 1 && strcmp(contVariableLocation,'initial') == 1
       contVariableIndexMarker = initialBCcount;
   elseif strcmp(allStates(i,1),contVariable) == 1 && strcmp(contVariableLocation,'terminal') == 1
       contVariableIndexMarker = terminalBCcount;
   end
   initialBCcount = initialBCcount + 1;
   terminalBCcount = terminalBCcount + 1;
end

if strcmp(contVariableLocation,'path') == 1
    contVariableIndexMarker = initialBCcount + terminalBCcount + numArcs + 1;
end

if strcmp(contVariableLocation,'initial') == 1
    contVariableBCIndex = contVariableIndexMarker;
elseif strcmp(contVariableLocation,'terminal') == 1
    contVariableBCIndex = initialBCcount + contVariableIndexMarker;
elseif strcmp(contVariableLocation,'path') ==1
    contVariableBCIndex = contVariableIndexMarker;
end
totalBCcount = initialBCcount + terminalBCcount + numArcs + numPathConstraints + (numArcs-1)*(2*numStates+numArcs);

% Calculate and select derivative values
error = zeros(totalBCcount,1);
error(contVariableBCIndex) = 1;

for i = 1 : numData 
    PSI1(:,:,i) = inv(MNPhi(:,:,i));
    derivative(:,i) = PSI1(:,:,i)*error;
end
Derivative = derivative';


%%%%%%%%%%%%%%%%%%% 
%% Curve Fitting %% 
%%%%%%%%%%%%%%%%%%%

% Determine range to view fits
fitRangeV = [contValues(1) contValues(numData)]; % fit range vector [initial final]
fitRange = fitRangeV(2) - fitRangeV(1); % fit range value
fitRangeScale = 1 + roundn(predict_range,-2); % How long post final point prediction vs data range
fitPoints = 100; % number of datapoints over solutions range
    
fitX = linspace(fitRangeV(1), fitRangeV(1) + fitRange*fitRangeScale, fitPoints*fitRangeScale+1)';
[fitx_points, drop] = size(fitX);
    
% Polyfit data constants
for i = 1 : totalBCcount
    regFit(:,i) = polyfit(contValues',M(:,i),numData-1);

end

% Calculate regular polyfit curves
for i = 1 : totalBCcount
    regFitY(:,i) = regFit(1,i)*ones(fitPoints*fitRangeScale+1,1);
end
     
for i = 1 : totalBCcount
    for j = 1 : numData-1
        regFitY(:,i) = regFitY(:,i).*fitX + regFit(j+1,i)*ones(fitPoints*fitRangeScale+1,1);
    end
end
    
%%%%% Quadrature Fit
for i = 1 : totalBCcount
    QuadData(:,i) = [contValues'; M(numData,i); Derivative(:,i)];
    Quadfun = str2func(quadfun);
    QuadY(:,i) = Quadfun(QuadData(:,i),fitX)';
end

for i = 1 : totalBCcount
    errormag(i) = QuadY(fitPoints+1,i) - M(numData,i);
    errorrat(i) = 1-(QuadY(fitPoints+1,i)/M(numData,i));
    QuadY(:,i) = QuadY(:,i) - errormag(i);
end
   
%%%%% Enhance Poly Fit
if doEnhancePolyfit == true
    for i = 1 : totalBCcount
        EpolyData(:,i) = [contValues'; M(:,i); Derivative(:,i)];
        Epolyfun = str2func(epolyfun);
        EpolyY(:,i) = Epolyfun(EpolyData(:,i),fitX)';
    end
    
    for i = 1 : totalBCcount
        errormag(i) = EpolyY(fitPoints+1,i) - M(numData,i);
        errorrat(i) = 1-(EpolyY(fitPoints+1,i)/M(numData,i));
        EpolyY(:,i) = EpolyY(:,i) - errormag(i);
    end
end
   
% Calculate STT1 prediction curves
for i = 1 : numData
    STT1X(:,i) = linspace(contValues(i), contValues(i) + fitRange*(fitRangeScale-1), fitPoints*(fitRangeScale-1)+1)';
    for j = 1 : totalBCcount
        STT1Y(:,i,j) = M(i,j) + Derivative(i,j) * (STT1X(:,i)-contValues(i));
    end
end
  


%%%%%%%%%%%%%%%%%%%%% 
%% Step Prediction %% 
%%%%%%%%%%%%%%%%%%%%%

% Calculate relative error of each fit vs STT1 prediction
for i = 1 : totalBCcount
    regFit_error(:,i) = abs((STT1Y(:,numData,i)-regFitY((fitPoints+1):end,i))/((M(numData,i)+M(1,i))/2));
    Quad_error(:,i) = abs((STT1Y(:,numData,i)-QuadY((fitPoints+1):end,i))/((M(numData,i)+M(1,i))/2));
    if doEnhancePolyfit == true
        Epoly_error(:,i) = abs((STT1Y(:,numData,i)-EpolyY((fitPoints+1):end,i))/((M(numData,i)+M(1,i))/2));
    end
end
  
% Determines when enhance poly becomes largest error   
max_error_tolerance = 1*numData; % max error keeps cap on exploding predictions
%max_error_tolerance = data.in.gpuSolve.Tol*max_error_scale/sqrt(predictionScale); % max error keeps cap on exploding predictions
min_error_tolerance = data.in.gpuSolve.Tol*min_error_scale; % min error keeps above noise level
%min_error_tolerance = data.in.gpuSolve.Tol*min_error_scale/predictionScale;

for i = 1 : totalBCcount
    for j = 1 : fitPoints*(fitRangeScale-1)+1
        if doEnhancePolyfit == true
            if abs(regFit_error(j,i)) > max_error_tolerance || abs(Epoly_error(j,i)) > max_error_tolerance
                find_error1(j,i) = sign(-1);
            elseif abs(regFit_error(j,i)) < min_error_tolerance || abs(Epoly_error(j,i)) < min_error_tolerance
                find_error1(j,i) = sign(1);
            else
                find_error1(j,i) = sign(regFit_error(j,i)-Epoly_error(j,i));
            end
            if abs(Quad_error(j,i)) > max_error_tolerance || abs(Epoly_error(j,i)) > max_error_tolerance
                find_error2(j,i) = sign(-1);
            elseif abs(Quad_error(j,i)) < min_error_tolerance || abs(Epoly_error(j,i)) < min_error_tolerance
                find_error2(j,i) = sign(1);
            else
                find_error2(j,i) = sign(Quad_error(j,i)-Epoly_error(j,i));
            end
        else
            if abs(regFit_error(j,i)) > max_error_tolerance || abs(Quad_error(j,i)) > max_error_tolerance
                find_error3(j,i) = sign(-1);
            elseif abs(regFit_error(j,i)) < min_error_tolerance || abs(Quad_error(j,i)) < min_error_tolerance
                find_error3(j,i) = sign(1);
            else
                find_error3(j,i) = sign(regFit_error(j,i)-Quad_error(j,i));
            end
        end
    end
    
    % Combine error results if using Enhanced Polyfit
    if doEnhancePolyfit == true
        find_errorS(:,i) = sign(find_error1(:,i) + find_error2(:,i));
    else
        find_errorS(:,i) = sign(find_error3(:,i));
    end
    [value(i),step_predict_index(i)] = min(find_errorS(2:end,i));
    
    % Determine index of max alloted step size for each variable
    if value(i) > -1
        step_predict_index(i) = fitx_points;
    else
        step_predict_index(i) = step_predict_index(i) + fitPoints;
    end
    step_predict(i) = fitX(step_predict_index(i))-fitX(fitPoints+1);
end

%find_error1
%find_error2

% Determines which predictions relate to BC of interest
for i = 1 : numStates
    if markerInitialFixTF(i,1) == 1
        step_predictImportant(i) = i + numStates;
    else
        step_predictImportant(i) = i;
    end
end

% Reduces prediction results to just those of interest
step_predictImportant(numStates+1) = 2*numStates + 1;
step_predictImportant = sort(step_predictImportant);
step_predict = step_predict(step_predictImportant)*scales(numData,contVariableIndex);
step_predict_index = step_predict_index(step_predictImportant);

% Select min step from step prediction of each variable
if step_predict(1) > 0
    min_step_predict = min(step_predict);
else
    min_step_predict = max(step_predict);
end

if strcmp(contVariableType,'other') == 1
    min_step_predict = min_step_predict/contValuesScales(numData);
end

if postData == true
    numberIndex = linspace(1,totalBCcount,totalBCcount);
    [step_predictImportant; step_predict_index; step_predict]'
    min_step_predict*[1 180/pi]
end


%%%%%%%%%%%%%%%%%% 
%% Fit Plotting %% 
%%%%%%%%%%%%%%%%%%

lineWidth = 1;
markerSize = 24;
fontSize = 16;

% Actual Plots
if doPlots == true
for i = 1 : totalBCcount
    if ismember(i,step_predictImportant) == 1
        
        figure(i)
        hold on
        grid on
        set(gca,'FontSize',fontSize)
        plot(contValues,M(:,i),'.black','MarkerSize',markerSize)
        plot(fitX,regFitY(:,i),'blue','LineWidth',lineWidth)
        plot(fitX,QuadY(:,i),'green','LineWidth',lineWidth)
        if doEnhancePolyfit == true
            plot(fitX,EpolyY(:,i),'cyan','LineWidth',lineWidth)
        end
        for j = 1 : numData
            plot(STT1X(:,j),STT1Y(:,j,i),'red','LineWidth',lineWidth)
        end
        legend('Prior Solution','Poly. Fit','Quad. Fit','E. Poly. Fit','1st Order Diff.','Location','SouthEast')
        xlabel('Nondimensional Step Size')
        ylabel('Longitude Costate Value')

        figure(i+totalBCcount)
        hold on
        grid on
        set(gca,'FontSize',fontSize)
        plot(fitX((fitPoints+1):end)-fitX(fitPoints+1),regFit_error(:,i),'blue','LineWidth',lineWidth)
        plot(fitX((fitPoints+1):end)-fitX(fitPoints+1),Quad_error(:,i),'green','LineWidth',lineWidth)
        if doEnhancePolyfit == true
            plot(fitX((fitPoints+1):end)-fitX(fitPoints+1),Epoly_error(:,i),'cyan','LineWidth',lineWidth)
        end
        %plot(0.349e-3*ones(1,100),linspace(0,0.16,100),'black','LineWidth',lineWidth)
        legend('Poly. Fit','Quad. Fit','E. Poly. Fit','Location','NorthWest')
        xlabel('Nondimensional Step Size')
        ylabel('Percentage Error from 1st Order Differential')

    end
end
end

return