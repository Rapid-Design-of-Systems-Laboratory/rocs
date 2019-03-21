function analyzeOptimization(inputfunc)
 
%clear all; % clear all causes a segmentation fault with Mathematica's mathlink.
close all; clc; 

%%%%%%%%%%%% 
%% Inputs %% 
%%%%%%%%%%%% 
 
[in] = inputfunc();
 
numData = length(in.outputFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Through Results Files %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for ctrData = 1 : 1 : numData 
     
    %%%%%%%%%%%%%%%%% 
    %% Assign Data %% 
    %%%%%%%%%%%%%%%%% 
     
    % Load data file 
    data = load(in.outputFile{ctrData}); 
     
    % Clear variables 
    clear scale constraintVal const traj VEH CONT ID setCONT 
    
		% Determine number of states and controls - change to handle singular
		% arcs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		numStates = data.out.oc.num.origStates;
		numControl = data.out.oc.num.controls;
    numOrigControl = length(data.in.oc.control(:,1));
    
    % Load constraints 
		if isfield(data.in,'constraintVal')
			constraintNames = fieldnames(data.in.constraintVal);
		else
			constraintNames = {};
		end
		
    for ctrConstraint = 1 : 1 : length(constraintNames)
      eval(['constraintVal.',constraintNames{ctrConstraint},'= data.in.constraintVal.(constraintNames{ctrConstraint}){1};']); 
    end 
     
    % Load constants 
		if isfield(data.in,'const')
			constNames = fieldnames(data.in.const);
		else
			constNames = {};
		end
		
    for ctrConst = 1 : 1 : length(constNames) 
      eval(['const.',constNames{ctrConst},'= data.in.const.(constNames{ctrConst}){1};']); 
    end
    
    % % Load trajectory parameters 
    % trajNames = fieldnames(data.in.traj); 
    % for ctrTraj = 1 : 1 : length(trajNames) 
    %   eval(['traj.',trajNames{trajNames},'= data.in.traj.(trajNames{trajNames});']); 
    % end
		
		% Load optimal control data
		oc = data.in.oc;
     
%    % Load vehicle parameters 
%    vehNames = fieldnames(data.out.VEH); 
%    for ctr = 1 : 1 : length(vehNames) 
%        eval(['VEH.',vehNames{ctr},'= data.out.VEH.(vehNames{ctr});']); 
%    end 
     
    % Load input continuation data 
    CONT = data.in.CONT; 
    
		% Load initial guess 
    if data.in.convertParametersToStates
        IG.(data.in.oc.independentVariable{1}) = data.out.IG.sol.x*data.out.IG.sol.y(2*numStates+1);
    else
        IG.(data.in.oc.independentVariable{1}) = data.out.IG.sol.x*data.out.IG.sol.parameters(1); 
    end
    
    for ctr = 1 : 1 : numControl
      if ctr <= numOrigControl
        IG.(data.in.oc.control{ctr,1}) = data.out.IG.sol.control(ctr,:);
      else
		  %data.out.IG.sol.control(ctr,:)
		  % BUG HERE
        IG.(['uEpsilon',int2str(ctr-numOrigControl)]) = data.out.IG.sol.control(1,:);
      end
    end
    IG.d2Hdu2 = data.out.IG.sol.d2Hdu2;
    
    for ctr = 1 : 1 : numStates 
        eval(['IG.',char(data.in.oc.state(ctr)),' = data.out.IG.sol.y(ctr,:);']); 
        eval(['IG.lam',upper(char(data.in.oc.state(ctr))),' = data.out.IG.sol.y(ctr+numStates,:);']); 
    end
 		
    % Load output continuation data 
    if in.plot.continuation.flag
			
      for ctrSet = 1 : 1 : length(data.out.setCONT) 
        for ctrCont = 1 : 1 : length(data.out.setCONT(ctrSet).CONT) 
        	
					% Convert to normal time 
          if data.in.convertParametersToStates 
            setCONT(ctrSet).CONT(ctrCont).t = convert2ActualTime( ... 
                data.out.setCONT(ctrSet).CONT(ctrCont).sol.x, ... 
                data.out.setCONT(ctrSet).CONT(ctrCont).sol.y(2*numStates+1:end,1));
          else
            if ~data.out.setCONT(ctrSet).CONT(ctrCont).sol.errorFlag
              setCONT(ctrSet).CONT(ctrCont).t = convert2ActualTime( ... 
                  data.out.setCONT(ctrSet).CONT(ctrCont).sol.x, ... 
                  data.out.setCONT(ctrSet).CONT(ctrCont).sol.parameters);
            else
              setCONT(ctrSet).CONT(ctrCont).t = NaN;
            end
          end
					% Assign control
					for ctr = 1 : 1 : numControl
            
            if ctr <= numOrigControl
              setCONT(ctrSet).CONT(ctrCont).(data.in.oc.control{ctr,1}) = data.out.setCONT(ctrSet).CONT(ctrCont).sol.control(ctr,:);
            else
              setCONT(ctrSet).CONT(ctrCont).(['uEpsilon',int2str(ctr-numOrigControl)]) = data.out.setCONT(ctrSet).CONT(ctrCont).sol.control(ctr,:);
            end
          	
					end
					
					% Assign sufficiency conditions
					setCONT(ctrSet).CONT(ctrCont).d2Hdu2 = data.out.setCONT(ctrSet).CONT(ctrCont).sol.d2Hdu2; 
                 
          % Assign states and costates 
          for ctr = 1 : 1 : numStates 
              
					  % Time history of states and costates
            eval(['setCONT(ctrSet).CONT(ctrCont).',char(data.in.oc.state(ctr)),' = data.out.setCONT(ctrSet).CONT(ctrCont).sol.y(ctr,:);']); 
            eval(['setCONT(ctrSet).CONT(ctrCont).lam',upper(char(data.in.oc.state(ctr))),' = data.out.setCONT(ctrSet).CONT(ctrCont).sol.y(ctr+numStates,:);']);
							
						% Initial and final states 
            setCONT(ctrSet).CONT(ctrCont).x0(ctr) = data.out.setCONT(ctrSet).CONT(ctrCont).sol.y(ctr,1); 
            setCONT(ctrSet).CONT(ctrCont).xF(ctr) = data.out.setCONT(ctrSet).CONT(ctrCont).sol.y(ctr,end); 
               
          end
          
        end 
      end
			         
      % Determine starting index for plotting continuation process 
      if in.plot.contStartCtr == 0 
          startCtr = length(setCONT(in.cont.Index(ctrData)).CONT); 
      elseif in.plot.contStartCtr < 1 
          startCtr = floor(in.plot.contStartCtr*length(setCONT(in.cont.Index(ctrData)).CONT)); 
      else 
          startCtr = in.plot.contStartCtr; 
      end 
       
      % Determine continuation data not plotted 
      if in.plot.contSkip == 0 
          contSkip = length(setCONT(in.cont.Index(ctrData)).CONT) - startCtr; 
      else 
          contSkip = in.plot.contSkip; 
      end
         
    end
    
		% Determine input variable set
		inputSet = {};
		index = 0;
		for ctrState = 1 : 1 : numStates
			index = index + 1;
			inputSet{index} = char(data.in.oc.state(ctrState));
		end
		for ctrCostate = 1 : 1 : numStates
			index = index + 1;
			inputSet{index} = ['lam',upper(char(data.in.oc.state(ctrCostate)))];
		end
		for ctrControl = 1 : 1 : numControl
			index = index + 1;
      if ctrControl <= numOrigControl
        inputSet{index} = char(data.in.oc.control{ctrControl,1});
      else
        inputSet{index} = ['uEpsilon',int2str(ctrControl-numOrigControl)];
      end
		end
    for ctrConstraint = 1 : 1 : length(constraintNames) 
			index = index + 1;
			inputSet{index} = constraintNames{ctrConstraint};
    end 
    for ctrConstant = 1 : 1 : length(constNames)
			index = index + 1;
			inputSet{index} = constNames{ctrConstant};
    end
		index = index + 1;
		inputSet{index} = oc.independentVariable{1};
		
		%%%%%%%%%%%
		%% Plots %%
		%%%%%%%%%%%
		
		indexFig = 0;
		
		hLegend = [];
		for ctrPlot = 1 : 1 : length(in.figure)
			
			% Determine number of columns in plot
			try length(in.figure(ctrPlot).plot(2)); % if works then know subplot
				numRow = 2;
				numCol = ceil(length(in.figure(ctrPlot).plot)/2);
			catch
				numRow = 1;
				numCol = 1;
			end
			
			% Create figure
			indexFig = indexFig + 1;
			plotHandle.figure(ctrPlot).fig = figure(indexFig);
			
			% Movie actions
			if in.figure(ctrPlot).movie.make
				set(gcf,'Position',in.figure(ctrPlot).movie.figPos); 
				set(gcf,'color',[1 1 1]);
		    aviobj = avifile(in.figure(ctrPlot).movie.name,'Compression',in.movie.compression,'FPS',in.movie.fps,'quality',in.movie.quality); % Open movie file 
			end
			
			% Legend information
			legendStr = 'legend(hLegend(indPlot).handles';
			
			% Loop through all subplots
			indPlot = 0;
			xFunc = {};
			yFunc = {};
			for ctrRow = 1 : 1 : numRow
				for ctrCol = 1 : 1 : numCol
          
					% Set settings for subplot
					indPlot = indPlot + 1;
					subplot(numRow,numCol,indPlot); grid on; hold on;
					hLegend(indPlot).handles = [];
          
          % Check if odd number of plots in subplot. If so, skip last plot.
          try
            in.figure(ctrPlot).plot(indPlot).x{1};
          catch
            continue;
          end
          
          % Determine number of plotting expressions
          numPlotExpressions = length(in.figure(ctrPlot).plot(indPlot).x);
					
					% Convert plotted expression into function for easy evaluation
          for ctrPlotExpression = 1 : 1 : numPlotExpressions
            xFunc{ctrRow,ctrCol,ctrPlotExpression} = matlabFunction(sym(in.figure(ctrPlot).plot(indPlot).x{ctrPlotExpression}{1}),'vars',inputSet);
            yFunc{ctrRow,ctrCol,ctrPlotExpression} = matlabFunction(sym(in.figure(ctrPlot).plot(indPlot).y{ctrPlotExpression}{1}),'vars',inputSet);
          end
				
					if in.plot.presentation
						presentation_subplot;
					end
				
					% Plot initial guess
					if in.plot.initialGuess.flag
						% Perform calculations for plotting
						numPoints = length(IG.(char(data.in.oc.state(1))));
						plotX = NaN(numPoints,numPlotExpressions);
						plotY = NaN(numPoints,numPlotExpressions);
						inputVals = cell(size(inputSet));
						for ctrVal = 1 : 1 : numPoints
				
							% Construct input cell to matlabFunction
							index = 0;
							for ctrState = 1 : 1 : numStates
								index = index + 1;
								inputVals{index} = IG.(char(data.in.oc.state(ctrState)))(ctrVal);
							end
							for ctrCostate = 1 : 1 : numStates
								index = index + 1;
								inputVals{index} = IG.(['lam',upper(char(data.in.oc.state(ctrCostate)))])(ctrVal);
							end
							for ctrControl = 1 : 1 : numControl
								index = index + 1;
                if ctrControl <= numOrigControl
                  inputVals{index} = IG.(char(data.in.oc.control{ctrControl,1}))(ctrVal);
                else
                  inputVals{index} = IG.(['uEpsilon',int2str(ctrControl-numOrigControl)])(ctrVal);
                end
							end
					    for ctrConstraint = 1 : 1 : length(constraintNames) 
								index = index + 1;
								inputVals{index} = constraintVal.(constraintNames{ctrConstraint});
					    end 
					    for ctrConstant = 1 : 1 : length(constNames)
								index = index + 1;
								inputVals{index} = const.(constNames{ctrConstant});
					    end
							index = index + 1;
							inputVals{index} = IG.(data.in.oc.independentVariable{1})(ctrVal);
						
              for ctrPlotExpression = 1 : 1 : numPlotExpressions
                plotX(ctrVal,ctrPlotExpression) = xFunc{ctrRow,ctrCol,ctrPlotExpression}(inputVals{:});
                plotY(ctrVal,ctrPlotExpression) = yFunc{ctrRow,ctrCol,ctrPlotExpression}(inputVals{:});
              end
              
            end
            
						% Plot data
            for ctrPlotExpression = 1 : 1 : numPlotExpressions
              plotHandle.plot(ctrPlot).initialGuess = ...
              plot(plotX(:,ctrPlotExpression),plotY(:,ctrPlotExpression),in.plot.initialGuess.style,'LineWidth',in.plot.linewidth);
            end
							
						% Legend information
						if indPlot == 1
							legendStr = [legendStr,',''Initial Guess'''];
						end
						hLegend(indPlot).handles = [hLegend(indPlot).handles,plotHandle.plot(ctrPlot).initialGuess];
						
					end
					
				end
			end
				 
			% Plot continuation data
			if data.in.cont.method(in.cont.Index(ctrData)) == 1 || data.in.cont.method(in.cont.Index(ctrData)) == 2 && in.plot.continuation.flag
		
				% Legend data
				legendStr = [legendStr,',''Indirect Solution'''];
		
        % Loop through continuation data and plot 
        for ctrCont = startCtr : contSkip : length(setCONT(in.cont.Index(ctrData)).CONT)
					
					indPlot = 0;
					for ctrRow = 1 : 1 : numRow
						for ctrCol = 1 : 1 : numCol
							
							% Set settings for subplot
							indPlot = indPlot + 1;
							subplot(numRow,numCol,indPlot); grid on; hold on;
              
              % Check if odd number of plots in subplot. If so, skip last plot.
              try
                in.figure(ctrPlot).plot(indPlot).x{1};
              catch
                continue;
              end
							
							% Perform calculations for plotting
							numPoints = length(setCONT(in.cont.Index(ctrData)).CONT(ctrCont).(char(data.in.oc.state(1))));
							plotX = NaN(numPoints,numPlotExpressions);
							plotY = NaN(numPoints,numPlotExpressions);
							inputVals = cell(size(inputSet));
							
							
							
							% HARDCODED!!! REPLACE WITH GENERIC LOGIC IN FUTURE
% %                             I = find(9*sin(setCONT(in.cont.Index(ctrData)).CONT(ctrCont).ang)+5 < 0);
%              I = find(setCONT(in.cont.Index(ctrData)).CONT(ctrCont).alfa < 0);
% %                             setCONT(in.cont.Index(ctrData)).CONT(ctrCont).ang(I) = -setCONT(in.cont.Index(ctrData)).CONT(ctrCont).ang(I);
%              setCONT(in.cont.Index(ctrData)).CONT(ctrCont).alfa(I) = -setCONT(in.cont.Index(ctrData)).CONT(ctrCont).alfa(I);
%
%                             for ctrI = 1 : 1 : length(I)
%                                 if setCONT(in.cont.Index(ctrData)).CONT(ctrCont).bank(I(ctrI)) < 0
%                                     setCONT(in.cont.Index(ctrData)).CONT(ctrCont).bank(I(ctrI)) = setCONT(in.cont.Index(ctrData)).CONT(ctrCont).bank(I(ctrI))+pi;
%                                 else
%                                     setCONT(in.cont.Index(ctrData)).CONT(ctrCont).bank(I(ctrI)) = setCONT(in.cont.Index(ctrData)).CONT(ctrCont).bank(I(ctrI))-pi;
%                                 end
%
%              end

								% % Perform custom calculations prior to plotting
% 								for ctrCustomCalc = 1 : 1 : length(in.customCalc)
%
% 									eval(['I = find(setCONT(in.cont.Index(ctrData)).CONT(ctrCont).',in.customCalc{ctrCustomCalc}.condition,');']);
%
% 									for ctrAction = 1 : 1 : length(in.customCalc{ctrCustomCalc}.action)
% 										eval(['setCONT(in.cont.Index(ctrData)).CONT(ctrCont).',in.customCalc{ctrCustomCalc}.action,';']);
% 									end
%
% 								end
							
							
							
							for ctrVal = 1 : 1 : numPoints
								
								% Construct input cell to matlabFunction
								index = 0;
								for ctrState = 1 : 1 : numStates
									index = index + 1;
									inputVals{index} = setCONT(in.cont.Index(ctrData)).CONT(ctrCont).(char(data.in.oc.state(ctrState)))(ctrVal);
								end
								for ctrCostate = 1 : 1 : numStates
									index = index + 1;
									inputVals{index} = setCONT(in.cont.Index(ctrData)).CONT(ctrCont).(['lam',upper(char(data.in.oc.state(ctrCostate)))])(ctrVal);
								end
								for ctrControl = 1 : 1 : numControl
									index = index + 1;
                  if ctrControl <= numOrigControl
                    inputVals{index} = setCONT(in.cont.Index(ctrData)).CONT(ctrCont).(char(data.in.oc.control{ctrControl,1}))(ctrVal);
                  else
                    inputVals{index} = setCONT(in.cont.Index(ctrData)).CONT(ctrCont).(['uEpsilon',int2str(ctrControl-numOrigControl)])(ctrVal);
                  end
								end
						    for ctrConstraint = 1 : 1 : length(constraintNames) 
									index = index + 1;
									inputVals{index} = constraintVal.(constraintNames{ctrConstraint});
						    end 
						    for ctrConstant = 1 : 1 : length(constNames)
									index = index + 1;
									inputVals{index} = const.(constNames{ctrConstant});
						    end
								index = index + 1;
								inputVals{index} = setCONT(in.cont.Index(ctrData)).CONT(ctrCont).(data.in.oc.independentVariable{1})(ctrVal);
								
                for ctrPlotExpression = 1 : 1 : numPlotExpressions
                  plotX(ctrVal,ctrPlotExpression) = xFunc{ctrRow,ctrCol,ctrPlotExpression}(inputVals{:});
                  plotY(ctrVal,ctrPlotExpression) = yFunc{ctrRow,ctrCol,ctrPlotExpression}(inputVals{:});
                end
						
							end
							
							% Plot data
							if in.plot.colorSegments.flag % color each segment
                Isegments = find(diff(data.out.setCONT(in.cont.Index(ctrData)).CONT(ctrCont).sol.x) == 0);
                IsegmentsBegin = [1 Isegments+1]; % +1 due to repeated point at segment juncture
                IsegmentsEnd = [Isegments numPoints];
                for ctrPlotExpression = 1 : 1 : numPlotExpressions
                  for ctrSegments = 1 : 1 : length(Isegments)+1
                    plotHandle.plot(ctrPlot).continuation(ctrData) = ...
                      plot(plotX(IsegmentsBegin(ctrSegments):IsegmentsEnd(ctrSegments),ctrPlotExpression), ...
                      plotY(IsegmentsBegin(ctrSegments):IsegmentsEnd(ctrSegments),ctrPlotExpression), ...
                      in.plot.colorSegments.colorSet{ctrSegments},'LineWidth',in.plot.linewidth);
                  end
                end
							else % do not color each segment
                if ~isfield(in.figure(ctrPlot).plot(indPlot),'style')
                  for ctrPlotExpression = 1 : 1 : numPlotExpressions
                    plotHandle.plot(ctrPlot).continuation(ctrData) = ...
                    plot(plotX(:,ctrPlotExpression),plotY(:,ctrPlotExpression),in.plot.continuation.style{ctrData},'LineWidth',in.plot.linewidth);
                  end
                else
                  for ctrPlotExpression = 1 : 1 : numPlotExpressions
                    plotHandle.plot(ctrPlot).continuation(ctrData) = ...
                    plot(plotX(:,ctrPlotExpression),plotY(:,ctrPlotExpression),in.figure(ctrPlot).plot(indPlot).style{ctrPlotExpression},'LineWidth',in.plot.linewidth);
                  end
                end
							end
						  
							if in.figure(ctrPlot).movie.make
								plotHandleToDelete(ctrRow,ctrCol) = plotHandle.plot(ctrPlot).continuation(ctrData);
							end
							
							% Plotting items
							if ctrCont == startCtr
								
								% Legend information
								hLegend(indPlot).handles = [hLegend(indPlot).handles,plotHandle.plot(ctrPlot).continuation(ctrData)];
								if indPlot == 1
									legendStr = [legendStr,',''Location'',''',in.figure(ctrPlot).plot(indPlot).legend.location,''');'];
								end
								
								if in.plot.legend
									eval(legendStr);
								end
								
								% Labels
								xlabel([in.figure(ctrPlot).plot(indPlot).x{1}{3},', ',in.figure(ctrPlot).plot(indPlot).x{1}{2}]);
								ylabel([in.figure(ctrPlot).plot(indPlot).y{1}{3},', ',in.figure(ctrPlot).plot(indPlot).y{1}{2}]);
								
								% Title
								if in.plot.title
									title([in.figure(ctrPlot).plot(indPlot).y{1}{3},' vs. ',in.figure(ctrPlot).plot(indPlot).x{1}{3}]);
								end
								
								% Font size
								if in.plot.presentation
									presentation_subplot;
								end
								
							end
							
							if isfield(in.figure(ctrPlot).plot(indPlot),'axis')
								axis(in.figure(ctrPlot).plot(indPlot).axis);
							end
							
						end
					end
					
					% Movie actions
					if in.figure(ctrPlot).movie.make
						
						% Add frame to movie
						aviobj = addframe(aviobj,getframe(gcf));
						
						% Delete continuation data
            delete(plotHandleToDelete);
						
					end
						
				end
				
			end
			
	    % Close movie file
			if in.figure(ctrPlot).movie.make
	    	aviobj = close(aviobj);
			end
		
		end
		
end
 
%%%%%%%%%%%%%%%%%%%%%% 
%% Plot Cost Trades %% 
%%%%%%%%%%%%%%%%%%%%%% 
 
if in.plot.continuation.flag && isfield(in,'trade') 
     
    % Initialize 
    cost = nan(length(setCONT(in.cont.Index(ctrData)).CONT),1); 
    tradeVal = nan(length(setCONT(in.cont.Index(ctrData)).CONT),1); 
     
    for ctr = startCtr : contSkip : length(setCONT(in.cont.Index(ctrData)).CONT) 
         
        % Assign state variables 
        for ctrStates = 1 : 1 : numStates 
            eval([char(data.in.oc.state(ctrStates)),' = setCONT(in.cont.Index(ctrData)).CONT(ctr).',char(data.in.oc.state(ctrStates)),';']); 
        end 
         
        % Convert Lagrangian to elementwise multiplication and division 
        lagrangeVectorForm = elementwiseConversion(char(data.in.oc.Lagrange)); 
         
        % Load constants for subsequent evaluations 
        constNames = fieldnames(data.in.const); 
        for ctrConst = 1 : 1 : length(constNames) 
            eval([constNames{ctrConst},'= setCONT(in.cont.Index(ctrData)).CONT(ctr).in.const.(constNames{ctrConst});']); 
        end 
         
        % Load vehicle parameters for subsequent evaluations 
        vehNames = fieldnames(data.out.VEH); 
        for ctrVEH = 1 : 1 : length(vehNames) 
            eval([vehNames{ctrVEH},'= setCONT(in.cont.Index(ctrData)).CONT(ctr).out.VEH.(vehNames{ctrVEH});']); 
        end 
         
        % Load initial and final states for subsequent evaluations 
        if isfield(setCONT(in.cont.Index(ctrData)).CONT(ctr),'x0') 
            x0 = setCONT(in.cont.Index(ctrData)).CONT(ctr).x0; 
        end 
        if isfield(setCONT(in.cont.Index(ctrData)).CONT(ctr),'xF') 
            xF = setCONT(in.cont.Index(ctrData)).CONT(ctr).xF; 
        end
				
				% Load control
				eval([char(data.in.oc.control),' = setCONT(in.cont.Index(ctrData)).CONT(ctrCont).(data.in.oc.control);']);
         
        % Compute Lagrangian cost 
        eval(['cost(ctr) = trapz(setCONT(in.cont.Index(ctrData)).CONT(ctr).time,(',lagrangeVectorForm,')*in.trade.scaleCost);']); 
         
        % Convert trade value to elementwise multiplication and division 
        tradeVectorForm = elementwiseConversion(char(in.trade.func{ctrData})); 
         
        % Compute trade value
        eval(['tradeVal(ctr) = ',tradeVectorForm,';']); 
         
    end
    
    % Sort trade value and cost so in ascending order for plotting 
    [tradeValSorted,I] = sort(tradeVal);
    costSorted = cost(I);
     
    % Plot trade
		if ~exist(figHandleCost)
    	figHandleCost = figure; grid on; hold on;
		else
			figure(figHandleCost);
		end
    plotHandle.trade.cost(ctrData) = plot(tradeValSorted,costSorted,in.trade.plotStyle{ctrData},'LineWidth',in.plot.linewidth);
		
		% Plotting items
    xlabel(in.trade.xLabel); 
    ylabel(in.trade.yLabel); 
    title(in.trade.title);
		if in.plot.presentation
    	presentation_plot;
		end
		
    if numData > 1 
        legend_str = 'legend(plotHandle.trade.cost'; 
        legend_str = [legend_str,',',in.trade.legend.str]; 
        legend_str = [legend_str,',''Location'',',in.trade.legend.location,');']; 
        eval(legend_str); 
    end
    
    % Compute cost modifier function (for extra analysis) 
    eval(['costModifier = ',in.trade.costModifier.func,';']); 
     
    % Sort modified cost value 
    costModifierSorted = costModifier(I); 
     
    % Plot modified cost function
		if ~exist(figHandleCostModified)
    	figure; grid on; hold on;
		else
			figure(figHandleCostModified);
		end
    plotHandle.trade.modifiedCost(ctrData) = plot(tradeValSorted,costModifierSorted,plotStyle{ctrData},'LineWidth',in.plot.linewidth); 
    
		% Plotting items
    xlabel(in.trade.xLabel); 
    ylabel(in.trade.costModifier.yLabel); 
    title(in.trade.title);
		if in.plot.presentation
    	presentation_plot;
		end
		
    if numData > 1
      legend_str = 'legend(plotHandle.trade.modifiedCost';
      legend_str = [legend_str,',',in.trade.legend.str];
      legend_str = [legend_str,',''Location'',',in.trade.legend.location,');'];
			if in.plot.legend
				eval(legend_str);
			end
    end
		
end
 
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
function [t] = convert2ActualTime(tNormalized,parameters) 
 
% Converts normalized time from BVP4C to actual time 
 
% Initialize actual time vector 
t = nan(size(tNormalized)); 
 
% Remove costate jumps from end of parameters vector and augment parameters 
% vector with zero for first arc starting time 
maxNormalizedTime = max(tNormalized); 
timeSet = abs(parameters(1:maxNormalizedTime)); 

% Determine number of arcs for conversion 
numArcs = length(timeSet); 

% Loop through arcs and perform conversion
timeSave = 0;
for ctrArc = 1 : 1 : numArcs
     
    I = find((tNormalized > ctrArc-1) & (tNormalized < ctrArc));
    
    % Protect for no intermediate points for arc (occurs for simple problems
    % like Brachistochrone)
    if isempty(I)
      I1 = find(tNormalized == ctrArc-1,1,'last'); % first point
      I2 = find(tNormalized == ctrArc,1,'first'); % last point
      I = [I1:1:I2];
    else
      I = [I(1)-1 I I(end)+1]; %#ok<*AGROW> % Include one of repeated points 
    end
    
    t(I) = (tNormalized(I)-(ctrArc-1))*timeSet(ctrArc) + timeSave;
		timeSave = timeSave + timeSet(ctrArc);
     
end 

% for ctrArc = 0 : 1 : numParameters-2 
%      
%     I = find((tNormalized > ctr) & (tNormalized < ctr+1)); 
%     I = [I(1)-1 I I(end)+1]; %#ok<*AGROW> % Include one of repeated points 
%     t(I) = (tNormalized(I)-ctr)*(parametersAug(ctr+2)) + parametersAug(ctr+1);
% 		timeSave =  
%      
% end 
 
return 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
function [outStr] = elementwiseConversion(inStr) 
 
% Initialize output string 
outStr = inStr; 
 
% Convert "*" to ".*" 
k_star = strfind(outStr,'*'); 
for ctr = 1 : 1 : length(k_star) 
     
    % Find index of "*" 
    k_star = strfind(outStr,'*'); 
    ind = k_star(ctr); 
     
    % Split string for inserting since strings not same length 
    outStr = [outStr(1:ind-1),'.*',outStr(ind+1:end)]; 
     
end 
 
% Convert "/" to "./" 
k_slash = strfind(outStr,'/'); 
for ctr = 1 : 1 : length(k_slash) 
     
    % Find index of "/" 
    k_slash = strfind(outStr,'/'); 
    ind = k_slash(ctr); 
     
    % Split string for inserting since strings not same length 
    outStr = [outStr(1:ind-1),'./',outStr(ind+1:end)]; 
     
end 
 
% Convert "^" to ".^" 
k_power = strfind(outStr,'^'); 
for ctr = 1 : 1 : length(k_power) 
     
    % Find index of "^" 
    k_power = strfind(outStr,'^'); 
    ind = k_power(ctr); 
     
    % Split string for inserting since strings not same length 
    outStr = [outStr(1:ind-1),'.^',outStr(ind+1:end)]; 
     
end 
 
return 

