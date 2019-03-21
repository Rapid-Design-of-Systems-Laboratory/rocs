function varargout = mikegui(varargin)
% MIKEGUI Mike's developing GUI for optimization.

% NASA JSC - DM42 / Michael J. Grant on July 2008

% Notes: 'HandleVisibility','callback' does not seem to allow children to be
% deleted by clf. Removed from all objects except menu items since they are 
% deleted (maybe try using multiple tabs).

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

% clear all; pack; close all; clc;
pos = [0.01 0.05 0.95 0.85];
pos_axes = [0.20 0.10 0.70 0.70];
vel_pos = [0.02 0.66 0.10 0.15];
plot_set = {'b.','r.','g.','m.','k.','b*','r*','g*','m*','k*', ...
  'b+','r+','g+','m+','k+','bo','ro','go','mo','ko','bx','rx','gx','mx', ...
  'kx','bs','rs','gs','ms','ks'};
symbol_set = {'ro','gx','k+','rs','rd','rv','rp','rh'};

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non-UI Declarations %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
try
  data = load('results.mat'); % Optimization data
  in = data.in;
  od = data.od;
  clear data;
catch
  error('results.mat not found');
end

%%%%%%%%%%%%%%%%%%%%
%% Error Checking %%
%%%%%%%%%%%%%%%%%%%%

% Send error if not enough plot sets available
if length(plot_set) < od.dim
  error(['Not enough plotting types. Increase number of plotting ', ...
    'types in Inputs section to match number of trade space dimensions.']);
  return;
end

% Send error if not enough symbols available for constraints
if isfield(in,'constr') && length(symbol_set) < length(in.constr(:,1))
  error(['Not enough symbol types for constraints. Increase number of ', ...
    'symbols in Inputs section to match number of constraints.']);
end

%%%%%%%%%%%%%%%%%%%
%% Preprocessing %%
%%%%%%%%%%%%%%%%%%%

% Determine Pareto front (not just archive at end of optimization)
pareto_front = build_pareto_front(od.fit(:,:,1:od.iter),[],in.maximization);

%%%%%%%%%%%%%%%%%%%%%
%% UI Declarations %%
%%%%%%%%%%%%%%%%%%%%%

% Main GUI figure
hMainFigure = figure('Name', mfilename, ... % 'MenuBar','none', 'Toolbar','none', 
  'NumberTitle','off', 'Color', get(0, 'defaultuicontrolbackgroundcolor'), ...
  'Units','normalized', 'Position',pos, 'HandleVisibility','callback');

% Navigation menu
hNavMenu = uimenu('Parent',hMainFigure, 'Label','Menu', ...
  'HandleVisibility','callback');

% PSO or MOPSO
if strcmp(in.method,'pso') || strcmp(in.method,'mopso')
  % Velocity menu
  hVelMenuitem = uimenu('Parent',hNavMenu, 'Label','Velocity', ...
    'Callback', @hVelMenuitemCallback, 'HandleVisibility','callback');
  % Inertia menu
  hInertiaMenuitem = uimenu('Parent',hNavMenu, 'Label','Inertia', ...
    'Callback', @hInertiaMenuitemCallback, 'HandleVisibility','callback');
  % Trade space / objective space
  hTradeMenuitem = uimenu('Parent',hNavMenu, 'Label','Trade Space / Objective Space', ...
    'Callback', @hTradeMenuitemCallback, 'HandleVisibility','callback');
end

% MOPSO only
if strcmp(in.method,'mopso')
  % Pareto front menu
  hParetoMenuitem = uimenu('Parent',hNavMenu, 'Label','Pareto Front', ...
    'Callback', @hParetoMenuitemCallback, 'HandleVisibility','callback');
end

% Summary menu
hSummaryMenuitem = uimenu('Parent',hNavMenu, 'Label','Summary', ...
  'Callback', @hSummaryMenuitemCallback, 'HandleVisibility','callback');

% Copy menu
hCopyMenuitem = uimenu('Parent',hNavMenu, 'Label','Copy Figure', ...
  'HandleVisibility','callback');

hCopyKeepitem = uimenu('Parent',hCopyMenuitem, 'Label','New Window', ...
  'Callback', @hCopyMenuitemCallback, 'HandleVisibility','callback');

hCopyRemoveitem = uimenu('Parent',hCopyMenuitem, 'Label','Copy Figure', ...
  'Callback', @hCopyMenuitemCallback, 'HandleVisibility','callback');

% Movie menu
hMovieMenuitem = uimenu('Parent',hNavMenu, 'Label','Create Movie', ...
  'HandleVisibility','callback', 'Callback', @hMovieMenuitemCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Callback and Utility Functions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %----------------------------------------------------------------------
  % Create movie menu
  function hMovieMenuitemCallback(hObject, eventdata)
    
    num_fig = 1;
    
    % Initialize plot information
    pos1 = [0.05 0.50 0.46 0.45];
    pos2 = [0.08 0.06 0.33 0.33];
    pos3 = [0.58 0.50 0.40 0.45];
    pos4 = [0.50 0.08 0.33 0.33];
    pos5 = [0.83 0.08 0.15 0.30];
    
    plot_list = [{'N/A'},in.ts_label,in.obj_label,{'Iteration','Inertia','Cd'}];
    
    clf(hMainFigure);
    % Turn on figure toolbar
    set(gcf,'toolbar','figure');
    
%     % Axis listbox label
%     hPlotText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
%       'Position',[0.03 0.73+0.1 0.18 0.02], ...
%       'String','Axis 1             Axis 2             Axis 3', 'Style','text', ...
%       'HorizontalAlignment','Left');

    % Add listboxes for plotting
    % Plot 1
    hPlot1a = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.02 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    hPlot1b = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.09 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    hPlot1c = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.16 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');

    % Plot 2
    hPlot2a = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.26 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    hPlot2b = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.33 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    hPlot2c = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.40 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    
    % Plot 3
    hPlot3a = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.50 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    hPlot3b = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.57 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    hPlot3c = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.64 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    
    % Plot 4
    hPlot4a = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.74 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    hPlot4b = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.81 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    hPlot4c = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.89 0.58+0.1 0.06 0.15], 'String',plot_list, 'Style','listbox');
    
    % Movie text box
    hPlotText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.16 0.40+0.1 0.06 0.15], 'String',{'Iter'}, 'Style','listbox','Max',2);
    
    % Make movie button
    hMakeMovieButton = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) 0.15 vel_pos(3) 0.05], ...
      'String','Create Movie', 'Callback', @hMakeMovieButtonCallback);
    
%     % Plot type static text label
%     hPlotText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
%       'Position',[vel_pos(1) vel_pos(2)+vel_pos(4)+0.001+0.1 vel_pos(3) 0.02], ...
%       'String','Plot Type:', 'Style','text', 'HorizontalAlignment','Left');
%     
%     % Plot type popup menu
%     hPlotOptionMenu = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
%       'Position',[vel_pos(1) vel_pos(2)+vel_pos(4)-0.03+0.1 vel_pos(3)+0.03 0.03], ...
%       'String',{'Scatter Plot','Parallel Coordinates'}, 'HorizontalAlignment','Left', ...
%       'Style','popupmenu', 'Callback', @hPlotOptionMenuCallback);
%     
% 
%     
%     hTimeText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
%       'Position',[vel_pos(1) .8 0.8 0.02], ...
%       'String',['Execution Time: ',time_string], 'Style','text', 'HorizontalAlignment','Left');
%     
%     hIterText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
%       'Position',[vel_pos(1) 0.75 0.8 0.02], ...
%       'String',['Iterations: ',num2str(od.iter),' (max ',num2str(in.max_iter),')'], ...
%       'Style','text', 'HorizontalAlignment','Left');
    
    %----------------------------------------------------------------------
    % Create movie
    function hMakeMovieButtonCallback(hObject, eventdata)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Output - Move to Summary Later %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Determine convergence
      if od.iter == in.max_iter

        fprintf('\nSolution not necessarily converged.\n');
        fprintf('Maximum number of iterations performed (%g iterations)\n',od.iter);
  
      else
  
        fprintf('\nSolution Converged\n');
        fprintf('Number of Iterations: %g\n',od.iter);
  
      end

      % Execution time
      fprintf('Time elapsed: %g\n',od.time);

      % Determine number of objectives in optimization
      if in.num_obj == 1

        % Single-objective optimization

        % Negate fitness if performing maximization
        if in.maximization
          od.best_fit = -od.best_fit;
          od.fit = -od.fit;
        end

        [Y,I] = min(od.best_fit);
        I = min(I);
        [Y,p] = min(od.fit(1,:,I));
        p = min(p);

        if in.maximization
          od.best_fit = -od.best_fit;
          od.fit = -od.fit;
        end

        % Data
        fprintf('\nOptimal Solution Attained in Iteration %d out of %d',I,od.iter);
        fprintf('\nOptimal Solution Value: %g',od.best_fit(I));
        fprintf('\nTrade Space Location:\n');
        fprintf('      %g\n',od.pop(:,p,I));
  
      else
  
        % Multi-objective optimization
  
        % Check to make sure more than one non-dominated solution
        if length(find(isnan(od.cd) == 0)) == 1
    
          % Only one non-dominated solution. Objectives may not be in competition
          fprintf('Warning! Only one non-dominated solution found!\n');
    
        end

      end

      fprintf('\n\n');
      
      %%%%%%%%%%%%%%%%%%
      %% Movie Figure %%
      %%%%%%%%%%%%%%%%%%
      
      %% New %%
      
      % Create new figure to create movie from
      fig = figure; %('Visible','off');
      set(gcf,'color','white');
      set(gcf,'Position',[1 31 1023 663]);
      make_movie = false;
      movie_name = 'optimization_movie.avi';
      pause_val = 0.01;
      
      % Get plot information
      Axis_Set1 = [get(hPlot1a,'Value') get(hPlot1b,'Value') get(hPlot1c,'Value')];
      Axis_Label1 = plot_list(Axis_Set1);
      Axis_Set2 = [get(hPlot2a,'Value') get(hPlot2b,'Value') get(hPlot2c,'Value')];
      Axis_Label2 = plot_list(Axis_Set2);
      Axis_Set3 = [get(hPlot3a,'Value') get(hPlot3b,'Value') get(hPlot3c,'Value')];
      Axis_Label3 = plot_list(Axis_Set3);
      Axis_Set4 = [get(hPlot4a,'Value') get(hPlot4b,'Value') get(hPlot4c,'Value')];
      Axis_Label4 = plot_list(Axis_Set4);
      
      
      function [data] = get_plot_data(param,iter)
        
        % This function is used to get values from a specific variable
        % (specified by param) at the current iteration (specified by iter)
        
        switch param
          
          case 'Iteration' % Obtain all of the iteration data
            data = [1 : 1 : od.iter];
            
          case 'Inertia' % Obtain all inertia data
            data = od.w(1:od.iter);
            
          case 'Cd' % Obtain mean normalized crowdedness distance
            data = od.cd_mean(1:od.iter);
            
          case in.ts_label
            I = find(strcmp(param,in.ts_label));
            data = od.pop(I,:,counter);
            
          case in.obj_label
            I = find(strcmp(param,in.obj_label));
            data = od.fit(I,:,counter);
            
          case 'pos_arch'
            data = od.pos_arch(:,:,counter);
            
          case 'fit_arch'
            data = od.fit_arch(:,:,counter);
            
          case 'N/A'
            data = zeros(1,in.num_ind);
            
          otherwise % An invalid param was used
            error('Invalid param value used');
          
        end
        
      end
      
      
      
      % Multi-objective particle swarm

      
      % Create movie object and start recording
      if make_movie
        aviobj = avifile(movie_name,'Compression','Indeo5','FPS',8);
      end
      
      for counter = 1 : 1 : od.iter
        
        % Plot 1 - Top left
        % Originally trade space plot
        subplot('position',pos1);
        % Plot swarm
        plot3(get_plot_data(Axis_Label1{1}),get_plot_data(Axis_Label1{2}),get_plot_data(Axis_Label1{3}),'r.');
        hold on;
        % Plot front
        data = get_plot_data('pos_arch');
        plot3(data(1,:),data(2,:),data(3,:),'b.');
        title('Evolution of Pareto Front and Swarm - Trade Space');
        xlabel(Axis_Label1{1});
        ylabel(Axis_Label1{2});
        zlabel(Axis_Label1{3});
        grid on;
        hold off;
        
        % Plot 2 - Bottom left - Current assumes 2-D plot with iteration bar
        % Originally inertia vs. iteration
        subplot('position',pos2);
        plot(get_plot_data(Axis_Label2{1}),get_plot_data(Axis_Label2{2}));
        hold on;
        grid on;
        axis_info = get(gca);
        plot([counter counter],axis_info.YLim,'r');
        title([Axis_Label2{2},' vs. ',Axis_Label2{1}]);
        xlabel(Axis_Label2{1});
        ylabel(Axis_Label2{2});
        grid on;
        hold off;
        
        % Plot 3 - Top right
        % Originally fitness plot
        subplot('position',pos3);    
        % Set plot size - Change to be over all iterations
        trade_min = min(od.fit_arch(:,:,counter),[],2);
        trade_max = max(od.fit_arch(:,:,counter),[],2);
        axis_arg = [];
        for ctr = 1 : 1 : in.num_obj
          axis_arg = [axis_arg trade_min(ctr) trade_max(ctr)];
        end
        % Plot front
        data = get_plot_data('fit_arch');
        data(3,:) = zeros(1,in.num_arch);
        plot3(data(1,:),data(2,:),data(3,:),'b.');
        hold on;
        axis(axis_arg);
        grid on;
        % Plot swarm
        plot3(get_plot_data(Axis_Label3{1}),get_plot_data(Axis_Label3{2}),get_plot_data(Axis_Label3{3}),'rx');
        title('Evolution of Pareto Front and Swarm - Fitness Space');
        xlabel(Axis_Label3{1});
        ylabel(Axis_Label3{2});
        if in.num_obj == 3
          zlabel(Axis_Label3{3});
        end
        hold off;
        % Change view if only plotting 2 vectors
        if in.num_obj == 2
          view(2);
        end
        
        % Plot 4 - Bottom right - Current assums 2-D plot with iteration bar
        % Originally average crowdedness distance vs. iteration
        % Original title: 'Mean Normalized Crowdedness Distance vs. Iteration'
        % Original ylabel: 'Mean Normalized Crowdedness Distance'
        subplot('position',pos4);
        plot(get_plot_data(Axis_Label4{1}),get_plot_data(Axis_Label4{2}));
        hold on;
        axis_info2 = get(gca);
        plot([counter counter],axis_info2.YLim,'r');
        grid on;
        title([Axis_Label4{2},' vs. ',Axis_Label4{1}]);
        xlabel(Axis_Label4{1});
        ylabel(Axis_Label4{2});
        hold off;
      
        % Text box
        subplot('position',pos5);
        if counter == 1
          set(gca,'Visible','off');
        end
        cla;
        % Text information
        text(0.1,0.8,['Iteration: ',num2str(counter)]);

        % Save frame if making movie and pause if necessary
        if make_movie
          aviobj = addframe(aviobj,getframe(gcf)); % Add figure to movie
        end
        pause(pause_val);
        
      end
      
      % If creating movie, close movie object
      if make_movie
        aviobj = close(aviobj);
      end

%         % Plot parallel coordinate system
% 
%         [m,n] = size(od.fit_arch(:,:,od.iter));
%         max_vals = max(od.fit_arch(:,:,od.iter),[],2);
%         max_vals = ones(m,1);
%         arch_norm = od.fit_arch(:,:,od.iter) ./ ...
%           max_vals(:,ones(size(od.fit_arch(:,:,od.iter),2),1));
% 
%         figure; hold on;
%         set(gcf,'color','white');
%         set(gcf,'Position',[1 31 1023 663]);
% 
%         for member = 1 : 1 : n
% 
%           plot([1:1:m],arch_norm(:,member));
% 
%         end
%         set(gca,'xtick',[0:1:in.num_obj]) 
%         xlabel('Objective');
%         ylabel('J');
%         title('Parallel Coordinate System of Pareto Front');
% 
%         % Big plot of Pareto Front
%         figure;
%         set(gcf,'color','white');
%         set(gcf,'Position',[1 31 1023 663]);
%         if in.num_obj == 2
%           plot(od.fit_arch(1,:,od.iter),od.fit_arch(2,:,od.iter),'b.');
%         elseif in.num_obj == 3
%           plot3(od.fit_arch(1,:,od.iter),od.fit_arch(2,:,od.iter), ...
%             od.fit_arch(3,:,od.iter),'b.');
%           zlabel('J3');
%         end
%         title('Pareto Front');
%         xlabel('J1');
%         ylabel('J2');
%         grid on;

      
      
      %% Old %%
      
      
      % Determine number of objectives in optimization
      if in.num_obj == 1

        % Single-objective optimization - generic plots

        % Contour of positions
        x = [];
        y = [];
        fitness = [];

        if od.dim ~= 1

          % Construct fitness contour
          for iter_counter = 1 : 1 : od.iter
            fitness = [fitness od.fit(:,:,iter_counter)];
            x = [x od.pop(1,:,iter_counter)];
            y = [y od.pop(2,:,iter_counter)];
          end

          x_min = min(x);
          x_max = max(x);
          y_min = min(y);
          y_max = max(y);

          step_x = (x_max - x_min)/100;
          step_y = (y_max - y_min)/100;
          xg = [x_min : step_x : x_max];
          yg = [y_min : step_y : y_max];

          [X,Y] = meshgrid(xg,yg);
          warning off; % Turn off warning for duplicate data
          fitness_plot = griddata(x,y,fitness,X,Y);
          warning on;

        end

        % Create subplot of current solution state
        subplot('position',pos1); hold on; set(gca,'box','on'); grid on;
        title('Position of Individuals in Population - Current');
        xlabel('x');
        xlim([in.pos_min(1) in.pos_max(1)]);
        if od.dim ~= 1
          ylabel('y');
          ylim([in.pos_min(2) in.pos_max(2)]);
          contour(X,Y,fitness_plot);
        end

        if strcmp(in.method,'pso')
          % PSO - Create subplot of inertia
          subplot('position',pos2); set(gca,'box','on'); grid on; hold on;
          title('Inertia vs. Iteration');
          xlabel('Iteration'); ylabel('Inertia');
          plot([1:1:od.iter],od.w(1:od.iter));
          axis_info = get(gca);
        elseif strcmp(in.method,'ga')
          % GA - Do nothing
        end

        % Create subplot of solution state history
        subplot('position',pos3); set(gca,'box','on'); hold on; grid on;
        title('Position of Individuals in Population - History');
        xlabel('x');
        xlim([in.pos_min(1) in.pos_max(1)]);
        if od.dim ~= 1
          ylabel('y');
          ylim([in.pos_min(2) in.pos_max(2)]);
          contour(X,Y,fitness_plot);
        end

        % Create subplot of best fitness curve
        subplot('position',pos4); set(gca,'box','on'); grid on; hold on;
        title('Best Fitness vs. Iteration');
        xlabel('Iteration'); ylabel('Best Fitness');
        plot([1:1:od.iter],od.best_fit(1:od.iter));
        axis_info2 = get(gca);

        % Create "subplot" of text information
        subplot('position',pos5); set(gca,'Visible','off'); hold on;

        % Create movie object and start recording
        if make_movie
      %     aviobj = avifile('optimization_movie.avi','Compression','Indeo5','FPS',8);
          aviobj = avifile('optimization_movie','FPS',8);
        end

        % Plot data and create movie frame by frame
        for counter = 1 : 1 : od.iter
          % Plot particle swarm info
          subplot('position',pos1);
          if od.dim == 1
            h1 = plot(od.pop(1,:,counter),zeros(1,length(od.pop(1,:,counter))),'.');
          else
            h1 = plot(od.pop(1,:,counter),od.pop(2,:,counter),'.');
          end
          subplot('position',pos4);
          h4 = plot([counter counter],axis_info2.YLim,'r');
          subplot('position',pos3);
          if od.dim == 1
            h2 = plot(od.pop(1,:,counter),zeros(1,length(od.pop(1,:,counter))),'.');
          else
            h2 = plot(od.pop(1,:,counter),od.pop(2,:,counter),'.');
          end
          subplot('position',pos5);
          h5 = text(0.1,0.8,['Iteration: ',num2str(counter)]);
          if strcmp(in.method,'pso')
            subplot('position',pos2);
            h6 = plot([counter counter],axis_info.YLim,'r');
          elseif strcmp(in.method,'ga')
            h6 = [];
          end
          if make_movie
            frame = getframe(gcf);
            aviobj = addframe(aviobj,frame); % Add figure to movie
          end
          if pausing; pause(0.1); end
          if counter ~= od.iter
            delete(h1,h4,h5,h6);
          end
        end

        if make_movie
          aviobj = close(aviobj);
        end

        if strcmp(in.method,'pso')

          % Fitness Coefficient of Variance vs. Iteration
          figure;
          plot([1:od.iter],od.cov(1:od.iter));
          title('Fitness COV vs. Iteration');
          xlabel('Iteration');
          ylabel('Coefficient of Variance');
          grid on;

          % Velocity vs. Iteration - Each Dimension
          for dim_counter = 1 : 1 : od.dim
            figure; hold on; grid on;
            % Formatting
            title(['Velocity of Each Individual - Dim ',num2str(dim_counter)]);
            xlabel('Iteration');
            ylabel('Velocity');
            % Plot velocity of each individual for each iteration
            for ind_counter = 1 : 1 : in.num_ind
              data(1:od.iter) = od.vel(dim_counter,ind_counter,[1:od.iter]);
              plot([1:od.iter],data,'.');
            end
      %       % Size window to bounds
      %       axis([1 od.iter in.vel_min(dim_counter) in.vel_max(dim_counter)]);
          end

        end

      else

        


      end

      % Close figure if making movie and no pausing
      if make_movie && pause_val == 0
        close(fig);
      end
      
    end

  end

  %----------------------------------------------------------------------
  % Copy figure menu
  function hCopyMenuitemCallback(hObject, eventdata)

    % Get axis for copying (assumes only one axis on plot!)
    hPlotAxes = gca;
    hLeg = findobj(gcf,'Tag','legend');
    
    % Copy object
    hplotcopy = copyobj(hPlotAxes,hMainFigure);
    
    % Paste copied object to temporary figure
    if hObject == hCopyRemoveitem
      htempfig = figure('Visible','off');
    elseif hObject == hCopyKeepitem
      htempfig = figure;
    end
    set(hplotcopy,'Parent',htempfig);
    set(hplotcopy,'Position','default');
    Legend_Set = get(hLeg,'String');
    Legend_Pos = get(hLeg,'Location');
    Legend_String = 'legend(''';
    for ctr = 1 : 1 : length(Legend_Set)
      Legend_String = [Legend_String,Legend_Set{ctr}];
      if ctr ~= length(Legend_Set)
        Legend_String = [Legend_String,''','''];
      else
        Legend_String = [Legend_String,''''];
      end
    end
    Legend_String = [Legend_String,',''Location'',''',Legend_Pos,''');'];
    
    eval(Legend_String);
    
    % Copy the figure to the clipboard
    if hObject == hCopyRemoveitem
      eval(['print -dmeta -r0 -f',num2str(htempfig)]);
      close(htempfig); % Close temporary figure
    end

  end

  %----------------------------------------------------------------------
  % Trade space / objective space analysis
  function hTradeMenuitemCallback(hObject, eventdata)
  
    clf(hMainFigure);
    % Turn on figure toolbar
    set(gcf,'toolbar','figure');
    
    % Trade space / objective space front plot axes
    hPlotAxes = axes('Parent', hMainFigure, 'Units', 'normalized', ...
      'Position',[0.25 0.10 0.70 0.70]);
    
    % Plot type static text label
    hPlotText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) vel_pos(2)+vel_pos(4)+0.001+0.1 vel_pos(3) 0.02], ...
      'String','Plot Type:', 'Style','text', 'HorizontalAlignment','Left');
    
    % Plot type popup menu
    hPlotOptionMenu = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) vel_pos(2)+vel_pos(4)-0.03+0.1 vel_pos(3)+0.03 0.03], ...
      'String',{'Scatter Plot','Parallel Coordinates'}, 'HorizontalAlignment','Left', ...
      'Style','popupmenu', 'Callback', @hPlotOptionMenuCallback);
    
    % Set up for default
    index = get(hPlotOptionMenu,'Value');
    
    switch index
      
      case 1 % Scatter plot
        
        % Axis listbox label
        hPlotText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
          'Position',[0.03 0.73+0.1 0.18 0.02], ...
          'String','Axis 1             Axis 2             Axis 3', 'Style','text', ...
          'HorizontalAlignment','Left');
        
        % Add listboxes for plotting
        % Trade space point list box
        hPlot1 = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
          'Position',[0.02 0.58+0.1 0.06 0.15], 'String',[in.ts_label,in.obj_label], 'Style','listbox');
        hPlot2 = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
          'Position',[0.09 0.58+0.1 0.06 0.15], 'String',[in.ts_label,in.obj_label], 'Style','listbox');
        hPlot3 = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
          'Position',[0.16 0.58+0.1 0.06 0.15], 'String',[in.ts_label,in.obj_label], 'Style','listbox');
        
      case 2 % Parallel coordinates plot

    end
    
    %     hPlotOptionMenuCallback(hPlotOptionMenu, eventdata);
    
    % Constraint listbox label
    hPlotText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.02 0.53+0.1 0.18 0.02], 'String','Constraints:', ...
      'Style','text', 'HorizontalAlignment','Left');

%     % Constraint listbox - CURRENTLY DOES NOTHING
%     hConstraintPlot = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
%       'Position',[0.02 0.38+0.1 0.16 0.15], 'String',in.constr(:,1), 'Style','listbox');
    
    % Constraint plotting listbox label
    hPlotText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.02 0.33+0.1 0.18 0.02], 'String','Constraints Plot Type:', ...
      'Style','text', 'HorizontalAlignment','Left');

    % Constraint plotting listbox
    hConstraintTypePlot = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[0.02 0.18+0.1 0.16 0.15], 'String',{'Leave','Remove','Color','Symbol'}, ...
      'Style','listbox');
    
    % Plot button
    hPlotTSButton = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) 0.15 vel_pos(3) 0.05], ...
      'String','Plot', 'Callback', @hPlotTSButtonCallback);

    %----------------------------------------------------------------------
    % Plot trade space analysis
    function hPlotTSButtonCallback(hObject, eventdata)
      
      % Update trade space / objective space plot function

      % Get trade space dimensions / objectives to be plotted
      Axis_Set = [get(hPlot1,'Value') get(hPlot2,'Value') get(hPlot3,'Value')];
      for ctr = 1 : 1 : length(Axis_Set)
        if Axis_Set(ctr) <= length(in.ts_label)
          Axis_Data(ctr,:) = reshape(od.pop(Axis_Set(ctr),:,1:od.iter),1,in.num_ind*od.iter);
        else
          Axis_Data(ctr,:) = reshape(od.fit_save(Axis_Set(ctr)-length(in.ts_label),:,1:od.iter),1,in.num_ind*od.iter);
        end
      end
      
      Axis_Label_Set = get(hPlot1,'String');
      Axis_Label = Axis_Label_Set(Axis_Set);
      
      % Enfore constraints on plot
      plot_option_constraint_set = get(hConstraintTypePlot,'String');
      plot_option_constraint = plot_option_constraint_set{get(hConstraintTypePlot,'Value')};
      
      switch plot_option_constraint
        
        case 'Leave'
          
          % Do nothing
          
        case 'Remove'
          
          % Remove points that broke constraint from plot
          % Determine individuals that broke constraint and set to NaN
          I_broke = [];
          for ctr = 1 : 1 : od.num_constr
            constr_set = reshape(od.constr(ctr,:,1:od.iter),1,in.num_ind*od.iter);
            I_satisfy = eval(['find(constr_set',in.constr{ctr,3},num2str(in.constr{ctr,2}),')']);
            I_broke = [I_broke setdiff([1:1:in.num_ind*od.iter],I_satisfy)];
          end
          I_broke = unique(I_broke);
          Axis_Data(:,I_broke) = NaN;
          
        case 'Color'
          
          % Recolor points if the violate a constraint
          I_broke = [];
          for ctr = 1 : 1 : od.num_constr
            constr_set = reshape(od.constr(ctr,:,1:od.iter),1,in.num_ind*od.iter);
            I_satisfy = eval(['find(constr_set',in.constr{ctr,3},num2str(in.constr{ctr,2}),')']);
            I_broke = [I_broke setdiff([1:1:in.num_ind*od.iter],I_satisfy)];
          end
          I_broke = unique(I_broke);
          Axis_Data_Constr = Axis_Data(:,I_broke);
          Axis_Data(:,I_broke) = NaN;
          
        case 'Symbol'
          
          % Change symbol of points that violated constraint
          % Recolor points if the violate a constraint
          I_broke = [];
          for ctr = 1 : 1 : od.num_constr
            constr_set = reshape(od.constr(ctr,:,1:od.iter),1,in.num_ind*od.iter);
            I_satisfy = eval(['find(constr_set',in.constr{ctr,3},num2str(in.constr{ctr,2}),')']);
            I_broke = setdiff([1:1:in.num_ind*od.iter],I_satisfy);
            Axis_Data_Constr(ctr).val = Axis_Data(:,I_broke);
            Axis_Data(:,I_broke) = NaN;
          end
          
      end
      
      % Set up axes
      set(gcf,'CurrentAxes',hPlotAxes);
      cla; % Clear axes
      
      % Plot and label
      plot3(Axis_Data(1,:),Axis_Data(2,:),Axis_Data(3,:),'b.');
      hold on;
      if strcmp(plot_option_constraint,'Leave')
        legend_string = 'legend(''All'')';
      elseif strcmp(plot_option_constraint,'Remove')
        legend_string = 'legend(''Constraint Satisfied'')';
      elseif strcmp(plot_option_constraint,'Color')
        plot3(Axis_Data_Constr(1,:),Axis_Data_Constr(2,:),Axis_Data_Constr(3,:),'r.');
        legend_string = 'legend(''Constraint Satisfied'',''Constraint Violated'',''Location'',''NorthEastOutside'');';
%         axis([-2 1.5 -1 1 0 20]);
      elseif strcmp(plot_option_constraint,'Symbol')
        legend_string = 'legend(''Constraint Satisfied'',';
        for ctr = 1 : 1 : od.num_constr
          plot3(Axis_Data_Constr(ctr).val(1,:),Axis_Data_Constr(ctr).val(2,:), ...
            Axis_Data_Constr(ctr).val(3,:),symbol_set{ctr});
          legend_string = [legend_string '''',in.constr{ctr},''''];
          if ctr ~= od.num_constr
            legend_string = [legend_string,','];
          else
            legend_string = [legend_string,');'];
          end
        end
      end
      grid on;
      xlabel(Axis_Label{1});
      ylabel(Axis_Label{2});
      zlabel(Axis_Label{3});
      eval(legend_string);

    end

  end

  %----------------------------------------------------------------------
  % Inertia analysis
  function hInertiaMenuitemCallback(hObject, eventdata)
    
    clf(hMainFigure);
    % Turn on figure toolbar
    set(gcf,'toolbar','figure');
    
    %%%%%%%%%%%%%%%%%%%%%
    %% UI Declarations %%
    %%%%%%%%%%%%%%%%%%%%%
    
    % Velocity plot axes
    hPlotAxes = axes('Parent', hMainFigure, 'Units', 'normalized', ...
      'Position',pos_axes);

    % Update inertia plot
    set(gcf,'CurrentAxes',hPlotAxes);
    cla; % Clear axes
    hold on;
    grid on;
    xlabel('Iteration');
    ylabel('Inertia');
    title('Inertia vs. Iteration');
    plot([1:1:od.iter],od.w(1:od.iter));
    hold off;
    
  end

  %----------------------------------------------------------------------
  % Summary sheet
  function hSummaryMenuitemCallback(hObject, eventdata)
    
    clf(hMainFigure);
    % Turn on figure toolbar
    set(gcf,'toolbar','figure');
    
    % Pareto front point static text label (develop time string)
    time_string = [];
    remainder = od.time;
    % Hours
    if remainder >= 3600
      % Optimization took hours
      time_string = [time_string num2str(floor(od.time/3600)),' hours '];
      remainder = mod(od.time,3600);
    end
    % Minutes
    if remainder >= 60
      time_string = [time_string num2str(floor(remainder/60)),' minutes '];
      remainder = mod(remainder,60);
    end
    % Seconds
    time_string = [time_string num2str(remainder),' seconds '];

    hTimeText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) .8 0.8 0.02], ...
      'String',['Execution Time: ',time_string], 'Style','text', 'HorizontalAlignment','Left');
    
    hIterText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) 0.75 0.8 0.02], ...
      'String',['Iterations: ',num2str(od.iter),' (max ',num2str(in.max_iter),')'], ...
      'Style','text', 'HorizontalAlignment','Left');

  end

  %----------------------------------------------------------------------
  % Pareto analysis
  function hParetoMenuitemCallback(hObject, eventdata)
    
    clf(hMainFigure);
    % Turn on figure toolbar
    set(gcf,'toolbar','figure');
    
    %%%%%%%%%%%%%%%%%%%%%
    %% UI Declarations %%
    %%%%%%%%%%%%%%%%%%%%%
    
    % Pareto front plot axes
    hPlotAxes = axes('Parent', hMainFigure, 'Units', 'normalized', ...
      'Position',pos_axes);
    
    % Pareto front point list box
    hParetoMenu = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',vel_pos, 'String',{}, 'Style','listbox');

    % Pareto front point static text label
    hParetoText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) vel_pos(2)+vel_pos(4)+0.001 vel_pos(3) 0.02], ...
      'String','Pareto Front Point:', 'Style','text', 'HorizontalAlignment','Left');

    % Trade space point static text label
    hTradeText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) vel_pos(2)-0.05 vel_pos(3) 0.04], ...
      'String','Trade Space Position:', 'Style','text', 'HorizontalAlignment','Left');

    % Trade space point list box
    hTradeMenu = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) vel_pos(2)-0.20 vel_pos(3) vel_pos(4)], 'String',{}, 'Style','listbox');

    % Get trade space location button
    hGetTradeButton = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) vel_pos(2)-0.28 vel_pos(3) 0.05], ...
      'String','Locate', 'Callback', @hParetoButtonCallback);
    
    
    
%     %% COLOR CODE - MANUAL %%
%     % Determine maximum and minimum values
%     max_val = max(od.pos_arch(1,:,od.iter));
%     min_val = min(od.pos_arch(1,:,od.iter));
%     range = max_val-min_val;
%     col = [(od.pos_arch(1,:,od.iter)-min_val)/range; zeros(1,in.num_arch); 1-(od.pos_arch(1,:,od.iter)-min_val)/range];
%     S = 100*ones(1,in.num_arch);
% %     for ctr = 1 : 1 : in.num_ind
% %       temp_max = max(od.pos_arch(1,ctr,od.iter));
% %       temp_min = min(od.pos_arch(1,ctr,od.iter));
% %       if temp_max > max_val
% %         max_val = temp_max;
% %       end
% %       if temp_min < min_val
% %         min_val = temp_min;
% %       end
% %     end
% %     range = max_val-min_val;
% %     
% %     for n=1:length(fit)
% %       if isnan(fit(n))~=1 & fit(n)>6000
% %           col=[(fit(n)-little)/range 0 1-(fit(n)-little)/range];
% %           plot3(popb1(n),popgamma(n),fit(n),'.','Color',col)
% %           hold on
% %       end
% %     end
%     
%     col = col';
    
    % Plot Pareto Front
    % Set up axes
    set(gcf,'CurrentAxes',hPlotAxes);
    cla; % Clear axes
    hold on;
    grid on;
    title('Pareto Front');
    % Plot Pareto front
    if in.num_obj == 2
      % 2-objective optimization
      plot(pareto_front(1,:),pareto_front(2,:),'b.');
%       plot(od.fit_arch(1,:,od.iter),od.fit_arch(2,:,od.iter),'b.');
%       scatter(od.fit_arch(1,:,od.iter),od.fit_arch(2,:,od.iter),S,col,'.');
      grid on;
      xlabel(in.obj_label{1});
      ylabel(in.obj_label{2});
    elseif in.num_obj == 3
      % 3-objective optimization
      plot3(pareto_front(1,:),pareto_front(2,:),pareto_front(3,:),'b.');
%       plot3(od.fit_arch(1,:,od.iter),od.fit_arch(2,:,od.iter), ...
%         od.fit_arch(3,:,od.iter),'b.');
%       scatter3(od.fit_arch(1,:,od.iter),od.fit_arch(2,:,od.iter), ...
%         od.fit_arch(3,:,od.iter),S,col,'.');
      grid on;
      xlabel(in.obj_label{1});
      ylabel(in.obj_label{2});
      zlabel(in.obj_label{3});
    else
      error('Too many objectives');
    end
    
    %----------------------------------------------------------------------
    function hParetoButtonCallback(hObject, eventdata)
      
      h = datacursormode;
      info = getCursorInfo(h);

      % Get location on Pareto Front
      pareto_string = cell(1,length(info.Position));
      for ctr = 1 : 1 : length(info.Position)
        pareto_string{ctr} = num2str(info.Position(ctr));
      end

      % Output location to listbox
      set(hParetoMenu,'String',pareto_string);

      % Locate corresponding trade space location
      for ctr = 1 : 1 : od.iter
        [TF,LOC] = ismember(info.Position,od.fit(:,:,ctr)','rows');
        I = find(TF == 1);
        if ~isempty(I)
          % Solution match found
          trade_space_pos = od.pop(:,LOC,ctr);
          break;
        end
      end

      % Output location to listbox
      trade_string = cell(1,length(info.Position));
      for ctr = 1 : 1 : length(trade_space_pos)
        trade_string{ctr} = num2str(trade_space_pos(ctr));
      end
      set(hTradeMenu,'String',trade_string);
      format long;
      trade_space_pos % MJG output to screen
      format short;
    end
    
  end
  
  %----------------------------------------------------------------------
  % Velocity analysis
  function hVelMenuitemCallback(hObject, eventdata)
    
    vel_plot_type = {'Normalized','Absolute'};
    
    clf(hMainFigure);
    % Turn on figure toolbar
    set(gcf,'toolbar','figure');
    
    %%%%%%%%%%%%%%%%%%%%%
    %% UI Declarations %%
    %%%%%%%%%%%%%%%%%%%%%
    
    % Velocity plot axes
    hPlotAxes = axes('Parent', hMainFigure, 'Units', 'normalized', ...
      'Position',pos_axes);

    % Velocity dimensions list box
    hVelMenu = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',vel_pos, 'String',in.ts_label, 'Style','listbox', 'Max',2);

    % Velocity dimensions static text label
    hVelText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) vel_pos(2)+vel_pos(4)+0.001 vel_pos(3) 0.02], ...
      'String','Dimensions:', 'Style','text', 'HorizontalAlignment','Left');

    % Plot type static text label
    hVelText = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) vel_pos(2)-0.05 vel_pos(3) 0.02], ...
      'String','Plot Type:', 'Style','text', 'HorizontalAlignment','Left');

    % Plot type popup (normalized or absolute velocities)
    hVelPlotMenu = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) vel_pos(2)-0.10 vel_pos(3) 0.05], ...
      'String',vel_plot_type, 'Style','popupmenu');

    % Plot velocity button
    hPlotButton = uicontrol('Parent', hMainFigure, 'Units','normalized', ...
      'Position',[vel_pos(1) vel_pos(2)-0.15 vel_pos(3) 0.05], ...
      'String','Plot Velocities', 'Callback', @hPlotButtonCallback);
    
    %----------------------------------------------------------------------
    function hPlotButtonCallback(hObject, eventdata)
      
      % This function updates the velocity plot function. Velocities in any
      % trade space dimension may be plotted simultaneously. The legend is
      % automatically generated and placed on the top right outside of the axes.
      
      % Get dimensions to be plotted
      plot_dim = get(hVelMenu,'Value');
      plot_type = get(hVelPlotMenu,'Value');
      
      % Construct matrix of plot data - iteration vector
      Axis_Set = [];
      for ctr = 1 : 1 : od.iter
        Axis_Set = [Axis_Set ctr*ones(1,in.num_ind)]; % Iteration number
      end
      
      % Construct matrix of plot data - trade space dimension vectors
      for ctr = 1 : 1 : length(plot_dim)
        Axis_Set(ctr+1,:) = reshape(od.vel(plot_dim(ctr),:,1:od.iter),1,in.num_ind*od.iter);
        % Normalize data if normalization selected - assumes in.vel_max and
        % in.vel_min are of same magnitude
        if strcmp(vel_plot_type(plot_type),'Normalized')
          Axis_Set(ctr+1,:) = Axis_Set(ctr+1,:)/in.vel_max(ctr);
        end
      end
      
      % Set up plot
      set(gcf,'CurrentAxes',hPlotAxes);
      cla; % Clear axes
      hold on;
      grid on;
      xlabel('Iteration');
      ylabel([vel_plot_type{plot_type},' Velocity']);
      title([vel_plot_type{plot_type},' Velocity vs. Iteration']);
      
      % Plot velocities for each iteration
      for ctr = 1 : 1 : length(plot_dim)
        plot(Axis_Set(1,:),Axis_Set(ctr+1,:),plot_set{plot_dim(ctr)});
      end
      hold off;

      % Construct legend command
      plot_dim_set = get(hVelMenu,'String');
      legend_string = 'legend('; 
      for ctr = 1 : 1 : length(plot_dim) % int2str(plot_dim(ctr_dim))
        legend_string = [legend_string,'''',plot_dim_set{plot_dim(ctr)},''''];
        if ctr ~= length(plot_dim)
          legend_string = [legend_string,','];
        else
          legend_string = [legend_string,',''Location'',''NorthEastOutside'');'];
        end
      end
      eval(legend_string);
      
    end
    
  end

  %----------------------------------------------------------------------


end % end of mikegui
