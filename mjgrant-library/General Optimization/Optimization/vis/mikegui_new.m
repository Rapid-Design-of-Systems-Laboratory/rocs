function varargout = mikegui(varargin)
% MIKEGUI Mike's developing GUI for optimization.

% NASA JSC - DM42 / Michael J. Grant on July 2008

% Notes: 'HandleVisibility','callback' does not seem to allow children to be
% deleted by clf. Removed from all objects except menu items since they are 
% deleted (maybe try using multiple tabs).

%%%%%%%%%%%%
%% Inputs %%
%%%%%%%%%%%%

clear all; pack; close all; clc;
pos = [0.01 0.05 0.4 0.2];
pos_axes = [0.20 0.10 0.70 0.70];
vel_pos = [0.02 0.66 0.10 0.15];
plot_set = {'b.','r.','g.','m.','k.','b*','r*','g*','m*','k*', ...
  'b+','r+','g+','m+','k+','bo','ro','go','mo','ko'};
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
  errordlg('results.mat not found');
end

%%%%%%%%%%%%%%%%%%%%
%% Error Checking %%
%%%%%%%%%%%%%%%%%%%%

% Send error if not enough plot sets available
if length(plot_set) < od.dim
  errordlg(['Not enough plotting types. Increase number of plotting ', ...
    'types in Inputs section to match number of trade space dimensions.']);
  return;
end

% Send error if not enough symbols available for constraints
if isfield(in,'constr') && length(symbol_set) < length(in.constr(:,1))
  errordlg(['Not enough symbol types for constraints. Increase number of ', ...
    'symbols in Inputs section to match number of constraints.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UI Declarations - Main Window %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main GUI figure
hMainFigure = figure('Name', mfilename, 'Toolbar','none', ... % 'MenuBar','none'
  'NumberTitle','off', 'Color', get(0, 'defaultuicontrolbackgroundcolor'), ...
  'Units','normalized', 'Position',pos, 'HandleVisibility','callback');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Callback and Utility Functions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



