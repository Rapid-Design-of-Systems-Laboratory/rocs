function [str2] = matlab2math_str(str)
%
% This function converts a string from Matlab typesetting to Mathematica for
% symbolic integration.
%

% Versions:
%   V0 - 08/19/08
%     Initial creation of script for various trig functions
%

% Author: GT SSDL / Michael Grant on Oct. 19, 2008

%%%%%%%%%%%%%%%%
%% Initialize %%
%%%%%%%%%%%%%%%%

% Initialize string to modify
str2 = str;

% Determine open and close parenthesis indices
[ind_op,ind_cp] = find_parenthesis_indices(str);

%%%%%%%%%%%%%%%%%%%%
%% Convert String %%
%%%%%%%%%%%%%%%%%%%%

% Convert "--(" to "+(" - only need to do this for Mathematica, not when going
% back to Matlab
k_mm = strfind(str2,'--(');
while ~isempty(k_mm)
  
  % Convert string of "--" - find index of start
  ind = k_mm(1);
  
%   % Determine closing parenthesis index
%   I = find(ind_op>ind,1);
%   
%   % Convert parenthesis closing the atan argument first since going to insert
%   % longer string
%   str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcSin
  % Add space so add equivalent length string so don't mess up open and close
  % parenthesis indices
  str2 = [str2(1:ind-1),' +(',str2(ind+3:end)];
  
  % Need to reset atan string index and open and close parenthesis positions
  k_mm = strfind(str2,'--(');
  
%   % Find open and close brackets
%   [ind_op,ind_cp] = find_parenthesis_indices(str2)
%   fprintf('\nhere\n');
  
end

% Convert asinh to arcsinh - search before asinh
k_asinh = strfind(str2,'asinh(');
while ~isempty(k_asinh)
  
  % Convert string of atanh - find index of start
  ind = k_asinh(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the atanh argument first since going to insert
  % longer string
  str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcSinh
  str2 = [str2(1:ind-1),'ArcSinh[',str2(ind+6:end)];
  
  % Need to reset asin string index and open and close parenthesis positions
  k_asinh = strfind(str2,'asinh(');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert asinh to arcsinh - search before asinh
k_acsch = strfind(str2,'acsch(');
while ~isempty(k_acsch)
  
  % Convert string of acsch - find index of start
  ind = k_acsch(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the atanh argument first since going to insert
  % longer string
  str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcSinh
  str2 = [str2(1:ind-1),'ArcCsch[',str2(ind+6:end)];
  
  % Need to reset asin string index and open and close parenthesis positions
  k_acsch = strfind(str2,'acsch(');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert atanh to arctanh - search before atanh
k_atanh = strfind(str2,'atanh(');
while ~isempty(k_atanh)
  
  % Convert string of atanh - find index of start
  ind = k_atanh(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the atanh argument first since going to insert
  % longer string
  str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcTanh
  str2 = [str2(1:ind-1),'ArcTanh[',str2(ind+6:end)];
  
  % Need to reset atan string index and open and close parenthesis positions
  k_atanh = strfind(str2,'atanh(');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert atan to arctan - search before tan
k_atan = strfind(str2,'atan(');
while ~isempty(k_atan)
  
  % Convert string of atan - find index of start
  ind = k_atan(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the atan argument first since going to insert
  % longer string
  str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcTan
  str2 = [str2(1:ind-1),'ArcTan[',str2(ind+5:end)];
  
  % Need to reset atan string index and open and close parenthesis positions
  k_atan = strfind(str2,'atan(');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end



% Convert conj
k_conj = strfind(str2,'conj(');
while ~isempty(k_conj)
  
  % Convert string of atan - find index of start
  ind = k_conj(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the atan argument first since going to insert
  % longer string
  str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcTan
  str2 = [str2(1:ind-1),'Conjugate[',str2(ind+5:end)];
  
  % Need to reset atan string index and open and close parenthesis positions
  k_conj = strfind(str2,'conj(');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert abs
k_abs = strfind(str2,'abs(');
for ctr = 1 : 1 : length(k_abs)
  
  % Convert string of cos
  ind = k_abs(ctr);
  str2(ind:ind+3) = 'Abs[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the cos argument
  str2(ind_cp(I)) = ']';
  
end

% Convert cosh
k_cosh = strfind(str2,'cosh(');
for ctr = 1 : 1 : length(k_cosh)
  
  % Convert string of cosh
  ind = k_cosh(ctr);
  str2(ind:ind+4) = 'Cosh[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the cot argument
  str2(ind_cp(I)) = ']';
  
end

% Convert sinh
k_sinh = strfind(str2,'sinh(');
for ctr = 1 : 1 : length(k_sinh)
  
  % Convert string of sinh
  ind = k_sinh(ctr);
  str2(ind:ind+4) = 'Sinh[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the cot argument
  str2(ind_cp(I)) = ']';
  
end

% Convert acot to ArcCot
k_acot = strfind(str2,'acot(');
while ~isempty(k_acot)
  
  % Convert string of atan - find index of start
  ind = k_acot(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the atan argument first since going to insert
  % longer string
  str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcSin
  str2 = [str2(1:ind-1),'ArcCot[',str2(ind+5:end)];
  
  % Need to reset atan string index and open and close parenthesis positions
  k_acot = strfind(str2,'acot(');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert asec to ArcSec
k_asec = strfind(str2,'asec(');
while ~isempty(k_asec)
  
  % Convert string of atan - find index of start
  ind = k_asec(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the atan argument first since going to insert
  % longer string
  str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcSin
  str2 = [str2(1:ind-1),'ArcSec[',str2(ind+5:end)];
  
  % Need to reset atan string index and open and close parenthesis positions
  k_asec = strfind(str2,'asec(');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert asin to ArcSin
k_asin = strfind(str2,'asin(');
while ~isempty(k_asin)
  
  % Convert string of atan - find index of start
  ind = k_asin(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the atan argument first since going to insert
  % longer string
  str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcSin
  str2 = [str2(1:ind-1),'ArcSin[',str2(ind+5:end)];
  
  % Need to reset atan string index and open and close parenthesis positions
  k_asin = strfind(str2,'asin(');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert acos to ArcCos
k_acos = strfind(str2,'acos(');
while ~isempty(k_acos)
  
  % Convert string of atan - find index of start
  ind = k_acos(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the atan argument first since going to insert
  % longer string
  str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcTan
  str2 = [str2(1:ind-1),'ArcCos[',str2(ind+5:end)];
  
  % Need to reset atan string index and open and close parenthesis positions
  k_acos = strfind(str2,'acos(');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert cos
k_cos = strfind(str2,'cos(');
for ctr = 1 : 1 : length(k_cos)
  
  % Convert string of cos
  ind = k_cos(ctr);
  str2(ind:ind+3) = 'Cos[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the cos argument
  str2(ind_cp(I)) = ']';
  
end

% Convert sin
k_sin = strfind(str2,'sin(');
for ctr = 1 : 1 : length(k_sin)
  
  % Convert string of sin
  ind = k_sin(ctr);
  str2(ind:ind+3) = 'Sin[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the sin argument
  str2(ind_cp(I)) = ']';
  
end

% Convert csc
k_csc = strfind(str2,'csc(');
for ctr = 1 : 1 : length(k_csc)
  
  % Convert string of csc
  ind = k_csc(ctr);
  str2(ind:ind+3) = 'Csc[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the csc argument
  str2(ind_cp(I)) = ']';
  
end

% Convert sec
k_sec = strfind(str2,'sec(');
for ctr = 1 : 1 : length(k_sec)
  
  % Convert string of sec
  ind = k_sec(ctr);
  str2(ind:ind+3) = 'Sec[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the sec argument
  str2(ind_cp(I)) = ']';
  
end

% Convert cot
k_cot = strfind(str2,'cot(');
for ctr = 1 : 1 : length(k_cot)
  
  % Convert string of cot
  ind = k_cot(ctr);
  str2(ind:ind+3) = 'Cot[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the cot argument
  str2(ind_cp(I)) = ']';
  
end

% Convert tan
k_tan = strfind(str2,'tan(');
for ctr = 1 : 1 : length(k_tan)
  
  % Convert string of tan
  ind = k_tan(ctr);
  str2(ind:ind+3) = 'Tan[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the tan argument
  str2(ind_cp(I)) = ']';
  
end

% Convert pi
k_pi = strfind(str2,'pi');
for ctr = 1 : 1 : length(k_pi)
  
  % Convert string of pi
  ind = k_pi(ctr);
  str2(ind:ind+1) = 'Pi';
  
end

% Convert sqrt
k_sqrt = strfind(str2,'sqrt(');
for ctr = 1 : 1 : length(k_sqrt)
  
  % Convert string of sqrt
  ind = k_sqrt(ctr);
  str2(ind:ind+4) = 'Sqrt[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the sqrt argument
  str2(ind_cp(I)) = ']';
  
end

% Convert exp
k_exp = strfind(str2,'exp(');
for ctr = 1 : 1 : length(k_exp)
  
  % Convert string of exp
  ind = k_exp(ctr);
  str2(ind:ind+3) = 'Exp[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the exp argument
  str2(ind_cp(I)) = ']';
  
end

% Convert log
k_log = strfind(str2,'log(');
for ctr = 1 : 1 : length(k_log)
  
  % Convert string of log
  ind = k_log(ctr);
  str2(ind:ind+3) = 'Log[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the log argument
  str2(ind_cp(I)) = ']';
  
end

% Convert Function - these are custom functions provided by user
k_Function = strfind(str2,'Function(');
for ctr = 1 : 1 : length(k_Function)
  
  % Convert string of Function
  ind = k_Function(ctr);
  str2(ind:ind+8) = 'Function[';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the log argument
  str2(ind_cp(I)) = ']';
  
end

% Convert imag to Im
k_imag = strfind(str2,'imag(');
while ~isempty(k_imag)
  
  % Convert string of imag - find index of start
  ind = k_imag(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the imag argument first since going to insert
  % longer string
  str2(ind_cp(I)) = ']';
  
  % Split string for inserting since strings not same length and insert ArcTan
  str2 = [str2(1:ind-1),'Im[',str2(ind+5:end)];
  
  % Need to reset imag string index and open and close parenthesis positions
  k_imag = strfind(str2,'imag(');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_op,ind_cp] = find_parenthesis_indices(str)

  % Find all parenthesis indices (used to determine ending of cos/sin
  % arguments, etc)
  k_op = strfind(str,'('); % open parenthesis
  k_cp = strfind(str,')'); % close parenthesis

  % Determine indices of matching open and close parenthesis
  k_p_sort = sort([k_op,k_cp]);

  ind_op = nan(1,length(k_op));
  ind_cp = nan(1,length(k_cp));
  ctr_op = 0; % Array counter
  for ctr = 1 : 1 : length(k_p_sort)

    if ismember(k_p_sort(ctr),k_op)
      % Add to op index
      ctr_op = ctr_op + 1;
      ind_op(ctr_op) = k_p_sort(ctr);
    else
      % Determine corresponding array location of closing parenthesis. Start at
      % end of open parenthesis array and work backwards until find last open
      % parenthesis.
      I = find(isnan(ind_op) == 0);
      while_ctr = 0;
      while 1
        if isnan(ind_cp(I(end-while_ctr)))
          ind_cp(I(end-while_ctr)) = k_p_sort(ctr);
          break;
        else
          while_ctr = while_ctr + 1;
        end
      end
    end

  end
  
return
  
