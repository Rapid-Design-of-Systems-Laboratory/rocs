function [str2] = math2matlab_str(str)
%
% This function converts a string from Mathematica typesetting to Matlab for
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

% Find open and close brackets
[ind_op,ind_cp] = find_parenthesis_indices(str);

%%%%%%%%%%%%%%%%%%%%
%% Convert String %%
%%%%%%%%%%%%%%%%%%%%

% Convert arctanh to atanh
k_atanh = strfind(str2,'ArcTanh[');
for ctr = 1 : 1: length(k_atanh)
% while ~isempty(k_atanh)
  
  % Convert string of atanh - find index of start
  ind = k_atanh(ctr);
  str2(ind:ind+7) = '  atanh(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the ArcTanh argument first since going to insert
  % shorter string
  str2(ind_cp(I)) = ')';
  
%   % Split string for inserting since strings not same length and insert ArcTanh
%   str2 = [str2(1:ind-1),'atanh(',str2(ind+8:end)];
%   
%   % Need to reset atanh string index and open and close parenthesis positions
%   k_atanh = strfind(str2,'ArcTanh[');
%   
%   % Find open and close brackets
%   [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert arccsch to acsch
k_acsch = strfind(str2,'ArcCsch[');
for ctr = 1 : 1: length(k_acsch)
% while ~isempty(k_atanh)
  
  % Convert string of acsch - find index of start
  ind = k_acsch(ctr);
  str2(ind:ind+7) = '  acsch(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the ArcCsch argument first since going to insert
  % shorter string
  str2(ind_cp(I)) = ')';
  
%   % Split string for inserting since strings not same length and insert ArcTanh
%   str2 = [str2(1:ind-1),'atanh(',str2(ind+8:end)];
%   
%   % Need to reset atanh string index and open and close parenthesis positions
%   k_atanh = strfind(str2,'ArcTanh[');
%   
%   % Find open and close brackets
%   [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert arcsinh to asinh
k_asinh = strfind(str2,'ArcSinh[');
for ctr = 1 : 1: length(k_asinh)
  
  % Convert string of atanh - find index of start
  ind = k_asinh(ctr);
  str2(ind:ind+7) = '  asinh(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the ArcSinh argument first since going to insert
  % shorter string
  str2(ind_cp(I)) = ')';
  
end

% Convert ArcCot to acot
k_acot = strfind(str2,'ArcCot[');
for ctr = 1 : 1 : length(k_acot)
% while ~isempty(k_asec)
  
  % Convert string of atan - find index of start
  ind = k_acot(ctr);
  str2(ind:ind+6) = '  acot(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the ArcTan argument first since going to insert
  % shorter string
  str2(ind_cp(I)) = ')';
  
%   % Split string for inserting since strings not same length and insert ArcTan
%   str2 = [str2(1:ind-1),'asec(',str2(ind+7:end)];
%   
%   % Need to reset atan string index and open and close parenthesis positions
%   k_asec = strfind(str2,'ArcSec[');
%   
%   % Find open and close brackets
%   [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert ArcSec to asec
k_asec = strfind(str2,'ArcSec[');
for ctr = 1 : 1 : length(k_asec)
% while ~isempty(k_asec)
  
  % Convert string of atan - find index of start
  ind = k_asec(ctr);
  str2(ind:ind+6) = '  asec(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the ArcTan argument first since going to insert
  % shorter string
  str2(ind_cp(I)) = ')';
  
%   % Split string for inserting since strings not same length and insert ArcTan
%   str2 = [str2(1:ind-1),'asec(',str2(ind+7:end)];
%   
%   % Need to reset atan string index and open and close parenthesis positions
%   k_asec = strfind(str2,'ArcSec[');
%   
%   % Find open and close brackets
%   [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert ArcSin to asin
k_asin = strfind(str2,'ArcSin[');
for ctr = 1 : 1 : length(k_asin)
% while ~isempty(k_asin)
  
  % Convert string of atan - find index of start
  ind = k_asin(ctr);
  str2(ind:ind+6) = '  asin(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the ArcTan argument first since going to insert
  % shorter string
  str2(ind_cp(I)) = ')';
  
%   % Split string for inserting since strings not same length and insert ArcTan
%   str2 = [str2(1:ind-1),'asin(',str2(ind+7:end)];
%   
%   % Need to reset atan string index and open and close parenthesis positions
%   k_asin = strfind(str2,'ArcSin[');
%   
%   % Find open and close brackets
%   [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert ArcCos to acos
k_acos = strfind(str2,'ArcCos[');
for ctr = 1 : 1 : length(k_acos)
% while ~isempty(k_acos)
  
  % Convert string of atan - find index of start
  ind = k_acos(ctr);
  str2(ind:ind+6) = '  acos(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the ArcTan argument first since going to insert
  % shorter string
  str2(ind_cp(I)) = ')';
  
%   % Split string for inserting since strings not same length and insert ArcTan
%   str2 = [str2(1:ind-1),'acos(',str2(ind+7:end)];
%   
%   % Need to reset atan string index and open and close parenthesis positions
%   k_acos = strfind(str2,'ArcCos[');
%   
%   % Find open and close brackets
%   [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert arctan to atan
k_atan = strfind(str2,'ArcTan[');
for ctr = 1 : 1 : length(k_atan)
% while ~isempty(k_atan)
  
  % Convert string of atan - find index of start
  ind = k_atan(ctr);
  str2(ind:ind+6) = '  atan(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the ArcTan argument first since going to insert
  % shorter string
  str2(ind_cp(I)) = ')';
  
%   % Split string for inserting since strings not same length and insert ArcTan
%   str2 = [str2(1:ind-1),'atan(',str2(ind+7:end)];
%   
%   % Need to reset atan string index and open and close parenthesis positions
%   k_atan = strfind(str2,'ArcTan[');
%   
%   % Find open and close brackets
%   [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert cosh
k_cosh = strfind(str2,'Cosh[');
for ctr = 1 : 1 : length(k_cosh)
  
  % Convert string of cosh
  ind = k_cosh(ctr);
  str2(ind:ind+4) = 'cosh(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the cot argument
  str2(ind_cp(I)) = ')';
  
end

% Convert sinh
k_sinh = strfind(str2,'Sinh[');
for ctr = 1 : 1 : length(k_sinh)
  
  % Convert string of sinh
  ind = k_sinh(ctr);
  str2(ind:ind+4) = 'sinh(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the cot argument
  str2(ind_cp(I)) = ')';
  
end

% Convert conj
k_conj = strfind(str2,'Conjugate[');
for ctr = 1 : 1 : length(k_conj)
% while ~isempty(k_conj)
  
  % Convert string of atanh - find index of start
  ind = k_conj(ctr);
  str2(ind:ind+9) = '     conj(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the ArcTanh argument first since going to insert
  % shorter string
  str2(ind_cp(I)) = ')';
  
%   % Split string for inserting since strings not same length and insert ArcTanh
%   str2 = [str2(1:ind-1),'conj(',str2(ind+10:end)];
%   
%   % Need to reset atanh string index and open and close parenthesis positions
%   k_conj = strfind(str2,'Conjugate[');
%   
%   % Find open and close brackets
%   [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert abs
k_abs = strfind(str2,'Abs[');
for ctr = 1 : 1 : length(k_abs)
  
  % Convert string of abs
  ind = k_abs(ctr);
  str2(ind:ind+3) = 'abs(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the cos argument
  str2(ind_cp(I)) = ')';
  
end

% Convert csc
k_csc = strfind(str2,'Csc[');
for ctr = 1 : 1 : length(k_csc)
  
  % Convert string of csc
  ind = k_csc(ctr);
  str2(ind:ind+3) = 'csc(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the csc argument
  str2(ind_cp(I)) = ')';
  
end

% Convert sec
k_sec = strfind(str2,'Sec[');
for ctr = 1 : 1 : length(k_sec)
  
  % Convert string of sec
  ind = k_sec(ctr);
  str2(ind:ind+3) = 'sec(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the sec argument
  str2(ind_cp(I)) = ')';
  
end

% Convert cot
k_cot = strfind(str2,'Cot[');
for ctr = 1 : 1 : length(k_cot)
  
  % Convert string of cot
  ind = k_cot(ctr);
  str2(ind:ind+3) = 'cot(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the cot argument
  str2(ind_cp(I)) = ')';
  
end

% Convert tan
k_tan = strfind(str2,'Tan[');
for ctr = 1 : 1 : length(k_tan)
  
  % Convert string of tan
  ind = k_tan(ctr);
  str2(ind:ind+3) = 'tan(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the tan argument
  str2(ind_cp(I)) = ')';
  
end

% Convert cos
k_cos = strfind(str2,'Cos[');
for ctr = 1 : 1 : length(k_cos)
  
  % Convert string of cos
  ind = k_cos(ctr);
  str2(ind:ind+3) = 'cos(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the cos argument
  str2(ind_cp(I)) = ')';
  
end

% Convert sin
k_sin = strfind(str2,'Sin[');
for ctr = 1 : 1 : length(k_sin)
  
  % Convert string of sin
  ind = k_sin(ctr);
  str2(ind:ind+3) = 'sin(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the sin argument
  str2(ind_cp(I)) = ')';
  
end

% Convert pi
k_pi = strfind(str2,'Pi');
for ctr = 1 : 1 : length(k_pi)
  
  % Convert string of pi
  ind = k_pi(ctr);
  str2(ind:ind+1) = 'pi';
  
end

% Convert sqrt
k_sqrt = strfind(str2,'Sqrt[');
for ctr = 1 : 1 : length(k_sqrt)
  
  % Convert string of sqrt
  ind = k_sqrt(ctr);
  str2(ind:ind+4) = 'sqrt(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the sqrt argument
  str2(ind_cp(I)) = ')';
  
end

% Convert atan to E to exp(1)
k_e = strfind(str2,'E^');
while ~isempty(k_e)
  
  % Convert string of e^ - find index of start
  ind = k_e(1);
  
  % Split string for inserting since strings not same length and insert exp(1)^
  str2 = [str2(1:ind-1),'exp(1)^',str2(ind+2:end)];
  
  % Need to reset atan string index and open and close parenthesis positions
  k_e = strfind(str2,'E^');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% Convert log
k_log = strfind(str2,'Log[');
for ctr = 1 : 1 : length(k_log)
  
  % Convert string of log
  ind = k_log(ctr);
  str2(ind:ind+3) = 'log(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the log argument
  str2(ind_cp(I)) = ')';
  
end

% Convert Function
k_Function = strfind(str2,'Function[');
for ctr = 1 : 1 : length(k_Function)
  
  % Convert string of Function
  ind = k_Function(ctr);
  str2(ind:ind+8) = 'Function(';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the log argument
  str2(ind_cp(I)) = ')';
  
end

% The following does not work since can change variable names with a capital "I"
% % Convert I to i
% k_i = strfind(str2,'I');
% while ~isempty(k_i)
%   
%   % Convert string of e^ - find index of start
%   ind = k_i(1);
%   
%   % Split string for inserting since strings not same length and insert 1i
%   str2 = [str2(1:ind-1),' i',str2(ind+1:end)];
%   
%   % Need to reset atan string index and open and close parenthesis positions
%   k_i = strfind(str2,'I');
%   
%   % Find open and close brackets
%   [ind_op,ind_cp] = find_parenthesis_indices(str2);
%   
% end

% Undo all of the 1inf expressions that result
k_1inf = strfind(str2,'1inf');
for ctr = 1 : 1 : length(k_1inf)
  
  % Convert string of 1inf
  ind = k_1inf(ctr);
  str2(ind:ind+3) = ' inf';
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the log argument
  str2(ind_cp(I)) = ')';
  
end

% Convert scientific notation "*^" to "*10^"
k_sci = strfind(str2,'*^');
while ~isempty(k_sci)
  
  % Convert string of "*^" - find index of start
  ind = k_sci(1);
  
  % Split string for inserting since strings not same length and insert "*10^"
  str2 = [str2(1:ind-1),'*10^',str2(ind+2:end)];
  
  % Need to reset sci string index
  k_sci = strfind(str2,'*^');
  
end

% Convert ".*" to ".0*"
k_dec = strfind(str2,'.*');
while ~isempty(k_dec)
  
  % Convert string of ".*" - find index of start
  ind = k_dec(1);
  
  % Split string for inserting since strings not same length and insert ".0*"
  str2 = [str2(1:ind-1),'.0*',str2(ind+2:end)];
  
  % Need to reset sci string index
  k_dec = strfind(str2,'.*');
  
end

% % Find any new line feeds and remove them
% ascii_vals0 = abs(str2);
% line_feed0 = find(ascii_vals0 == 10);
% for ctr = 1 : 1 : length(line_feed0)
%   ascii_vals = abs(str2);
%   line_feed = find(ascii_vals == 10);
%   ind = line_feed(1);
%   str2 = [str2(1:ind-1),str2(ind+1:end)];
% end

% Find any new line feeds and convert them to a space
ascii_vals = abs(str2);
line_feed = find(ascii_vals == 10);
str2(line_feed) = ' '; %#ok<FNDSB>

% Convert '>' in string to a space
num_gt = length(strfind(str2,'>'));
for ctr = 1 : 1 : num_gt
  
  k_gt = strfind(str2,'>');
  
  % Convert string of gt
  ind = k_gt(1);
  str2(ind) = ' ';
  
end

% Convert '\' in string to a space
num_gt = length(strfind(str2,'\'));
for ctr = 1 : 1 : num_gt
  
  k_gt = strfind(str2,'\');
  
  % Convert string of gt
  ind = k_gt(1);
  str2(ind) = ' ';
  
end

% Convert Im to imag
k_imag = strfind(str2,'Im[');
while ~isempty(k_imag)
  
  % Convert string of Im - find index of start
  ind = k_imag(1);
  
  % Determine closing parenthesis index
  I = find(ind_op>ind,1);
  
  % Convert parenthesis closing the Im argument
  str2(ind_cp(I)) = ')';
	
  % Split string for inserting since strings not same length and insert imag(
  str2 = [str2(1:ind-1),'imag(',str2(ind+3:end)];
  
  % Need to reset Im string index and open and close parenthesis positions
  k_imag = strfind(str2,'Im[');
  
  % Find open and close brackets
  [ind_op,ind_cp] = find_parenthesis_indices(str2);
  
end

% % Covert all multiple spaces to single spaces
% num_dblspace = length(strfind(str,'  '));
% while num_dblspace ~= 0
%   
%   for ctr = 1 : 1 : num_dblspace
% 
%     k_dblspace = strfind(str2,'  ');
% 
%     % Convert string of gt
%     ind = k_dblspace(1);
%     str2 = [str2(1:ind),str2(ind+2:end)];
% 
%   end
%   
%   num_dblspace = length(strfind(str2,'  '));
%   
% end

% Eliminate all white space. White space may break up variable names which
% Mathematica later interprets as two separate variables that are multiplied
% together.
str2(str2==' ') = ''; % Remove all whitespace

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_op,ind_cp] = find_parenthesis_indices(str)

% Find all bracket indices (used to determine ending of cos/sin arguments, etc)
k_op = strfind(str,'['); % open bracket
k_cp = strfind(str,']'); % close bracket

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
  

