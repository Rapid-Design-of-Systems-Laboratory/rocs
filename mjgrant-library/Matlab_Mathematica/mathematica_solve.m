function [integrand_int_strSet] = mathematica_solve(syms_str_in, ...
  syms_val_in,var,assumptions)
  
% Convert to Mathematica string
integrand_math = matlab2math_str(char(syms_str_in));
val_math = matlab2math_str(char(syms_val_in));

% Solve equation for variable of interest. Interestingly, solving
% Exp[FakeVariable]*(expression) == 0 can result in a very different (and
% compact result without conditional statements than simply solving expression
% == 0. Setting NumberMarks to false removes ` character that appears when
% Mathematica is communicating precision information.
if strcmp(val_math,'0')
  integrand_int_math = math(['InputForm[Simplify[Solve[Exp[FakeVariable]*(', ...
    char(integrand_math),') ==',char(val_math),',',char(var), ...
    '],Assumptions -> {',assumptions,'}],NumberMarks -> False]']);
else
    integrand_int_math = math(['InputForm[Simplify[Solve[', ...
    char(integrand_math),' ==',char(val_math),',',char(var), ...
    '],Assumptions -> {',assumptions,'}],NumberMarks -> False]']);
end

% if strcmp(val_math,'0')
%   integrand_int_math = math(['InputForm[Solve[Exp[FakeVariable]*(', ...
%     char(integrand_math),') ==',char(val_math),',',char(var), ...
%     '],NumberMarks -> False]']);
% else
%   integrand_int_math = math(['InputForm[Solve[', ...
%     char(integrand_math),' ==',char(val_math),',',char(var), ...
%     '],NumberMarks -> False]']);
% end

% Convert back to Matlab expression
integrand_int_str = math2matlab_str(integrand_int_math);

% Protect for custom functions that have commas. If Solve[ in solution, then know
% Mathematica didn't finish correctly.
if strfind(integrand_int_str,'Solve[')
	integrand_int_strSet{1} = integrand_int_str;
	return;
end

% Get all solutions from Mathematica string
I1 = strfind(integrand_int_str,'{');
I2 = strfind(integrand_int_str,'}');
I_remove = union(I1,I2);
I = (1:1:length(integrand_int_str));
I_keep = setdiff(I,I_remove);
integrand_int_str = integrand_int_str(I_keep);
I_divider = [0,strfind(integrand_int_str,','),length(integrand_int_str)+1];
integrand_int_strSet = cell(length(I_divider)-1,1);

% For R2010b, need to not use any symbolic conversions. Results change as
% Matlab is upgraded.
for ctr = 1 : 1 : length(I_divider)-1
  vString = integrand_int_str(I_divider(ctr)+1:I_divider(ctr+1)-1);
  % Remove "expr -" that results from Mathematica after ">" removed from "expr
  % ->". After trim, know is first length(expr)+2 characters
  vStringTrim = strtrim(vString);
  integrand_int_strSet{ctr} = vStringTrim(length(char(var))+2:end);
end

% % To try and prevent Matlab segmentation faults when run for second time.
% math('quit');

return
