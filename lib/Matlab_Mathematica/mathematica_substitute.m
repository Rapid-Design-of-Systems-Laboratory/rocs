function [syms_str_out] = mathematica_substitute(syms_str_in,var1,var2)

numInt = length(syms_str_in);

% Determine if input sym or str. Added this functionality to only work in
% strings since symbolic engine changes expressions as Matlab is updated.
varData = whos('syms_str_in');

if strcmp(varData.class,'sym')
  syms_str_out = sym(nan(numInt,1));
elseif strcmp(varData.class,'cell')
  syms_str_out = cell(numInt,1);
end

for ctr = 1 : numInt

  % Assign varibles in Mathematica environment
  for ctr2 = 1 : 1 : length(var1)
    math(['InputForm[', matlab2math_str(char(var1(ctr2))) ...
      ,' = ',matlab2math_str(char(var2(ctr2))),']']);
  end
  
  % Convert to Mathematica string
  integrand_math = matlab2math_str(char(syms_str_in(ctr)));
  
  % Integrate in Mathematica
  integrand_int_math = math(['InputForm[Evaluate[', ...
    integrand_math,']]'])
  
%   % Convert back to Matlab symbolic variables - old way
  integrand_int_str = math2matlab_str(integrand_int_math)
%   syms_str_out(ctr1,ctr2) = sym(integrand_int_str);
  
  % Assign output
  if strcmp(varData.class,'sym')
    syms_str_out(ctr) = sym(integrand_int_str);
  elseif strcmp(varData.class,'cell')
    syms_str_out{ctr} = integrand_int_str;
  end
  
  % Clear variables in Matlab environment
  for ctr2 = 1 : 1 : length(var1)
    math(['InputForm[Clear[', matlab2math_str(char(var1(ctr2))),']]']);
  end
  
end

return
  
