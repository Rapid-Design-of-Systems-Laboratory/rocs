function [syms_str_out] = mathematica_subtract(syms_str_in1,syms_str_in2)

numInt = length(syms_str_in1);

% Determine if input sym or str. Added this functionality to only work in
% strings since symbolic engine changes expressions as Matlab is updated.
varData = whos('syms_str_in1');

if strcmp(varData.class,'sym')
  syms_str_out = sym(nan(numInt,1));
elseif strcmp(varData.class,'cell')
  syms_str_out = cell(numInt,1);
end

for ctr = 1 : numInt

  % Convert to Mathematica string
  integrand_math1 = matlab2math_str(char(syms_str_in1(ctr)));
  integrand_math2 = matlab2math_str(char(syms_str_in2(ctr)));
  
  % Integrate in Mathematica
  integrand_int_math = math(['InputForm[Evaluate[(', ...
    integrand_math1,')-(',integrand_math2,')]]']);
  
  % Convert back to Matlab symbolic variables
  integrand_int_str = math2matlab_str(integrand_int_math);
%   syms_str_out(ctr1,ctr2) = sym(integrand_int_str);

  % Assign output
  if strcmp(varData.class,'sym')
    syms_str_out(ctr) = sym(integrand_int_str);
  elseif strcmp(varData.class,'cell')
    syms_str_out{ctr} = integrand_int_str;
  end
  
end

return

