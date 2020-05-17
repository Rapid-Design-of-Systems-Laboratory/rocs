function [syms_str_out] = mathematica_expand(syms_str_in,var)

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

  % Convert to Mathematica string
  integrand_math = matlab2math_str(char(syms_str_in(ctr)));

  % Expand in Mathematica
  if nargin == 1
    integrand_int_math = math(['InputForm[Expand[',integrand_math,']]']);
  else
    integrand_int_math = math(['InputForm[Expand[',integrand_math,',',char(var),']]']);
  end
  
  % Convert back to Matlab symbolic variables
  integrand_int_str = math2matlab_str(integrand_int_math);

  % Assign output
  if strcmp(varData.class,'sym')
    syms_str_out(ctr) = sym(integrand_int_str);
  elseif strcmp(varData.class,'cell')
    syms_str_out{ctr} = integrand_int_str;
  end
  
end

return
  
