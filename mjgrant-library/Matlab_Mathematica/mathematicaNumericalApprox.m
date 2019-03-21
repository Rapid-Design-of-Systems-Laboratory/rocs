function [syms_str_out] = mathematicaNumericalApprox(syms_str_in)
  
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

  % Solve equation for variable of interest
  integrand_int_math = math(['InputForm[N[', ...
    char(integrand_math),']]']);

  % Convert back to Matlab expression
  integrand_int_str = math2matlab_str(integrand_int_math);

  % Convert back to Matlab symbolic variables - old way
  integrand_int_str = math2matlab_str(integrand_int_math)

  % Assign output
  if strcmp(varData.class,'sym')
    syms_str_out(ctr) = sym(integrand_int_str);
  elseif strcmp(varData.class,'cell')
    syms_str_out{ctr} = integrand_int_str;
  end
  
end



return

