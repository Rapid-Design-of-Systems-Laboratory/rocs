function [syms_str_out] = mathematica_simplify_str(syms_str_in,assumptions)

[m,n] = size(syms_str_in);
syms_str_out = sym(nan(m,n));

for ctr1 = 1 : 1 : m
for ctr2 = 1 : 1 : n

  % Convert to Mathematica string
  integrand_math = matlab2math_str(char(syms_str_in(ctr1,ctr2)));
  
  % Integrate in Mathematica
  integrand_int_math = math(['InputForm[Simplify[', ...
    integrand_math,',Assumptions -> {',assumptions,'}]]']);
  
  % Convert back to Matlab symbolic variables
  integrand_int_str = math2matlab_str(integrand_int_math);

  if m == 1 && n == 1
    syms_str_out = sym(integrand_int_str);
  else
    syms_str_out(ctr1,ctr2) = sym(integrand_int_str);
  end
  
%   math('quit');
%   pause(1);
  
end
end

return
  
