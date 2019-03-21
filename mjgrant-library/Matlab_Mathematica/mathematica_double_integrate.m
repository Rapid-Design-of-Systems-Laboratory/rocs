function [syms_str_out] = mathematica_double_integrate( ...
  syms_str_in,assumptions,v1,v2)

syms_str_out = sym(nan(size(syms_str_in)));

for ctr = 1 : 1 : length(syms_str_in)

  % Convert to Mathematica string
  integrand_math = matlab2math_str(char(syms_str_in(ctr)));
  v1_math = matlab2math_str(char(v1));
  v2_math = matlab2math_str(char(v2));
  
  % Integrate in Mathematica
  integrand_int_math = math(['InputForm[Integrate[', ...
    integrand_math,',{v,',v1_math,',',v2_math,'},{u,u1,u2}', ...
    ',Assumptions -> {',assumptions,'}]]']);
  
  % Convert back to Matlab symbolic variables
  integrand_int_str = math2matlab_str(integrand_int_math);
  syms_str_out(ctr) = sym(integrand_int_str);
  
%   math('quit');
%   pause(1);
  
end

return
  
  