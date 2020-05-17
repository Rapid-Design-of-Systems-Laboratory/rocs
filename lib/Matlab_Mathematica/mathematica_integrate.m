function [syms_str_out] = mathematica_integrate(syms_str_in,assumptions,int_var,simplification)

% Default to partial simplification if not specified
if nargin < 4
  simplification = 'partial';
end

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
  integrand_math = matlab2math_str(char(syms_str_in(ctr)))
  
  % Seems that integration of summation of many trig functions can result in
  % integral blowing up into a huge expression. But, if integrate each term
  % separately, get a nice expression.
  
  % Parse string and separate + and appropriate - (subtraction of terms) with
  % spaces in front and back to denote separation in term
  % Get pairs of parenthesis
  [ind_op,ind_cp] = find_parenthesis_indices(integrand_math);
  
  % Determine all + and - of terms if outside of all parenthesis pairs
  I_plus = strfind(integrand_math,'+');
  I_minus = strfind(integrand_math,'-');
  I_plus_minus = union(I_plus,I_minus);
  I_term_divider = [];
  for index = 1 : 1 : length(I_plus_minus)
    outside = true;
    for index2 = 1 : 1 : length(ind_op)
      if (I_plus_minus(index) > ind_op(index2)) && ...
          (I_plus_minus(index) < ind_cp(index2))
        outside = false;
      end
    end
    if outside
      I_term_divider = [I_term_divider,I_plus_minus(index)];
    end
  end
%   I_term_divider = union(I_plus,I_minus); % Automatically sorted
  I_terms = [0,I_term_divider,length(integrand_math)+1];
  
  integrand_int_math = [];
  for term_ctr = 1 : 1 : length(I_term_divider)+1
    fprintf('Term ctr is %d out of %d',term_ctr,length(I_term_divider)+1);
    integrand_term_math = ...
      integrand_math(I_terms(term_ctr)+1:I_terms(term_ctr+1)-1)
    
    if ~isempty(integrand_term_math) % Protect for leading negative sign
      
      % Integrate in Mathematica
      if strcmp(simplification,'none')
        leadingChar = 'InputForm[Integrate[';
        endingChar = ']]';
      elseif strcmp(simplification,'partial')
        leadingChar = 'InputForm[Simplify[Integrate[';
        endingChar = ']]]';
      elseif strcmp(simplification,'full')
        leadingChar = 'InputForm[FullSimplify[Integrate[';
        endingChar = ']]]';
      else
        error('Incorrect simplification input.');
      end
      integrand_term_int_math = math([leadingChar, ...
        integrand_term_math,',',int_var,',Assumptions -> {',assumptions,'}',endingChar]);
      
      % Add parenthesis around result in case multiple terms result from
      % integration (protects if minus sign leads the multiple terms resulting
      % from the integration)
      integrand_term_int_math = ['(',integrand_term_int_math,')']
      
      integrand_int_math = [integrand_int_math,integrand_term_int_math];
      
    end
    
    if term_ctr ~= length(I_term_divider)+1
      integrand_int_math = [integrand_int_math, ...
        integrand_math(I_term_divider(term_ctr))];
    end
    
  end
  
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
  
