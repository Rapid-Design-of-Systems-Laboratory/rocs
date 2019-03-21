function [combos] = eomTermsSTT(m)
%
% This function identifies the combination of terms associated with the state
% transition tensor (STT) first order ordinary differential equations derived 
% from the Taylor Series approximation.
%
% Input:
%   m = maximum order of state transition tensor
%
% Output:
%   combos{p,pJac} = cell array of each term in rate equations associated with 
%     order p and Jacobian term pJac. Note, this is a lower diagonal matrix.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Terms for Taylor Series Approximation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
combos = cell(m,1);

% Loop through each order of Taylor Series approximation
for p = m : -1 : 1
  
  % Loop through each Jacobian term
  for pJac = m : -1 : 1
  
    % Construct matrix of possible Taylor Series combinations
    vecCombos = repmat((1 : 1 : m),1,pJac);
    permSet = unique(nchoosek(vecCombos,pJac),'rows');

    % Determine if sum of terms matches order p
    I = find(sum(permSet,2) == p);

    % Save each Taylor Series combination that matches order p
    for comboCtr = 1 : 1 : length(I)

      % Save Taylor Series combination
      combos{p,pJac}.term{comboCtr,1} = permSet(I(comboCtr),:);

      % Determine multiplier to match coefficient (1/factorial(p)) of term
      multTerms = prod(factorial(permSet(I(comboCtr),:)));
      
      % Determine multiplier that comes from Jacobian term
      multJac = factorial(pJac);
      
      % Determine target for Jacobian term
      multTarget = factorial(p)/multJac;

      % Save multiplier for Taylor Series combination term
      combos{p,pJac}.multiplier{comboCtr,1} = multTarget/multTerms;
      
    end
    
  end
  
end

return

