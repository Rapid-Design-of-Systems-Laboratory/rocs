function [str,var,reduceRatio] = equationReduction(strReduce,str)

% Remove whitespace
strTemp = char(strReduce);
I = strfind(strTemp,' ');
strTemp(I) = '';
var = strTemp;
sizeVar = length(var);

% Find string
sumReduce = 0;
for ctr = 1 : 1 : length(str)
  
  % Determine first element of each string
  K = strfind(var,str{ctr});
%   length(K)
  
  % Determine how much reduce string by
  sumReduce = sumReduce + length(str{ctr})*length(K) - length(K)*(1+length(int2str(ctr)));
  
  % Perform substitution
  while ~isempty(K)

    % Convert string of atan - find index of start
    ind = K(1);

    % Split string for inserting since strings not same length and insert ArcSin
    var = [var(1:ind-1),'q',int2str(ctr),var(ind+length(str{ctr}):end)];

    % Need to reset atan string index and open and close parenthesis positions
    K = strfind(var,str{ctr});

  end
  
end

reduceRatio = sumReduce/sizeVar;

return