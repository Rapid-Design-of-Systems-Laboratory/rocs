function y = smoothFcn(z,noise)
% Objective function
if nargin < 2
    noise = 0;
end
LB = [-5 -5];      %Lower bound 
UB = [5 5];        %Upper bound

y = zeros(1,size(z,1));
for i = 1:size(z,1)
    x = z(i,:);
    if any(x<LB) || any(x>UB)
        y(i) = Inf;
    else
    y(i) = x(1)^3 - x(2)^2 + ...
        100*x(2)/(10+x(1)) + noise*randn;
    end
end

