function [val] = ncr(top,bottom)

val = factorial(top)/(factorial(bottom)*factorial(top-bottom));

return

