function output = hangGliderEndpoint(input)

xf = input.phase(1).finalstate(1);

output.objective = -xf;