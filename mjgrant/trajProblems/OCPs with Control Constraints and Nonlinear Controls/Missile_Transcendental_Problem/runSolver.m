if ~exist('./autocode','dir')
	mkdir('autocode');
end
if ~exist('./data','dir')
	mkdir('data');
end

% runCombinedProcess(@new_thrust_aoaapprox);
runCombinedProcess(@new_thrust_aoadot);