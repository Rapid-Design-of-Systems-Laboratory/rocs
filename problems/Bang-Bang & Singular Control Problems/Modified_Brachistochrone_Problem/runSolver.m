if ~exist('./autocode','dir')
	mkdir('autocode');
end
if ~exist('./data','dir')
	mkdir('data');
end

% runCombinedProcess(@mod_brach_approx);
runCombinedProcess(@mod_brach_trig);