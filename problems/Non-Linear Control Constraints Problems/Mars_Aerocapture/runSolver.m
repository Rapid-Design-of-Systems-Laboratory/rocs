if ~exist('./autocode','dir')
	mkdir('autocode');
end
if ~exist('./data','dir')
	mkdir('data');
end

runCombinedProcess(@Aerocapture_bounded);
% runCombinedProcess(@Aerocapture_unbounded);