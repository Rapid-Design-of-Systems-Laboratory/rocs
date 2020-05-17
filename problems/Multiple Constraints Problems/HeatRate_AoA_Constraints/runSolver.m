if ~exist('./autocode','dir')
	mkdir('autocode');
end
if ~exist('./data','dir')
	mkdir('data');
end

runCombinedProcess(@main_heatrate_aoa);
