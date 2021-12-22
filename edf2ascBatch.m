function edf2ascBatch(folder)

saveDf = cd(folder);
fileList = cellstr(ls('*exp_*.edf'));

for i = 1:numel(fileList)
    
    [stat,~] = system(sprintf('edf2asc %s',fileList{i}));
    
%     if stat
%         warning('Can not convert file %s ',fileList{i});
%     end 

end

cd(saveDf);

end