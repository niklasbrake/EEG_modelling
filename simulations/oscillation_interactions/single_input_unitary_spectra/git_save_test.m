
filename = mfilename('fullpath');
filename = 'E:\Research_Projects\004_Propofol';
gitRepo = fileparts(filename);
parentDirectory = fileparts(gitRepo);
while ~strcmp(gitRepo,parentDirectory);
    parentDirectory = gitRepo;
    gitRepo = fileparts(gitRepo);
    F = dir(gitRepo);
    gitFound = sum(cellfun(@(x)strcmp(x,'.git'),{F(:).name}));
    if(gitFound==1)
        break;
    elseif(gitFound>1)
        error('Multiple git repositories discovered.');
    end
end


if(strcmp(gitRepo,parentDirectory))
    r = 'No accountability found. Code run outside git environment.';
    warning(r)
    infoString = r;
else
    command = sprintf('git -C %s rev-parse HEAD',gitRepo);
    [~,hash] = system(command);
    hash = strip(hash);

    m = length(gitRepo);
    file = filename(m+2:end);

    infoString = strip(sprintf('Code accountable: %s\nGit repository: %s\nLast hash: %s\n',file,gitRepo,hash))
end


