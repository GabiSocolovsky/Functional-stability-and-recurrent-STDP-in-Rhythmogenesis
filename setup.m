function setup()
    projectRoot = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(projectRoot,'src')));
    addpath(genpath(fullfile(projectRoot,'scripts')));
    % data/ typically doesn’t need to be on path, but you can add it if needed:
    % addpath(genpath(fullfile(projectRoot,'data')));
    fprintf('Paths added. You can now run scripts in scripts/.\n');
end