function run_all()
%RUN_ALL  Run every script headless and regenerate results/figures.
%   matlab -batch "addpath('scripts'); run_all"
%   (A function: the scripts start with "clear" and are run in the base
%   workspace so they cannot wipe this loop's variables.)
here = fileparts(mfilename('fullpath'));
files = [dir(fullfile(here, 's0*.m')); dir(fullfile(here, 'demo_*.m'))];
names = erase({files.name}, '.m');                  % every lab script present, in order
for k = 1:numel(names)
    fprintf('\n##### %s\n', names{k});
    t0 = tic;
    evalin('base', sprintf('run(''%s'');', fullfile(here, [names{k} '.m'])));
    close all force;
    fprintf('##### %s done in %.1f s\n', names{k}, toc(t0));
end
end
