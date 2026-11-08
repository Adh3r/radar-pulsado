% STARTUP  Add the radar library (src/ and its subfolders) to the MATLAB path.
%   Every script calls it automatically; safe to call repeatedly.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));
