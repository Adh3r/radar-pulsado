function root = repo_root()
%REPO_ROOT  Absolute path of the repository root.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % src/util -> root
end
