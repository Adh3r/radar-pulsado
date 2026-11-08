function file = save_figure(fig, name)
%SAVE_FIGURE  Export a figure to results/figures/<name>.png (light theme, 150 dpi).
%   file = save_figure(fig, name)
outdir = fullfile(repo_root(), 'results', 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end
try theme(fig, 'light'); catch, end          % R2025a+ figure themes
set(fig, 'Color', 'w');
file = fullfile(outdir, [name '.png']);
exportgraphics(fig, file, 'Resolution', 150, 'BackgroundColor', 'white');
end
