function paper_save(h, out_dir, name)
% Save a paper panel as vector PDF (for Illustrator) and 300 dpi PNG.
% The page equals the figure size, so the axes keep their cm positions
% (no tight cropping): panels of one kind align when placed at the origin.

if ~isfolder(out_dir); mkdir(out_dir); end
pos = get(h, 'Position'); % cm
set(h, 'PaperUnits', 'centimeters', 'PaperSize', pos(3:4), 'PaperPosition', [0 0 pos(3:4)]);
print(h, fullfile(out_dir, [name '.pdf']), '-dpdf', '-painters');
print(h, fullfile(out_dir, [name '.png']), '-dpng', '-r300');
fprintf('  saved %s\n', fullfile(out_dir, name));
end
