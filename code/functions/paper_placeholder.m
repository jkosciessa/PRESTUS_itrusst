function paper_placeholder(out_dir, name, fig_size, axes_pos, y_limits, y_label, x_labels, note)
% Empty panel with the axes box and y-range of its neighbouring panel, as
% placeholder for panels provided by collaborators (e.g. BabelBrain).

h = paper_figure(fig_size);
ax = axes(h, 'Units', 'centimeters', 'Position', axes_pos);
n = numel(x_labels);
set(ax, 'XLim', [0.4 max(n, 1)+0.6], 'YLim', y_limits, 'XTick', 1:n, ...
    'XTickLabel', x_labels, 'TickLabelInterpreter', 'tex', 'Box', 'on', 'XTickLabelRotation', 0, ...
    'XColor', [0.6 0.6 0.6], 'YColor', [0.6 0.6 0.6]);
ylabel(ax, y_label, 'Interpreter', 'tex');
text(ax, mean(xlim(ax)), mean(ylim(ax)), note, 'HorizontalAlignment', 'center', ...
    'Color', [0.6 0.6 0.6], 'Interpreter', 'none');
paper_save(h, out_dir, name);
close(h);
end
