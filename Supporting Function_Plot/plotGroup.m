function plotGroup(t, data, names, figTitle)
    figure; 
    n = size(data,2);
    for i = 1:n
        subplot(n,1,i);
        plot(t, data(:,i), 'LineWidth', 1.5);
        ylabel(names{i}, 'Interpreter', 'latex', 'FontSize', 14);
        grid on;
    end
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    sgtitle(figTitle, 'Interpreter', 'latex', 'FontSize', 14);
end
