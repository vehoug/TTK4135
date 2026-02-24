function plot_mpc_results(y_t, u_t, N, title_str)
    t_x = 0:N;
    t_u = 0:N-1;
    
    figure('Color', 'w', 'Name', title_str);
    
    subplot(2, 1, 1);
    plot(t_x, y_t, '-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.2);
    grid on;
    ylabel('$y_t$');
    title([title_str, ': State Trajectory']);
    
    subplot(2, 1, 2);
    stairs(t_u, u_t', '-o', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
    grid on;
    xlabel('$k$');
    ylabel('$u_t$');
    title([title_str, ': Control Inputs']);
    
    linkaxes(get(gcf, 'Children'), 'x');
    xlim([0 N]);
end
