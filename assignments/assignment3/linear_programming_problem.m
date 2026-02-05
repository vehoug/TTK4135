%% (b) Contour Plot with Inequality Constraints

x1_l = -1.5; x1_h = 6.5;
x2_l = -1.5; x2_h = 6.5;
res = 0.1;
[x1, x2] = meshgrid(x1_l:res:x1_h, x2_l:res:x2_h);

f = -3*x1 - 2*x2;
levels = (-32:2:6)';

figure; hold on; grid on;
[C, h] = contour(x1, x2, f, levels, 'Color', .7*[1 1 1], 'DisplayName', '$-3x_1 -2x_2$');
clabel(C, h, 'FontSize', 8);

yline(0, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');

x_vals = linspace(x1_l, x1_h, 100);

line1 = 8 - 2*x_vals;
plot(x_vals, line1, 'r', 'LineWidth', 2, 'DisplayName', '$2x_1 + x_2 \leq 8$');

line2 = (15 - x_vals)/3;
plot(x_vals, line2, 'b', 'LineWidth', 2, 'DisplayName', '$x_1 + 3x_2 \leq 15$');

vx = [0, 4, 1.8, 0];
vy = [0, 0, 4.4, 5];
patch(vx, vy, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Feasible Region');

xlabel('$x_1$ (Product A in 1000 kg)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$ (Product B in 1000 kg)', 'Interpreter', 'latex', 'FontSize', 12);
title('LP Reactor Problem Visualization', 'Interpreter', 'latex', 'FontSize', 14);

legend('show', 'Interpreter', 'latex', 'Location', 'northeast');

axis([x1_l x1_h x2_l x2_h]);
hold off;

%% (c) Simplex Method Solution

clc;

c = [-3 -2  0  0]';
A = [ 2  1  1  0;
      1  3  0  1];

b  = [8 15]';
x0 = [0  0  8  15]';

[x, fval, iterates] = simplex(c, A, b, x0, 'report');

iter_x1_x2 = iterates(1:2, :);

iter_1 = iter_x1_x2(:,1);
iter_2 = iter_x1_x2(:,2);
iter_3 = iter_x1_x2(:,3);

%% (d) Visualize Simplex Iterations and Optimal Solution

hold on; 

plot(iter_x1_x2(1,:), iter_x1_x2(2,:), 'k--o', ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'y', ...
    'DisplayName', 'Simplex Path');

for i = 1:size(iter_x1_x2, 2)
    text(iter_x1_x2(1,i) + 0.2, iter_x1_x2(2,i) + 0.2, ...
        sprintf('Iteration %d', i), ...
        'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'bold');
end

plot(x(1), x(2), 'kp', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'm', 'DisplayName', 'Optimal Solution');

hold off;