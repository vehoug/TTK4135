%% Quadratic Programming Example from Nocedal & Wright (2006)

G = [ 2  0;
      0  2;];

c = [-2 -5]';

A = [ 1 -2;
     -1 -2;
     -1  2;
      1  0;
      0  1;];

b = [-2 -6 -2  0  0]';

x0 = [ 2  0]';
W0 = [ 0  0  1  0  0]';

[x,fval,iterates,W,lam,foo] = qp_163(G, c, A, b, x0, W0, 'report')

%% Visualizing Objective Function, Contraints, Iterations and Solution

x1_l = -1; x1_h = 5;
x2_l = -1; x2_h = 4;
res = 0.05;
[X1, X2] = meshgrid(x1_l:res:x1_h, x2_l:res:x2_h);

Z = X1.^2 + X2.^2 - 2*X1 - 5*X2;

levels = -10:1:10; 

figure('Color', 'w', 'Position', [100, 100, 800, 600]); 
hold on; grid on; axis equal;

[C, h] = contour(X1, X2, Z, levels, 'Color', 0.7*[1 1 1], ...
    'DisplayName', '$f(x) = x_1^2 + x_2^2 - 2x_1 - 5x_2$');
clabel(C, h, 'FontSize', 8, 'Color', 'k');

xline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
yline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');

x_vals = linspace(x1_l, x1_h, 100);

line1 = 0.5*x_vals + 1;
plot(x_vals, line1, 'r', 'LineWidth', 1.5, 'DisplayName', '$x_1 - 2x_2 \geq -2$');

line2 = -0.5*x_vals + 3;
plot(x_vals, line2, 'b', 'LineWidth', 1.5, 'DisplayName', '$-x_1 - 2x_2 \geq -6$');

line3 = 0.5*x_vals - 1;
plot(x_vals, line3, 'm', 'LineWidth', 1.5, 'DisplayName', '$-x_1 + 2x_2 \geq -2$');

vx = [0, 2, 4, 2, 0];
vy = [0, 0, 1, 2, 1];

patch(vx, vy, 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
    'DisplayName', 'Feasible Region');

plot(iterates(1,:), iterates(2,:), 'k--o', ...
    'LineWidth', 2, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'y', ...
    'DisplayName', 'Active Set Path');

for i = 1:size(iterates, 2)
    text(iterates(1,i) + 0.1, iterates(2,i) + 0.15, ...
        sprintf('k=%d', i-1), ...
        'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'bold');
end

plot(x(1), x(2), 'kp', 'MarkerSize', 16, ...
    'MarkerFaceColor', 'c', 'DisplayName', 'Optimal Solution $x^*$');

xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 14);
title('QP Active Set Method (Nocedal \& Wright Ex 16.3)', ...
    'Interpreter', 'latex', 'FontSize', 16);

legend('show', 'Interpreter', 'latex', 'Location', 'northeastoutside');
axis([x1_l x1_h x2_l x2_h]);
hold off;