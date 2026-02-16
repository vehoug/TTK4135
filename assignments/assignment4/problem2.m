%% Constraint and Contour Visualization

x1_l = -1.5; x1_h = 6.5;
x2_l = -1.5; x2_h = 6.5;
res = 0.01;
[ x1 , x2 ] = meshgrid ( x1_l : res : x2_h , x1_l : res : x2_h );

f = -(3 -0.4* x1 ).* x1 - (2 -0.2* x2 ).* x2 ;

levels = ( -12:2:8)';

figure('Color', 'w', 'Position', [100, 100, 800, 600]); 
hold on; grid on; axis equal;

[C , h ] = contour ( x1 , x2 , f , levels , 'Color', .7*[1 1 1], ...
    'DisplayName', '$q(x) = -(3-0.4x_1)x_1 - (2-0.2x_2)x_2$');
clabel(C, h, 'FontSize', 8, 'Color', 'k');

xline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
yline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');

x_vals = linspace(x1_l, x1_h, 100);

line1 = 8 - 2*x_vals;
plot(x_vals, line1, 'r', 'LineWidth', 2, 'DisplayName', '$2x_1 + x_2 \leq 8$ ($R_I$)');

line2 = (15 - x_vals)/3;
plot(x_vals, line2, 'b', 'LineWidth', 2, 'DisplayName', '$x_1 + 3x_2 \leq 15$ ($R_{II}$)');

vx = [0, 4, 1.8, 0];
vy = [0, 0, 4.4, 5];
patch(vx, vy, 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'Feasible Region');

xlabel('$x_1$ (Product A tonnes)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$ (Product B tonnes)', 'Interpreter', 'latex', 'FontSize', 12);
title('Quadratic Programming: Production Planning', 'Interpreter', 'latex', 'FontSize', 14);
legend('show', 'Interpreter', 'latex', 'Location', 'northeast');
axis([x1_l x1_h x2_l x2_h]);

%% QP Solution

G = [0.8 0; 0 0.4];
c = [ -3 -2]';
A = [ -2 , -1;
-1 , -3;
1 , 0;
0 , 1];
b = [ -8 , -15 ,0 ,0]';
x0 = [0 ,0]';
W0 = [0 0 0 0 0]';
[x , fval , iterates ,W , lam , foo ] = qp_163 (G ,c ,A ,b , x0 , W0 , 'report')

fprintf('Optimal Production:\n');
fprintf('Product A: %.4f tonnes\n', x(1));
fprintf('Product B: %.4f tonnes\n', x(2));

qp_path_x1 = iterates(1, :);
qp_path_x2 = iterates(2, :);

plot(qp_path_x1, qp_path_x2, 'k--d', 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'c', ...
    'DisplayName', 'QP Solver Path');

for i = 1:length(qp_path_x1)
    text(qp_path_x1(i) + 0.15, qp_path_x2(i) + 0.15, ...
        sprintf('iteration %d', i-1), 'Interpreter', 'latex', ...
        'FontSize', 9, 'Color', [0.2 0.2 0.2]);
end

plot(x(1), x(2), 'kp', 'MarkerFaceColor', 'm', 'MarkerSize', 14, ...
    'DisplayName', sprintf('QP Optimum: (%.2f, %.2f)', x(1), x(2)));

legend('show', 'Interpreter', 'latex', 'Location', 'northeast');
hold off;