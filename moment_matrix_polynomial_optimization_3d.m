%% Basic Setup (generate sample points)
% Number of sample points
clear; clc; close all;
N = 100;

% Define the hyperbolic paraboloid 4x^2 - y^2 + 8x - 6y - 4z - 5 = 0
% The standard form is given bt z = (x + 1)^2 - 1/4 * (y - 3)^2

% Generate random points on the surface
x_arr = rand(N, 1) * 2 - 2;
y_arr = rand(N, 1) * 2 + 2;
x_arr
y_arr
z_arr = (x_arr + 1) .^ 2 - 1.0 / 4 * (y_arr - 3) .^ 2;

% Optionally, add slight noise to simulate real-world data
%noise_level = 0.05;  % Adjust noise level as needed
%x_noisy = x + noise_level * randn(size(x));
%y_noisy = y + noise_level * randn(size(y));

% Combine into a matrix of sample points
%data = [x_noisy', y_noisy'];
data = [x_arr, y_arr, z_arr];

figure;
% plot3(x_arr, y_arr, z_arr, 'b-', 'LineWidth', 2); hold on;
plot3(x_arr, y_arr, z_arr, 'ro', 'MarkerSize', 4, 'LineWidth', 1);
%plot(x_noisy, y_noisy, 'ro', 'MarkerSize', 4, 'LineWidth', 1);
axis equal;
xlim([-2 0]);
ylim([2 4]);
xlabel('x');
ylabel('y');
zlabel('z');
title('Sample Points on the Hyperbolic Paraboloid');
legend('Ground Truth Circle', 'Sample Points');
grid on;
hold off;

%% Define the monomial basis up to degree 3
% Columns: [x, y, z, x^2, y^2, z^2, x*y, x*z, y*z] moment matrix
M = [data(:,1), data(:,2), data(:,3) ...
         data(:,1).^2, data(:,2).^2, data(:,3).^2 ...
         data(:,1).*data(:,2), data(:,1).*data(:,3), data(:,2).*data(:,3)];
% want to find nullspace of moment matrix M v \in Null(M) Mv = 0

%% Solve the least square system
coeffs = lsqr(M, ones(N, 1), 1e-5, 20);     % canonical choice of constant term to be 1
  
% Display the polynomial equation
fprintf('\nRecovered Polynomial Equation (Degree 2, 3D):\n');
fprintf('1 + %.4f*x + %.4f*y + %.4f*z + %.4f*x^2 + %.4f*y^2 + %.4f*z^2 + %.4f*xy + %.4f*yz + %.4f*xz = 0\n', ...
    coeffs(1), coeffs(2), coeffs(3), ...
    coeffs(4), coeffs(5), coeffs(6), ...
    coeffs(7), coeffs(8), coeffs(9));

%% Compare with Ground Truth:

% Plot the recovered circle
% Solve for y in terms of x or x in terms of y, or use implicit plotting

% Create a fine grid for plotting the recovered circle
fplot_handle = @(x_val, y_val, z_val) ...
    1 + coeffs(1)*x_val + coeffs(2)*y_val + coeffs(3)*z_val + ...
    coeffs(4)*x_val.^2 + coeffs(5)*y_val.^2 + coeffs(6)*z_val.^2 + ...
    coeffs(7)*x_val.*y_val + coeffs(8)*x_val.*z_val + coeffs(9)*y_val.*z_val;

% Generate points where the recovered polynomial is approximately zero
% Using implicit equation solving

% Define a grid
x_range = linspace(min(data(:,1))-1, max(data(:,1))+1, 400);
y_range = linspace(min(data(:,2))-1, max(data(:,2))+1, 400);
z_range = linspace(min(data(:,3))-1, max(data(:,3))+1, 400);
[X_grid, Y_grid, Z_grid] = meshgrid(x_range, y_range, z_range);
grid3d = coeffs(1)*X_grid + coeffs(2)*Y_grid + coeffs(3)*Z_grid + ...
    coeffs(4)*X_grid.^2 + coeffs(5)*Y_grid.^2 + coeffs(6)*Z_grid.^2 + ...
    coeffs(7)*X_grid.*Y_grid + coeffs(8)*X_grid.*Z_grid + coeffs(9)*Y_grid.*Z_grid;

%% Plot the data points and recovered circle
figure;
% plot3(x_arr, y_arr, z_arr, 'ro', 'MarkerSize', 4, 'LineWidth', 1); hold on;
% contour3(X_grid, Y_grid, Z_grid, grid3d, [0 0], 'b-', 'LineWidth', 2);

% Plot ground truth
[x_gt, y_gt] = meshgrid(-2:0.1:0, 2:0.1:4);
z_gt = (x_gt + 1) .^ 2 - 1.0 / 4 * (y_gt - 3) .^ 2;
mesh3d = mesh(x_gt, y_gt, z_gt, "FaceAlpha", 0.5, "FaceColor", "#78290f", "EdgeColor", "none")

% Plot least mean square solution
lms_solution = @(x, y, z) -1 + ...
    coeffs(1) * x + coeffs(2) * y + coeffs(3) * z + ...
    coeffs(4) * x.^2 + coeffs(5) * y.^2 + coeffs(6) * z.^2 + ...
    coeffs(7) * x.*y + coeffs(8) * z.*y + coeffs(9) * x.*z;

domain = [-2 0 2 4 0 5]; 
hold on;
fimplicit3(lms_solution, domain, "EdgeColor", "none", "FaceAlpha", 0.5, "FaceColor", "#219ebc");

axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('Recovered Surface vs. Ground Truth Surface');
% legend('Noisy Sample Points', 'Recovered Circle', 'Ground Truth Circle', 'Location', 'BestOutside');
grid on;
hold off;