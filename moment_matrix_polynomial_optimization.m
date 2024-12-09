%% Basic Setup (generate sample points)
% Number of sample points
N = 100;

% Define the circle equation: x^2 + y^2 + 2x - 3y + 1 = 0
% Convert to standard form to find center and radius
% (x + 1)^2 + (y - 1.5)^2 = (1 + 1^2 + 1.5^2) = 1 + 1 + 2.25 = 4.25
% So, center = (-1, 1.5), radius = sqrt(4.25) ≈ 2.0616

% Circle parameters
h = -1;            % x-coordinate of center
k = 1.5;           % y-coordinate of center
r = sqrt(4.25);    % Radius

% Generate N angles between 0 and 2*pi
theta = linspace(0, 2*pi, N);

% Generate (x, y) points on the circle
x = h + r * cos(theta);
y = k + r * sin(theta);

% Optionally, add slight noise to simulate real-world data
%noise_level = 0.05;  % Adjust noise level as needed
%x_noisy = x + noise_level * randn(size(x));
%y_noisy = y + noise_level * randn(size(y));

% Combine into a matrix of sample points
%data = [x_noisy', y_noisy'];
data = [x', y'];

figure;
plot(x, y, 'b-', 'LineWidth', 2); hold on;
plot(x, y, 'ro', 'MarkerSize', 4, 'LineWidth', 1);
%plot(x_noisy, y_noisy, 'ro', 'MarkerSize', 4, 'LineWidth', 1);
axis equal;
xlabel('x');
ylabel('y');
title('Sample Points on the Circle');
legend('Ground Truth Circle', 'Sample Points');
grid on;
hold off;

%% Define the monomial basis up to degree 3
% Columns: [1, x, y, x^2, x*y, y^2] moment matrix
M = [ones(N,1), data(:,1), data(:,2), ...
         data(:,1).^2, data(:,1).*data(:,2), data(:,2).^2];  % want to find nullspace of moment matrix M v \in Null(M) Mv = 0

%% Perform Singular Value Decomposition
[U, S, V] = svd(M, 'econ');  % 'econ' for economy size decomposition

% Display singular values
singular_values = diag(S);
disp('Singular Values:');
disp(singular_values);

%% Extract nullspace of M
% Define tolerance for considering singular values as zero
tol = 1e-8;

% Identify indices of singular values below the tolerance
null_space_indices = find(singular_values < tol);

% Extract corresponding nullspace vectors
null_space_vectors = V(:, null_space_indices);

% Display nullspace vectors
disp('Nullspace Vectors:');
disp(null_space_vectors);

%% Form the polynomial equation coefficients
% Assuming one nullspace vector corresponds to the circle
if size(null_space_vectors, 2) >=1
    null_vec = null_space_vectors(:,1);
    
    % Normalize the nullspace vector
    null_vec = null_vec / null_vec(end);  % Normalize by y^2 coefficient
    
    % Display the nullspace vector
    disp('First Nullspace Vector (Degree 2):');
    disp(null_vec);
    
    % Form the polynomial equation coefficients
    % [1, x, y, x^2, x*y, y^2]
    coefficients = null_vec;
    
    % Display the polynomial equation
    fprintf('\nRecovered Polynomial Equation (Degree 2):\n');
    fprintf('%.4f + %.4f*x + %.4f*y + %.4f*x^2 + %.4f*x*y + %.4f*y^2 = 0\n', ...
            coefficients(1), coefficients(2), coefficients(3), ...
            coefficients(4), coefficients(5), coefficients(6));
end

%% Compare with Ground Truth: x^2 + y^2 + 2x -3y +1=0

% Plot the recovered circle
% Solve for y in terms of x or x in terms of y, or use implicit plotting

% Create a fine grid for plotting the recovered circle
fplot_handle = @(x_val, y_val) coefficients(1) + coefficients(2)*x_val + ...
                                coefficients(3)*y_val + coefficients(4)*x_val.^2 + ...
                                coefficients(5)*x_val.*y_val + coefficients(6)*y_val.^2;

% Generate points where the recovered polynomial is approximately zero
% Using implicit equation solving

% Define a grid
x_range = linspace(min(data(:,1))-1, max(data(:,1))+1, 400);
y_range = linspace(min(data(:,2))-1, max(data(:,2))+1, 400);
[X_grid, Y_grid] = meshgrid(x_range, y_range);
Z = coefficients(1) + coefficients(2)*X_grid + ...
    coefficients(3)*Y_grid + coefficients(4)*X_grid.^2 + ...
    coefficients(5)*X_grid.*Y_grid + coefficients(6)*Y_grid.^2;

%% Plot the data points and recovered circle
figure;
plot(x, y, 'ro', 'MarkerSize', 4, 'LineWidth', 1); hold on;
contour(X_grid, Y_grid, Z, [0 0], 'b-', 'LineWidth', 2);
% Plot ground truth circle
theta_fine = linspace(0, 2*pi, 400);
x_gt = h + r * cos(theta_fine);
y_gt = k + r * sin(theta_fine);
plot(x_gt, y_gt, 'g--', 'LineWidth', 2);
axis equal;
xlabel('x');
ylabel('y');
title('Recovered Circle vs. Ground Truth');
legend('Noisy Sample Points', 'Recovered Circle', 'Ground Truth Circle', 'Location', 'BestOutside');
grid on;
hold off;