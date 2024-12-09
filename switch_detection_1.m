% Switching Point Detection in Multiple Trajectories
% Author: Kaito Iwasaki
% Date: 12/1/2024 (last updated)

clear; clc; close all;

%% 1. Define System Parameters
% Precompute constants for efficiency
exp_pi_over_2 = exp(pi/2);

% define system matrices A1 and A2 for switching modes dx/dt = A_1x on xy > 0
% and dx/dt = A_2x on xy < 0
A1 = [-1, 1/exp_pi_over_2;
      -exp_pi_over_2, -1];

A2 = [-1, exp_pi_over_2;
      -1/exp_pi_over_2, -1];

%% 2. Define Initial Conditions
% Define N (large, can change) different initial conditions spread across various quadrants
initialConditions = [
    1,    0.5;
    -1,   -0.5;
    0.5,  1;
    -0.5, -1;
    1,   -1;
    -1,    1;
    2,    1;
    -2,   -1;
    1.5, -1.5;
    -1.5, 1.5;
    0.8,  0.8;
    -0.8, -0.8;
    1,   -0.2;
    -1,    0.2;
    0.5,  1.2   
    -0.5, -1.2;
    1.2, -0.8;
    -1.2, 0.8;
    2.5,  1.5;
    -2.5, -1.5  
    0.3,  0.3;
    -0.3, -0.3;
    0.4,  0.5;
    -0.4, 0.5;
    0.2, -0.4; 
    0.1, -0.2 
    -0.1, 0.2 
];

numTrajectories = size(initialConditions, 1); % # of trajectories (= N) (just to keep record)

%% 3. Define Time Span for Simulation
tStart = 0;
tEnd = 20;  % Total time

%% 4. Initialize Storage for Trajectories
% preallocate cell arrays for efficiency
trajectories = cell(numTrajectories, 1);

%% 5. Generate Trajectories
for traj = 1:numTrajectories
    fprintf('Simulating Trajectory %d of %d...\n', traj, numTrajectories);
    
    % intialize
    x0 = initialConditions(traj, 1);
    y0 = initialConditions(traj, 2);
    initialCondition = [x0; y0];
    
    % Define the ODE function handle
    odeFunc = @(t, x) switchingODE(x, A1, A2);
    
    % Simulate the traj
    [t, x] = ode45(odeFunc, [tStart, tEnd], initialCondition, odeset('RelTol',1e-8, 'AbsTol',1e-10));
    
    %store traj data
    trajectories{traj} = [t, x];
end

%% 6. Detect Switching Points Using Difference Quotients
% define threshold multiplier (e.g., 3 times the median)
thresholdMultiplier = 5.5; % this method is sensitive to the choice of this threshold multiplier

% initialize cell array to store detected switching points for all trajectories
detectedSwitchingPointsAll = cell(numTrajectories, 1);

for traj = 1:numTrajectories
    trajData = trajectories{traj};
    t = trajData(:,1);  % time vector
    x = trajData(:,2);  % x(t) vector
    y = trajData(:,3);  % y(t) vector
    
    % compute difference dx = x(t_i+1) - x(t_i) for each time step dt = t_i+1 - t_i
    dt = diff(t);
    dx = diff(x);
    dy = diff(y);
    
    % prevent division by zero by replacing zero dt with a very small number
    dt(dt == 0) = eps;
    
    % calculate difference quotients to approximate dx/dt
    dq_x = abs(dx ./ dt);
    dq_y = abs(dy ./ dt);
    
    % determine dynamic thresholds based on median difference quotients
    % this part may be improved better
    median_dq_x = median(dq_x);
    median_dq_y = median(dq_y);
    threshold_x = thresholdMultiplier * median_dq_x;
    threshold_y = thresholdMultiplier * median_dq_y;
    
    % identify indices where difference quotients exceed thresholds
    switchIdx_x = find(dq_x > threshold_x);
    switchIdx_y = find(dq_y > threshold_y);
    
    % combine any overlaps
    switchIdx = unique([switchIdx_x; switchIdx_y]);
    
    % record detected switching points (using the earlier time point)
    detectedSwitchingPoints = [x(switchIdx), y(switchIdx)];
    
    % store detected switching points for this trajectory
    detectedSwitchingPointsAll{traj} = detectedSwitchingPoints;
end

%% 7. Plot All Trajectories and Highlight Detected Switching Points
figure('Name', 'Trajectories with Detected Switching Points', 'NumberTitle', 'off');
hold on; grid on; axis equal;

colors = lines(numTrajectories);

for traj = 1:numTrajectories
    trajData = trajectories{traj};
    plot(trajData(:,2), trajData(:,3), 'Color', colors(traj,:), 'LineWidth', 1.0);
end

for traj = 1:numTrajectories
    detectedSwitchingPoints = detectedSwitchingPointsAll{traj};
    if ~isempty(detectedSwitchingPoints)
        plot(detectedSwitchingPoints(:,1), detectedSwitchingPoints(:,2), ...
             'o', 'Color', 'k', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    end
end

xlabel('x(t)');
ylabel('y(t)');
title('Trajectories with Detected Switching Points');
hold off;

%% 8. Plot All Detected Switching Points Separately
% (Optional): Visualize all switching points across trajectories
figure('Name', 'All Detected Switching Points', 'NumberTitle', 'off');
hold on; grid on; axis equal;

for traj = 1:numTrajectories
    switchingPoints = detectedSwitchingPointsAll{traj};
    if ~isempty(switchingPoints)
        plot(switchingPoints(:,1), switchingPoints(:,2), '.', ...
             'Color', 'k', 'MarkerSize', 10, 'LineWidth', 1.5);
    end
end

xlabel('x(t)');
ylabel('y(t)');
title('All Detected Switching Points');
hold off;

%% 9. Summary of Detected Switching Points
% Display the number of detected switching points for each trajectory
fprintf('\nSummary of Detected Switching Points:\n');
for traj = 1:numTrajectories
    numSwitches = size(detectedSwitchingPointsAll{traj}, 1);
    fprintf('Trajectory %d: %d switching points detected.\n', traj, numSwitches);
end

%% 10. Moment Matrix 
switchingData = [];

for traj = 1:numTrajectories
    switchingData = cat(1,switchingData,detectedSwitchingPointsAll{traj});
end

[rownum,colnum]=size(switchingData);
N = rownum; % # of sample data points

% Columns: [1, x, y, x^2, x*y, y^2] moment matrix of degree 2 (modify degree as needed)
M = [ones(N,1), switchingData(:,1), switchingData(:,2), ...
         switchingData(:,1).^2, switchingData(:,1).*switchingData(:,2), switchingData(:,2).^2];

%% 11. Perform Singular Value Decomposition
[U, S, V] = svd(M, 'econ');  % 'econ' for economy size decomposition

% Display singular values
singular_values = diag(S);
disp('Singular Values:');
disp(singular_values);

%% 12. Extract nullspace of M
% Define tolerance for considering singular values as zero
tol = 0.5;

% Identify indices of singular values below the tolerance
null_space_indices = find(singular_values < tol);

% Extract corresponding nullspace vectors
null_space_vectors = V(:, null_space_indices);

% Display nullspace vectors
disp('Nullspace Vectors:');
disp(null_space_vectors);

%% 13. Form the polynomial equation coefficients
% Assuming one nullspace vector corresponds to the circle
if size(null_space_vectors, 2) >=1
    null_vec = null_space_vectors(:,1);
    
    % Normalize the nullspace vector
    %null_vec = null_vec / null_vec(end);  % Normalize by y^2 coefficient
    
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

%% 14. Compare with Ground Truth: xy=0

% Plot the recovered circle
% Solve for y in terms of x or x in terms of y, or use implicit plotting

% Create a fine grid for plotting the recovered switching surface
fplot_handle = @(x_val, y_val) coefficients(1) + coefficients(2)*x_val + ...
                                coefficients(3)*y_val + coefficients(4)*x_val.^2 + ...
                                coefficients(5)*x_val.*y_val + coefficients(6)*y_val.^2;

% Generate points where the recovered polynomial is approximately zero using implicit equation solving

% Define a grid
x_range = linspace(min(switchingData(:,1))-1, max(switchingData(:,1))+1, 400);
y_range = linspace(min(switchingData(:,2))-1, max(switchingData(:,2))+1, 400);
[X_grid, Y_grid] = meshgrid(x_range, y_range);
Z = coefficients(1) + coefficients(2)*X_grid + ...
    coefficients(3)*Y_grid + coefficients(4)*X_grid.^2 + ...
    coefficients(5)*X_grid.*Y_grid + coefficients(6)*Y_grid.^2;

%% Plot the data points and recovered surface
figure;
plot(switchingData(:,1), switchingData(:,2), 'ro', 'MarkerSize', 4, 'LineWidth', 1); hold on;
%plot(x, y, 'ro', 'MarkerSize', 4, 'LineWidth', 1); hold on;
contour(X_grid, Y_grid, Z, [0 0], 'b', 'LineWidth', 1);

% Plot ground truth surface
[x, y] = meshgrid(-4:0.1:4);
f = x.*y; 
contour(x, y, f, [36 36], 'r-');
axis equal;
xlabel('x');
ylabel('y');
title('Recovered Switching Surface vs. Ground Truth');
legend('Sample Points', 'Recovered Surface', 'Ground Truth', 'Location', 'BestOutside');
grid on;
hold off;

%% 10. Local Function Definitions

% Define the ODE with switching
function dxdt = switchingODE(x, A1, A2)
    % switching condition to be used later
    product = x(1) * x(2);
    
    if product < 0
        % Use A1 dynamics
        dxdt = A1 * x;
    else
        % Use A2 dynamics
        dxdt = A2 * x;
    end
end