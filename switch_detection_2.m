% Switching Point Detection in Multiple Trajectories
% Author: Kaito Iwasaki
% Date: 12/1/2024 (last updated)

clear; clc; close all;

%% 1. Define Initial Conditions
% Define N (large, can change) different initial conditions spread across various quadrants
initialConditions = [
    -10, 10;
    -10, 9.5;
    -10,  9;
    -10, 8.5;
    -10, 8;
    -10, 7.5;
    -10, 7;
    -10, 6.5;
    -10, 6;
    -10, 5.5;
    -10, 5;
    -10, 4.5;
    -10, 4;
    -10, 3.5;
    -10, 3;    
    -10, 2.5;
    -10, 2;
    -10, 1.5;
    -10, 1;
    -10, 0.5;  
    -10, 0;
    -10, -0.5;
    -10, -1;
    -10, -1.5;
    -10, -2; 
    -10, -2.5;
    -10, -3;
    -10, -3.5;
    -10, -4;
    -10, -4.5;
    -10, -5;
    -10, -5.5;
    -10, -6;
    -10, -6.5;
    -10, -7;
    -10, -7.5;
    -10, -8;
    -10, -8.5;
    -10, -9;
    -10, -9.5;
    -10, -10;
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
    odeFunc = @(t, x) switchingODEonAlgVariety(x);
    
    % Simulate the traj
    [t, x] = ode45(odeFunc, tStart:0.01:tEnd, initialCondition, odeset('RelTol',1e-8, 'AbsTol',1e-10));
    
    %store traj data
    trajectories{traj} = [t, x];
end

%% 6. Detect Switching Points Using Difference Quotients
% define threshold multiplier (e.g., 3 times the median)

% initialize cell array to store detected switching points for all trajectories
detectedSwitchingPointsAll = cell(numTrajectories, 1);
initStepLength = 16;                   % Init number of steps taken into account when testing
subdivTolerance = 0.2;                 % If most derivatives deviates this much from average, subdivide

for traj = 1:numTrajectories
    trajData = trajectories{traj};
    t = trajData(:, 1);  % time vector
    pos = trajData(:, 2:3);

    size(t)

    timeDiff = diff(t);
    posDiff = diff(pos);
    
    fdPosGradient = posDiff ./ timeDiff; % Proxy gradient by finite difference

    curStartStepIndex = 1;

    switches = [];

    while( curStartStepIndex < size(t, 1) )

        curLength = initStepLength;

        while(true)
            
            if (curStartStepIndex + curLength < size(t, 1))
                curTestSteps = fdPosGradient(curStartStepIndex:1:curStartStepIndex + curLength, :);
            else
                curTestSteps = fdPosGradient(curStartStepIndex:1:end, :);
                curLength = size(curTestSteps, 1);
            end

            curDiffs = diff(curTestSteps);   % Finite difference of gradient

            diff_disp = curDiffs';

            % Test for change in norm
            curDiffNorms = (vecnorm(curDiffs'))';
            medianNorm = median(curDiffNorms);
            normUpperBound = medianNorm + subdivTolerance;
            normLowerBound = medianNorm - subdivTolerance;
            normOutliers = find(curDiffNorms > medianNorm + subdivTolerance | curDiffNorms < medianNorm - subdivTolerance);

            norm_disp = curDiffNorms';

            % Test for change in direction
            curDots = zeros(size(curDiffs, 1) - 1, 1);
            for i = 1:(size(curDiffs, 1) - 1)
                curDots(i) = dot(curDiffs(i), curDiffs(i+1));
            end
            curDotsNorms = (vecnorm(curDots'))';
            medianDotNorm = median(curDotsNorms);
            dotUpperBound = medianDotNorm + subdivTolerance * medianNorm;
            dotLowerBound = medianDotNorm - subdivTolerance * medianNorm;
            dotOutliers = find(curDotsNorms > medianDotNorm + subdivTolerance * medianNorm | curDotsNorms < medianDotNorm - subdivTolerance * medianNorm);

            dot_disp = curDots';

            normOutlierCnt = size(normOutliers, 1);
            dotOutlierCnt = size(dotOutliers, 1);

            if (normOutlierCnt >= 2 || dotOutlierCnt >= 2) % Subdivision
                curLength = curLength / 2;
            elseif (normOutlierCnt == 1) % Detected
                idx = normOutliers(1);
                switches(end + 1, :) = pos(curStartStepIndex + idx, :);
                break
            elseif (dotOutlierCnt == 1)  % Detected
                idx = dotOutliers(1);
                switches(end + 1, :) = pos(curStartStepIndex + idx, :);
                break
            else
                break
            end

            % disp('iter');
            % disp('---------------------------------------------')
        
        end

        % disp('next batch');
        % disp('=============================================')

        curStartStepIndex = curStartStepIndex + curLength;
    end

    disp('traj');

    detectedSwitchingPointsAll{traj} = switches;
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

xin = linspace(-10, 10, 500);
yin = linspace(-10, 15, 500);

% Create a meshgrid for evaluating the function
[X, Y] = meshgrid(xin, yin);

% Define the equation
Z = 2*X.^3 + 2*X.^2.*Y + X.*Y.^2 - Y.^3 + X.^2 - 8*X.*Y + 15*Y + 2*X + 20*Y + 8;

% Plot the contour for Z = 0
contour(X, Y, Z, [0 0], 'LineWidth', 2);

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

% Plot the contour for Z = 0
contour(X, Y, Z, [0 0], 'LineWidth', 2);

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

% Columns: [1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3] moment matrix of degree 3 (modify degree as needed)
M = [ones(N,1), switchingData(:,1), switchingData(:,2), ...
         switchingData(:,1).^2, switchingData(:,1).*switchingData(:,2), switchingData(:,2).^2, ...
         switchingData(:,1).^3, (switchingData(:,1).^2).*(switchingData(:,2)), (switchingData(:,1)).*(switchingData(:,2).^2), switchingData(:,2).^3];

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
    % [1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3]
    coefficients = null_vec;
    
    % Display the polynomial equation
    fprintf('\nRecovered Polynomial Equation (Degree 3, 2D):\n');
    fprintf('%.4f + %.4f*x + %.4f*y + %.4f*x^2 + %.4f*x*y + %.4f*y^2 + %.4f*x^3 + %.4f*x^2*y + %.4f*x*y^2 + %.4f*y^3 = 0\n', ...
            coefficients(1), coefficients(2), coefficients(3), ...
            coefficients(4), coefficients(5), coefficients(6), ...
            coefficients(7), coefficients(8), coefficients(9), coefficients(10));
end

%% 14. Compare with Ground Truth: xy=0

% Plot the recovered circle
% Solve for y in terms of x or x in terms of y, or use implicit plotting

% Create a fine grid for plotting the recovered switching surface
fplot_handle = @(x_val, y_val) coefficients(1) + coefficients(2)*x_val + ...
                                coefficients(3)*y_val + coefficients(4)*x_val.^2 + ...
                                coefficients(5)*x_val.*y_val + coefficients(6)*y_val.^2 + ...
                                coefficients(7)*x_val.^3 + coefficients(8)*x_val.^2.*y_val + ...
                                coefficients(9)*x_val.*y_val.^2 + coefficients(10)*y_val.^3;

% Generate points where the recovered polynomial is approximately zero using implicit equation solving

% Define a grid
x_range = linspace(min(switchingData(:,1))-1, max(switchingData(:,1))+1, 400);
y_range = linspace(min(switchingData(:,2))-1, max(switchingData(:,2))+1, 400);
[X_grid, Y_grid] = meshgrid(x_range, y_range);
Z = coefficients(1) + coefficients(2)*X_grid + ...
    coefficients(3)*Y_grid + coefficients(4)*X_grid.^2 + ...
    coefficients(5)*X_grid.*Y_grid + coefficients(6)*Y_grid.^2 + ...
    coefficients(7)*X_grid.^3 + coefficients(8)*X_grid.^2.*Y_grid + ...
    coefficients(9)*X_grid.*Y_grid.^2 + coefficients(10)*Y_grid.^3;

%% Plot the data points and recovered surface
figure;
plot(switchingData(:,1), switchingData(:,2), 'ro', 'MarkerSize', 4, 'LineWidth', 1); hold on;
%plot(x, y, 'ro', 'MarkerSize', 4, 'LineWidth', 1); hold on;
contour(X_grid, Y_grid, Z, [0 0], 'b', 'LineWidth', 1);

% Plot ground truth surface
[x, y] = meshgrid(-10:0.1:10);
f = 2*x.^3 + 2*x.^2.*y + x.*y.^2 - y.^3 + x.^2 - 8*x.*y + 15*y + 2*x + 20*y + 8; ; 
contour(x, y, f, [0 0], 'r-');
axis equal;
xlabel('x');
ylabel('y');
title('Recovered Switching Surface vs. Ground Truth');
legend('Sample Points', 'Recovered Surface', 'Ground Truth', 'Location', 'BestOutside');
grid on;
hold off;


% Define the ODE with switching
function dxdt = switchingODEonAlgVariety(x)

    % switching condition to be used later
    variety = 2*x(1)^3 + 2*x(1)^2*x(2) + x(1)*x(2)^2 - x(2)^3 + x(1)^2 - 8*x(1)*x(2) + 15*x(2) + 2*x(1) + 20*x(2) + 8;
    
    if variety < 0
        % dx/dt = 1, dy/tdt = 0
        dxdt = [1;0];
    else
        % dx/dt = 1, dy/dt = 1
        dxdt = [1;1];
    end
end