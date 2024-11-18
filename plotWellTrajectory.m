function [M, I, A, N, E, T, MIANETTable] = plotWellTrajectory(surfaceLocation, payZonePositions, payZoneShapes, kopDepth, interval, varargin)
    % Check if 'singleColor' is specified in the optional arguments
    singleColor = any(strcmp(varargin, 'singleColor'));

    % Check if 'plot2D' is specified in the optional arguments
    plot2D = any(strcmp(varargin, 'plot2D'));

    % Check if 'saveExcel' is specified in the optional arguments
    saveExcel = any(strcmp(varargin, 'saveExcel'));

    % Check if 'showLabels' is specified in the optional arguments
    showLabels = any(strcmp(varargin, 'showLabels'));

    % Initialize variables for the well path
    [M, I, A, N, E, T] = initializePath(surfaceLocation);

    % Move to KOP and calculate path, storing the values separately for plotting
    [M_kop, I_kop, A_kop, N_kop, E_kop, T_kop] = moveToKOPAndCalculatePath(M, I, A, N, E, T, kopDepth, interval, payZonePositions);

    % Setup 3D plot
    setup3DPlot();
    
    % Start with a color index
    colorIndex = 1;

    % Plot the initial segment from Surface Location to KOP
    colorIndex = plotSegmentWithColor(N_kop, E_kop, T_kop, colorIndex, singleColor);  % Initial segment from surface to KOP
    
    % Set the path variables from KOP for further segments
    M = M_kop;
    I = I_kop;
    A = A_kop;
    N = N_kop;
    E = E_kop;
    T = T_kop;
    
      if license('test', 'Optimization_Toolbox')
        % Run optimization if the toolbox is available
        payZoneOrder = optimizePayZoneOrder(payZonePositions);
    else
        % Use the given pay zone order and log a message
        warning('Optimization Toolbox not available. Proceeding with default pay zone order.');
        payZoneOrder = 1:length(payZonePositions); % Default order
      end

    % Loop through each pay zone and calculate the direction transitions
    for idx = 1:length(payZoneOrder)
        i = payZoneOrder(idx);
        % Get target positions for the current pay zone
        [N_target, E_target, T_target] = getPayZoneTarget(payZonePositions{i});
    
         % Calculate transition path to next pay zone
        [M_values, N_values, E_values, T_values] = calculateCurvedPath(N(end), E(end), T(end), N_target, E_target, T_target, interval);
    
        % Adjust M_values to continue from the previous measured depth value
        M_values = M(end) + M_values;
        % Calculate inclination and azimuth for the current transition
        [I_end, A_end] = calculateDirection(N, E, T, N_values, E_values, T_values);
    
       % Ensure all new values are column vectors before concatenating
        M = [M; M_values(2:end)']; % Transpose to ensure it's a column vector
        N = [N; N_values(2:end)']; % Transpose to ensure it's a column vector
        E = [E; E_values(2:end)']; % Transpose to ensure it's a column vector
        T = [T; T_values(2:end)']; % Transpose to ensure it's a column vector
        
        % Use the length of any of the appended arrays (since they are all the same length) to create I and A
        I = [I; repmat(I_end, length(N_values(2:end)), 1)];
        A = [A; repmat(A_end, length(N_values(2:end)), 1)];
        
        % Plot the current segment
        colorIndex = plotSegmentWithColor(N_values, E_values, T_values, colorIndex, singleColor);
        
        % Plot the pay zone shape
        plotPayZoneShape(payZoneShapes{i}, payZonePositions{i}, N_target, E_target, T_target);
        % Plot the label if 'showLabels' is specified
        if showLabels
            % Set a reasonable offset distance (e.g., 10 units above the pay zone)
            offset = -250;  % Adjust this value based on your plot scale
            % Plot the label with offset
            validShapes = {'Cylindrical', 'Cuboid', 'Spherical'};
            if any(strcmp(payZoneShapes{i}, validShapes))
                plotPayZoneLabelWithOffset(payZoneShapes{i}, payZonePositions{i}, N_target, E_target, T_target, sprintf('Pay Zone %d', i), offset);
            else
                warning('Unknown pay zone shape: %s', payZoneShapes{i});
            end
        end
    end
    % Call plotWellPathWithCoating to add the coating
    coatingRadius = 5;  % Set the radius of the coating
    

    % Create the MIANET Table
    MIANETTable = createMIANETTable(M, I, A);

     % Check if the user wants to save the table to Excel
    if saveExcel
        saveMIANETTableToExcel(MIANETTable);
    end
    % Label and finalize the plot
    xlabel('Easting');
    ylabel('Northing');
    zlabel('True Vertical Depth');
    title('3D Drilling Path through Multiple Pay Zones');
    view(3);
    set(gca, 'ZDir', 'reverse');
    axis tight;

    % Plot the 2D path if 'plot2D' keyword is passed
    if plot2D
        plot2DPath(N, E, payZoneShapes, payZonePositions, singleColor);
    end
end
 
%%HELPER FUNCTIONS

%Function to save MIANET Table to excel
function saveMIANETTableToExcel(MIANETTable, filename)
    % This function saves the given MIANETTable to an Excel file
    %
    % Inputs:
    %   MIANETTable - A table containing Measured Depth, Inclination, Azimuth,
    %                 Northing, Easting, and True Vertical Depth
    %   filename    - Optional string specifying the filename (including .xlsx extension)
    %                 If not provided, a default name "MIANETTable.xlsx" will be used.

    % Check if the filename was provided
    if nargin < 2
        filename = 'MIANETTable.xlsx'; % Default filename
    end

    % Write the table to an Excel file
    try
        writetable(MIANETTable, filename);
        fprintf('The MIANETTable has been successfully saved to "%s".\n', filename);
    catch ME
        fprintf('An error occurred while trying to save the file:\n%s\n', ME.message);
    end
end

% Optimization Function for Pay Zone Order
function optimizedOrder = optimizePayZoneOrder(payZonePositions)
    % Calculate distances between pay zones and solve TSP

    % Number of pay zones
    numPayZones = length(payZonePositions);

    % Calculate distances from the surface to each pay zone
    surfaceLocation = [0, 0, 0]; % Assuming surface location is [0, 0, 0]
    distancesFromSurface = zeros(1, numPayZones);
    for i = 1:numPayZones
        [N_i, E_i, T_i] = getPayZoneTarget(payZonePositions{i});
        distancesFromSurface(i) = sqrt((N_i - surfaceLocation(1))^2 + ...
                                       (E_i - surfaceLocation(2))^2 + ...
                                       (T_i - surfaceLocation(3))^2);
    end

    % Find the first pay zone to visit (closest to the surface)
    [~, closestPayZone] = min(distancesFromSurface);
    optimizedOrder = closestPayZone;

    % Calculate distances between pay zones
    distances = zeros(numPayZones);
    for i = 1:numPayZones
        for j = 1:numPayZones
            if i ~= j
                [N_i, E_i, T_i] = getPayZoneTarget(payZonePositions{i});
                [N_j, E_j, T_j] = getPayZoneTarget(payZonePositions{j});
                distances(i, j) = sqrt((N_j - N_i)^2 + (E_j - E_i)^2 + (T_j - T_i)^2);
            else
                distances(i, j) = Inf; % Set to a large value to avoid looping back to the same zone
            end
        end
    end

    % Solve TSP starting from the closest pay zone to the surface
    currentZone = closestPayZone;
    remainingZones = setdiff(1:numPayZones, currentZone);

    while ~isempty(remainingZones)
        % Find the nearest unvisited pay zone
        [~, idx] = min(distances(currentZone, remainingZones));
        nextZone = remainingZones(idx);
        optimizedOrder = [optimizedOrder, nextZone];
        currentZone = nextZone;
        remainingZones = setdiff(remainingZones, currentZone);
    end
end

%MIANETTable creation function
function MIANETTable = createMIANETTable(M, I, A)
 
    % Ensureing all arrays have the same length
    if ~(length(M) == length(I) && length(M) == length(A))
        error('M, I, and A must have the same length.');
    end
    % Using conversion function to get NET values
    [N, E, T] = convertMIAToNET(M, I, A);
    
    % Create the table
    MIANETTable = table(M, I, A, N, E, T, ...
        'VariableNames', {'MeasuredDepth', 'Inclination', 'Azimuth', 'Northing', 'Easting', 'TrueVerticalDepth'});
end

%Payzonelabel Function
function plotPayZoneLabelWithOffset(shape, positions, N_target, E_target, T_target, labelText, offset)
    % This function plots the label on top of a pay zone shape, with a specified offset.
    % Inputs:
    %   shape - The type of pay zone ('Cylindrical', 'Cuboid', 'Spherical')
    %   positions - The positions defining the pay zone dimensions
    %   N_target, E_target, T_target - Target coordinates of the pay zone
    %   labelText - String containing the name of the pay zone
    %   offset - Distance to add or subtract to/from the label's height (positive values go above)

    % Handling different shapes
    switch shape
        case 'Cylindrical'
            % Calculate the top of the cylinder and add the offset
            height = max(positions(5), positions(6)) + offset;
            text(N_target, E_target, height, labelText, 'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        
        case 'Cuboid'
            % Calculate the highest point of the cuboid and add the offset
            T_max = max(positions(5), positions(6)) + offset;
            text(N_target, E_target, T_max, labelText, 'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        
        case 'Spherical'
            % Calculate the top of the sphere and add the offset
            radius = max(abs(positions(1) - positions(2)) / 2, abs(positions(3) - positions(4)) / 2);
            verticalStretch = 5; % As defined in the plotPayZoneShape function
            height = T_target + radius * verticalStretch + offset;
            text(N_target, E_target, height, labelText, 'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        
        otherwise
            % Handle unknown shapes gracefully
            fprintf('Warning: Unknown pay zone shape: %s. No label was added.\n', shape);
    end
end

%function to convert MIA data to NET values using minimum curvature method
function [N, E, T] = convertMIAToNET(M, I, A)
    % Initialize arrays for Northing, Easting, and True Vertical Depth
    N = zeros(size(M));
    E = zeros(size(M));
    T = zeros(size(M));
    
    % Loop through intervals to calculate NET values
    for i = 2:length(M)
        % Calculate measured depth increment
        deltaM = M(i) - M(i-1);

        % Convert inclinations and azimuths to radians
        I1 = deg2rad(I(i-1));  % Previous inclination
        I2 = deg2rad(I(i));    % Current inclination
        A1 = deg2rad(A(i-1));  % Previous azimuth
        A2 = deg2rad(A(i));    % Current azimuth

        % Calculate the dogleg severity (angle between two vectors)
        dogleg = acos(cos(I2) * cos(I1) + ...
            sin(I2) * sin(I1) * cos(A2 - A1));

        % Calculate the ratio factor (RF)
        if dogleg < 1e-6
            RF = 1;  % Straight-line approximation for small dogleg
        else
            RF = 2 / dogleg * sin(dogleg / 2);
        end

        % Compute directional increments using minimum curvature
        dN = deltaM * RF * (sin(I2) * cos(A2) + sin(I1) * cos(A1)) / 2;
        dE = deltaM * RF * (sin(I2) * sin(A2) + sin(I1) * sin(A1)) / 2;
        dT = deltaM * RF * (cos(I2) + cos(I1)) / 2;

        % Update cumulative values for Northing, Easting, and TVD
        N(i) = N(i-1) + dN;
        E(i) = E(i-1) + dE;
        T(i) = T(i-1) + dT;
    end
end

% function to initailze path for surface location
function [M, I, A, N, E, T] = initializePath(surfaceLocation)
    M = [0];  % Initial measured depth is 0
    I = [0];  % Initial inclination is 0
    A = [0];  % Initial azimuth is 0
    N = [surfaceLocation(1)];  % Start at surface location (Northing)
    E = [surfaceLocation(2)];  % Start at surface location (Easting)
    T = [surfaceLocation(3)];  % Start at surface location (True Vertical Depth)
end

%function to movepath from surface location to KOP
function [M, I, A, N, E, T] = moveToKOPAndCalculatePath(M, I, A, N, E, T, kopDepth, interval, payZonePositions)
    % Move vertically to the Kick-Off Point (KOP)
    T_current = T(end);
    while T_current < kopDepth
        T_current = min(T_current + interval, kopDepth); % Ensure we stop at KOP depth
        N_new = N(end);  % No change in Northing
        E_new = E(end);  % No change in Easting

        % Calculate distance increment
        distance_increment = sqrt((N_new - N(end))^2 + (E_new - E(end))^2 + (T_current - T(end))^2);
        M_new = M(end) + distance_increment;

        % Append the new values
        M = [M; M_new];
        N = [N; N_new];
        E = [E; E_new];
        T = [T; T_current];
        I = [I; 0];  % Inclination is 0 while moving vertically
        A = [A; 0];  % Azimuth is 0 while moving vertically
    end
    
end


%Function to plot 3D path
function setup3DPlot()
    figure;
    hold on;
    grid on;
    xlabel('Easting');
    ylabel('Northing');
    zlabel('True Vertical Depth');
    title('3D Drilling Path through Multiple Pay Zones');
    view(3);
    set(gca, 'ZDir', 'reverse');
    axis tight;
end

%Function to plot different paths with different colors
function colorIndex = plotSegmentWithColor(N_values, E_values, T_values, colorIndex, singleColor)
    % Define colors to cycle through for each segment
    colors = {'b', 'r', 'g'};  % Blue, Red, Green

    % Determine the color to use
    if singleColor
        colorToUse = colors{1};  % Use the first color for all segments
    else
        colorToUse = colors{colorIndex};  % Cycle through colors
    end

    % Plot the current segment in 3D with the selected color
    plot3(N_values, E_values, T_values, '-o', 'MarkerSize', 5, ...
          'MarkerFaceColor', colorToUse, 'Color', colorToUse, 'LineWidth', 1.5);

    % Add the coating around the current segment with the same color
    coatingRadius = 5;  % Adjust the coating radius as needed
    plotWellPathWithCoating(N_values, E_values, T_values, coatingRadius, colorToUse);

    % Update the color index to cycle to the next color (if not using single color)
    if ~singleColor
        colorIndex = mod(colorIndex, length(colors)) + 1;
    end
end

%Function to plot the well path with a cynlinder like coating
function plotWellPathWithCoating(N_values, E_values, T_values, coatingRadius, color)
    % This function plots the well path with a cylindrical coating around it
    % Inputs:
    % - N_values, E_values, T_values: Coordinates of the well path
    % - coatingRadius: Radius of the coating around the well path
    % - color: Color of the coating (same as the path segment color)

    % Number of points along the well path
    numPoints = length(N_values);

    % Generate circle coordinates around each path point
    numCirclePoints = 20;  % Number of points around each circle for smoothness
    theta = linspace(0, 2 * pi, numCirclePoints);

    % Circle coordinates in local reference frame (around each path point)
    x_circle = coatingRadius * cos(theta);
    y_circle = coatingRadius * sin(theta);

    % Loop through each segment of the well path to construct the coating
    for i = 1:numPoints-1
        % Local frame around the current path segment
        N_center = [N_values(i), N_values(i+1)];
        E_center = [E_values(i), E_values(i+1)];
        T_center = [T_values(i), T_values(i+1)];

        % Calculate the direction vector of the segment
        dirVec = [N_center(2) - N_center(1), E_center(2) - E_center(1), T_center(2) - T_center(1)];
        dirVec = dirVec / norm(dirVec);  % Normalize the direction vector

        % Define two perpendicular vectors to create a local coordinate frame
        if all(dirVec == [0, 0, 1])
            perpVec1 = [1, 0, 0];
        else
            perpVec1 = cross(dirVec, [0, 0, 1]);
            perpVec1 = perpVec1 / norm(perpVec1);
        end
        perpVec2 = cross(dirVec, perpVec1);

        % Compute the circle coordinates around each point of the segment
        N_tube = N_center' + perpVec1(1) * x_circle + perpVec2(1) * y_circle;
        E_tube = E_center' + perpVec1(2) * x_circle + perpVec2(2) * y_circle;
        T_tube = T_center' + perpVec1(3) * x_circle + perpVec2(3) * y_circle;

        % Plot the surface between points with the specified color
        surf(N_tube, E_tube, T_tube, 'FaceColor', color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
end

% Function to plot 2D path
function plot2DPath(N, E, payZoneShapes, payZonePositions,singleColor)
    figure;

    % Plot the well path using Northing and Easting
    plot(E, N, '-o', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
    hold on;

    % Loop through each pay zone and plot it on the 2D plane
    for i = 1:length(payZonePositions)
        positions = payZonePositions{i};
        shape = payZoneShapes{i};

        % Calculate center coordinates in 2D (Northing and Easting)
        N_center = mean([positions(1), positions(2)]);
        E_center = mean([positions(3), positions(4)]);

        if strcmp(shape, 'Cylindrical')
            % Plot a circle to represent a cylindrical zone
            radius = abs(positions(1) - positions(2)) / 2;
            theta = linspace(0, 2 * pi, 100);
            x = E_center + radius * cos(theta);
            y = N_center + radius * sin(theta);
            plot(x, y, 'r--', 'LineWidth', 1.5);  % Dashed red circle

        elseif strcmp(shape, 'Cuboid')
            % Ensure that E_min/E_max and N_min/N_max are properly ordered
            N_min = min(positions(2), positions(1));
            N_max = max(positions(2), positions(1));
            E_min = min(positions(4), positions(3));
            E_max = max(positions(4), positions(3));

            % Plot a rectangle to represent a cuboid zone
            rectangle('Position', [E_min, N_min, E_max - E_min, N_max - N_min], ...
                      'EdgeColor', 'r', 'LineStyle', '--', 'LineWidth', 1.5);  % Dashed red rectangle

        elseif strcmp(shape, 'Spherical')
            % Plot a circle to represent the projected area of the spherical zone
            radius = max(abs(positions(1) - positions(2)) / 2, abs(positions(3) - positions(4)) / 2);
            theta = linspace(0, 2 * pi, 100);
            x = E_center + radius * cos(theta);
            y = N_center + radius * sin(theta);
            plot(x, y, 'b--', 'LineWidth', 1.5);  % Dashed blue circle
        end
    end

    % Set up the plot labels and title
    xlabel('Easting');
    ylabel('Northing');
    title('2D Well Path (Northing vs. Easting) with Target Zones');
    grid on;
    axis equal;  % Set equal scaling for both axes for accurate trajectory representation
    axis tight;  % Ensures the plot limits fit the data
    legend('Well Path', 'Target Zones');
    hold off;
end

function [N_target, E_target, T_target] = getPayZoneTarget(payZonePosition)
    N_target = mean([payZonePosition(1), payZonePosition(2)]);
    E_target = mean([payZonePosition(3), payZonePosition(4)]);
    T_target = mean([payZonePosition(5), payZonePosition(6)]);
end

%Function to calculate path
function [M_values, N_values, E_values, T_values] = calculateCurvedPath(N_current, E_current, T_current, N_target, E_target, T_target, interval)
    % Define intermediate control points for smooth curves
    N_mid1 = (3 * N_current + N_target) / 3.8;
    E_mid1 = (3 * E_current + E_target) / 3.8;
    T_mid1 = (2 * T_current + T_target) / 2.8;

    N_mid2 = (N_current + 3 * N_target) / 3.8;
    E_mid2 = (E_current + 3 * E_target) / 3.8;
    T_mid2 = (T_current + 2 * T_target) / 2.8;

    % Define control points for cubic spline interpolation
    control_points_N = [N_current, N_mid1, N_mid2, N_target];
    control_points_E = [E_current, E_mid1, E_mid2, E_target];
    control_points_T = [T_current, T_mid1, T_mid2, T_target];

    % Generate a smooth path using cubic spline interpolation
    t = linspace(0, 1, interval);  % Dynamic interval ensures smooth transitions
    N_values = spline([0, 0.33, 0.66, 1], control_points_N, t);
    E_values = spline([0, 0.33, 0.66, 1], control_points_E, t);
    T_values = spline([0, 0.33, 0.66, 1], control_points_T, t);

    % Apply smoothing to the path points
    N_values = smoothdata(N_values, 'sgolay', 5);  % Smooth Northing
    E_values = smoothdata(E_values, 'sgolay', 5);  % Smooth Easting
    T_values = smoothdata(T_values, 'sgolay', 5);  % Smooth True Vertical Depth

    % Initialize measured depth (M) values array
    M_values = zeros(size(N_values));
    M_values(1) = 0; % Starting point, relative distance is zero

    % Calculate cumulative measured depth along the path
    for i = 2:length(N_values)
        % Calculate the distance between consecutive points in 3D space
        dN = N_values(i) - N_values(i - 1);
        dE = E_values(i) - E_values(i - 1);
        dT = T_values(i) - T_values(i - 1);
        distance_increment = sqrt(dN^2 + dE^2 + dT^2);

        % Update measured depth
        M_values(i) = M_values(i - 1) + distance_increment;
        
    end
end

function [I_end, A_end, directionType] = calculateDirection(N, E, T, N_values, E_values, T_values)
    % Initial point in transition
    N_current = N_values(1);
    E_current = E_values(1);
    T_current = T_values(1);

    I_start = atan2d(sqrt((N_current - N(end))^2 + (E_current - E(end))^2), T_current - T(end));
    A_start = atan2d(E_current - E(end), N_current - N(end));
    if A_start < 0
        A_start = A_start + 360;
    end

    % Final point in transition
    N_end = N_values(end);
    E_end = E_values(end);
    T_end = T_values(end);

    I_end = atan2d(sqrt((N_end - N_current)^2 + (E_end - E_current)^2), T_end - T_current);
    A_end = atan2d(E_end - E_current, N_end - N_current);
    if A_end < 0
        A_end = A_end + 360;
    end

    % Determine direction type
    if I_end > I_start && A_end == A_start
        directionType = "Build Only";
    elseif I_end < I_start && A_end == A_start
        directionType = "Drop Only";
    elseif I_end == I_start && A_end ~= A_start
        directionType = "Turn Only";
    elseif I_end > I_start && A_end ~= A_start
        directionType = "Build and Turn";
    elseif I_end < I_start && A_end ~= A_start
        directionType = "Drop and Turn";
    else
        directionType = "No Build, No Drop, No Turn";
    end
end

function plotPayZoneShape(shape, positions, N_target, E_target, T_target,i)
    if strcmp(shape, 'Cylindrical')
        radius = abs(positions(1) - positions(2)) / 2;
        height = abs(positions(5) - positions(6));
        [X, Y, Z] = cylinder(radius, 50);
        Z = Z * height + min(positions(5), positions(6));
        X = X + N_target;
        Y = Y + E_target;
        surf(X, Y, Z, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'k');

    elseif strcmp(shape, 'Cuboid')
        N_min = positions(2); N_max = positions(1);
        E_min = positions(4); E_max = positions(3);
        T_min = positions(6); T_max = positions(5);
        vertices = [N_min, E_min, T_min; N_max, E_min, T_min; N_max, E_max, T_min; N_min, E_max, T_min; 
                    N_min, E_min, T_max; N_max, E_min, T_max; N_max, E_max, T_max; N_min, E_max, T_max];
        faces = [1, 2, 3, 4; 5, 6, 7, 8; 1, 5, 8, 4; 2, 6, 7, 3; 1, 2, 6, 5; 4, 3, 7, 8];
        patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
        
    elseif strcmp(shape, 'Spherical')
        radius = max(abs(positions(1) - positions(2)) / 2, abs(positions(3) - positions(4)) / 2);
        verticalStretch = 5;
        [X, Y, Z] = sphere(50);
        X = X * radius + N_target;
        Y = Y * radius + E_target;
        Z = Z * radius * verticalStretch + T_target;
        surf(X, Y, Z, 'FaceColor', 'blue', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
    end
     
end

%%For more information contact the Software developer for Group04 :Damisola
%%Deboh-ajiga:08136050347