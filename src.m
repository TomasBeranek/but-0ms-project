% File    src.m
% Author  Tomas Beranek <xberan46@stud.fit.vutbr.cz>
% Brief   Visualization of kinematic string
% Date    19.4.2022
% Up2date sources can be found at https://github.com/TomasBeranek/but-0ms-project

file = "manipulator.xlsx";

[num, txt, raw] = xlsread(file);

% load number of segments and number of strings from .xlsx file
lengthsNum = raw{2,3};
instancesNum = raw{3,3};

assert(~isnan(lengthsNum), "ERROR: Incorrect format of number of lengths!");
assert(~isnan(instancesNum), "ERROR: Incorrect format of number of instances!");

% load lengths and angles from .xlsx file
try
    l = cell2mat(raw(6, 3:lengthsNum+2));
catch E
    if (E.identifier == 'MATLAB:cell2mat:MixedDataTypes')
        error("ERROR: Incorrect format of lengths!");
    end
end

try
    a = cell2mat(raw(8:instancesNum+7, 3:lengthsNum+2));
catch E
    if (E.identifier == 'MATLAB:cell2mat:MixedDataTypes')
        error("ERROR: Incorrect format of instances!");
    end
end

% to catch missing values and
assert(~sum(isnan(l)), "ERROR: Incorrect format of lengths!");
assert(~sum(isnan(a), 'all'), "ERROR: Incorrect format of instances!");

% check that length is a positive non-zero number
l(l <= 0) = nan;
assert(~sum(isnan(l)), "ERROR: Atleast one segment has zero length!");

% calculate joint coordinates and display them in a graph
[x,y] = showKinematicString(l,a);

% clear contents
raw(instancesNum+8:end, 1:end) = {nan};

% save results
raw(instancesNum+9,3) = {'polohy počátků ramen'};
raw(instancesNum+10,1) = {'x'};
raw(instancesNum*2+10)  = {'y'};
raw(instancesNum+10:instancesNum*2+9, 2:(lengthsNum+2)) = num2cell(x');
raw(instancesNum*2+10:instancesNum*3+9,2:(lengthsNum+2)) = num2cell(y');

% write everything back the original .xslx file
writecell(raw, file, 'UseExcel', false);

function [x] = calculateCoordinate(l, a, coordinateType)
    stringsTotal = size(a,1);
    positionsTotal = length(l);

    % prepare angles for summming
    sumPrep = flipud(repelem(a', 1, positionsTotal));
    % triangular-like mask for summming
    sumMask = repmat(fliplr(tril(ones(positionsTotal, positionsTotal))), 1, stringsTotal);
    % mask out unnecessary elements and sum the rest
    anglesSum = sum(sumPrep .* sumMask);

    % calculate intermediate results dependent on the typ of coordinate (X or Y)
    if coordinateType == 'X'
        intermRes = (cosd(anglesSum) .* (repmat(l, 1, stringsTotal))).';
    elseif coordinateType == 'Y'
        intermRes = (sind(anglesSum) .* (repmat(l, 1, stringsTotal))).';
    end

    % prepare intermediate results for summing
    intermSumPrep = flipud(repmat(intermRes, 1, positionsTotal));
    % triangular-like mask for summing
    intermSumMask = repmat(fliplr(tril(ones(positionsTotal, positionsTotal))), stringsTotal, 1);
    % apply triangular-like mask
    intermSum = intermSumPrep .* intermSumMask;
    % matrix for custom summing
    customSumDot = repelem(fliplr(diag(diag(ones(stringsTotal, stringsTotal)))), 1, positionsTotal);
    % dot product to get final results as a sum of intermediate ones
    x = customSumDot * intermSum;
end

function [x, y] = showKinematicString(l, a)
    x = calculateCoordinate(l, a, 'X');
    y = calculateCoordinate(l, a, 'Y');

    % add starting point (0,0) for each string
    x = [zeros(size(a,1), 1), x]';
    y = [zeros(size(a,1), 1), y]';
    
    hold on
    % plot base
    plot(0,0, '^', MarkerSize=10, MarkerEdgeColor='b', MarkerFaceColor='b');

    % calculate axes for segment
    segmentAxisLen = 0.2;

    xAxis = [x(1,:); repelem(x(2:end-1,:), 2, 1); x(end,:)];
    xAxis = reshape(xAxis, 2, numel(xAxis)/2);
    xOrig = xAxis;
    xAxis = xAxis(2,:) - xAxis(1,:);


    yAxis = [y(1,:); repelem(y(2:end-1,:), 2, 1); y(end,:)];
    yAxis = reshape(yAxis, 2, numel(yAxis)/2);
    yOrig = yAxis;
    yAxis = yAxis(2,:) - yAxis(1,:);
    
    l = vecnorm([xAxis; yAxis]);

    xAxis = (xAxis ./ l) * segmentAxisLen;
    yAxis = (yAxis ./ l) * segmentAxisLen;

    segEndX = xOrig(2,:);
    segEndY = yOrig(2,:);

    axisEndX = segEndX + xAxis;
    axisEndY = segEndY + yAxis;
    
    % plot segment axes
    plot([segEndX; axisEndX], [segEndY; axisEndY], '-', LineWidth=2, Color='g');

    % plot movement interpolation
    endPointsX = x(end,:);
    
    minX = min(endPointsX, [], 'all');
    maxX = max(endPointsX, [], 'all');

    range = minX:0.1:maxX;

    interpolatedPointsY = interp1(x(end,:), y(end,:), range, 'spline');
    plot(range,interpolatedPointsY,'-', LineWidth=2, Color='m');

    % plot segments
    plot(x,y,'-o', LineWidth=2 ,MarkerSize=10, MarkerEdgeColor='b', Color='black', MarkerIndices = 2 : length(x) - 1, MarkerFaceColor='b');

    % plot endpoints
    plot(x(end,:), y(end,:), 'o', MarkerSize=10, MarkerEdgeColor='r', MarkerFaceColor='r');
    hold off

    title('Kinematic Strings');
    xlabel('X Axis');
    ylabel('Y Axis');
end
