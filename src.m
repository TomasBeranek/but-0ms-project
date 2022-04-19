% File    project2.m
% Author  Tomas Beranek <xberan46@stud.fit.vutbr.cz>
% Brief   Visualization of kinematic string
% Date    19.4.2022
% Up2date sources can be found at 

% values from the first assignment
l = [2, 1, 0.5];
a = [[70, -100, -20]; [80, -120, -30]; [90, -140, -50]; [110, -160, -30]];

[num, txt, raw] = xlsread("cv9.xlsx");

% load number of segments and number of strings from .xlsx file
lengthsNum = raw{2,1};
instancesNum = raw{2,2};

% load lengths and angles from .xlsx file
l = cell2mat(raw(4, 1:lengthsNum));
a = cell2mat(raw(6:5+instancesNum, 1:lengthsNum));

% calculate joint coordinates and display them in a graph
[x,y] = showKinematicString(l,a);

% safe results
raw(6+instancesNum,1) = {'Polohy kloubů X'};
raw(7+instancesNum:6+2*instancesNum, 1:(lengthsNum + 1)) = num2cell(x');
raw(7+2*instancesNum,1) = {'Polohy kloubů Y'};
raw(8+2*instancesNum:7+3*instancesNum, 1:(lengthsNum+1)) = num2cell(y');

% write everything back the original .xslx file
writecell(raw, "cv9.xlsx", 'UseExcel', false);

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

    plot(x,y,'-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', 'r', 'Color', 'black');
    title('Kinematic Strings');
    xlabel('X Axis');
    ylabel('Y Axis');
end