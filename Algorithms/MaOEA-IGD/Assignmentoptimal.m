function [assignment, cost] = Assignmentoptimal(distMatrix)
%ASSIGNMENTOPTIMAL Compute optimal assignment by Munkres algorithm
% ASSIGNMENTOPTIMAL(DISTMATRIX) computes the optimal assignment (minimum
% overall costs) for the given rectangular distance or cost matrix, for
% example the assignment of tracks (in rows) to observations (in
% columns). The result is a column vector containing the assigned column
% number in each row (or 0 if no assignment could be done).
%
% [ASSIGNMENT, COST] = ASSIGNMENTOPTIMAL(DISTMATRIX) returns the
% assignment vector and the overall cost.
%
% The distance matrix may contain infinite values (forbidden
% assignments). Internally, the infinite values are set to a very large
% finite number, so that the Munkres algorithm itself works on
% finite-number matrices. Before returning the assignment, all
% assignments with infinite distance are deleted (i.e. set to zero).
%
%
%
% Markus Buehren
% Last modified 05.07.2011

% save original distMatrix for cost computation
originalDistMatrix = distMatrix;

% check for negative elements
if any(distMatrix(:) < 0)
error('All matrix elements have to be non-negative.');
end

% get matrix dimensions
[nOfRows, nOfColumns] = size(distMatrix);

% check for infinite values
finiteIndex = isfinite(distMatrix);
infiniteIndex = find(~finiteIndex);
if ~isempty(infiniteIndex)
% set infinite values to large finite value
maxFiniteValue = max(max(distMatrix(finiteIndex)));
if maxFiniteValue > 0
infValue = abs(10 * maxFiniteValue * nOfRows * nOfColumns);
else
infValue = 10;
end
if isempty(infValue)
% all elements are infinite
assignment = zeros(nOfRows, 1);
cost = 0;
return
end
distMatrix(infiniteIndex) = infValue;
end

% memory allocation
coveredColumns = zeros(1, nOfColumns);
coveredRows = zeros(nOfRows, 1);
starMatrix = zeros(nOfRows, nOfColumns);
primeMatrix = zeros(nOfRows, nOfColumns);

% preliminary steps
if nOfRows <= nOfColumns
minDim = nOfRows;

% find the smallest element of each row
minVector = min(distMatrix, [], 2);

% subtract the smallest element of each row from the row
distMatrix = distMatrix - repmat(minVector, 1, nOfColumns);

% Steps 1 and 2
for row = 1:nOfRows
for col = find(distMatrix(row,:)==0)
if ~coveredColumns(col)%~any(starMatrix(:,col))
starMatrix(row, col) = 1;
coveredColumns(col) = 1;
break
end
end
end

else % nOfRows > nOfColumns
minDim = nOfColumns;

% find the smallest element of each column
minVector = min(distMatrix);

% subtract the smallest element of each column from the column
distMatrix = distMatrix - repmat(minVector, nOfRows, 1);

% Steps 1 and 2
for col = 1:nOfColumns
for row = find(distMatrix(:,col)==0)'
if ~coveredRows(row)
starMatrix(row, col) = 1;
coveredColumns(col) = 1;
coveredRows(row) = 1;
break
end
end
end
coveredRows(:) = 0; % was used auxiliary above
end

if sum(coveredColumns) == minDim
% algorithm finished
assignment = buildassignmentvector__(starMatrix);
else
% move to step 3
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim); %#ok
end

% compute cost and remove invalid assignments
[assignment, cost] = computeassignmentcost__(assignment, originalDistMatrix, nOfRows);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function assignment = buildassignmentvector__(starMatrix)

[maxValue, assignment] = max(starMatrix, [], 2);
assignment(maxValue == 0) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, cost] = computeassignmentcost__(assignment, distMatrix, nOfRows)

rowIndex = find(assignment);
costVector = distMatrix(rowIndex + nOfRows * (assignment(rowIndex)-1));
finiteIndex = isfinite(costVector);
cost = sum(costVector(finiteIndex));
assignment(rowIndex(~finiteIndex)) = 0;
end
% Step 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step2__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim)

% cover every column containing a starred zero
maxValue = max(starMatrix);
coveredColumns(maxValue == 1) = 1;

if sum(coveredColumns) == minDim
% algorithm finished
assignment = buildassignmentvector__(starMatrix);
else
% move to step 3
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);
end
end
% Step 3: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim)

zerosFound = 1;
while zerosFound

zerosFound = 0;
for col = find(~coveredColumns)
for row = find(~coveredRows')
if distMatrix(row,col) == 0

primeMatrix(row, col) = 1;
starCol = find(starMatrix(row,:));
if isempty(starCol)
% move to step 4
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step4__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, row, col, minDim);
return
else
coveredRows(row) = 1;
coveredColumns(starCol) = 0;
zerosFound = 1;
break % go on in next column
end
end
end
end
end

% move to step 5
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step5__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);
end
% Step 4: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step4__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, row, col, minDim)

newStarMatrix = starMatrix;
newStarMatrix(row,col) = 1;

starCol = col;
starRow = find(starMatrix(:, starCol));

while ~isempty(starRow)

% unstar the starred zero
newStarMatrix(starRow, starCol) = 0;

% find primed zero in row
primeRow = starRow;
primeCol = find(primeMatrix(primeRow, :));

% star the primed zero
newStarMatrix(primeRow, primeCol) = 1;

% find starred zero in column
starCol = primeCol;
starRow = find(starMatrix(:, starCol));

end
starMatrix = newStarMatrix;

primeMatrix(:) = 0;
coveredRows(:) = 0;

% move to step 2
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step2__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);
end
% Step 5: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step5__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim)

% find smallest uncovered element
uncoveredRowsIndex = find(~coveredRows');
uncoveredColumnsIndex = find(~coveredColumns);
[s, index1] = min(distMatrix(uncoveredRowsIndex,uncoveredColumnsIndex));
[s, index2] = min(s); %#ok
h = distMatrix(uncoveredRowsIndex(index1(index2)), uncoveredColumnsIndex(index2));

% add h to each covered row
index = find(coveredRows);
distMatrix(index, :) = distMatrix(index, :) + h;

% subtract h from each uncovered column
distMatrix(:, uncoveredColumnsIndex) = distMatrix(:, uncoveredColumnsIndex) - h;

% move to step 3
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);
end