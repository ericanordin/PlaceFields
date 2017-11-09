function blankedMatrix = blankCentre(unblankedMatrix, matrixLength, centreRadius)
%blankCentre converts the centre of the place field to 0s for correlation
%maximization based on the annulus
%   After the centre of the matrix is converted to 0s, the rotation
%   maximization considers only the area around the blank centre for
%   maximization

% Input:
% unblankedMatrix = matrix before centre blanking
% matrixLength = side length of matrix (assumes 2D square)
% centreRadius = radius of the circle marking the centre of the place cell

% Output:
% blankedMatrix = matrix after centre blanking

if ~exist('centreRadius', 'var')
    centreRadius = 10; %Default = 10 pixels
end

indexRange = 1:matrixLength;
indexRange = indexRange-ceil(matrixLength/2);

[X,Y] = meshgrid(indexRange, indexRange);
%Produces indeces relative to 0

circle_i = sqrt(X.^2 + Y.^2)<centreRadius; 
%centreRadius pixels on either side of the centre point
%Index of 0s and 1s with 1s in the centre circle

blankedMatrix = unblankedMatrix;

blankedMatrix(circle_i) = 0; %Create zeros in central circle where circle_i == 1
%Only values corresponding to circle_i == 1 are impacted

end

