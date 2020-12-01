%% Function for the calculation of the semivariogram

function [value] = semiVariogram(distance,X)
position = 1;
difference = zeros(1,length(X)-distance);

%% Calculating the difference
while position+distance<=length(X)
    difference(position) = X(position) - X(position+distance);
    position = position+1;
end

%% The semivariogram value for the distance
% squaring the differences
square = difference .* difference;

% Sum the squares
squareSum = sum(square);

% Calculating the value
value = 1/(2 * length(difference))*squareSum;

end
