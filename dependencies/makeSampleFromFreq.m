function X = makeSampleFromFreq(V, f)

% Generates a matrix X of rows of V such that v(i, :) occurs f(i) times, 
% in increasing order of i, such that the number of rows of X is sum(f)
%
% USAGE:  X = makeVectorFromFreq(V, f)
%
% INPUTS:  V - column vector or matrix of values
%          f - column vector of required frequencies 
% OUTPUTS: X - column vector or matrix of samples of rows of v with the
% specified frequencies 
% 
% EXAMPLE: X = makeVectorFromFreq([1 2 3 4 1 2 3 4]', [0 0 1 0 0 2 0 2]')
%    returns X = [3 2 2 4 4]'
%
% NOTE: makeVectorFromFreq is the inverse operation of 'tabulate'
%    e.g. tbl = tabulate(v);
%    v = makeVectorFromFreq(tbl(:, 1), tbl(:, 2)) 
% should recover the same starting vector v

nElements = sum(f);
[nRows, nCols] = size(V);
X = zeros(nElements, nCols);
iSample = 1;
for iRow = 1:nRows
    X(iSample:iSample+f(iRow)-1, :) = repmat(V(iRow, :), f(iRow), 1) ;
    iSample = iSample+f(iRow);
end


