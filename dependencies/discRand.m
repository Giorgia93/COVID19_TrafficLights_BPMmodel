function x = discRand(v, nRows, mCols)
% Generate random numbers from an arbitrary discrete probability distribution with probability mass function described by the vector v
% USAGE: x = discRand(v, nRows, mCols)
% INPUTS: v - 1 x n vector of probability weights (these do NOT need to sum
% to 1 as this normalisation is done by the fucntion) of the random
% variable taking values [0, 1, 2, ...]
%         nRows - number of rows in the output matrix
%	  mCols - number of columns in the output matirxc
% OUTPUTS: x - matrix of random numbers generated


cdf = cumsum(v);
cdf = cdf/cdf(end);
x = 0:length(cdf);
y = [0, cdf];

keepFlag = [y(1:end-1) < y(2:end), 1];
x = x(keepFlag == 1);
y = y(keepFlag == 1);

r = rand(nRows, mCols);

try
    x = interp1(y, x, r, 'previous');
catch
   v
   v(end)
   v > 0
   cdf
   keepFlag
   x
   y
   
end



