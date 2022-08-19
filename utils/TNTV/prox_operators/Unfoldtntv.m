
function [X] = Unfoldtntv( X, dim, i )
X = reshape(shiftdim(X,i-1), dim(i), []);