function [bi, nBoxes] = discretize(x, nMin)

% discretize through binary division
%
% bi = discretize(x, nMin)
%
% Convert a multidimensional real-valued data set x
% into a one-dimensional series of box indices bi.
% Indices are allocated only for occupied boxes.
%
% x:            data (samples x variables)
% nMin:         minimum number of points per box
% bi:           box index sequence
%
% the syntax [bi, nBoxes] = discretize(x, nMin)
% is used for internal recursive calls
%
% Copyright (C) 2009 Carsten Allefeld


persistent nic

if nargout < 2
    fprintf('discretize:\n')
    nic = 0;
end

n = size(x, 1);
bi = ones(n, 1);
nBoxes = 1;

% check whether to subdivide
if n < 2 * nMin, return, end

% dominant direction
[v, l] = eig(double(cov(x)));
[~, ind] = max(diag(l));
pc = x * v(:, ind);

% seek large gap, or do median cut
spc = sort(pc);
dpc = diff(spc);
i = find(dpc / (spc(end) - spc(1)) > 0.5);

% divide at gap, or median
if isempty(i) || (i < nMin) || (n - i < nMin)
    i = floor((n + 1) / 2);     % median cut
else
    nic = nic + 1;              % irregular cut
end
threshold = spc(i);

% determine parts, and discretize them
ind = (pc > threshold);
[bi0, nBoxes0] = discretize(x(~ind, :), nMin);
[bi1, nBoxes1] = discretize(x(ind, :), nMin);

% put sub-discretizations together
bi(~ind) = bi0;
bi(ind) = bi1 + nBoxes0;
nBoxes = nBoxes0 + nBoxes1;

% statistics
if nargout < 2
    fprintf('  %d irregular cuts occurred\n', nic)
    bon = hist(bi, 1 : nBoxes);
    fprintf('  %d boxes occupied, occupation numbers: min %d, median %d, max %d.\n', ...
        nBoxes, min(bon), median(bon), max(bon))
end


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

