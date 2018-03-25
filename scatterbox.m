function scatterbox(x, bi)

% multicolor scatterplot of classified data
%
% scatterbox(x, bi)
%
% x:    data (points x variables)
% bi:   corresponding box indices
%
% Generates a scatterplot of the data points, using (as far as possible)
% different colors for points belonging to different classes ("boxes").
% One-dimensional data are blown up to two dimensions to make the plot
% better readable, using randomly generated coordinates. For more than
% three-dimensional data, only the first three variables are used.
%
% Copyright (C) 2009 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

x = x(~isnan(bi), :);
bi = bi(~isnan(bi), :);

nDim = size(x, 2);
newplot

colMax = 15;
if max(bi) <= colMax
    % small number of boxes: unique colors
    colors = hsv(max(bi));
    ci = bi;
else
    % else, randomly assigned colors
    colors = hsv(colMax);
    bi2ci = floor(rand(1, max(bi)) * colMax) + 1;
    ci = bi2ci(bi);
end

% blow up 1d -> 2d
if nDim == 1
    x(:, 2) = rand(size(x)) * 0.01 * (max(x) - min(x));
    axis equal
end

% plot, for each color separately
for i = unique(ci(:)')
    if nDim <= 2
        line(x(ci == i, 1), x(ci == i, 2), 'Color', colors(i, :), ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 1)
    else
        line(x(ci == i, 1), x(ci == i, 2), x(ci == i, 3), 'Color', colors(i, :), ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 1)
    end
end

% set 3d-view
if nDim >= 3
    view(3)
end

