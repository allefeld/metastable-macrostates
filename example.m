% example: 4-phase system with hierarchical structure
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

clear
fprintf('----------\n')


% parameters
nPoints = 1e6;
nMin = 200;
a = 0.01;
b = [0.03, 0.05];

% generate data
xc = rand(1, 2) - 0.5;
x = zeros(nPoints, 2);
for i = 1 : nPoints
    xc = xc + a * (xc - 2 * xc .^ 3) + b .* randn(1, 2);
    x(i, :) = xc;
end
        
% discretize
bi = discretize(x, nMin);
fprintf('\n')


% determine almost invariant sets
result = almostInvariant(bi);
nEval = size(result.l, 1);
ci = result.bi2ci(bi);
q = max(ci);


% plot microstates in data space
figure
scatterbox(x, bi)
axis equal
title('microstates', 'FontSize', 12)
xlabel('x_1')
ylabel('x_2')

% plot eigenvalue spectrum
figure
subplot(3, 1, 1)
plot(1 : nEval, abs(result.l), '.-')
ylabel('|\lambda_k|')
xlim([0.5, nEval + 0.5])
subplot(3, 1, 2)
plot(2 : nEval, result.T(2 : nEval), '.-')
ylabel('T(k)')
xlim([0.5, nEval + 0.5])
subplot(3, 1, 3)
plot(1 : nEval - 1, result.F, '.-')
ylabel('F(k)')
xlabel('k')
xlim([0.5, nEval + 0.5])

% plot macrostates in eigenvector space
figure
scatterbox(result.o, result.bi2ci)
xlabel('o_1')
ylabel('o_2')
zlabel('o_3')
axis equal
axis vis3d
rotate3d on
title('eigenvector space', 'FontSize', 12)

% plot macrostates in data space
figure
scatterbox(x, ci)
axis equal
title('macrostates', 'FontSize', 12)

% plot macrostate sequence
figure
plot(ci, '.', 'MarkerSize', 1)
ylim([0.5, q + 0.5])
xlabel('time')
ylabel('macrostate')
set(gca, 'YTick', 1 : q)

