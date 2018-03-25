function [chi, ci] = pccap(X, lambda, useScaling)

% PCCA+ algorithm of Deuflhard & Weber (2005)
%
% [chi, ci] = pccap(X)
% [chi, ci] = pccap(X, lambda, useScaling = false)
%
% X:        strongest eigenvectors (columns)
% lambda:   corresponding eigenvalues
% chi:      characteristic functions of almost invariant sets
% ci:       indices of ais with largest chi
%
% The Perron index k is derived from the number of eigenvectors given.
% If eigenvalues are not specified, the "scaling" objective function (I1)
% is used for optimization instead of "metastability" (I2).
% Coefficients of the dominant eigenvector are assumed to be equal to 1,
% the first (largest) eigenvalue is assumed to be 1.
%
% Reference: P. Deuflhard & M. Weber, "Robust Perron cluster analysis in
% conformation dynamics", Linear Algebra Appl. 398 (2005), 161-184
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


% parameters
[N, k] = size(X);
if nargin < 3, useScaling = (nargin == 1); end

fprintf('  pccap:\n')
fprintf('    Perron index: %d\n', k)
if useScaling
    fprintf('    using "scaling" objective function\n')
else
    fprintf('    using "metastability" objective function\n')
end

% fix some stuff
X(:, 1) = 1;

% "Subalgorithm A"
% initial simplex guess
l = simplexGuess(X(:, 2 : end));
% corresponding transform & reduced matrix
Abar = inv(X(l, :));
AP = Abar;
AP(1, :) = [];
AP(:, 1) = [];
A0 = full(AP, X);

% objective function
if useScaling
    objective = @(ap) (-scaling(ap, X));
else
    objective = @(ap) (-metastability(ap, X, lambda));
end

% optimization
[AP, ~, conv] = fminsearch(objective, AP, optimset('MaxIter', 10000, 'MaxFunEvals', 100000));
if ~conv, error('optimization did not converge!'), end
A = full(AP, X);

% characteristic functions
chi = X * A;

% additional information
p0 = A(1, :);
I1 = scaling(AP, X);
if nargin > 1
    w = (A' * diag(lambda) * A) ./ (A(1, :)' * ones(1, k));
    I2 = metastability(AP, X, lambda);
end

% attribution to almost invariant sets
[~, ci] = max(chi, [], 2);

% diagnostic output
fprintf('    scaling: I1 = %.3f\n', I1)
fprintf('    delta = %.3f\n', 1 - I1 / k)
if nargin > 1
    I2 = metastability(AP, X, lambda);
    fprintf('    metastability: I2 = %.3f < %.3f\n', I2, sum(lambda))
end

% visual output
if nargout == 0
    figure
    plot(X(:, 2), X(:, 3), '.')
    hold all
    vert = inv(Abar);
    plot(vert([1 : end, 1], 2), vert([1 : end, 1], 3), 'o-')
    vert = inv(A0);
    plot(vert([1 : end, 1], 2), vert([1 : end, 1], 3), 'o-')
    vert = inv(A);
    plot(vert([1 : end, 1], 2), vert([1 : end, 1], 3), 'ko-')

    figure
    subplot(3, 1, 1)
    plot(X * Abar, '.-')
    axis([0.5, N + 0.5, -0.1, 1.1])
    subplot(3, 1, 2)
    plot(X * A0, '.-')
    axis([0.5, N + 0.5, -0.1, 1.1])
    subplot(3, 1, 3)
    plot(X * A, '.-')
    axis([0.5, N + 0.5, -0.1, 1.1])
end


function I1 = scaling(AP, X)
% "scaling" objective function, Eq. 4.19
I1 = sum(max(X * full(AP, X)));


function I2 = metastability(AP, X, lambda)
% "metastability" objective function, Eq. 4.22
A = full(AP, X);
I2 = sum(sum(diag(lambda) * A .^ 2) ./ A(1, :));


function A = full(AP, X)
% "Subalgorithm B"
A = blkdiag(0, AP);
A(2 : end, 1) = -sum(A(2 : end, 2 : end), 2);
A(1, :) = -min(X(:, 2 : end) * A(2 : end, :));
A = A / sum(A(1, :));


function l = simplexGuess(X)
% determine nodes of a simplex approximating the data
% part of "Subalgorithm A"
[N, dim] = size(X);
% point farthest from origin
[~, l] = max(sum(X .^ 2, 2));
% translate it to the origin
X0 = X - ones(N, 1) * X(l, :);
for i = 1 : dim
    % point farthest from origin
    [~, ln] = max(sum(X0 .^ 2, 2));
    l = [l, ln];
    % remove corresponding subspace
    v = X0(ln, :);
    v = v / sqrt(sum(v .^ 2));
    X0 = X0 - X0 * v' * v;
end
