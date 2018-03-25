function result = almostInvariant(bi, param)

% search for almost invariant sets
% 
% result = almostInvariant(bi)
% result = almostInvariant(bi, param)
%
% bi:             box indices (microstate sequence)
%
% param is a struct with the optional fields
%   segOff, segLen: for temporally discontinuous data, specify the offsets
%                   and lengths of continuous segments
%   q:              force number of almost invariant sets to search for
%
% result is a struct with the fields
%   bi2ci:          mapping from microstates to macrostates
%   l:              eigenvalue spectrum
%   T:              timescale spectrum
%   F:              timescale separation factors
%   o:              positions of microstates in eigenvector space
%   chi:            almost characteristic functions across microstates
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


fprintf('almostInvariant:\n\n')

% parameters, heuristic setting
nEval = 20;

if nargin < 2
    param = struct;
end

if ~exist('param.segOff', 'var')
    param.segOff = 0;
    param.segLen = numel(bi);
end


% determine transition matrix
[R, vpi] = transitionMatrix(bi, param.segOff, param.segLen, true);
fprintf('\n')

% compute eigenvalues and left eigenvectors
[A, l] = eigs(R', nEval, 'lr', setfield([], 'disp', 0)); %#ok<SFLD>
l = diag(l);
[~, ind] = sort(abs(l), 'descend');
l = l(ind);
A = A(:, ind)';

% normalize left eigenvectors
A = diag(1 ./ sqrt( sum(A .^ 2 * diag(sparse(vpi)), 2) )) * A;

% compute timescales & separation factors
T = [Inf; -1 ./ log(abs(l(2 : nEval)))];
F = T(1 : nEval - 1) ./ T(2 : nEval);

% if number of clusters not specified, select
if exist('param.q', 'var')
    q = param.q;
else
    [~, ind] = max(F(2 : end));
    q = ind + 1;
end

% determine timescale for number of clusters q
Tq = T(q + 1);
fprintf('  q      = %d\n', q)
fprintf('  F(q)   = %.2f\n', F(q))
fprintf('  T(q+1) = %.1f\n\n', Tq)

% compute positions in eigenvector space
o = A(2 : q, :)';

% determine macrostate as microstate clusters via PCCAP+
[chi, bi2ci] = pccap(A(1 : q, :)');

% package results
result = struct('bi2ci', bi2ci, ...
    'T', T, 'l', l, 'F', F, ...
    'o', o, 'chi', chi);
