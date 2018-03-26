function [TO, p] = transitionMatrix(bi, segOff, segLen, reversibilize)

% estimate transition matrix of a Markov process
% from state index time series
%
% [TO, p] = transitionMatrix(bi, segOff, segLen, reversibilize = false)
%
% bi:            time series of state/box indices
% segOff:        offsets of data segments (discontinuities), zero-based
% segLen:        length of data segments
% reversibilize: whether the process should be reversibilized
% TO:            transition operator matrix (sparse)
% p:             stationary measure of T0
%
% It is assumed that all indices from 1 to max(bi) actually occur.
% The transition matrix is given in column-stochastic (operator) form:
% TO(i, j) is an estimate of Pr{ bi(t + 1) = i ¦ bi(t) = j }.
%
% Copyright (C) 2009 Carsten Allefeld

    fprintf('  transitionMatrix:\n')

    bi = bi(:);
    segOff = segOff(:);
    nBoxes = max(bi);
    if nargin < 4, reversibilize = false; end
    fprintf('    %d states, %d transitions\n', nBoxes, sum(segLen) - size(segLen, 1))
    
    % count transitions
    c = sparse(nBoxes, nBoxes);
    for si = 1 : size(segOff, 1)
        c = c + sparse(bi(segOff(si) + (2 : segLen(si))), ...
            bi(segOff(si) + (1 : segLen(si) - 1)), 1, nBoxes, nBoxes);
    end

    % reversibilize by counting backwards transitions
    if reversibilize
        c = c + c';
    end

    % make column-stochastic
    if ~all(sum(c) > 0), warning('there are boxes with no detectable outflow!'), end
    TO = c * diag(sparse(1 ./ sum(c)));
    
    % estimate the invariant measure
    if reversibilize
        p = sum(c)' / sum(sum(c));
    else
        [v, d] = eigs(TO, 1, 'lr', setfield([], 'disp', 0));
        p = v / sum(v);
    end

%     % make reversibilization exact – necessary?
%     if reversibilize
%         if ~all(p > 0), error('there are boxes with no detectable inflow!'), end
%         % compute time-reversed transition matrix
%         TOit = (TO * diag(sparse(p)))' * diag(sparse(1 ./ p));
%         % compute reversibilized transition matrix
%         TO = (TO + TOit) / 2;
%     end


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

