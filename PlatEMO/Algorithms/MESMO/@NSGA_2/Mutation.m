function [ x_mutated ] = Mutation( obj, x_crossed )
% MUTATION 
%   Mutate parameters values using gaussian mutation

mutation_flag = rand( obj.n_pop, obj.n_var ); % rand for mutation selection

% Estimate mutation value
scale_mutation_new = obj.scale_mutation - ...
    obj.shrink_mutation * obj.scale_mutation * obj.n_gen / obj.max_gen;

mutation_val = scale_mutation_new .* (obj.ub - obj.lb);

% Mutation
x_mutated = x_crossed + ...
    ( mutation_flag < obj.fraction_mutation ) .* ...
    bsxfun( @times, randn( obj.n_pop, obj.n_var ), mutation_val );

% Keep parameters into bounds
A = bsxfun( @lt, x_mutated, obj.lb );
B = bsxfun( @gt, x_mutated, obj.ub );

x_mutated = x_mutated .* (~A & ~B) + ...
    bsxfun( @times, A, obj.lb ) + bsxfun( @times, B, obj.ub );

end





% ==========================================================================
%
%    This file is part of SBDOT.
%
%    SBDOT is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    SBDOT is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with SBDOT.  If not, see <http://www.gnu.org/licenses/>.
%
%    Use of SBDOT is free for personal, non-profit, pure academic research
%    and educational purposes. Restrictions apply on commercial or funded
%    research use. Please read the IMPORTANT_LICENCE_NOTICE file.
%
% ==========================================================================


