function [ int_Front_population,int_rank_pop ] = Fast_non_dominated_sort( obj, y, pop_number )
% FAST_NON_DOMINATED_SORT
%   Sort data based on the front rank

y_comparaison = permute( y, [3 2 1] );

domination = false( size(y,1), size(y,1) );

domination(:,:) = permute( all( bsxfun( @le, y, y_comparaison ),2) ...
    & any( bsxfun( @lt, y, y_comparaison ), 2 )...
    ,[1 3 2]);

% Build domination matrix (1 if row is undominated by col)
if obj.constraint_logical
    
    feasible = all( obj.g <= 0, 2 );
    distance_feas = sqrt( sum( ( (obj.g > 0) .* obj.g ) .^ 2, 2)) ;
    feasible_comparaison = permute( feasible, [3 2 1] );
    distance_comparaison = permute( distance_feas, [3 2 1] );
    
    % Both solutions are feasible
    both_feas1 = permute( bsxfun( @and, feasible ,feasible_comparaison ), [1 3 2] );
    
    % One solution is feasible
    one_feas = permute( bsxfun( @and, feasible, ~feasible_comparaison ), [1 3 2] );
    
    % Both solution are unfeasible
    both_violate1 = permute( bsxfun( @and, ~feasible, ~feasible_comparaison ), [1 3 2] );
    both_violate2 = permute( bsxfun( @lt, distance_feas, distance_comparaison ), [1 3 2] );
    
    both_violate = ( both_violate1 & both_violate2 );
    
    both_feas =( domination & both_feas1);
    
    domination_final = both_feas | one_feas| both_violate;
    
else
    
    domination_final = domination;
    
end

% Sort by front
front = 1;
int_rank_pop = zeros(obj.n_pop,1);
int_Front_population = cell(1,1);

while ~isempty( domination_final )
    
    order = sum( domination_final, 1 )' == 0;
    int_Front_population{ front, : } = pop_number( order );
    int_rank_pop( pop_number( order ), : ) = ones( sum(order), 1 ) * front;
    pop_number(order) = [];
    domination_final = domination_final( ~order, ~order );
    front = front + 1;
    
end

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


