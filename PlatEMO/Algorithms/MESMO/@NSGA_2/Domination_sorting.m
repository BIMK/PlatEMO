function [] = Domination_sorting( obj )
% DOMINATION_SORTING 
% Sort population based on domination operator

[ int_Front_population, int_rank_pop ] = obj.Fast_non_dominated_sort( obj.y, 1:size(obj.y,1) );

[Front_population] = obj.Crowding_distance_sorting( int_Front_population, int_rank_pop );

obj.x = obj.x( Front_population, : );
obj.y = obj.y( Front_population, : );

if obj.constraint_logical
    obj.g = obj.g( Front_population, : );
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


