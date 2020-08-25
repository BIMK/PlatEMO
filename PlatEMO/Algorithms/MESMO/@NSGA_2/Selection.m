function [ x_selected ] = Selection( obj )
% SELECTION 
% Binary tournament selection

randnum = ceil( obj.n_pop .* rand( 2 * obj.n_pop ,1 ) );

pop1 = randnum( 1:2:2*obj.n_pop ); % first pop
pop2 = randnum( 2:2:2*obj.n_pop ); % second pop

% Binary tournament using the crowded comparison operator
rank_cond = ( obj.rank_pop(pop1,:) < obj.rank_pop(pop2,:) );
crowd_cond = ( obj.rank_pop(pop1,:) == obj.rank_pop(pop2,:) &....
    obj.crowding_distance(pop1,:) > obj.crowding_distance(pop2,:) );

selection_pop_bin = ( rank_cond | crowd_cond );

x_selected = [ pop1(selection_pop_bin,:) ; pop2(~selection_pop_bin,:) ];

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


