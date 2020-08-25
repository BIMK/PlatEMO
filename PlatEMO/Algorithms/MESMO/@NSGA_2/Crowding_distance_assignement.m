function [ crowding_distance ] = Crowding_distance_assignement( obj, I, y )
% CROWDING_DISTANCE_ASSIGNEMENT
% 	Calculate crowding distance

n = length(I);
max_y = max( y (I,:) );
min_y = min( y(I,:) );

crowding_distance_temp = Inf *ones( n, size(y,2) );
[ y_sort, ind_sort ] = sort( y(I,:), 1 );
crowding_distance_temp( 2:n-1 , : ) = bsxfun( @rdivide,...
    y_sort(3:n,:) - y_sort(1:n-2,:), max_y - min_y );

for i=1:size(y,2)
    
    crowding_distance( ind_sort(:,i), i ) = crowding_distance_temp( :, i );
    
end

crowding_distance = sum( crowding_distance, 2 );

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


