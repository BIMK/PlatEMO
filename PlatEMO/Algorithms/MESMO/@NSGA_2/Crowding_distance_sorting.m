function [ Front_population ] = Crowding_distance_sorting( obj, int_Front_population, int_rank_pop )
% CROWDING_DISTANCE_SORTING
%   Sort population based on front number and the crowding distance

N=0;
iter=1;
crowding_distance_temp=[];
Front_population=[];

while obj.n_pop > ( N + length( int_Front_population{iter,:} ) )
    
    int_crowding_distance = obj.Crowding_distance_assignement( int_Front_population{iter,:}, obj.y );
    Front_population = [ Front_population ; int_Front_population{iter,:}' ];
    crowding_distance_temp = [ crowding_distance_temp ; int_crowding_distance ];
    N = N + length( int_Front_population{iter,:} );
    iter = iter+1;
    
end

if obj.n_pop == ( N + length( int_Front_population{iter,:} ) )
    
    int_crowding_distance = obj.Crowding_distance_assignement ( int_Front_population{iter,:} , obj.y );
    Front_population = [ Front_population ; int_Front_population{iter,:}' ];
    obj.crowding_distance = [ crowding_distance_temp ; int_crowding_distance ];
    
else
    
    int_crowding_distance = obj.Crowding_distance_assignement( int_Front_population{iter,:}, obj.y );
    sort_matrix = flipud( sortrows( [ int_crowding_distance, int_Front_population{iter,:}' ], 1 ) );
    Front_population = [ Front_population ; sort_matrix( 1:obj.n_pop-N,2 ) ];
    obj.crowding_distance = [ crowding_distance_temp ; sort_matrix( 1:obj.n_pop-N , 1) ];
    
end

obj.rank_pop = int_rank_pop( Front_population );

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


