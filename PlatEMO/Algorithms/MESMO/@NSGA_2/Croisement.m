function [ x_crossed ] = Croisement( obj, x_selected )
% CROISEMENT
% Crossing parameters values

param_pop1 = x_selected( 1:2:obj.n_pop, : ); % parent 1
param_pop2 = x_selected( 2:2:obj.n_pop, : ); % parent 2

croisement_flag = rand( obj.n_pop/2, obj.n_var ); % rand for crossing selection
randNum = rand( obj.n_pop/2,obj.n_var ); % rand to change the parents

param_pop1_croisement = param_pop1 + ...
    ( croisement_flag < obj.fraction_croisement ) .* randNum.* ...
    obj.ratio_croisement .* ( param_pop1 - param_pop2 );
param_pop2_croisement = param_pop2 -...
    ( croisement_flag < obj.fraction_croisement ) .* randNum .* ...
    obj.ratio_croisement .* ( param_pop1 - param_pop2 );

x_crossed = [ param_pop1_croisement ; param_pop2_croisement ]; % childrens

% Keep parameters into bounds
A = bsxfun( @lt, x_crossed, obj.lb );
B = bsxfun( @gt, x_crossed, obj.ub );

x_crossed = x_crossed .* (~A & ~B) + ...
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


