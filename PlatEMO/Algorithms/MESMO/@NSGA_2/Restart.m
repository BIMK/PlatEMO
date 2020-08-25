function [] = Restart( obj, iter_sup )
% RESTART
% Restart an optimization sequence with 50 generations or an user defined nummber
%
% Syntax :
%   objEGO.restart()
%   objEGO.restart(iter_sup)

if nargin==0
    
    obj.max_gen = 10 + obj.max_gen ;
    
else
    
    obj.max_gen = obj.max_gen + iter_sup;
    
end

obj.Domination_sorting();
obj.Opt();

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


