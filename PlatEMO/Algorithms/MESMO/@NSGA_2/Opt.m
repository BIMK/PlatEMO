function [] = Opt( obj )
% OPT 
%   Launch NSGA-II algorithm

if obj.display
    fprintf(sprintf('NSGA-II current generation : %3.0i / %3.0i',obj.n_gen,obj.max_gen))
end
while obj.n_gen <= obj.max_gen
    
    % Selection
    x_selected = obj.Selection();
    % Crossing
    x_crossed = obj.Croisement( obj.x(x_selected,:) );
    % Mutation
    x_mutation = obj.Mutation( x_crossed );
    
    % Function evaluation
    [ y_new, g_new ]=feval( obj.function_name, x_mutation );
    
    % Combine previous and new population
    obj.x = [ obj.x ; x_mutation ];
    obj.y = [ obj.y ; y_new ];
    if obj.constraint_logical
        obj.g = [ obj.g; g_new ];
    end
    
    % Sorting fronts and extract only n_pop points
    obj.Domination_sorting();
    
    % Display
    if obj.display
        fprintf(repmat('\b',1,9));
        fprintf(sprintf('%3.0i / %3.0i',obj.n_gen,obj.max_gen))
    end
    
    % Next generation
    obj.n_gen=obj.n_gen+1;
    
    % Store history
    obj.hist(obj.n_gen).x = obj.x;
    obj.hist(obj.n_gen).y = obj.y;
    if obj.constraint_logical
        obj.hist(obj.n_gen).g = obj.g;
    end
end

if obj.display, fprintf('\n');end

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


