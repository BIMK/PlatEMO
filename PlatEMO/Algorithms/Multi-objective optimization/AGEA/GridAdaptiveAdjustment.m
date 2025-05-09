function [Grid, Grid_index, div] = GridAdaptiveAdjustment(Grid, zmin, gmax, div, N) 
    [Grid, Grid_index] = EnvironmentalSelection_AGEA(Grid, zmin, gmax, div); 
    if length(Grid) > N && div > 2
        div = div - 1;
        [Grid_temp, Grid_index_temp] = EnvironmentalSelection_AGEA(Grid, zmin, gmax, div);
        if length(Grid_temp) > N  && div >= 2
            Grid = Grid_temp;
            Grid_index = Grid_index_temp;
        else
            div = div + 1;
        end
    end
    if length(Grid) < N && div < 200
        div = div + 1;
    end
end