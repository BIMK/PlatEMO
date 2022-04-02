function obj = Infill_Maximal_Distance(x, sample_x)
obj = - min(pdist2(x, sample_x),[],2);
end