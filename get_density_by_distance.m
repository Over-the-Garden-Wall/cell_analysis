function density = get_density_by_distance(cell_num, varargin)
    
    density = get_density_all(cell_num, varargin{:});    
    density = squeeze(sum(sum(density,3),1));        

end