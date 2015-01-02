function [fileID]=create_hdf5_file_many_chunks(file, path, total_size, chunk_size, init_size, varargin)

if(nargin>5)
	type=varargin{1};
else
	type='float';
end

if(isequal(type,'float'))
	hdf5_type='H5T_NATIVE_FLOAT';
	matlab_type = 'single';
elseif(isequal(type,'int'))
	hdf5_type='H5T_NATIVE_INT';
	matlab_type = 'int32';
end

% if file exists, check to see if the parameters match
if exist(file,'file'),
	h5info = hdf5info(file);
	current_total_size = h5info.GroupHierarchy.Datasets.Dims;
	current_data_type = h5info.GroupHierarchy.Datasets.Datatype.Class;
	if isequal(current_total_size,total_size) && ...
		(isequal(current_data_type,'H5T_IEEE_F32LE') && isequal(type,'float')),
		% parameters match!
		return;
	end
else,
	delete(file)
	create_hdf5_file(file, path, total_size, chunk_size, chunk_size, varargin{:});
	points = mk_split_points(init_size,chunk_size);
	for k = 1:size(points,3),
		block = zeros([points(:,2,k)-points(:,1,k)+1]',matlab_type);
		write_hdf5_file(file, path, points(:,1,k), points(:,2,k), block)
	end
end
