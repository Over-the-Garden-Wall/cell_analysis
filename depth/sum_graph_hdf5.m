% sum an hdf5 graph into an in-out graph
function [h5outputfile] = sum_graph_hdf5(h5inputfile, h5inputpath, h5outputfile, h5outputpath, begin_coords, end_coords)

blksz = [300 300 300];
		
% open input file
[ndims, input_size, max_size]=get_hdf5_size(h5inputfile, h5inputpath);
input_dataID=H5F.open(h5inputfile,'H5F_ACC_RDONLY','H5P_DEFAULT');
input_datasetID=H5D.open(input_dataID, h5inputpath);
input_dataspaceID=H5D.get_space(input_datasetID);

% create and then open output file
create_hdf5_file(h5outputfile, h5outputpath, input_size(1:3), blksz, blksz, 'float');
output_dataID=H5F.open(h5outputfile,'H5F_ACC_RDWR','H5P_DEFAULT');
output_datasetID=H5D.open(output_dataID, h5outputpath);	
output_dataspaceID=H5D.get_space(output_datasetID);

if(isempty(begin_coords))
	begin_coords=ones([1, ndims]);
end

if(isempty(end_coords))
	end_coords=input_size';
end

sz=end_coords-begin_coords+1;

%% connect components of the blocks
for iblk = begin_coords(1):blksz(1):end_coords(1),
	for jblk = begin_coords(2):blksz(2):end_coords(2),
		for kblk = begin_coords(3):blksz(3):end_coords(3),
						
			% get this block of the component file
			block_begin_coords=[iblk jblk kblk];
			block_end_coords=[min(end_coords(1),iblk+blksz(1)-1) min(end_coords(2),jblk+blksz(2)-1) min(end_coords(3),kblk+blksz(3)-1)];				
			g_block=get_hdf5(input_datasetID, input_dataspaceID, [block_begin_coords 1], [block_end_coords 3]);

			inout_block=single(sum(g_block,4)./input_size(4));
			
			% write out new components
			write_hdf5(output_datasetID, output_dataspaceID, block_begin_coords, block_end_coords, inout_block);					
			
			fprintf(1, '.');
		end
	end
end

H5D.close(input_datasetID);
H5F.close(input_dataID);

H5D.close(output_datasetID);
H5F.close(output_dataID);

fprintf(1, '\n');

