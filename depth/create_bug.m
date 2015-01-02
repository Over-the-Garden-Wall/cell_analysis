% CREATE AN HDF5 file with 4 dimensional single array 
fileID=H5F.create('test.h5', 'H5F_ACC_EXCL', 'H5P_DEFAULT', 'H5P_DEFAULT');
cparms=H5P.create('H5P_DATASET_CREATE');
H5P.set_chunk(cparms, [3 10 10 10]);
dataspaceID=H5S.create_simple(4, [3 100 100 100], [-1 -1 -1 -1]);
datasetID=H5D.create(fileID, '/main', 'H5T_NATIVE_FLOAT', dataspaceID, cparms);
memspace=H5D.get_space(datasetID);
H5D.write(datasetID, 'H5T_NATIVE_FLOAT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', single(zeros([100 100 100 3]))); 

[numdims dims maxdims]=H5S.get_simple_extent_dims(dataspaceID);
dims

% READ a 50x50x2x3 block - works correctly
start_index=[0 45 50 50];
length=[3 2 50 50];
H5S.select_hyperslab(dataspaceID,'H5S_SELECT_SET',start_index, [1 1 1 1], length, [1 1 1 1]);
output_dataspace=H5S.create_simple(4, length, length);
block=H5D.read(datasetID,'H5ML_DEFAULT',output_dataspace,dataspaceID,'H5P_DEFAULT');
size(block)

% READ a 50x50x1x3 block - works INCORRECTLY, returns a 3x1x50x50 block
length=[3 1 50 50];
H5S.select_hyperslab(dataspaceID,'H5S_SELECT_SET',start_index, [1 1 1 1], length, [1 1 1 1]);
output_dataspace=H5S.create_simple(4, length, length);
block=H5D.read(datasetID,'H5ML_DEFAULT',output_dataspace,dataspaceID,'H5P_DEFAULT');
size(block)