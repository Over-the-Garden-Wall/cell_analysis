addpath /omelette/2/jk/e2198/bin

lstdir=dir('sac*');
for i=1:numel(lstdir)
    if ~(lstdir(i).isdir)
        continue;
    end
	cd(lstdir(i).name);
	d=dir('*.coord');
	findSegmentsListFromSkeletonCoordinate('/omniData/e2198',d.name);
	cd ../
end
