function writeFile(fileID, data, dim, mData)
	pid = 1; %preference constraint, useless here
	for i = 1:size(mData,2)
		fprintf(fileID, '%f ', mData(dim,i));
		fprintf(fileID,'pid:%d ',pid);
		pid = pid + 1;
		for j = 1:size(data,1)
			fprintf(fileID,'%d:%f ',j-1,data(j,i));
		end
		fprintf(fileID,'\n');
	end
end