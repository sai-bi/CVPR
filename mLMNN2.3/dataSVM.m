% This program is used for turning feature.mat into a file
% in svm-light format.
% @Bi Sai, 07-31-2013

% Specifications for svm-light format:
% Each line is a training instance:
% <line> .=. <target> <feature>:<value> <feature>:<value> ... <feature>:<value> # <info>
% <target> .=. +1 | -1 | 0 | <float>   (Class of the instance)
% <feature> .=. <integer> | "qid"      (feature id)
% <value> .=. <float>                  (feature value)
% <info> .=. <string>                  (extra info)

% for example, the line
% -1 1:0.43 3:0.12 9284:0.2 # abcdef
% specifies a negative example for which feature number 1 has the value 0.43, 
% feature number 3 has the value 0.12, feature number 9284 has the value 0.2, 
% and all the other features have value 0

% load('finalFeature.mat');
load('data/viper_features.mat');
load('uxLdim2.mat');


feature = ux;

fileID = fopen('trainHSV.txt','w+');

perm = randperm(632);
save('perm.mat','perm');
training1 = feature(:,perm(1:316));
training2 = feature(:,perm(1:316)+632);

test1 = feature(:,perm(317:632));
test2 = feature(:,perm(317:632)+632);

feature = uxLdim * feature;
feature = feature(1,:);

target1 = feature(:,perm(1:316));
target2 = feature(:,perm(1:316)+632);

target3 = feature(:,perm(317:632));
target4 = feature(:,perm(317:632)+632);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write training1 to trainHSV
qid = 1;
for i = 1:270
   
	fprintf(fileID,'%d ',target1(i));
    fprintf(fileID,'qid:%d ',qid);
    qid = qid+1;
	for j = 1:size(training1,1)
		fprintf(fileID, '%d:%f ',j-1,training1(j,i));
	end
	fprintf(fileID,'\n');
end
% write training2 to trainHSV
for i = 1:270
	fprintf(fileID,'%d ',target2(i));
    fprintf(fileID,'qid:%d ',qid);
    qid = qid+1;
	for j = 1:size(training2,1)
		fprintf(fileID, '%d:%f ',j-1,training2(j,i));
	end
	fprintf(fileID,'\n');
end

fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen('valHSV.txt','w+');
for i = 271:316
	fprintf(fileID,'%d ',target1(i));
    fprintf(fileID,'qid:%d ',qid);
    qid = qid+1;
	for j = 1:size(training1,1)
		fprintf(fileID, '%d:%f ',j-1,training1(j,i));
	end
	fprintf(fileID,'\n');
end

for i = 271:316
	fprintf(fileID,'%d ',target2(i));
    fprintf(fileID,'qid:%d ',qid);
    qid = qid+1;
	for j = 1:size(training2,1)
		fprintf(fileID, '%d:%f ',j-1,training2(j,i));
	end
	fprintf(fileID,'\n');
end
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen('testHSV.txt','w+');
for i = 1:316
	fprintf(fileID,'%d ',target3(i));
    fprintf(fileID,'qid:%d ',qid);
    qid = qid+1;
	for j = 1:size(test1,1)
		fprintf(fileID, '%d:%f ',j-1,test1(j,i));
	end
	fprintf(fileID,'\n');
end

for i = 1:316
	fprintf(fileID,'%d ',target4(i));
    fprintf(fileID,'qid:%d ',qid);
    qid = qid+1;
	for j = 1:size(test2,1)
		fprintf(fileID, '%d:%f ',j-1,test2(j,i));
	end
	fprintf(fileID,'\n');
end
fclose(fileID);


















