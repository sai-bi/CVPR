load('permuation.mat');
load('HSVPart.mat');
load('YUVPart.mat');

HSVSample = LdimHSV * HSVsample;
YUVSample = LdimYUV * YUVSample;
sample = [HSVSample;YUVSample];

HSVYUVTest1 = [LdimHSV * HSVTest1;LdimYUV * YUVTest1];
HSVYUVTest2 = [LdimHSV * HSVTest2;LdimYUV * YUVTest2];

%[LdimHSVYUV,det1] = multiLMNN(sample,[1:316 1:316],1,'maxiter',450,'checkup',0,'out);

save('HSVYUV.mat','HSVYUVTest1','HSVYUVTest2','LdimHSVYUV');


