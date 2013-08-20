function [L,Det]=lmnn2(x,y,varargin)
%
% function [L,Det]=lmnn(x,y,Kg,'Parameter1',Value1,'Parameter2',Value2,...);
%
% Input:
%
% x = input matrix (each column is an input vector) 
% y = labels 
% (*optional*)  L = initial transformation matrix (e.g eye(size(x,1)))
% (*optional*) Kg = attract Kg nearest similar labeled vectos
% 
% Important Parameters:
% diagonal = (default false) If set to true, a diagonal matrix is learned  
% stepsize = (default 1e-09)
% outdim = (default: size(x,1)) output dimensionality
% maxiter = maximum number of iterations (default: 1000)
% validation = (def 0) fraction of training data to be used as
%               validation set (best output is stored in Det.bestL)
% validationstep = (def 50) every "valcount" steps do validation
% quiet = {0,1} surpress output (default=0)  
%
%
% Specific parameters (for experts only):
% correction = (def 15) how many steps between each update 
%              The number of impostors are fixed for until next "correction"
% factor = (def 1.1) multiplicative factor by which the
%         "correction" gab increases
% obj = (def 1) if 1, solver solves in L, if 0, solver solves in L'*L
% thresho = (def 1e-9) cut off for change in objective function (if
%           improvement is less, stop)
% thresha = (def 1e-22) cut off for stepsize, if stepsize is
%           smaller stop
% scale = (def. 0) if 1, all data gets re-scaled s.t. average
%         distance to closest neighbor is 1
%
%
% Output:
%
% L = linear transformation xnew=L*x
%    
% Det.obj = objective function over time
% Det.nimp = number of impostors over time
% Det.pars = all parameters used in run
% Det.time = time needed for computation
% Det.iter = number of iterations
% Det.verify = verify (results of validation - if used)
%  
% Version 2.3
% copyright by Kilian Q. Weinbergerr (2005)
% Washington University in St. Louis
% contact kilian@wustl.edu
%


fprintf('LMNN stable version 2.2\n');
if(nargin==0)
 help lmnn;
 return;
end;

if(~isempty(varargin) & isnumeric(varargin{1}))
  % check if neighborhood or L have been passed on
  Kg=varargin{1};
  fprintf('Setting neighborhood to k=%i\n',Kg);
  if(length(varargin)>1 & ~ischar(varargin{2}))
   L=varargin{2};
   fprintf('Setting initial transformation!\n');
  end;
  
  % skip Kgand L parameters
  newvarargin={};copy=0;j=1;
  for i=1:length(varargin)
    if(ischar(varargin{i})) copy=1;end;
    if(copy)newvarargin{j}=varargin{i};j=j+1;end;
  end;
  varargin=newvarargin;  
  clear('newvarargin','copy');
else
 fprintf('Neigborhood size not specified. Setting k=3\n');
 Kg=3;
end;


pars.outdim=size(x,1);
if(exist('L','var')~=1)
 fprintf(['Initial starting point not specified.\nStarting with PCA.\n']);  
 %L=eye(size(x,1));
 L=pca(x)';
end;
tic

% checks
D=size(L,2);
x=x(1:D,:);
if(size(x,1)>length(L)) error('x and L must have matching dimensions!\n');end;



 % set parameters
 pars.diagonal=0;
 pars.stepsize=1e-07;
 pars.minstepsize=0;
 pars.tempid=-1;
 pars.maxiter=1000;
 pars.factor=1.1;
 pars.correction=15;
 pars.thresho=1e-7;
 pars.thresha=1e-22;
 pars.ifraction=1;
 pars.scale=0;
 pars.obj=1;
 pars.quiet=0;
 pars.classsplit=0;
 pars.validation=0;
 pars.validationstep=25;
 pars.earlystopping=0;
 pars.valrand=1;
 
 
 pars.aggressive=0;
 pars.stepgrowth=1.01;
 pars.weight1=0.5;
 pars.maximp=100000;
 pars.maximp0=1000000;
 pars.treesize=50;
 pars.checkup=2; %0=notree 1=tree 2=choose
 
 pars.targetlabels=[];

 
pars=extractpars(varargin,pars);
if pars.diagonal, pars.obj=3;L=eye(size(L)); end;

L=L(1:pars.outdim,:);


% verification dataset
%i=randperm(size(x,2));
if(pars.validation<0 | pars.validation>1)
    error('validation parameter should be >0 and <1. Thanks.');
end;
earlycounter=0;
[itr,ite]=makesplits(y,1-pars.validation,1,pars.classsplit,Kg+1,pars.valrand);
xv=x(:,ite);
yv=y(:,ite);
x=x(:,itr);
y=y(itr);

if(~isempty(pars.targetlabels))
    pars.targetlabels=pars.targetlabels(itr);
end;

verify=[];besterr=inf;
clear('xo','yo');
lowesterr=inf;
verify=zeros(1,pars.maxiter);
bestL=L;

if(~pars.quiet)
pars
end;


% Initializationip
[D,N]=size(x);
fprintf('%i input vectors with %i dimensions\n',N,D);

[gen,NN]=getGenLS(x,y,Kg,pars);
obj=zeros(1,pars.maxiter);
nimp=zeros(1,pars.maxiter);

if(~pars.quiet) fprintf('Total number of genuine pairs: %i\n',size(gen,2));end;


%tic;[imp]= checkup(L,x,y,NN(end,:),pars);toc
 [imp]= checkup(L,x,y,NN(end,:),pars,1);
%all(all(imp==imp0))
%keyboard;
%t2=toc;
%keyboard;

if(size(imp,2)>pars.maximp0)
 ip=randperm(size(imp,2));
 ip=ip(1:pars.maximp0);
 imp=imp(:,ip);
end;




if(~pars.quiet)fprintf('Total number of imposture pairs: %i\n',size(imp,2));end;
 dfG=vec(SOD(x,gen(1,:),gen(2,:)));


if(pars.scale)
 Lx=L*x;
 sc=sqrt(mean(sum( ((Lx-Lx(:,NN(end,:)))).^2,1)));
 L=L./(sc+2);
end;

df=zeros(D^2,1);
correction=pars.correction;
ifraction=pars.ifraction;
stepsize=pars.stepsize;
lastcor=1;



for nnid=1:Kg; a1{nnid}=[];a2{nnid}=[];end;
df=zeros(size(dfG));

% Main Loop
for iter=1:pars.maxiter
 % save old position
 Lold=L;dfold=df;
 for nnid=1:Kg; a1old{nnid}=a1{nnid};a2old{nnid}=a2{nnid};end;
if(iter>1)L=step(L,mat((dfG.*pars.weight1+df.*(1-pars.weight1))),stepsize,pars);end;


if(~pars.quiet)fprintf('%i.',iter);end;

Lx=L*x;
%Lx2=sum(Lx.^2);
totalactive=0;

 

g0=cdist(Lx,imp(1,:),imp(2,:));


kk=1;
Ni=zeros(Kg,N);
for nnid=kk:Kg
 Ni(nnid,:)=(sum((Lx-Lx(:,NN(nnid,:))).^2,1)+1);
end;


 g1=Ni(:,imp(1,:)); 
 g2=Ni(:,imp(2,:)); 
 act1=[];act2=[];

if(pars.validation>0 & (mod(iter,pars.validationstep)==0 | iter==1))
 %verify=[verify Ltest2in(Lx,y,L*xv,yv,Ni,Kg,pars)];
 try
    verify(iter)=knnclassifytree([],Lx,y,L*xv,yv,Kg,'train',0);
 catch
     lasterr()
     keyboard;
 end;
 fprintf('kNN validation error: %2.2f ',verify(iter)*100);
 
 if(verify(iter)<=besterr) 
     fprintf('<= %2.2f   :-) %i/%i\n',besterr*100,earlycounter,pars.earlystopping);besterr=verify(iter);bestL=L;Det.bestiter=iter;
     earlycounter=0;
 else
     fprintf('> %2.2f   :-( %i/%i\n',besterr*100,earlycounter,pars.earlystopping);earlycounter=earlycounter+1;
 end;
 if(pars.earlystopping>0 & earlycounter>pars.earlystopping)
       fprintf('Validation error is no longer improving!\n');break;
 end;
end;
%clear('Lx','Lx2');
 
% objv=dfG'*vec((L'*L));
for nnid=Kg:-1:kk
  act1=find(g0<g1(nnid,:)); 
  act2=find(g0<g2(nnid,:)); 

 active=[act1 act2];

 if(~isempty(a1{nnid}) | ~isempty(a2{nnid}))
try
  [plus1,minus1]=sd(act1(:)',a1{nnid}(:)');
  [plus2,minus2]=sd(act2(:)',a2{nnid}(:)');
catch 
 lasterr
 keyboard;
end;
 else
  plus1=act1;plus2=act2;
  minus1=[];minus2=[];
 end;



% [isminus2,i]=sort(imp(1,minus2));minus2=minus2(i);
 MINUS1a=[imp(1,minus1) imp(2,minus2)]; MINUS1b=[imp(1,[plus1 plus2])];
 MINUS2a=[NN(nnid,imp(1,minus1)) NN(nnid,imp(2,minus2))]; MINUS2b=[imp(2,[plus1 plus2])];

 [isplus2,i]= sort(imp(2,plus2));plus2=plus2(i);
 PLUS1a=[imp(1,plus1) isplus2]; PLUS1b=[imp(1,[minus1 minus2])];
 PLUS2a=[NN(nnid,imp(1,plus1)) NN(nnid,isplus2)]; PLUS2b=[imp(2,[minus1 minus2])];

 loss1=max(g1(nnid,:)-g0,0);
 loss2=max(g2(nnid,:)-g0,0);
 

% ;
 [PLUS ,pweight]=count([PLUS1a;PLUS2a]);
 [MINUS,mweight]=count([MINUS1a;MINUS2a]);


df2=SODW(x,PLUS(1,:),PLUS(2,:),pweight)-SODW(x,MINUS(1,:),MINUS(2,:),mweight);

 df4=SOD(x,PLUS1b,PLUS2b)-SOD(x,MINUS1b,MINUS2b);
 df=df+vec(df2+df4);

 a1{nnid}=act1;a2{nnid}=act2;
 totalactive=totalactive+length(active);

end;



if(any(any(isnan(df))))
  fprintf('Gradient has NaN value!\n');
  keyboard;
end;

%obj(iter)=objv;
obj(iter)=(dfG.*pars.weight1+df.*(1-pars.weight1))'*vec(L'*L)+totalactive.*(1-pars.weight1);

if(isnan(obj(iter)))
 fprintf('Obj is NAN!\n');
 keyboard;
end;

nimp(iter)=totalactive;
delta=obj(iter)-obj(max(iter-1,1));
if(~pars.quiet)fprintf(['  Obj:%2.2f Nimp:%i Delta:%2.4f max(G):' ...
						 ' %2.4f' ...
		    '             \n   '],obj(iter),nimp(iter),delta,max(max(abs(df))));
end;



if(iter>1 & delta>0 & correction~=pars.correction) 
 stepsize=stepsize*0.5;
 fprintf('***correcting stepsize***\n');
 if(stepsize<pars.minstepsize) stepsize=pars.minstepsize;end;
 if(~pars.aggressive)
  L=Lold;
  df=dfold;
  for nnid=1:Kg; a1{nnid}=a1old{nnid};a2{nnid}=a2old{nnid};end;
  obj(iter)=obj(iter-1);
 end;
% correction=1;
 hitwall=1;
else 
  if(correction~=pars.correction)stepsize=stepsize*pars.stepgrowth;end;
 hitwall=0;
end;

if(iter>10)
 if (max(abs(diff(obj(iter-3:iter))))<pars.thresho*obj(iter)  | stepsize<pars.thresha)
  if(pars.correction-correction>=5) 
     correction=1;
  else
    switch(pars.obj)
     case 0
      if(~pars.quiet)fprintf('Stepsize too small. No more progress!\n');end;
      break;
     case 1
      pars.obj=0;
      pars.correction=15;
      pars.stepsize=1e-9;
      correction=0;
      for nnid=1:Kg; a1{nnid}=[];a2{nnid}=[];end;
      df=zeros(size(dfG));
      if(~pars.quiet | 1) 
        fprintf('\nVerifying solution! %i\n',obj(iter)); 
      end;
     case 3
      if(~pars.quiet)fprintf('Stepsize too small. No more progress!\n');end;
      break;
     end;
  end;
 end;
end;

correction=correction-1;
if(correction==0)
   if(pars.quiet)fprintf('\n');end;
   [Vio]=checkup(L,x,y,NN(nnid,:),pars);

   Vio=setdiff(Vio',imp','rows')';
   if(pars.maximp<inf)
     i=randperm(size(Vio,2));
     Vio=Vio(:,i(1:min(pars.maximp,size(Vio,2))));
   end;
  
   ol=size(imp,2);
    [imp i1 i2]=unique([imp Vio].','rows');
    imp=imp.';
    if(size(imp,2)~=ol)
      for nnid=1:Kg;
	   a1{nnid}=i2(a1{nnid});
	   a2{nnid}=i2(a2{nnid});
      end;
    end;

   
   fprintf('Iteration %i: Obj: %f  Added %i active constraints. #Active constraints: %i\n\n',iter,obj(iter),size(imp,2)-ol,size(imp,2));     
   if(ifraction<1)
     i=1:size(imp,2);
     imp=imp(:,i(1:ceil(length(i)*ifraction)));
     if(~pars.quiet)fprintf('Only use %2.2f of them.\n',ifraction);end;
     ifraction=ifraction+pars.ifraction;
   end;
  
  % put next correction a little more into the future if no new impostors were added
   if(size(imp,2)-ol<=0) 
      pars.correction=min(pars.correction*2+2,300);
      correction=pars.correction-1;
   else
     pars.correction=round(pars.correction*pars.factor);
     correction=pars.correction;
   end;
   lastcor=iter;
end;

end;


if iter==pars.maxiter, fprintf('MAXIMUM Number of iterations reached. Terminating without convergence.\n');end;
% Output
Det.obj=obj(1:iter);
Det.nimp=nimp(1:iter);
Det.pars=pars;
Det.time=toc;
Det.iter=iter;

Det.verify=verify;

if(pars.validation>0)
 Det.minL=L;
 L=bestL;  
 Det.verify=verify;
end;




function [err,yy,Value]=Ltest2in(Lx,y,LxT,yTest,Ni,Kg,pars);
% function [err,yy,Value]=Ltest2(L,x,y,xTest,yTest,Kg,varargin);
%


% Initializationip
[D,N]=size(Lx);

Lx2=sum(Lx.^2,1);

MM=min(y);
y=y-MM+1;
un=unique(y);
Value=zeros(length(un),length(yTest));


B=500;
NTe=size(LxT,2);
for n=1:B:NTe
  nn=n:n+min(B-1,NTe-n);
  DD=distance(Lx,LxT(:,nn));  
 for i=1:length(un)
 % Main Loopfor iter=1:pars.maxiter 
  testlabel=un(i);
  
  enemy=find(y~=testlabel);
  friend=find(y==testlabel);

  Df=mink(DD(friend,:),Kg);
%  Value(i,nn)=sumiflessv2(DD,Ni(:,enemy),enemy)+sumiflessh2(DD,Df+1,enemy);
  Value(i,nn)=sumiflessv2(DD,Ni(:,enemy),enemy)+sumiflessh2(DD,Df,enemy)+sum(Df,1);  
 end;

end;

 fprintf('\n');
 [temp,yy]=min(Value);

 yy=un(yy)+MM-1;
err=sum(yy~=yTest)./length(yTest);
fprintf('KNN error:%2.2f%%\n',err*100);






function L=step(L,G,stepsize,pars);

% do step in gradient direction
if(size(L,1)~=size(L,2)) pars.obj=1;end;
switch(pars.obj)
  case 0    % updating Q
     Q=L'*L;
     Q=Q-stepsize.*G;
   case 1   % updating L
     G=2.*(L*G);
     L=L-stepsize.*G;     
     return;
  case 2    % multiplicative update - this is broken
     Q=L'*L;
     Q=Q-stepsize.*G+stepsize^2/4.*G*inv(Q)*G;
     return;
  case 3
     Q=L'*L;
	 Q=Q-stepsize.*G;
	 Q=diag(Q);
 	 L=diag(sqrt(max(Q,0)));
     return;
  otherwise
   error('Objective function has to be 0,1,2\n');
end;

 
% decompose Q
[L,dd]=eig(Q);
dd=real(diag(dd));
L=real(L);
% reassemble Q (ignore negative eigenvalues)
j=find(dd<1e-10);
if(~isempty(j)) 
    if(~pars.quiet)fprintf('[%i]',length(j));end;
end;
dd(j)=0;
[temp,ii]=sort(-dd);
L=L(:,ii);
dd=dd(ii);
% Q=L*diag(dd)*L';
L=(L*diag(sqrt(dd)))';

%for i=1:size(L,1)
% if(L(i,1)~=0) L(i,:)=L(i,:)./sign(L(i,1));end;
%end;






function [gen,NN]=getGenLS(x,y,Kg,pars);
if(~pars.quiet);fprintf('Computing nearest neighbors ...\n');end; %#ok<SEPEX>
[D,N]=size(x);

if(~isempty(pars.targetlabels))
    y=pars.targetlabels;
end;
un=unique(y);
Gnn=zeros(Kg,N);
for c=un
 fprintf('%i nearest genuine neighbors for class %i:',Kg,c);
 i=find(y==c);
 nn=LSKnn(x(:,i),x(:,i),2:Kg+1);
 Gnn(:,i)=i(nn);
 fprintf('\r');
end;
fprintf('\n');
NN=Gnn;
gen1=vec(Gnn(1:Kg,:)')';
gen2=vec(repmat(1:N,Kg,1)')';
gen=[gen1;gen2];



function imp=checkup(L,x,y,NN,pars,~)
persistent treetime notreetime;
if(nargin==6)
    treetime=-1;
    notreetime=-1;
end;
fprintf('Updating working set.\n');
t1=toc;
if(pars.checkup==1 | (pars.checkup==2 & treetime<notreetime))
  imp=checkupmtree(L,x,y,NN,pars);treetime=toc-t1;
else
  imp=checkupnotree(L,x,y,NN,pars);notreetime=toc-t1;
end;
 

function imp=checkupmtree(L,x,y,NN,pars)
if(~pars.quiet);fprintf('[Tree] Computing nearest neighbors ...\n');end;
[D,N]=size(x);

mL=max(L');
L=L(find(mL),:);
Lx=L*x;
Ni=sum((Lx-Lx(:,NN)).^2,1)+1;
un=unique(y);

% build up ball trees
for c=1:length(un)
 classindex{c}=find(y==un(c));
 forest{c}.tree=buildmtreemex(Lx(:,classindex{c}),pars.treesize);
end;
imp=[];
for c=1:length(un)-1
if(~pars.quiet)fprintf('All impostors for class %i    \r',c);end;
 for c2=c+1:length(un)
     try
      limps=findNimex(forest{c2}.tree,Lx(:,classindex{c2}),Lx(:,classindex{c}),Ni(classindex{c2}),Ni(classindex{c}));         
     catch
         fprintf('The bizarre error happened!\n');
         fprintf('Check class index, c2 etc\n');
         fprintf('Line 629 in lmnn2\n');
         keyboard;
     end;
      
%    keyboard;
 	if(size(limps,2)>pars.maximp)
  	 ip=randperm(size(limps,2));
  	 ip=ip(1:pars.maximp);
  	 limps=limps(:,ip);
  	end;
	limps=[classindex{c}(limps(1,:));classindex{c2}(limps(2,:))];
    imp=[imp limps];
   end;
end;
try
 imp=unique(sort(imp)','rows')';
catch
 fprintf('Sorry, probably ran out of memory!');
 keyboard;  
end;



function imp=checkupnotree(L,x,y,NN,pars)
if(~pars.quiet) fprintf('Computing nearest neighbors ...\n');end;
[D,N]=size(x);

Lx=L*x;
Ni=sum((Lx-Lx(:,NN)).^2,1)+1;

un=unique(y);
imp=[];

%parfor c=un(1:end-1)
for c=un(1:end-1)
 if(~pars.quiet)fprintf('All nearest impostor neighbors for class %i :',c);end;
 i=find(y==c);
 index=find(y>c);
 limps=LSImps2(Lx(:,index),Lx(:,i),Ni(index),Ni(i),pars);
 if(size(limps,2)>pars.maximp)
  ip=randperm(size(limps,2));
  ip=ip(1:pars.maximp);
  limps=limps(:,ip);
 end;
 imp=[imp [i(limps(2,:));index(limps(1,:))]];

 if(~pars.quiet)fprintf('\r');end;
end;

try
 imp=unique(sort(imp)','rows')';
catch
 fprintf('Sorry, probably ran out of memory!');
 keyboard;  
end;



function limps=LSImps2(X1,X2,Thresh1,Thresh2,pars);
B=2000;
[D,N2]=size(X2);
N1=size(X1,2);
limps=[];
for i=1:B:N2
  BB=min(B,N2-i);
  try
%  newlimps=findimps3Dac(X1,X2(:,i:i+BB), Thresh1,Thresh2(i:i+BB));  % This line 
  newlimps=findimps3Dm(X1,X2(:,i:i+BB), Thresh1,Thresh2(i:i+BB));
  if(~isempty(newlimps) & newlimps(end)==0)    
    [minv,endpoint]=min(min(newlimps));
    newlimps=newlimps(:,1:endpoint-1);
  end;
  newlimps=unique(newlimps','rows')';
  catch
    lasterr
    keyboard;
    end;
  newlimps(2,:)=newlimps(2,:)+i-1;
  limps=[limps newlimps];
  if(~pars.quiet)fprintf('(%i%%) ',round((i+BB)/N2*100)); end;
end;
if(~pars.quiet)fprintf(' [%i] ',size(limps,2));end;




function NN=LSKnn(X1,X2,ks,pars);
B=750;
[D,N]=size(X2);
NN=zeros(length(ks),N);
DD=zeros(length(ks),N);

for i=1:B:N
  BB=min(B,N-i);
  fprintf('.');
  Dist=distance(X1,X2(:,i:i+BB));
  fprintf('.');
  [dist,nn]=mink(Dist,max(ks));
  clear('Dist');
  fprintf('.'); 
  NN(:,i:i+BB)=nn(ks,:);
  clear('nn','dist');
  fprintf('(%i%%) ',round((i+BB)/N*100)); 
end;
  



