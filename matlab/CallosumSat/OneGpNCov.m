function OneGpNCov(xlsfiles, numCov, threshold, loop)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function OneGpNCov(xlsfiles, numCov, threshold, loop)
% e.g. OneGpNCov({'xlsfile'}, 2, 0.05, 10e3)
%
% Input:
%       xlsfile     xls file of one group with covariates (response observations) in first n column 
%       numCov      number of covs
%       threshold   threshold of P, default is 0.05
%       loop        number of randomisation trials, default is 10e3
%
% This function does one group multiple regression with covariates (response observations) in first n column  
% This a generalisation of OneGpCov.m
% Jian Chen 23/04/07


if nargin < 3
    threshold = 0.05;
    loop = 10e3;
end

if nargin < 4
    loop = 10e3;
end


[pathstr, name, ext, versn] = fileparts(xlsfiles{1});
    
if strcmp(ext, '.xls')
    Mat = xlsread(xlsfiles{1});
else if strcmp(ext, '.txt')
    Mat = load(xlsfiles{1});
    else
        disp('Cheack input file format!');
        return
    end
end

sumMat = sum(Mat, 2);
posnan = find(isnan(sumMat));
Mat(posnan, :) =[];

n1=size(Mat, 1);

%Center covariates at overall mean
for i = 1:numCov 
    Mat(:, i) =Mat(:,i) - mean(Mat(:,i));
end

slices = size(Mat, 2) - numCov;  
n = size(Mat, 1); %NUMBER OF SAMPLES
m = size(Mat, 2); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IS THERE A SIGNIFICANT EFFECT OF COV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CREATE DESIGN MATRIX
DesMat = [ones(n,1) Mat(:,1:numCov)];

%FULL MODEL HAS numCov+1 COEFFICIENTS
B=zeros(size(DesMat,2),slices);

%REDUCED MODEL HAS ONE COEFFICIENT
BRedGp=zeros(size(DesMat,2)-numCov,slices);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE FULL AND REDUCED MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FULL MODEL
InvF=pinv(DesMat);

%REDUCED MODEL 
InvRGp=pinv(DesMat(:,[1 :end-numCov]));

for count=numCov+1:m
  %SSEF = Y'Y-b'X'Y OF FULL MODEL  
  B(:,count-numCov)=InvF*Mat(:,count);
 
  SSEF(count-numCov)=(Mat(:,count))'*Mat(:,count)-(B(:,count-numCov)'*DesMat'*Mat(:,count));
  SST(count-numCov) = sum((Mat(:, count) - repmat(mean(Mat(:, count)), n, 1)).^2);
  %SSEFSt(count-numCov)=SSEF(count-numCov);
  
  %SSERGp = Y'Y-b'X'Y OF REDUCED MODEL
  BRedGp(:,count-numCov)=InvRGp*Mat(:,count);
  SSERGp=(Mat(:,count))'*Mat(:,count)-(BRedGp(:,count-numCov)'*DesMat(:,[1:end-numCov])'*Mat(:,count));
  
  %F* = MSR/MSE (chapter 16, p706)
  FGp(count-numCov, 1)=(SSERGp-SSEF(count-numCov))/SSEF(count-numCov);
end;



%Calculate effect size, (Chapter 2, p74, Kutner & Neter, 5 ed.)
Rsqr = (1-SSEF./SST);
ESize = Rsqr./(1-Rsqr);

for count=1:length(ESize)
    disp(['At slice ' num2str(count) ' Effect size is: ' num2str(ESize(count))]);
end; 
	

clear SSEF;
MatR=zeros(size(Mat));
Fstore=zeros(slices,loop);

tic
for lcount=1:loop
  randn('state',sum(100*clock));
  [Rand,Sort]=sort(randn(n,1));
  MatR=[Mat(Sort,:)];
  BRand=zeros(size(DesMat,2),slices);
  BRandRedGp=zeros(size(DesMat,2)-numCov,slices);

  %FOR UNEQUAL SAMPLE SIZES COMPARE FULL AND REDUCED MODELS
  for countR=numCov+1:m
    BRand(:,countR-numCov)=InvF*MatR(:,countR);
    SSEF=(MatR(:,countR))'*MatR(:,countR)-(BRand(:,countR-numCov)'*DesMat'*MatR(:,countR));
    %SSEFSt(countR-numCov)=SSEF;
    
    BRandRedGp(:,countR-numCov)=InvRGp*MatR(:,countR);
    SSERGp=(MatR(:,countR))'*MatR(:,countR)-(BRandRedGp(:,countR-numCov)'*DesMat(:,[1:end-numCov])'*MatR(:,countR));
    FRGp(countR-numCov, lcount)=(SSERGp-SSEF)/SSEF;
  end;

end; 
toc

%OMNIBUS TEST FOR EFFECT OF GROUP
disp('OMNIBUS TEST FOR EFFECT OF COVARIATE P');
P=length(find(max(FRGp)>max(FGp)))/loop;
disp(P);


%%%%%%%%%%%%%%%%%%%%%%
     
%STEP DOWN TEST FOR F

%%%%%%%%%%%%%%%%%%%%%%
if P < threshold
    
	%%%%%%%%%%%%%%%%%%%%%%
         
	%STEP DOWN TEST FOR F
	
	%%%%%%%%%%%%%%%%%%%%%%
	disp('STEP DOWN TEST FOR F')
	
	[FSort,Ind]=sort(FGp);
	FSort=flipud(FSort);
	Ind=flipud(Ind);
	FRGpSort=(flipud(sort(FRGp')))';
	Index=zeros(size(FGp));
	
	for count=1:length(FSort)
       Max=max(FRGp(find(Index==0),:));
       P=length(find(Max>FSort(count)))/loop;
       if P<threshold
         Index(count)=1;
         
         disp(['At slice ' num2str(Ind(count)) ' P is: ' num2str(P) ]);
         disp(['At slice ' num2str(Ind(count)) ' B is: ' num2str(B(2:end, Ind(count))')]);
     
       end;
	end; 
	
	Nodes=Ind(find(Index==1));
	 
else
    disp ('No significance is found for omnibus test');
end

return


