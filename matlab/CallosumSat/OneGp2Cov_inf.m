function OneGp2Cov_inf(xlsfiles, threshold, loop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function OneGp2Cov_ide(xlsfile, threshold, loop)
% e.g. OneGp2Cov_inf({'xlsfile'}, 0.05, 10e3)
% This program test if the regression model yields same responds function,
% i.e. no interaction between ordinal cov and quatitative cov
%
% Kutner & Neter, 5 ed, chapter 8, p326
%
% Input:
%       xlsfile     xls file of group 1 with 2 covariates, first column is qualititative variable, 
%                   second column is quantitative variable
%       threshold   threshold of P, default is 0.05
%       loop        number of randomisation trials, default is 10e3
%
% Jian Chen 01/09/08


if nargin < 2
    threshold = 0.05;
    loop = 10e3;
end

if nargin < 3
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

numCov = 2; 

%First cov is qualitative
%Center second quantitative cov at overall mean
Mat(:, numCov) = Mat(:, numCov) - mean(Mat(:, numCov));

slices = size(Mat, 2) -numCov;  
n = size(Mat, 1); %NUMBER OF SAMPLES
m = size(Mat, 2); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%IS THERE A SIGNIFICANT EFFECT OF GROUP?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE DESIGN MATRIX FOR REGRESSION MODEL WITH INTERACTION

DesMat = [ones(n,1) Mat(:,1:numCov), Mat(:,1).*Mat(:,numCov)];

%FULL MODEL HAS 4 COEFFICIENTS
B=zeros(size(DesMat,2),slices);

%REDUCED MODEL HAS THREE COEFFICIENT (NO INTERACTION)
%(chapter 8, p326)
BRedGp=zeros(size(DesMat,2)-1,slices);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE FULL AND REDUCED MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FULL MODEL
InvF=pinv(DesMat);

%REDUCED MODEL (Testing for Effect of Interaction)
InvRGp=pinv(DesMat(:,[1:end-1]));

for count=numCov+1:m
  %SSEF = Y'Y-b'X'Y OF FULL MODEL  
  B(:,count-numCov) = InvF*Mat(:, count);
 
  SSEF(count-numCov) = (Mat(:,count))'*Mat(:,count)-(B(:,count-numCov)'*DesMat'*Mat(:,count));% error sum squares
  SST(count-numCov) = sum((Mat(:, count) - repmat(mean(Mat(:, count)), n, 1)).^2);
  SSEFSt(count-numCov)=SSEF(count-numCov);
  
  %SSERGp = Y'Y-b'X'Y OF REDUCED MODEL
  BRedGp(:,count-numCov)=InvRGp*Mat(:,count);
  SSERGp=(Mat(:,count))'*Mat(:,count)-(BRedGp(:,count-numCov)'*DesMat(:,[1:end-1])'*Mat(:,count));
  
  %F* = MSR/MSE (chapter 16, p706)
  FGp(count-numCov, 1)=(SSERGp-SSEF(count-numCov))/SSEF(count-numCov);
end;

%Calculate effect size
%Rsqr = 1-(SSEF./SST);
%ESize = Rsqr./(1-Rsqr);


MatR=zeros(size(Mat));
Fstore=zeros(slices,loop);

tic
for loopcount =1:loop
  randn('state',sum(100*clock));
  [Rand,Sort]=sort(randn(n,1));
  MatR=[Mat(Sort,:)];
  BRand=zeros(size(DesMat,2),slices);
  BRandRedGp=zeros(size(DesMat,2)-1,slices);

  %FOR UNEQUAL SAMPLE SIZES COMPARE FULL AND REDUCED MODELS
  for countR=numCov+1:m
    BRand(:,countR-numCov)=InvF*MatR(:,countR);
    SSEF=(MatR(:,countR))'*MatR(:,countR)-(BRand(:,countR-numCov)'*DesMat'*MatR(:,countR));
    SSEFSt(countR-numCov)=SSEF;
    
    BRandRedGp(:,countR-numCov)=InvRGp*MatR(:,countR);
    SSERGp=(MatR(:,countR))'*MatR(:,countR)-(BRandRedGp(:,countR-numCov)'*DesMat(:,[1:end-1])'*MatR(:,countR));
    FRGp(countR-numCov, loopcount)=(SSERGp-SSEF)/SSEF;
  end;

end; 
toc

%OMNIBUS TEST FOR EFFECT OF GROUP
disp('OMNIBUS TEST FOR EFFECT OF INTERACTION P');
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
       Max= max (FRGp(find(Index==0),:));
       P=length(find(Max > FSort(count)))/loop
       if P < threshold
         Index(count)=1;
         disp(['At slice ' num2str(Ind(count)) ' P is: ' num2str(P) ]);   
       end;
	end; 
	
	Nodes=Ind(find(Index==1));
	 
else
    disp ('No significance is found for omnibus test');
end

return


