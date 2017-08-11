function OneGpCov(xlsfiles, threshold, loop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function OneGpCov(xlsfile, threshold, loop)
% e.g. OneGpCov({'xlsfile'}, 0.05, 10e3)
%
% Input:
%       xlsfile     xls file of group 1 with covariate (response observations) in first column              
%       threshold   threshold of P, default is 0.05
%       loop        number of randomisation trials, default is 10e3
%
% This function does one group simple regression with covariate (response observations) in first column  
%
% Jian Chen 23/04/07


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

%Center ageonset at overall mean
Mat(:, 1) =Mat(:,1) - mean(Mat(:,1));
slices = size(Mat, 2) - 1;  
n = size(Mat, 1); %NUMBER OF SAMPLES
m = size(Mat, 2); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IS THERE A SIGNIFICANT EFFECT OF COV, ie. AGE?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CREATE DESIGN MATRIX OF FULL MODEL
DesMat = [ones(n,1) Mat(:,1)];

%FULL MODEL HAS 2 COEFFICIENTS
B=zeros(size(DesMat,2),slices);

%REDUCED MODEL HAS ONE COEFFICIENT (NO COV EFFECT)
BRedGp=zeros(size(DesMat,2)-1,slices);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE FULL AND REDUCED MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FULL MODEL
InvF=pinv(DesMat);

%REDUCED MODEL 
InvRGp=pinv(DesMat(:,[1 :end-1]));

for count=2:m
  %SSEF = Y'Y-b'X'Y OF FULL MODEL  
  B(:,count-1) = InvF*Mat(:,count);
 
  SSEF(count-1) = (Mat(:,count))'*Mat(:,count)-(B(:,count-1)'*DesMat'*Mat(:,count));% error sum squares
  SST(count-1) = sum((Mat(:, count) - repmat(mean(Mat(:, count)), n, 1)).^2);
  %SSEFSt(count-1)=SSEF(count-1);
  
  %SSERGp = Y'Y-b'X'Y OF REDUCED MODEL
  BRedGp(:,count-1)=InvRGp*Mat(:,count);
  SSERGp=(Mat(:,count))'*Mat(:,count)-(BRedGp(:,count-1)'*DesMat(:,[1:end-1])'*Mat(:,count));
  
  %F* = MSR/MSE (chapter 16, p706)
  FGp(count-1, 1)=(SSERGp-SSEF(count-1))/SSEF(count-1);
end;


%Calculate effect size
Rsqr = 1 - (SSEF./SST);
ESize = Rsqr./(1-Rsqr);

for count=1:length(ESize)
    disp(['At slice ' num2str(count) ' Effect size is: ' num2str(ESize(count))]);
end; 



clear SSEF;

MatR=zeros(size(Mat));
Fstore=zeros(slices,loop);

tic
for count=1:loop
  randn('state',sum(100*clock));
  [Rand,Sort]=sort(randn(n,1));
  MatR=[Mat(Sort,:)];
  BRand=zeros(size(DesMat,2),slices);
  BRandRedGp=zeros(size(DesMat,2)-1,slices);

  %FOR UNEQUAL SAMPLE SIZES COMPARE FULL AND REDUCED MODELS
  for countR=2:m
    BRand(:,countR-1)=InvF*MatR(:,countR);
    SSEF=(MatR(:,countR))'*MatR(:,countR)-(BRand(:,countR-1)'*DesMat'*MatR(:,countR));
    %SSEFSt(countR-1)=SSEF;
    
    BRandRedGp(:,countR-1)=InvRGp*MatR(:,countR);
    SSERGp=(MatR(:,countR))'*MatR(:,countR)-(BRandRedGp(:,countR-1)'*DesMat(:,[1:end-1])'*MatR(:,countR));
    FRGp(countR-1,count)=(SSERGp-SSEF)/SSEF;
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
	
	[FSort,Ind]=sort(FGp, 'descend');
	FRGpSort= (sort(FRGp', 'descend'))';
    
	Index=zeros(size(FGp));
	
	for count=1:length(FSort)
       Max =max(FRGp(find(Index==0),:));
       P=length(find(Max > FSort(count)))/loop;
       if P < threshold
         Index(count)=1;
         
         disp(['At slice ' num2str(Ind(count)) ' P is: ' num2str(P) ]);
         disp(['At slice ' num2str(Ind(count)) ' B is: ' num2str(B(2, Ind(count)))]);
     
       end;
	end; 
	
	Nodes=Ind(find(Index==1));
	 
else
    disp ('No significance is found for omnibus test');
end



%%%%%%%%%%%%%%%%%%%%%%%

%No multiple comparison

%%%%%%%%%%%%%%%%%%%%%%%

disp ('If no multiple comparison, then: ');
for i = 1 : slices
    p(i) = length(find(FRGp(i, :) > FGp(i)))/loop;
    if p(i) < threshold 
        disp(['At slice ' num2str(i) ' p is: ' num2str(p(i)) ]);
    end
end

return


