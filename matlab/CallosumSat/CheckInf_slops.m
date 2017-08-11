
function CheckInf_slops(xlsfiles, threshold, loop)
%function CheckInf_slop(xlsfiles,  threshold, loop)
% e.g. CheckInf_slop({'xlsfile1', 'xlsfile2'}, 0.05, 10e3)
%
% Kutner & Neter, 5 ed, chapter 22, p932
%
% Input:
%       xlsfiles    group's xls files with first column is a cov, i.e. age
%       threshold   default is 0.05 
%       loop        default is 10e3
%
% 
% Jian Chen 23/04/07

if nargin < 2
    threshold = 0.05;
    loop = 10e3;
end

if nargin < 3
    loop = 10e3;
end

numMatFiles = length(xlsfiles);

Mat =[];
grouplab = [];

for num = 1:numMatFiles
        
    [pathstr, name, ext, versn] = fileparts(xlsfiles{num});
    
    if strcmp(ext, '.xls')
        matdata = xlsread(xlsfiles{num});
    else if strcmp(ext, '.txt')
        matdata = load(xlsfiles{num});
        else
            disp('Cheack input file format!');
            return
        end
    end
    
    sumMat = sum(matdata, 2);
    posnan = find(isnan(sumMat));
    matdata(posnan, :) =[];
    
    groupn(num) = size(matdata, 1);
    try
        Mat = [Mat; matdata];
    catch
        disp('Error: You are xls file should have same number of columns, please check!');
        return
    end
    grouplab =[grouplab; num*ones(groupn(num),1)];
end


%ADD CALSS LAB IN FRIST COLUMN
Mat=[grouplab Mat];

%Center covariate at overall mean
Mat(:, 2) =Mat(:,2) - mean(Mat(:,2));

slices = size(Mat, 2) - 2;  
n = size(Mat, 1); %NUMBER OF SAMPLES
m = size(Mat, 2); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MATRIX DESIGN ACCORDING TO NETER'S
%5TH EDITION, PAGE 932
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DesMat=[];

if numMatFiles == 2
    Des=[1 1 
         2 -1];
end

if numMatFiles == 3
   Des=[1 1 0 
        2 0 1  
        3 -1 -1];
end

if numMatFiles == 4
    Des=[1 1 0 0
        2 0 1 0 
        3 0 0 1
        4 -1 -1 -1];
end


for count=1:size(Mat,1);
           Pos=find(Des(:,1)==Mat(count,1));
           DesMat=[DesMat; 1 Des(Pos,2:end)];
end
       

%APPEND COVARIATE AT LAST COLUMN OF DesMat
DesMat = [DesMat Mat(:,2)];
 

%APPEND CROSS-PRODUCTION INTERACTION TERMS, Chapter 22, p 932
groups = size(Des, 2)-1;

for intfer = 1:groups
    DesMat = [DesMat Mat(:,2).*DesMat(:, 1+intfer)];
end


%FULL MODEL HAS 4 COEFFICIENTS
B=zeros(size(DesMat, 2), slices);

%REDUCED MODEL HAS 3 COEFFICIENT 
BRedGp=zeros(size(DesMat, 2)-groups, slices);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FOR UNEQUAL SAMPLE SIZES COMPARE FULL AND REDUCED MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FULL MODEL
InvF=pinv(DesMat);

%REDUCED MODEL HAS NO INTERFACTION COEFFICIENT 
InvRGp=pinv(DesMat(:,[1:end-groups]));

for count=3:m
  %SSEF = Y'Y-b'X'Y OF FULL MODEL  
  B(:,count-2)=InvF*Mat(:,count);
  SSEF=(Mat(:,count))'*Mat(:,count)-(B(:,count-2)'*DesMat'*Mat(:,count));
  SSEFSt(count-2)=SSEF;

  %SSERGp = Y'Y-b'X'Y OF REDUCED MODEL
  BRedGp(:,count-2)=InvRGp*Mat(:,count);
  SSERGp=(Mat(:,count))'*Mat(:,count)-(BRedGp(:,count-2)'*DesMat(:,[1:end-groups])'*Mat(:,count));

  %F* = MSR/MSE (chapter 16, p706)
  FGp(count-2, 1)=(SSERGp-SSEF)/SSEF;
end;


MatR=zeros(size(Mat));
Fstore=zeros(slices,loop);

tic
for count=1:loop
  %count

  randn('state',sum(100*clock));
  [Rand,Sort]=sort(randn(n,1));
  MatR=[Mat(Sort,:)];
  BRand=zeros(size(DesMat, 2),slices);
  BRandRedGp=zeros(size(DesMat, 2)-groups, slices);


  %FOR UNEQUAL SAMPLE SIZES COMPARE FULL AND REDUCED MODELS
  for countR=3:m
    BRand(:,countR-2)=InvF*MatR(:,countR);
    SSEF=(MatR(:,countR))'*MatR(:,countR)-(BRand(:,countR-2)'*DesMat'*MatR(:,countR));
    SSEFSt(countR-2)=SSEF;
    BRandRedGp(:,countR-2)=InvRGp*MatR(:,countR);
    SSERGp=(MatR(:,countR))'*MatR(:,countR)-(BRandRedGp(:,countR-2)'*DesMat(:,[1:end-groups])'*MatR(:,countR));
    FRGp(countR-2,count)=(SSERGp-SSEF)/SSEF;
  end;

end; 
toc

%OMNIBUS TEST FOR EFFECT OF COVARATE
disp('%OMNIBUS TEST FOR EFFECT OF COVARATE');
P=length(find(max(FRGp)>max(FGp)))/loop;
disp(P);



%%%%%%%%%%%%%%%%%%%%%

%STEP DOWN TEST FOR F

%%%%%%%%%%%%%%%%%%%%%
if P < threshold
    disp('THERE IS A SIGNIFICANCE FOR INTERACTION');
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
          disp(['At slice ' num2str(Ind(count)) ' P is: ' num2str(P)]);
       end;
    end; 
    Nodes=Ind(find(Index==1));

else
    disp('THERE IS NO SIGNIFICANCE DIFFERENCE BETWEEN THE SLOPS (INTERACTION)');
end
