
function CheckInf_2factors(xlsfiles, threshold, loop)
%function CheckInf_2factors(xlsfiles,  threshold, loop)
% e.g. CheckInf_2factors({'ctl-drd3-glycarrier.xls', 'ctl-drd3-serhomo.xls', 'scz-drd3-glycarrier.xls', 'scz-drd3-serhomo.xls'}, 0.2, 10e5)
%
% Kutner & Neter, 5 ed, chapter 22, p932
% Input:
%       xlsfiles    group's xls files 
%       threshold   default is 0.05 
%       loop        default is 10e3
%
% 
% Jian Chen 23/05/09

if nargin < 2
    threshold = 0.05;
    loop = 10e3;
end

if nargin < 3
    loop = 10e3;
end

numMatFiles = length(xlsfiles);

Mat =[];

gID = [11 12 21 22];
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
    grouplab =[grouplab; gID(num)*ones(groupn(num),1)];
end


%ADD CALSS LAB IN FRIST COLUMN
Mat=[grouplab Mat];

slices = size(Mat, 2) - 1;  
n = size(Mat, 1); %NUMBER OF SAMPLES
m = size(Mat, 2); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MATRIX DESIGN ACCORDING TO NETER'S
%5TH EDITION, PAGE 932
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DesMat=[];

Des=[11 1  1  1
     12 1 -1 -1
     21 -1 1 -1
     22 -1 -1 1];


for count=1:size(Mat,1);
           Pos=find(Des(:,1)==Mat(count,1));
           DesMat=[DesMat; 1 Des(Pos,2:end)];
end
       

%APPEND COVARIATE AT LAST COLUMN OF DesMat
%DesMat = [DesMat Mat(:,2)];
 

%FULL MODEL HAS 4 COEFFICIENTS
B=zeros(size(DesMat, 2), slices);

%REDUCED MODEL HAS 3 COEFFICIENT 
BRedGp=zeros(size(DesMat, 2)-1, slices);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FOR UNEQUAL SAMPLE SIZES COMPARE FULL AND REDUCED MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FULL MODEL
InvF=pinv(DesMat);

%REDUCED MODEL HAS NO INTERFACTION COEFFICIENT 
InvRGp=pinv(DesMat(:,[1:end-1]));

skip = 1;

for count=skip+1:m
  %SSEF = Y'Y-b'X'Y OF FULL MODEL  
  B(:,count-skip)=InvF*Mat(:,count);
  SSEF=(Mat(:,count))'*Mat(:,count)-(B(:,count-skip)'*DesMat'*Mat(:,count));
  SSEFSt(count-skip)=SSEF;

  %SSERGp = Y'Y-b'X'Y OF REDUCED MODEL
  BRedGp(:,count-skip)=InvRGp*Mat(:,count);
  SSERGp=(Mat(:,count))'*Mat(:,count)-(BRedGp(:,count-skip)'*DesMat(:,[1:end-1])'*Mat(:,count));

  %F* = MSR/MSE (chapter 16, p706)
  FGp(count-skip, 1)=(SSERGp-SSEF)/SSEF;
end;


MatR=zeros(size(Mat));
Fstore=zeros(slices,loop);

tic
for lp =1:loop

  randn('state',sum(100*clock));
  [Rand,Sort]=sort(randn(n,1));
  MatR=[Mat(Sort,:)];
  BRand=zeros(size(DesMat, 2),slices);
  BRandRedGp=zeros(size(DesMat, 2)-1, slices);


  %FOR UNEQUAL SAMPLE SIZES COMPARE FULL AND REDUCED MODELS
  for countR=skip+1:m
    BRand(:,countR-skip)=InvF*MatR(:,countR);
    SSEF=(MatR(:,countR))'*MatR(:,countR)-(BRand(:,countR-skip)'*DesMat'*MatR(:,countR));
    SSEFSt(countR-skip)=SSEF;
    BRandRedGp(:,countR-skip)=InvRGp*MatR(:,countR);
    SSERGp=(MatR(:,countR))'*MatR(:,countR)-(BRandRedGp(:,countR-skip)'*DesMat(:,[1:end-1])'*MatR(:,countR));
    FRGp(countR-skip,lp)=(SSERGp-SSEF)/SSEF;
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
