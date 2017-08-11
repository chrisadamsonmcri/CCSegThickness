function ANCOVAGps(xlsfiles, threshold, loop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function ANCOVAGps(xlsfiles, threshold, loop)
% e.g. ANCOVAGps({'xlsfile1', 'xlsfile2'},  0.05, 10e3)
%
% Input:
%       xlsfiles    xls files of groups with covariate in first column
%       threshold   default is 0.05 
%       loop        default is 10e3
%
%This function does ANCOVA with covaraite in first column
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%IS THERE A SIGNIFICANT EFFECT OF GROUP?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CREATE DESIGN MATRIX FOR FACTOR LEVEL REGRESSION MODEL WITHOUT INTERACTIONS
DesMat=[];

if numMatFiles == 2
    Des=[1 1 
         2 -1];
     
    ConMat=[1 -1; 
        -1 1];
end

if numMatFiles == 3
   Des=[1 1 0 
        2 0 1  
        3 -1 -1];
    
   ConMat = [1 -1 0;
         -1 1 0;
          0 1 -1;
          0 -1 1;
          1 0 -1;
          -1 0 1];
end

if numMatFiles == 4
    Des=[1 1 0 0
        2 0 1 0 
        3 0 0 1
        4 -1 -1 -1];
    
        
    ConMat=[1 -1 0 0;
        -1 1 0 0;
         0 1 -1 0;
         0 -1 1 0;
         1 0 -1 0;
         -1 0 1 0;
         1 0 0 -1;
         -1 0 0 1;
         0 1 0 -1;
         0 -1 0 1;
         0 0 1 -1;
         0 0 -1 1];
end


for count=1:size(Mat,1);
           Pos=find(Des(:,1)==Mat(count,1));
           DesMat=[DesMat; 1 Des(Pos,2:end)];
end
       

%APPEND COVARIATE AT LAST COLUMN OF DesMat
DesMat = [DesMat Mat(:,2)];


%FULL MODEL HAS ALL COEFFICIENTS
B=zeros(size(DesMat, 2), slices);

%REDUCED MODEL HAS TWO COEFFICIENT (NO COV INTERACTION)
reduce = size(DesMat,2)-(numMatFiles-1);
BRedGp=zeros(reduce,slices);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FOR UNEQUAL SAMPLE SIZES COMPARE FULL AND REDUCED MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FULL MODEL
InvF=pinv(DesMat);

%REDUCED MODEL (Testing for Effect of group)
InvRGp=pinv(DesMat(:,[1 end]));

for count=3:m
  %SSEF = Y'Y-b'X'Y OF FULL MODEL  
  B(:,count-2)=InvF*Mat(:,count);
  SSEF=(Mat(:,count))'*Mat(:,count)-(B(:,count-2)'*DesMat'*Mat(:,count));
  SSEFSt(count-2)=SSEF;

  %SSERGp = Y'Y-b'X'Y OF REDUCED MODEL
  BRedGp(:,count-2)=InvRGp*Mat(:,count);
  SSERGp=(Mat(:,count))'*Mat(:,count)-(BRedGp(:,count-2)'*DesMat(:,[1 end])'*Mat(:,count));

  %F* = MSR/MSE (chapter 16, p706)
  FGp(count-2, 1)=(SSERGp-SSEF)/SSEF;
end;



%%%%%%%%%%
%CONTRASTS
%%%%%%%%%%

if numMatFiles == 2
    Mn1 = B(2, :);
    Mn2 = -B(2, :);
    Mns = [Mn1; Mn2];
end

if numMatFiles == 3
    Mn1 = B(2, :);
    Mn2 = B(3, :);
    Mn3 = -B(2, :)-B(3, :);
    Mns = [Mn1; Mn2; Mn3];
end

if numMatFiles == 4;
    Mn1 = B(2, :);
    Mn2 = B(3, :);
    Mn3 = B(4, :);
    Mn4 = -B(2, :)-B(3, :)-B(4, :);
    Mns = [Mn1; Mn2; Mn3; Mn4];
end


nis = []; %inverse of 1/n s
for num = 1:numMatFiles
    nis = [nis; 1/groupn(num)];
end


%INFERENCE FOR DIFFERENCE BETWEEN FACTOR LEVEL MEANS (chapter 17, p741)
%e.g L=ConMat*[Mn1; Mn2];
%sL2 = MSE*sum(ci2/ni), t* = L/sL (chapter 17, p742)
% e.g. S2=(ConMat.^2)*[1/n1; 1/n2];

L=ConMat*Mns;
S2=(ConMat.^2)*nis;


%BECAUSE WE ONLY COMPARE RELATIVE t*, NOT WORRY ABOUT SQUARE AND DEGREE OF FREEDOM
t=L./(S2*SSEFSt); 


%%%%%%%%%%%%%%
%RANDOMIZATION
%%%%%%%%%%%%%%

MatR=zeros(size(Mat));
Fstore=zeros(slices,loop);
tst=zeros(size(t,1),slices,loop);

tic
for count=1:loop
    randn('state',sum(100*clock));
    [Rand,Sort]=sort(randn(n,1));
    MatR=[Mat(Sort,:)];
    BRand=zeros(size(DesMat, 2),slices);
    BRandRedGp=zeros(reduce,slices);
    BRedGp=zeros(reduce,slices);

    for countR=3:m
    BRand(:,countR-2)=InvF*MatR(:,countR);
    SSEF=(MatR(:,countR))'*MatR(:,countR)-(BRand(:,countR-2)'*DesMat'*MatR(:,countR));
    SSEFSt(countR-2)=SSEF;

    BRandRedGp(:,countR-2)=InvRGp*MatR(:,countR);
    SSERGp=(MatR(:,countR))'*MatR(:,countR)-(BRandRedGp(:,countR-2)'*DesMat(:,[1 end])'*MatR(:,countR));
    FRGp(countR-2,count)=(SSERGp-SSEF)/SSEF;
    end;
  

    %%%%%%%%%%
    %CONTRASTS
    %%%%%%%%%%

    if numMatFiles == 2
        Mn1 = B(2, :);
        Mn2 = -B(2, :);
        Mns = [Mn1; Mn2];
    end

    if numMatFiles == 3
      Mn1 = B(2, :);
      Mn2 = B(3, :);
      Mn3 = -B(2, :)-B(3, :);
      Mns = [Mn1; Mn2; Mn3];
    end

    if numMatFiles == 4;
      Mn1 = B(2, :);
        Mn2 = B(3, :);
        Mn3 = B(4, :);
        Mn4 = -B(2, :)-B(3, :)-B(4, :);
        Mns =[Mn1; Mn2; Mn3; Mn4];
    end

    L=ConMat*Mns;
    S2=(ConMat.^2)*nis;
    tst(:,:,count)=L./(S2*SSEFSt);
end; 
toc

%OMNIBUS NEW FOR EFFECT OF GROUP
disp('OMNIBUS NEW FOR EFFECT OF GROUP');
P=length(find(max(FRGp)>max(FGp)))/loop;

disp(['P ' num2str(P)]); 


if P < threshold 
    %CONTRASTS
    disp(sprintf('p value of t-test\n'));
      
    tSort=sort(tst,3);
    for count=1:size(tst, 1)
        twk=t(count,:);
        tSortwk=squeeze(tSort(count,:,:));
        tP(count)=length(find(max(tSortwk)>max(twk)))/loop;
        disp (['contrast ' num2str(ConMat(count, :)) '      p: ' num2str(tP(count))]);
        disp(sprintf('\n'));
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %STEP DOWN NEW FOR F

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(sprintf('\n'));
    disp('STEP DOWN FOR F');
    
    [FSort,Ind]=sort(FGp);
    FSort=flipud(FSort);
    Ind=flipud(Ind);
    FRGpSort=(flipud(sort(FRGp')))';
    Index=zeros(size(FGp));

    for count=1:length(FSort)
       Max=max(FRGp(find(Index==0),:));
       P=length(find(Max>FSort(count)))/loop;
       if P < threshold
         Index(count)=1;
         disp(['At slice ' num2str(Ind(count)) ' P is: ' num2str(P)]);
       end;
    end; 

    %Nodes=Ind(find(Index==1));
    disp(sprintf('\n'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %STEP DOWN NEW FOR t (CONTRASTS)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('STEP DOWN FOR t');

    sigtP = find(tP < threshold); %find significant tp

    if(~isempty(sigtP))
       for i = 1 : length(sigtP)
            disp(sprintf('\n'));
            disp(['Contrast: ' num2str(ConMat(sigtP(i), :))]);
                      
            twk=t(sigtP(i),:);
            tSortwk=squeeze(tSort(sigtP(i),:,:));
            Index=zeros(size(twk));
            Switch=1;

            while Switch==1
                if all(Index) == 1
                    break;
                end

                Pos = find(twk==max(twk));
                tPslice =length(find(max(tSortwk)>max(twk)))/loop;
                
                if tPslice < threshold
                    Index(Pos) = 1;
                    twk(Pos) = 0;
                    tSortwk(Pos,:)=zeros(size(tSortwk(Pos,:))); 
                    disp(['Node ' num2str(Pos) '    tP  ' num2str(tPslice)]); 
                else Switch=0;
                end;
            end;

        end

    else
        disp('No significance is found for step dwon t test');
    end

else
    disp ('No significance is found for omnibus test');
end

return



