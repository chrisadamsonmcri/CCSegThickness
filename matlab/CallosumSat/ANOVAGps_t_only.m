function ANOVAGps_t_only(xlsfiles, covflag,  threshold, loop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function ANOVAGps(xlsfiles, covflag, threshold, loop)
% e.g. ANOVAGps({'xlsfile1', 'xlsfile2'}, 0, 0.05, 10e3)
%
% Input:
%       xlsfiles    group's xls files
%       covflag     indicate whether there is a covariate at first column
%                   (1 with covariate; 0 without)
%       threshold   default is 0.05 
%       loop        default is 10e3
%
% This function does ANOVA, if your input files have covariate in first
% column, assgin covflag = 1, otherwise assgin covflag = 0
%
% This function is simliar to ANOVAGps but without testing omnibus, the aim
% is to look at t-test value only, even the analysis fail omnibus test
%
% Jian Chen 23/04/07

if nargin < 3
    threshold = 0.05;
    loop = 10e3;
end

if nargin < 4
    loop = 10e3;
end

numMatFiles = length(xlsfiles);

Mat =[];
grouplab = [];

numGroups = numMatFiles;
for num = 1:numGroups
    
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

if covflag == 1
    Mat(:,1) =[];
end


%ADD CALSS LAB IN FRIST COLUMN
Mat=[grouplab Mat];


slices = size(Mat, 2) - 1;  
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
       
    
%FULL MODEL HAS ALL COEFFICIENTS
B=zeros(size(DesMat, 2), slices);

%REDUCED MODEL HAS ONE COEFFICIENT (NO MEAN DIFFERENCE)
BRedGp=zeros(1,slices);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FOR UNEQUAL SAMPLE SIZES COMPARE FULL AND REDUCED MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FULL MODEL
InvF=pinv(DesMat);

%REDUCED MODEL (Testing for Effect of group)
InvRGp=pinv(DesMat(:,[1]));

for count=2:m
  %SSEF = Y'Y-b'X'Y OF FULL MODEL  
  B(:,count-1)=InvF*Mat(:,count);
  SSEF=(Mat(:,count))'*Mat(:,count)-(B(:,count-1)'*DesMat'*Mat(:,count));
  SSEFSt(count-1)=SSEF;

  %SSERGp = Y'Y-b'X'Y OF REDUCED MODEL
  BRedGp(:,count-1)=InvRGp*Mat(:,count);
  SSERGp=(Mat(:,count))'*Mat(:,count)-(BRedGp(:,count-1)'*DesMat(:,[1])'*Mat(:,count));

  %F* = MSR/MSE (chapter 16, p706)
  FGp(count-1,1)=(SSERGp-SSEF)/SSEF;
end;




%%%%%%%%%%
%CONTRASTS
%%%%%%%%%%
Mns = [];
nis = [];
for num = 1:numMatFiles
    Mn = mean(Mat(find(Mat(:,1)==num),2:m));
    Mns = [Mns; Mn];
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
  %count

  randn('state',sum(100*clock));
  [Rand,Sort]=sort(randn(n,1));
  MatR=[Mat(Sort,:)];
  BRand=zeros(size(DesMat, 2),slices);
  BRandRedGp=zeros(1,slices);

  %FOR UNEQUAL SAMPLE SIZES COMPARE FULL AND REDUCED MODELS
  for countR = 2:m
    BRand(:,countR-1)=InvF*MatR(:,countR);
    SSEF=(MatR(:,countR))'*MatR(:,countR)-(BRand(:,countR-1)'*DesMat'*MatR(:,countR));
    SSEFSt(countR-1)=SSEF;

    BRandRedGp(:,countR-1)=InvRGp*MatR(:,countR);
    SSERGp=(MatR(:,countR))'*MatR(:,countR)-(BRandRedGp(:,countR-1)'*DesMat(:,[1])'*MatR(:,countR));
    FRGp(countR-1,count)=(SSERGp-SSEF)/SSEF;
  end;

      
    %%%%%%%%%%
    %CONTRASTS
    %%%%%%%%%%
    Mns = [];
    for num = 1:numMatFiles
        Mn= mean(MatR(find(Mat(:,1)==num),2:m));
        Mns = [Mns; Mn];
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


% if P < threshold 
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

%     %STEP DOWN FOR F
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp(sprintf('\n'));
%     disp('STEP DOWN FOR F');
%     
%     [FSort,Ind]=sort(FGp);
%     FSort=flipud(FSort);
%     Ind=flipud(Ind);
%     FRGpSort=(flipud(sort(FRGp')))';
%     Index=zeros(size(FGp));
% 
%     for count=1:length(FSort)
%        Max=max(FRGp(find(Index==0),:));
%        P=length(find(Max>FSort(count)))/loop;
%        if P < threshold
%          Index(count)=1;
%          disp(['At slice ' num2str(Ind(count)) ' P is: ' num2str(P)]);
%        end;
%     end; 

    %Nodes=Ind(find(Index==1));
    disp(sprintf('\n'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %STEP DOWN FOR t (CONTRASTS)

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
% 
% else
    disp ('No significance is found for omnibus test');
end





