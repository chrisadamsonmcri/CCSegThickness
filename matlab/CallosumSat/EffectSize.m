function EffectSize(xlsfiles, covflag )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function EffectSize(xlsfiles, covflag )
% e.g. EffectSize({'xlsfile1', 'xlsfile2'}, covflag)
%
% Input:
%       xlsfiles    group's xls files
%       covflag     if 0 no cov at first column and 1 first column is a cov
%
% Jian Chen 21/02/08
%
% This function calculates effective size of TWO GROUPS at each node using Cohen's d
% d = (mean1-mean2)/sqrt((SD1.^2 + SD2.^2)/2)



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

if covflag == 1
    Mat(:,1) =[];
end


%ADD CALSS LAB IN FRIST COLUMN
Mat=[grouplab Mat];

slices = size(Mat, 2) - 1;  
n = size(Mat, 1); %NUMBER OF SAMPLES
m = size(Mat, 2); 


Mns = [];
Sds = [];

for num = 1:numMatFiles
    Mn = mean(Mat(find(Mat(:,1)==num),2:m));
    Mns = [Mns; Mn];
    Sd = std(Mat(find(Mat(:,1)==num),2:m));
    Sds = [Sds; Sd];
end

if numMatFiles == 2
    effsize = abs(Mns(1, :)-Mns(2, :))./sqrt((Sds(1,:).^2 + Sds(2,:).^2)/2);
    
end

disp (['node    effsize']);
for sliceNum = 1: slices;
    disp([num2str(sliceNum)   '     ' num2str(effsize(sliceNum))]);
end
