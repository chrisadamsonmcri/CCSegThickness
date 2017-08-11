function [cts groupn] = quadratic_reg(xlsfiles, model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function cts = quadratic_reg(xlsfiles, model)
%This function is doing quadratic regression by using matlab function: regstats
%
%Inputs:
%         xlsfiles, eg. {'ctl.xls', 'ptn.xls'}
%         model, The optional input MODEL specifies how the design matrix is created 
%                from DATA
%                model can be any of the following strings:
%                   'linear'        Constant and linear terms (the default)
%                   'interaction'   Constant, linear, and interaction terms
%                   'quadratic'     Constant, linear, interaction, and squared terms
%                   'purequadratic' Constant, linear, and squared terms
%
%Outputs:
%         cts, a struct of t statistics of the regression 
%         groupn, a vector of number of samples in each group


if nargin < 2
    model = 'quadratic';
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

%EXTACT PREDICTORS FROM FIRST COLUMN OF THE INPUT MATRIX
predictors = Mat(:,1);
Mat(:,1) =[];


%ADD CALSS LAB IN FRIST COLUMN
Mat=[grouplab Mat];


slices = size(Mat, 2) - 1;  
n = size(Mat, 1); %NUMBER OF SAMPLES
m = size(Mat, 2); 


numbslices = size(Mat, 2) - 1;  

startrow = 1;
for i = 1: numGroups
    
    for j = 1:numbslices
        cts(i, j) = regstats(Mat(startrow:startrow+groupn(i)-1, j+1),  predictors(startrow:startrow+groupn(i)-1), model, {'tstat'});
    end
    startrow = startrow + groupn(i);
end
