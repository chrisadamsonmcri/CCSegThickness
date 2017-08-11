function [t g1LinearP g2LinearP g1QuadP g2QuadP g1GTg2_linear g2GTg1_linear] = twoTailedOfBeta(cts, groupn, a)
%function twoTailedOfBeta(cts)
%
%twoTailedOfBeta does a two-tailed t-test on regression parameters from the
%output of function quadratic_reg
%
%Input: 
%       cts, a struct of t statistics of the regression 
%       groupn, a vector of number of samples in each group
%       a, alpha threshold of p


for i = 1: size(cts,2)
    % group 1 linear regression significance
    g1LinearP(i,1) = cts(1,i).tstat.pval(2);
    % group 1 quadratic regression significance
    g1QuadP(i,1) = cts(1,i).tstat.pval(3);
    
    % group 1 linear regression significance
    g2LinearP(i,1) = cts(2,i).tstat.pval(2);
    % group 1 quadratic regression significance
    g2QuadP(i,1) = cts(2,i).tstat.pval(3);
end

for i = 1:size(cts,2) %number of nodes
      for j = 1: 3    %number of parameters: intercept, linear, and quadratic
          
      beta1 = cts(1,i).tstat.beta(j); %parameter of group 1 
      beta2 = cts(2,i).tstat.beta(j); %parameter of group 2
      
      se1 = cts(1,i).tstat.se(j);  %standard error of group 1
      se2 = cts(2,i).tstat.se(j);  %standard error of group 2
      
      % t = (mean(X1) - mean(X2))/sqrt(se(X1).^2 + se(X2).^2)
      % Reference: Glantz, SA, 2002, Primer of Biostatistics, 5th ed.
      % McGraw-Hill, New York, p79.
      t(i, j) = (beta2 - beta1)/sqrt((se2.^2 + se1.^2)); 
      end 
end

df = sum(groupn)-2;

load ttable
col = find(ttable.a == a);
row = find(ttable.table(:,1)==df);

tthreshold = ttable.table(row, col+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Look for linear parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %find group 2 > group 1
g2GTg1_linear = find (t(:,2) > tthreshold);
if ~isempty(g2GTg1_linear)
    disp ('group1 great group 2 for linear parameter at:');
    for i = 1: length (g2GTg1_linear)
        disp(['Node ' num2str(g2GTg1_linear(i)) '    tP  ' num2str(t(g2GTg1_linear(i),2))]); 
    end
end


%find group 1 > group 2
g1GTg2_linear = find (-t(:,2) > tthreshold);
if ~isempty(g1GTg2_linear)
    disp ('group2 great group 1 for linear parameter at:');
    for i = 1: length (g1GTg2_linear)
        disp(['Node ' num2str(g1GTg2_linear(i)) '    tP  ' num2str(t(g1GTg2_linear(i),2))]); 
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Look for quadratic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%find group 2 > group 1
g2GTg1_quadratic= find (t(:,3) > tthreshold);
if ~isempty(g2GTg1_quadratic)
    disp ('group1 great group 2 for quadratic parameter at:');
    for i = 1: length (g2GTg1_quadratic)
        disp(['Node ' num2str(g2GTg1_quadratic(i)) '    tP  ' num2str(t(g2GTg1_quadratic(i),3))]); 
    end
end


%find group 1 > group 2
g2GTg1_quadratic = find (-t(:,3) > tthreshold);
if ~isempty(g2GTg1_quadratic)
    disp ('group2 great group 1 for quadratic parameter at:');
    for i = 1: length (g2GTg1_quadratic)
        disp(['Node ' num2str(g2GTg1_quadratic(i)) '    tP  ' num2str(t(g2GTg1_quadratic(i),3))]); 
    end
end
return