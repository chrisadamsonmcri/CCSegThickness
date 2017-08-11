function anova_by_std (statsOfGroups)
%function anova_by_std (statsOfGroups)
%statsOfGroups: 
%               a struct with 3 fields
%               number of sample in group
%               mean of the group
%               stdandard deviation of the group
%Reference:
%http://ccnmtl.columbia.edu/projects/qmss/anova/the_oneway_anova.html
numG = length(statsOfGroups.means);
if numG ==2
    
xmeans = statsOfGroups.means;
xn = statsOfGroups.nums;
xstd = statsOfGroups.xstds;

%Calculate the Variation Between Groups
Bmean = sum(xmeans)/length(xmeans);

BSS = xn(1)*(xmeans(1) - Bmean).^2 + xn(2)*(xmeans(2) - Bmean).^2;
% BSS =sum(xn.*(xmeans-repmat(mean(xmeans),length(xmeans),1)).^2);

%This sum of squares has a number of degrees of freedom equal to the number of groups minus 1
dfB = numG-1;

%We divide the BSS figure by the number of degrees of freedom 
%to get our estimate of the variation between groups, referred to as "Between Mean Squares" as:
BMS = BSS/dfB;


%Calculate the Variation Within Groups
%This is a sum referred to as the "within sum of squares" or WSS.  In formula terms, this is expressed as:
WSS =  (xn(1)-1)*xstd(1).^2 + (xn(2)-1)*xstd(2).^2;
%WSS = sum((xn-1).*xstd.^2);

%we need to adjust the WSS to transform it into an estimate of population variance, 
%an adjustment that involves a value for the number of degrees of freedom within. 
%To calculate this, we take a value equal to the number of cases in the total sample (N), 
%minus the number of groups (k). In formula terms dfW = (N-k),
dfW = sum(xn) - numG; 

%Then we can calculate the a value for "Within Mean Squares" as
WMS = WSS/dfW;

%Calculate the F test statistic: F = BMS /WMS
F = BMS/WMS;

TSS = BSS + WSS;
end


