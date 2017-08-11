function [P_at_Node,omniP]=Callosum_PairedT_Jian(Gp1FileName,Gp2FileName, ncov)

%This program will do a randomisation for paired t tests
%It requires data to be in two arrays one for each group
%The pairing is for columns in each array i.e.
%Column 1 in Array1 should be paired to Column 1 in Array2 
%Ptrue are the p values at each node
%omniP is the omnibus p value 
%USAGE: [P_at_Node,omniP]=Callosum_PairedT_Jian(Gp1,Gp2, ncov) 
%
%This was made change based on Callosum_Paired.m to make it works
%Jian Chen 15/06/10


Gp1=load(Gp1FileName);
Gp2=load(Gp2FileName);

if ncov ~=0
Gp1(:,1:ncov) =[];
Gp2(:,1:ncov) =[];
end

Diff=(Gp1-Gp2);
T=mean(Diff)./std(Diff);
[Tsort,iT]=sort(T);

%RANDOMIZATION
RandStore=[];
IndStore=[];
RandList=[1:size(Gp1,2)];

for count=1:100000
  rand('state',sum(100*clock));
  Rand=rand(2,size(Gp1,2)); 
  Pos=find(diff(Rand)>=0);
  Index=sum((Pos).^4);
  if count==1
    IndStore=[IndStore;Index];
    RandStore=[RandStore; Rand(1,:)>=Rand(2,:)];
  
   % if we havent already had this randomisation
   elseif min(min(abs(IndStore-Index)))>0
    IndStore=[IndStore;Index];
    RandStore=[RandStore; Rand(1,:)>=Rand(2,:)];
  end;  
end

%RANDOMISATION DISTRIBUTION
tmax=[];
Ci=zeros(size(T));
for count=1:size(RandStore,1)
  RandGp1=Gp1;
  RandGp2=Gp2;
  RandGp1(:,find(RandStore(count,:)==0))=Gp2(:,find(RandStore(count,:)==0));
  RandGp2(:,find(RandStore(count,:)==0))=Gp1(:,find(RandStore(count,:)==0));
  Diff=(RandGp1-RandGp2);
  t=mean(Diff)./std(Diff);
  tmax=[tmax;max(max(t))];

  %STEP DOWN TEST
  vi(1)=t(iT(1));
  for count2=2:length(t)
    vi(count2)=max([vi(count2-1) t(iT(count2))]);
  end;
  Ci=Ci+(vi>=Tsort);
end;
P=Ci./size(RandStore,1);
P=fliplr(P);
for count=2:length(P)
  P(count)=max([P(count-1) P(count)]);
end;
P=fliplr(P);
Ptrue(iT)=P;

P_at_Node = [iT' Ptrue'];


%OMNIBUS TEST
omniP=sum(tmax>=max(T))./size(RandStore,1);
