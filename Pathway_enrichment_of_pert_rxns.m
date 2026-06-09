clear
clc

% Load the mouse model
load('iEK1011_2.0.mat')

% Load the flux values
load('Final_Eflux_Control_Flux.mat')
load('Final_Eflux_Test_Flux.mat')

% Load the perturbed reactions (1.5 fold)
load('Up_rxns_1_5fold.mat')
load('Down_rxns_1_5fold.mat')

pert_rxns=union(Up_rxns_1_5fold(2:end,1), Down_rxns_1_5fold(2:end,1));

temp1=Up_rxns_1_5fold(2:end,:);
temp2=Down_rxns_1_5fold(2:end,:);

tot_pathways=union(temp1(:,5), temp2(:,5));
tot_pathways(1) = []; % empty string

% Finding up and down reactions in each pathway
temp3=zeros(size(tot_pathways,1),2);
for i=1:size(tot_pathways,1)
    find(temp1(:,5)==tot_pathways(i,1)) % Upregulated reactions
    temp3(i,1)=length(ans);
    
    clear ans
    
    find(temp2(:,5)==tot_pathways(i,1)) % Downregulated reactions
    temp3(i,2)=length(ans);
    
    clear ans
    
end

% Finding total reactions in each pathway
tot_rxns=[];
for i=1:size(tot_pathways,1)
    find(string(iEK1011_2_0_COBRA.subSystems(:,1))==tot_pathways(i,1));
    tot_rxns(i,1)=length(ans);
end

% Dividing number of up and down regulated reactions by number of total rxns in that pathway
temp4=[];
temp4(:,1)=temp3(:,1)./tot_rxns(:,1);
temp4(:,2)=temp3(:,2)./tot_rxns(:,1);
temp4(:,3)=tot_rxns(:,1)./tot_rxns(:,1);


%% For plotting
% Adding total number of reactions in brackets with pathways
combined=cell(size(tot_pathways,1),1);
for i=1:size(tot_pathways,1)
    combined{i}=[tot_pathways{i,1} ' (' num2str(tot_rxns(i,1)) ')'];
end
combined=string(combined);

tot_pathways=[["Pathways (Total number of reactions)", "Up-regulated reactions", "Down-regulated reactions", "tot_rxns"];...
    combined,temp4];

save tot_pathways tot_pathways



