clear
clc

load('mean_mat_Eflux_Control.mat')
load('mean_mat_Eflux_Test_final.mat')
load('Eflux_models.mat')

Final_Eflux_Control_Flux=mean_mat_Eflux_Control(:,10);
Final_Eflux_Test_Flux=mean_mat_Eflux_Test_final(:,25);

%% Perturbed reactions (Removed reactions that have opposite flux directions in the two cases and are perturbed)

FC=Final_Eflux_Test_Flux(:,1)./Final_Eflux_Control_Flux(:,1);


% Upregulated reactions (1.5 fold change)
idx=find(FC>(1.5));
temp=[string(Eflux_models.Control.rxns(idx)),Final_Eflux_Control_Flux(idx,1),Final_Eflux_Test_Flux(idx,1),...
    FC(idx),string(Eflux_models.Control.subSystems(idx))];

ans1=double(temp(:,2))>0 & double(temp(:,3))>0;
ans2=double(temp(:,2))<0 & double(temp(:,3))<0;
ind= find(ans1 | ans2);

x=zeros(size(temp,1),1);
x(ans1,1)=1;
x(ans2,1)=-1;

Up_rxns_1_5fold=[["Reactions","Control_flux_values","Test_flux_values","Fold_change","Subsystems","Direction"];...
    [temp(ind,:),x(ind,1)]];

clear idx temp ind ans1 ans2 x

% save Up_rxns_1_5fold Up_rxns_1_5fold


% Downregulated  reactions (1.5 fold change)

idx=find(FC<(1/(1.5)));
temp=[string(Eflux_models.Control.rxns(idx)),Final_Eflux_Control_Flux(idx,1),Final_Eflux_Test_Flux(idx,1),...
    FC(idx),string(Eflux_models.Control.subSystems(idx))];

ans1=double(temp(:,2))>0 & double(temp(:,3))>0;
ans2=double(temp(:,2))<0 & double(temp(:,3))<0;
ind= find(ans1 | ans2);

x=zeros(size(temp,1),1);
x(ans1,1)=1;
x(ans2,1)=-1;

Down_rxns_1_5fold=[["Reactions","Control_flux_values","Test_flux_values","Fold_change","Subsystems","Direction"];...
    [temp(ind,:),x(ind,1)]];

clear idx temp ind ans1 ans2 x

% save Down_rxns_1_5fold Down_rxns_1_5fold
