clear
clc

% initCobraToolbox(false)

% changeCobraSolver('gurobi', 'LP')
% gurobi_setup

%%%% Load the model %%%%
load('Eflux_models.mat');

% model=Eflux_models.Control;
model=Eflux_models.Test;

n_rxns=length(model.rxnNames);
inisample=n_rxns;

% upsample=10*inisample; % Control
upsample=25*n_rxns; % Test

% mean_mat=zeros(n_rxns,10); % Control
mean_mat=zeros(n_rxns,25); % Test


iter=0;
h = waitbar(0,'Have patience...');

for i=inisample:n_rxns:upsample
    storedata = gpSampler(model,i);
    iter=iter+1;
    
    M1=mean(storedata.points,2);
    mean_mat(:,iter)=M1;
    
    clear storedata M1
    waitbar(i/upsample,h) % sprintf('%12.9f',valueofpi)
    
end
 
mean_mat_Eflux_Control=mean_mat;
save mean_mat_Eflux_Control mean_mat_Eflux_Control
 
mean_mat_Eflux_Test_final=mean_mat;
save mean_mat_Eflux_Test_final mean_mat_Eflux_Test_final

