clear
clc

%% Load the processed data
load('final_processed_control.mat')
load('final_processed_test.mat')

%% Load the metabolic model of Mtb
load('iEK1011_2.0.mat')

[p,q,~]=intersect(final_processed_control(:,1),string(iEK1011_2_0_COBRA.genes));
% Out of 1012 unique genes in iEK1011, 1000 genes are present in the data.
% The genes in final_processed_control and final_processed_test are same

[~,q1,~]=intersect(final_processed_test(:,1),string(iEK1011_2_0_COBRA.genes));
% Out of 1012 unique genes in iEK1011, 1000 genes are present in the data.
% The genes in final_processed_control and final_processed_test are same

Metabolic_genes_data=[p,final_processed_control(q,2),final_processed_test(q1,2)];
Metabolic_genes_data=[["Metabolic Genes","Control (Positive sulfur)","Test (Negative sulfur)"];Metabolic_genes_data];

clear p q q1

%% Building complete metabolic genes data for E-flux method
Final_metabolic_genes_data=zeros(size(iEK1011_2_0_COBRA.genes,1),2);
tmp=Metabolic_genes_data(2:end,:);

for i=1:size(iEK1011_2_0_COBRA.genes,1)
    idx=find(tmp(:,1)==string(iEK1011_2_0_COBRA.genes(i,1)));
    if ~isempty(idx)
        Final_metabolic_genes_data(i,1:2)=double(tmp(idx,2:3));
    end
    clear idx i
end
Final_metabolic_genes_data=[["Metabolic Genes","Control (Positive sulfur)","Test (Negative sulfur)"];...
    [string(iEK1011_2_0_COBRA.genes),Final_metabolic_genes_data]];

% save Final_metabolic_genes_data Final_metabolic_genes_data


%% Building reaction expression array
num_DATA= Final_metabolic_genes_data(2:end,2:3);
rxn_exp_Eflux=zeros(size(num_DATA,1),size(num_DATA,2));

for m=1:size(num_DATA,2)
    x=num_DATA(:,m);
    
    for i=1:length(iEK1011_2_0_COBRA.rxns)
        m
        i
        if ~isempty(iEK1011_2_0_COBRA.grRules{i})
            ans1=strsplit(iEK1011_2_0_COBRA.grRules{i},'or');
            A=zeros(length(ans1),1);
            for j=1:length(ans1)
                ans2=strfind(ans1{j},'and');
                if isempty(ans2)
                    ans3=ans1{j}(setdiff(1:length(ans1{j}),union(union(strfind(ans1{j},' '),strfind(ans1{j},')')),strfind(ans1{j},'('))));
                    [a,b]=intersect(iEK1011_2_0_COBRA.genes,ans3);
                    A(j,1)=x(b);
                else
                    ans4=strsplit(ans1{j},'and');
                    B=zeros(length(ans4),1);
                    for k=1:length(ans4)
                        ans5=ans4{k}(setdiff(1:length(ans4{k}),union(union(strfind(ans4{k},'('),strfind(ans4{k},')')),strfind(ans4{k},' '))));
                        [a,b]=intersect(iEK1011_2_0_COBRA.genes,ans5);
                        B(k,1)=x(b);
                    end
                    if length(find(B))==0
                        A(j,1)=0;
                    else
                        A(j,1)=min(B(find(B)));
                    end
                end
            end
            rxn_exp_Eflux(i,m)=sum(A);
        end
    end
end

%  save rxn_exp_Eflux rxn_exp_Eflux

m=max(max(rxn_exp_Eflux)); % 59.6437
A=rxn_exp_Eflux/m;

% save A A


%% Evaluating flux values -  FBA
initCobraToolbox(false)
changeCobraSolver('gurobi')
% changeCobraSolver('glpk')

% flux_value=[]; % [Control, Test]
% Biomass_value=[];
model1=iEK1011_2_0_COBRA;
model1.c(find(model1.c))=0;

% % Different biomass reactions
% % 1233: 'BiomassGrowth_sMtb'
% % 1234: 'R_BIOMASS_2'; -- default 1 in c
% % 1235: 'BIOMASS_2 Including Universal Cofactors (Xavier et al., 2017)'

% model1.c(1233)=1;
% model1.c(1234)=1;
% model1.c(1235)=1;

category={'Control', 'Test'};
for i=1:size(A,2)
    model=model1;
    model.ub(find(A(:,i)))=A(find(A(:,i)),i);
    model.ub(find(A(:,i)==0))=1;
    
    model.lb(intersect(find(model.lb),find(A(:,i))))=-A(intersect(find(model.lb),find(A(:,i))),i);
    model.lb(setdiff(find(model.lb),find(A(:,i))))=-1;
    
    %     solution = optimizeCbModel(model);
    
    %         Control_model=model; % for i=1
    %         Test_model=model; % for i=2
    
    %     flux_value(:,i)=solution.x;
    %     Biomass_value(i,1)=solution.f;
    
    str=['Eflux_models.',category{i},'=','model'];
    eval(str);
    
    clear model solution
end

% flux_value=[["Reactions name", "Reactions ID", "Control", "Test"];...
%     [string(iEK1011_2_0_COBRA.rxnNames), string(iEK1011_2_0_COBRA.rxns), flux_value]];

% save Eflux_models Eflux_models
