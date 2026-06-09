clear
clc

%% Import the data
[~, ~, data] = xlsread('Data_Location\deg_genes.xlsx','deg_genes');

data = string(data);
data(ismissing(data)) = '';

genes=data(2:end,2);
control=double(data(2:end,12:14)); % positive sulfur
test=double(data(2:end,9:11)); % negative sulfur

%% Pre-processing of control

% Step1: Removing genes with 0 expression in all samples
ind=find(sum(control,2)~=0);
processed_control=control(ind,:);

g1=genes(ind,1); % Filtered genes

% Step2: Removing genes present in <50% samples
ind1=[];
for i=1:size(processed_control,1)
    ind1(i,1)=length(find(processed_control(i,:)));
end

find(ind1>=(size(processed_control,2)/2));
processed_control=processed_control(ans,:);
g1=g1(ans,1);

clear ind ind1

% Step3: Missing value imputation (Replace zero values with mean values of that gene in that condition)
y=[];
for i=1:size(processed_control,1)
    x=processed_control(i,:);
    x(x==0)=mean(x(find(x)));
    y(i,:)=x;
end
processed_control=y;

clear x y ans i

% Step4: Outlier identification and fixing
IQR=iqr(processed_control,2); % Finding inter quartile range of each gene
Q=quantile(processed_control,3,2); % Each row contains first, second, third quartile of each gene

for i=1:size(processed_control,1)
    idx1=find(processed_control(i,:)>max(Q(i,:))+1.5*IQR(i));
    idx2=find(processed_control(i,:)<min(Q(i,:))-1.5*IQR(i));
    outlier_idx=[idx1'; idx2'];
    non_outlier_idx=setdiff([1:size(processed_control,2)],outlier_idx);
    max_value=max(processed_control(i,non_outlier_idx));
    min_value=min(processed_control(i,non_outlier_idx));
    for j=1:size(outlier_idx,1)
        if processed_control(i,outlier_idx(j))>max(Q(i,:))+1.5*IQR(i)
            processed_control(i,outlier_idx(j))=max_value;
        elseif processed_control(i,outlier_idx(j))<min(Q(i,:))-1.5*IQR(i)
            processed_control(i,outlier_idx(j))=min_value;
        end
    end
    i
    
    clear idx1 idx2 outlier_idx non_outlier_idx max_value min_value
end

clear IQR Q i j

% Step5: Scaling genes by dividing with max value of that gene across all samples (gene-wise) -> converts  data 0 to 1
data=[];
for i=1:size(processed_control,1)
    x=double(processed_control(i,:));
    data(i,:)=x/max(x);
    
    clear x i 
end

% Step6: Combining data across samples
final_processed_control=[];
for i=1:size(data,1)
    x=mean(double(data(i,:)),2);
    final_processed_control(i,1)=x;
    
    clear x i
end
final_processed_control=[g1,final_processed_control];
clearvars -except genes final_processed_control test


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Pre-processing of test

% Step1: Removing genes with 0 expression in all samples
ind=find(sum(test,2)~=0);
processed_test=test(ind,:);

g1=genes(ind,1); % Filtered genes

% Step2: Removing genes present in <50% samples
ind1=[];
for i=1:size(processed_test,1)
    ind1(i,1)=length(find(processed_test(i,:)));
end

find(ind1>=(size(processed_test,2)/2));
processed_test=processed_test(ans,:);
g1=g1(ans,1);

clear ind ind1

% Step3: Missing value imputation (Replace zero values with mean values of that gene in that condition)
y=[];
for i=1:size(processed_test,1)
    x=processed_test(i,:);
    x(x==0)=mean(x(find(x)));
    y(i,:)=x;
end
processed_test=y;

clear x y ans i

% Step4: Outlier identification and fixing
IQR=iqr(processed_test,2); % Finding inter quartile range of each gene
Q=quantile(processed_test,3,2); % Each row contains first, second, third quartile of each gene

for i=1:size(processed_test,1)
    idx1=find(processed_test(i,:)>max(Q(i,:))+1.5*IQR(i));
    idx2=find(processed_test(i,:)<min(Q(i,:))-1.5*IQR(i));
    outlier_idx=[idx1'; idx2'];
    non_outlier_idx=setdiff([1:size(processed_test,2)],outlier_idx);
    max_value=max(processed_test(i,non_outlier_idx));
    min_value=min(processed_test(i,non_outlier_idx));
    for j=1:size(outlier_idx,1)
        if processed_test(i,outlier_idx(j))>max(Q(i,:))+1.5*IQR(i)
            processed_test(i,outlier_idx(j))=max_value;
        elseif processed_test(i,outlier_idx(j))<min(Q(i,:))-1.5*IQR(i)
            processed_test(i,outlier_idx(j))=min_value;
        end
    end
    i
    
    clear idx1 idx2 outlier_idx non_outlier_idx max_value min_value
end

clear IQR Q i j

% Step5: Scaling genes by dividing with max value of that gene across all samples (gene-wise) -> converts  data 0 to 1
data=[];
for i=1:size(processed_test,1)
    x=double(processed_test(i,:));
    data(i,:)=x/max(x);
    
    clear x i 
end

% Step6: Combining data across samples
final_processed_test=[];
for i=1:size(data,1)
    x=mean(double(data(i,:)),2);
    final_processed_test(i,1)=x;
    
    clear x i
end
final_processed_test=[g1,final_processed_test];
clearvars -except genes final_processed_control final_processed_test

save final_processed_control final_processed_control
save final_processed_test final_processed_test


