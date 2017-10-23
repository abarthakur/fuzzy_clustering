%read data
[features,labels,no_types,divisions]=read_iris_data();
% features= features(:,[3 4]);

%normalize data
%all values together
% features=(features-min(features(:)))/(max(features(:))- min(features(:)));
%normalise each feature

for i = 1 : size(features,2)
    features(:,i)= (features(:,i)-min(features(:,i)))/(max(features(:,i))- min(features(:,i)));
end
% features=features(:,[1 2]);
% disp (features1-features);

%feed data to this helper func
run_cluster( features,labels,no_types,divisions)



function [features,labels,no_types,divisions]=read_iris_data() 
    %read iris data set
    fprintf("Working on IRIS data set\n");
    datafile='./DATA/iris.csv';
    %read dataset into table
    data=readtable(datafile);
    %retrieve data into arrays. Documentation - access data from a table.
    features=data{:,2:end-1}; %first column is id, last column is label
    labels_in_data= data{:,end};
    divisions=[0,50,100,150];
    [labels,no_types]=convert_labels(labels_in_data);
end


function [labels,no_types] = convert_labels(labels_in_data)
    no_points=size(labels_in_data,1);
    labels=zeros(no_points,1);
    label_count=1;
    label_map=containers.Map;
    for i = 1:size(labels_in_data,1)
        label=char(labels_in_data(i));
        if isKey(label_map,label)~= 1
            label_map(label)=label_count;
            label_count=label_count+1;
        end
        labels(i)=label_map(label);
    end
    no_types=label_count-1;
end