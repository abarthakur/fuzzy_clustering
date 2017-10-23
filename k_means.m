function [ centers,labels,no_iterations,others] = k_means( points,no_clusters,epsilon,MAX_ITER,init_option,params )
%K_MEANS K Means Algorithm Implementation
    
    disp 'k_means called';
    
    no_points=size(points,1);
    feature_dim=size(points,2);

%     initialize centers 

    if init_option==1
        %initialize by Forgy method : choose k random centers
        random_permutation=randperm(no_points);
        new_centers=points(random_permutation(1:no_clusters),:);
        
    elseif init_option==2
        %initialize by random partition method : choose a random partition
        labels=randi([1 no_clusters],no_points,1);
        new_centers=calculate_centers(points,labels,feature_dim,no_points,no_clusters);
    elseif init_option==3
        new_centers=params(1).init_centers;
    end
    disp("initial centers :");
    disp(new_centers);
    
    old_centers=0;    
    no_iterations=0;
    
    save_centers=zeros(MAX_ITER,no_clusters,feature_dim);
    save_iter=1;
    save_centers(save_iter,:,:)=new_centers;
    save_iter=save_iter+1;
    
    while norm(new_centers-old_centers) > epsilon
        no_iterations=no_iterations+1;
%         disp(no_iterations);
%         disp(new_centers);
        
        %partition points
        for i = 1:no_points 
            point= points(i,:);
            min_dist = Inf ;
            for j = 1 : no_clusters
                c=new_centers(j,:);
                if  norm(point-c) < min_dist 
                    labels(i)=j;
                    min_dist=norm(point-c);
                end
            end
        end
        %update center
        old_centers=new_centers;
        new_centers=calculate_centers(points,labels,feature_dim,no_points,no_clusters);   
        
        save_centers(save_iter,:,:)=new_centers;
        save_iter=save_iter+1;
        
        if no_iterations==MAX_ITER
            break
        end
        
    end
    centers=new_centers;
    others(1).save_centers=save_centers(1:(save_iter-1),:,:);
        
    
end

function centers= calculate_centers(points,labels,feature_dim,no_points,no_clusters)
% Calculates the new centers
    centers=zeros(no_clusters,feature_dim);
    member_counts=zeros(no_clusters,1);
    for i = 1: no_points
        point = points(i,:);
        member_counts(labels(i))=member_counts(labels(i))+1;
        centers(labels(i),:)=centers(labels(i),:)+point;
    end
%     disp(centers);
%     disp(member_counts);
    for i = 1: no_clusters 
        if member_counts(i)~=0
            centers(i,:)= centers(i,:)/member_counts(i);
        else
            fprintf("Some cluster was assigned 0 points in an iteration");            
        end
    end
%     disp(centers);
    
end





