function [ centers,labels,no_iterations,others] = it2_fcm(points,no_clusters ,m1,m2 ,m,epsilon,MAX_ITER,init_option,params)
%IT2_FCM Interval Type 2 Fuzzy C Means with various m Algorithm Implementation
    disp 'IT2 FCM called'
    no_points= size(points,1);
    feature_dim=size(points,2);
    
    
    if init_option==1
        %initialize centers to results of FCM
        [new_centers,~,~,~,fcm_others]=fc_means( points,no_clusters ,2.0,epsilon,MAX_ITER,1,0 );
        others(1).fcm_save_centers=fcm_others(1).save_centers;
%         disp(new_centers);
    elseif init_option==2
        %initialize centers randomly from among the points
        random_permutation=randperm(no_points);
        new_centers=points(random_permutation(1:no_clusters),:);
    elseif init_option==3
        [new_centers,~,~,k_others]=k_means( points,no_clusters,epsilon,MAX_ITER,true );
        others(1).k_save_centers=k_others(1).save_centers;
    elseif init_option==4
        new_centers=params(1).init_centers;
    end
    

    %initialize centers randomly but within cluster
%     new_centers=zeros(no_clusters,feature_dim);
%     pts_per_cluster=no_points/no_clusters;
%     for j=1:no_clusters
%         random_permutation=randperm(pts_per_cluster);
%         disp(pts_per_cluster*(j-1) +random_permutation(1))
%         new_centers(j,:)= points(pts_per_cluster*(j-1) +random_permutation(1));
%     end

    disp('initial centers');
    disp(new_centers);

    dist=zeros(no_points,no_clusters);
    u_final=zeros(no_points,no_clusters);

    old_centers=0;
    
    %get sorted orders for KM
    
    [~,sorted_order] = sort(points,1);

    save_centers=zeros(MAX_ITER,no_clusters,feature_dim);
    save_iter=1;
    save_centers(save_iter,:,:)=new_centers;
    save_iter=save_iter+1;
    
    km_eps= params(1).km_eps;
    km_max_iter=params(1).km_max_iter;
    
    no_iterations=0;
    while norm(old_centers-new_centers)>epsilon
        no_iterations=no_iterations+1;
%         disp(no_iterations)

        %calculate distances
        for j = 1 : no_clusters 
            t1=points-new_centers(j,:);
            t1=t1.^2;
            t2= sum (t1,2);
            %sum along 2nd index, that is get n x 1 matrix
            t2 = sqrt(t2);
            dist(:,j)= t2;
        end

        
        %estimate upper & lower memberships
        u_mat1=calculate_fcm_memberships(points,new_centers,dist,m1);
        u_mat2=calculate_fcm_memberships(points,new_centers,dist,m2);
        u_upper=max(u_mat1,u_mat2);
        u_lower=min(u_mat1,u_mat2);
        u_upper_m=max(u_mat1.^m1,u_mat2.^m2);
        u_lower_m=min(u_mat1.^m1,u_mat2.^m2);

       
        old_centers=new_centers;
        %calculate center for every cluster
        for j = 1 : no_clusters
            [new_centers(j,:),u_final(:,j),~]=all_dimensions_km(points,sorted_order,u_lower_m(:,j),u_upper_m(:,j),1.0,km_eps,km_max_iter);
            %defuzzify to get crisp center
        end
%         disp(new_centers);
        
        save_centers(save_iter,:,:)=new_centers;
        save_iter=save_iter+1;
    
        if (no_iterations ==MAX_ITER)
            break
        end
    end
    centers=new_centers;
%     disp 'final centers';
%     disp(centers);
    
    others(1).save_centers=save_centers(1:(save_iter-1),:,:);

%hard partition
    
    u_final=(u_upper+u_lower)/2;


    labels=zeros(no_points,1);
    for i = 1 : no_points
        max_membership=-1;
        for j = 1 : no_clusters 
            if ( u_final (i,j) > max_membership)
                max_membership=u_final(i,j);
                labels(i) = j;
            end
        end
    end
end

function u_mat = calculate_fcm_memberships(points,new_centers,dist, m)

    no_points=size(points,1);
    no_clusters=size(new_centers,1);
    
    u_mat=zeros(no_points,no_clusters);
    
    for i = 1: no_points
        %calculate memberships for jth cluster
        for j = 1:no_clusters
            if dist(i,j)==0 %point is center of cluster j, therefore full membership
                u_mat(i,j)=1;
                continue;
            end
            denom=0;
            done=0;
            for k = 1: no_clusters
                if (dist(i,k)==0) %it is a center of different cluster k, 0 membership to jth cluster
                    u_mat(i,j)=0;
                    done=1;
                    break;
                end                       
                temp = dist (i,j)/dist(i,k);
                denom=denom+ temp^(2/(m-1));
            end
            if done~=1 %point was a center, hence already assigned
                u_mat(i,j)=1.0/denom;
            end
        end
    end 

end
