function [ centers,u_mat,labels,no_iterations,others] = fc_means( points,no_clusters ,m_val,epsilon,MAX_ITER,init_option,params )
%FC_MEANS Fuzzy C Means Algorithm Implementation

    disp 'fc_means called';

    no_points=size(points,1);
    feature_dim=size(points,2);
    if init_option==1
        %initialize memberships randomly
        u_mat = rand(no_points,no_clusters);
        %normalize for sum to 1
        u_mat= u_mat./sum(u_mat,2);
    %     disp(u_mat);
        new_centers=calculate_centers(points,u_mat,m_val);
    elseif init_option==2
        new_centers=params(1).init_centers;
        u_mat=zeros(no_points,no_clusters);
    end
        
%     disp(new_centers);
    old_centers=0;

    no_iterations=0;
    dist=zeros(no_points,no_clusters);
    
    save_centers=zeros(MAX_ITER,no_clusters,feature_dim);
    save_iter=1;
    
    while norm(new_centers-old_centers) > epsilon
      
        no_iterations=no_iterations+1;
%         disp(no_iterations);
        %calculate distances
                %calculate distances
        for j = 1 : no_clusters 
            t1=points-new_centers(j,:);
            t1=t1.^2;
            t2= sum (t1,2);
            %sum along 2nd index, that is get n x 1 matrix
            t2 = sqrt(t2);
            dist(:,j)= t2;
        end
        %update memberships
        u_mat=calculate_fcm_memberships(points,new_centers,dist,m_val);
          
        old_centers=new_centers;
        new_centers=calculate_centers(points,u_mat,m_val);   
        save_centers(save_iter,:,:)=new_centers;
        save_iter=save_iter+1;
        
        if no_iterations==MAX_ITER
            break
        end
    end
    centers=new_centers;
    
    others(1).save_centers=save_centers(1:(save_iter-1),:,:);

    %compare with MATLAB implementation of fcm
%     [centers,u_t]=fcm(points,no_clusters);
%     disp(centers);
%     u_mat=transpose(u_t);

    %hard partition
    labels=zeros(no_points,1);
    for i = 1 : no_points
        max_membership=0;
        for j = 1 : no_clusters 
            if ( u_mat (i,j) > max_membership)
                max_membership=u_mat(i,j);
                labels(i) = j;
            end
        end
    end
    
end

function centers= calculate_centers(points,u_mat,m_val)
% Calculates the new centers

    t1= u_mat.^m_val;
    t2= transpose(t1);

    % u_mat -> n x c , u_mat ' -> c x n , points -> n x d ,       
    %u_mat' * features -> c x d 
    centers=(t2*points);
    %If A is a matrix, then sum(A) returns a row vector containing the sum of each column.
    %u_mat -> n x c , membership_sums -> 1 x c
    membership_sums= sum (t1);
    %./ documentation : A is an M-by-N matrix and B is a scalar or 1-by-N row vector
    centers = transpose(transpose(centers)./membership_sums);
%     disp(centers);
    
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

