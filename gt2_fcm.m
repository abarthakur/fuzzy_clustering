function [centers,labels,no_iterations,others] = gt2_fcm(points,no_clusters ,m_set,no_planes,epsilon,MAX_ITER,init_option,params)
%GT2 FCM
    disp('GT2 FCM called');
    
    no_points=size(points,1);
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

    disp("initial centers:");
    disp(new_centers);
    
    
    

    
        
    %determine the intervals for alpha cuts of the fuzzifiers
    % m_set is an array of (m,primary_membership)
%     m_set
    [~,I]=sort(m_set(:,2));
%     I
    sorted_m=m_set(transpose(I),:);
%     sorted_m
    alpha=0;
    step=1.0/(no_planes-1);
    k=1;
    m_lower=zeros(no_planes,1);
    m_upper=zeros(no_planes,1);
    for i=1:no_planes
        while(sorted_m(k,2)<alpha && k < size(sorted_m,1))
            k=k+1;
        end
%         disp(sorted_m(k:end,1));
        m_max=max(sorted_m(k:end,1));
        m_min=min(sorted_m(k:end,1));
        m_lower(i)=m_min;
        m_upper(i)=m_max;
        alpha=alpha+step;
    end
%     disp(m_lower);
%     disp(m_upper);

    km_eps= params(1).km_eps;
    km_max_iter=params(1).km_max_iter;

    [~,sorted_order] = sort(points,1);
    
    dist=zeros(no_points,no_clusters);
    
    u_lower=zeros(no_planes,no_points,no_clusters);
    u_upper=zeros(no_planes, no_points,no_clusters);
    
    u_lower_m=zeros(no_planes,no_points,no_clusters);
    u_upper_m=zeros(no_planes, no_points,no_clusters);
    
    left_centers=zeros(no_planes, no_clusters,feature_dim);
    right_centers=zeros(no_planes, no_clusters,feature_dim);
    
    old_centers=0;
    no_iterations=0;
    
    save_centers=zeros(MAX_ITER,no_clusters,feature_dim);
    save_iter=1;
    save_centers(save_iter,:,:)=new_centers;
    save_iter=save_iter+1;
    
    while (norm(new_centers-old_centers)>epsilon)

        no_iterations=no_iterations+1;
        
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
        for i = 1 : no_planes
            m1=m_lower(i);
%             disp(m1);
            m2=m_upper(i);
            u_mat1=calculate_fcm_memberships(points,new_centers,dist,m1);
            u_mat2=calculate_fcm_memberships(points,new_centers,dist,m2);
            u_upper(i,:,:)=max(u_mat1,u_mat2);
            u_lower(i,:,:)=min(u_mat1,u_mat2);

            u_upper_m(i,:,:)=max(u_mat1.^m1,u_mat2.^m2);
            u_lower_m(i,:,:)=min(u_mat1.^m1,u_mat2.^m2);

        end
        %update centers

        for i = 1 : no_planes
            for j= 1: no_clusters

                [~,~,other_vals]=all_dimensions_km(points,sorted_order,u_lower_m(i,:,j).',u_upper_m(i,:,j).',1.0,km_eps,km_max_iter);
                v_left=other_vals(1).v_left;
                v_right=other_vals(1).v_right;
                left_centers(i,j,:)=v_left;
                right_centers(i,j,:)=v_right;
 
            end
        end
        %defuzzify
        old_centers=new_centers;
        new_centers=new_centers*0;
        alpha=0;
        step=1.0/(no_planes-1);
        denom=0;
        for i=1:no_planes
            lc=reshape(left_centers(i,:,:),[no_clusters,feature_dim]);
            rc=reshape(right_centers(i,:,:),[no_clusters,feature_dim]);
            new_centers=new_centers+lc.*alpha +rc.*alpha;
            denom=denom+2*alpha;
            alpha=alpha+step;
        end
        new_centers= new_centers./denom;
        
%         disp(new_centers);

        save_centers(save_iter,:,:)=new_centers;
        save_iter=save_iter+1;

        if (no_iterations>=MAX_ITER)
            break;
        end
    end

    centers=new_centers;
    
    others(1).save_centers=save_centers(1:(save_iter-1),:,:);
    

    u_final=zeros(no_points,no_clusters);
  
    alpha=0;
    step=1.0/(no_planes-1);

    for i=1 : no_planes
        ul=reshape( u_lower(i,:,:),[no_points,no_clusters]);
        ur=reshape( u_upper(i,:,:),[no_points,no_clusters]);
        u_final=u_final+ul.*alpha + ur.*alpha;
        denom=denom+2*alpha;
        alpha=alpha+step;
    end
    u_final= u_final./denom;


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