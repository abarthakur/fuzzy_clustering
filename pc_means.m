function [ centers,u_mat,labels,no_iterations,others] = pc_means( points,no_clusters ,m_val ,K_val ,alpha, epsilon, MAX_ITER, init_option, params)
%PC_Means Implementation with variations
    
    disp 'pc_means called';
    
    no_points= size(points,1);
 
    %initialization options
    
    if init_option==1
        %initialize centers and memberships from fc_means
        [fcm_centers,u_mat,~,~,fcm_others]=fc_means( points,no_clusters ,params(1).fcm_m_val,epsilon,MAX_ITER ,1,0);
        others(1).fcm_save_centers=fcm_others(1).save_centers;
%         disp(fcm_centers);
    elseif init_option==2
        %initialize memberships randomly
        disp 'random initialization'
        u_mat = rand(no_points,no_clusters);
        %umm, not really fcm centers. Will refactor one day.
        fcm_centers= calculate_centers(points,u_mat,m_val);
    elseif init_option==3
        fcm_centers=params(1).init_centers;
        u_mat=params(1).init_u;
    end
    
    % objective functions <-> update step options    
    if params(1).pcm_variation=="log"
        disp 'log pcm running'
        update_memberships= @log_memberships;
        %overwrite value of m to be used for calculating centers and eta
        %values
        m_val=1.0;
    else
        update_memberships= @usual_memberships;
    end
    
    
    %options related to calculation of eta
    
    eta_option=params(1).eta_option;
    %option 1 : run with given values of eta
    if eta_option == 1 
        eta=params(1).eta;

        disp('first eta values');
        disp(eta);

        [centers,u_mat,no_iterations,temp_others]=iterate(points,fcm_centers,no_clusters,m_val,epsilon,MAX_ITER,eta,0,update_memberships);
        others(1).save_centers=temp_others(1).save_centers;
        
    %option 2 : calculate eta on every iteration by 1st definition
    elseif eta_option == 2
        [centers,u_mat,no_iterations,temp_others]=iterate(points,fcm_centers,no_clusters,m_val,epsilon,MAX_ITER,[],K_val,update_memberships);
        others(1).save_centers=temp_others(1).save_centers;
        
    %option 3 : calculate eta once or twice
    elseif eta_option==3 || eta_option==4
        %first calculation of eta
%         centers=calculate_centers(points,u_mat,2.0);
        disp('initial centers');
        disp(fcm_centers);
        eta=calculate_eta_1(points,u_mat,fcm_centers,no_clusters,m_val,K_val);    
        disp('first eta values');
        disp(transpose(eta));
        [centers,u_mat,n1,it1_others]=iterate(points,fcm_centers,no_clusters,m_val,epsilon,MAX_ITER,eta,0,update_memberships);
        
        others(1).save_centers=it1_others(1).save_centers;
        others(1).eta=eta;
        
        n2=0;
        if eta_option==4
            disp(centers)
            eta=calculate_eta_2(points,u_mat,centers,no_clusters,alpha);
            disp('second eta values');
            disp(transpose(eta));   
            [centers,u_mat,n2,it2_others]=iterate(points,fcm_centers,no_clusters,m_val,epsilon,MAX_ITER,eta,0,update_memberships);
            others(1).save_centers=cat(1,it1_others(1).save_centers,it2_others(1).save_centers);
        end
        
        no_iterations=n1+n2;
    end
    
    
    
    
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
    
%     disp('pcm output u vals');
%     disp(u_mat);
    
end


function [centers,u_mat,no_iterations,others]=iterate(points,fcm_centers,no_clusters,m_val,epsilon,MAX_ITER,eta,K_val,update_memberships)
    
    no_points= size(points,1);
    feature_dim=size(points,2);

    no_iterations=0;
%     u_mat_old=0;
%     U_save=zeros(no_points,no_clusters*100);
%     while norm(u_mat_old- u_mat) > epsilon
    new_centers=fcm_centers;
    old_centers=0;
    u_mat=zeros(no_points,no_clusters);
    
    save_centers=zeros(MAX_ITER,no_clusters,feature_dim);
    save_iter=1;
    save_centers(save_iter,:,:)=new_centers;
    save_iter=save_iter+1;
    
    while(norm(new_centers-old_centers)>epsilon)
  
%             U_save(:,no_iterations*no_clusters+1:no_iterations*no_clusters+2)=u_mat;
            no_iterations=no_iterations+1;
            if (K_val ~=0) %calculate eta each time
                if no_iterations==1
                    eta=rand(1,no_clusters);
                else
                    eta= calculate_eta_1(points,u_mat,new_centers,no_clusters,m_val,K_val);
                end
            end
            %update memberships
%             u_mat_old=u_mat;
            u_mat=update_memberships(points,u_mat,new_centers,no_clusters,m_val,eta);
            old_centers=new_centers;
            new_centers=calculate_centers(points,u_mat,m_val);
            
            save_centers(save_iter,:,:)=new_centers;
            save_iter=save_iter+1;
            
            if no_iterations==MAX_ITER
                break
            end
    end
%     csvwrite('save.csv',U_save);
    centers=new_centers;    
    others(1).save_centers=save_centers(1:(save_iter-1),:,:);
end

function u_mat = log_memberships(points,u_mat,new_centers,no_clusters,~,eta)
    
    no_points= size(points,1);
    for i = 1:no_points 
        xi= points(i,:);
        for j = 1 : no_clusters 
            cj=new_centers(j,:);
            dij= norm(xi-cj);
%                 disp(dij^2);
            u_mat(i,j)=exp(-(dij^2)/eta(j));    
        end
    end    
end


function u_mat = usual_memberships(points,u_mat,new_centers,no_clusters,m_val,eta)
    
    no_points= size(points,1);
    for i = 1:no_points 
        xi= points(i,:);
        for j = 1 : no_clusters 
            cj=new_centers(j,:);
            dij= norm(xi-cj);
%                 disp(dij^2);
            denom = 1+ (dij^2/eta(j))^(1/(m_val-1));
            u_mat(i,j)=1/denom;    
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

function eta=calculate_eta_1(points,u_mat,centers,no_clusters,m,K_val)
    u_mat=u_mat.^m;
    eta=zeros(no_clusters,1);
    for j=1:no_clusters
        num=0;
        denom=0;
        for i = 1: size(points,1)
            num=num+u_mat(i,j)*(norm(points(i)-centers(j))^2);
            denom=denom+u_mat(i,j);
        end
        eta(j)=K_val* (num/denom);
    end
end


function eta=calculate_eta_2(points,u_mat,centers,no_clusters,alpha)
    eta=zeros(no_clusters,1);
    for j=1:no_clusters
        num=0;
        denom=0;
        for i = 1: size(points,1)
            if (u_mat(i,j) > alpha)
                num=num+norm(points(i)-centers(j))^2;
                denom=denom+1;            
            end
        end
        eta(j)=(num/denom);
    end
end
            
            

            
            
            

