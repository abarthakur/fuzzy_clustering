function [v_avg,u_avg,other_vals]=all_dimensions_km(points,sorted_order,u_lower,u_upper,m,epsilon,MAX_ITER)
    
    no_points=size(points,1);
    feature_dim = size(points,2);
    %u_left,u_right -> n x d , u_upper,u_lower -> n x 1
    u_left=zeros(no_points,feature_dim);
    u_right=zeros(no_points,feature_dim);
    %initialize u_left and u_right for each feature
    for d = 1 : feature_dim
        u_left(:,d)=(u_lower+u_upper)/2;
        u_right(:,d)=(u_lower+u_upper)/2;
    end
    
    %initialize centroid as weighted average
    t1= points.*(u_left.^m); %m_val
    t1 = sum (t1,1);
    init_center=t1./sum(u_left.^m,1); %m_val
    
    %calculate v_left
    new_center=init_center;
    old_center=zeros(1,feature_dim);
    for d = 1 : feature_dim
        %run 1 dimensional KM for this feature
        no_iterations=0;
        while norm(old_center(d)-new_center(d)) > epsilon
            no_iterations=no_iterations+1;
            %find switch point
            k=find_switch(new_center(d),points(:,d),sorted_order(:,d));
%             fprintf("dimension %d. v = %f, switch = %d\n",d,new_center(d),k);
            
            %set u vals , switching at k
            first_k_idx=sorted_order(1:k,d);
            other_idx=sorted_order(k+1:end,d);
            
            u_left(first_k_idx,d)=u_upper(first_k_idx);
            u_left(other_idx,d)=u_lower(other_idx);
            
            %calculate new center
            old_center(d)=new_center(d);
            t2=points(:,d).*u_left(:,d).^m;%m_val
            t2=sum(t2,1);            
            new_center(d)= t2./sum(u_left(:,d).^m,1);%m_val 
            if (no_iterations ==MAX_ITER)
                break
            end 
        end
    end
    v_left=new_center;
    
    %calculate v_right
    
    new_center=init_center;
    old_center=zeros(1,feature_dim);
    
    for d = 1 : feature_dim
        %run 1 dimensional KM for this feature
        no_iterations=0;
        while norm(old_center(d)-new_center(d)) > epsilon
            no_iterations=no_iterations+1;
            %find switch point
            k=find_switch(new_center(d),points(:,d),sorted_order(:,d));
            
            %set u vals , switching at k
            first_k_idx=sorted_order(1:k,d);
            other_idx=sorted_order(k+1:end,d);
            
            u_right(first_k_idx ,d)=u_lower(first_k_idx);
            u_right(other_idx,d)=u_upper(other_idx);
            
            %calculate new center
            old_center(d)=new_center(d);
            t2=points(:,d).*u_right(:,d).^m;%m_val
            t2=sum(t2,1);            
            new_center(d)= t2./sum(u_right(:,d).^m,1);%m_val
            if (no_iterations ==100)
                break
            end 
        end
    end
    v_right=new_center;
    
    u_left_reduced= sum(u_left,2)./feature_dim;
    u_right_reduced= sum(u_right,2)./feature_dim;
    u_avg= (u_left_reduced+u_right_reduced)/2;
    
    v_avg= (v_left+v_right)/2;
    
    other_vals(1).v_left = v_left;
    other_vals(1).v_right = v_right;
    
end

function k= find_switch(v,x_vals,sorted_order)
    no_points=size(x_vals);
    for i=1:no_points-1
        idx=sorted_order(i);
        idxplus1=sorted_order(i+1);
        if ( x_vals(idx) <= v && v <= x_vals(idxplus1))
            k=i;
            break;
        end
    end 
end