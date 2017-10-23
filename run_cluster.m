function run_cluster( features,labels,no_clusters,divisions)
%RUN_CLUSTER Runs the clustering algorithm and prints plots
    no_points=size(features,1);
    feature_dim=size(features,2);
    fprintf("run_cluster called\n\n");
    
    %parameters for various clustering algorithms

    %%K-Means
    
    epsilon= 10 ^(-3);
    MAX_ITER=200;

    random_permutation=randperm(no_points);
    new_centers=features(random_permutation(1:no_clusters),:);

    init_option=2;
    params(1).init_centers=new_centers;

    [ centers,output_labels,no_iterations,others] = k_means( features,no_clusters,epsilon,MAX_ITER,init_option,params );
    disp(centers);
    disp(no_iterations);
    estimate_error(centers,output_labels,labels,no_points,divisions);

    %%FCM
    epsilon= 10 ^(-3);
    MAX_ITER=200;

    m_val=2.0;

    random_permutation=randperm(no_points);
    new_centers=features(random_permutation(1:no_clusters),:);

    init_option=2;
    params(1).init_centers=new_centers;

    [ centers,u_mat,output_labels,no_iterations] = fc_means( features,no_clusters ,m_val,epsilon,MAX_ITER,init_option,params );
    disp(centers);
    disp(no_iterations);
    estimate_error(centers,output_labels,labels,no_points,divisions);
    

    %%PCM

    epsilon= 10 ^(-3);
    MAX_ITER=200;

    pcm_m_val=1.3;
    K_val=1.0;
    alpha=0.2;

    params(1).pcm_variation="log";
    params(1).eta_option=4;
    params(1).eta= [0.0212  , 0.0129  ,  0.0266];
    params(1).fcm_m_val=2.0;

    random_permutation=randperm(no_points);
    new_centers=features(random_permutation(1:no_clusters),:);
    params(1).init_centers=new_centers;
    params(1).init_u=calculate_fcm_memberships(features,new_centers,params(1).fcm_m_val);

    init_option=1;

    [ centers,u_mat,output_labels,no_iterations,others] = pc_means( features,no_clusters ,pcm_m_val , K_val ,alpha ,epsilon, MAX_ITER, init_option, params);
    disp(centers);
    disp(no_iterations);
    estimate_error(centers,output_labels,labels,no_points,divisions);


    %%IT2 FCM
    
    epsilon= 10 ^(-3);
    MAX_ITER=200;
    
    m1=2.0;
    m2=3.0;
    IT2_m_val=(m1+m2)/2;

    random_permutation=randperm(no_points);
    new_centers=features(random_permutation(1:no_clusters),:);

    init_option=4;
    params(1).init_centers=new_centers;

    params(1).km_eps = 10 ^(-3);
    params(1).km_max_iter = 100;

    [ centers,output_labels,no_iterations,others] = it2_fcm(features,no_clusters ,m1,m2 ,IT2_m_val,epsilon,MAX_ITER,init_option,params);
    disp(centers);
    disp(no_iterations);
    estimate_error(centers,output_labels,labels,no_points,divisions);


    %%GT2 FCM
    epsilon= 10 ^(-3);
    MAX_ITER=200;
    

    no_planes=5;
    N=17;
    m1=2.0;
    m2=4.0;
    m_set=generate_fuzzifier(m1,m2,4,no_planes,N);

    random_permutation=randperm(no_points);
    new_centers=features(random_permutation(1:no_clusters),:);

    init_option=3;
    params(1).init_centers=new_centers;

    params(1).km_eps = 10 ^(-3);
    params(1).km_max_iter = 100;


    [centers,output_labels,no_iterations,others] = gt2_fcm(features,no_clusters ,m_set,no_planes,epsilon,MAX_ITER,init_option,params);
    disp(centers);
    disp(no_iterations);
    estimate_error(centers,output_labels,labels,no_points,divisions);

   
    %tSNE visualisation of multi-dimensional data set

%     Y = tsne(features,'Algorithm','barneshut','NumPCAComponents',size(features,2));
%     figure;
%     gscatter(Y(:,1),Y(:,2),labels);
% %     hold on;    
%     figure;
%     gscatter(Y(:,1),Y(:,2),output_labels);

end


%N must be an odd no for this function to work properly for all types except type=4,5
%For type=4,5 , N+1 must be divisible by 9
function m_set=generate_fuzzifier(m1,m2,type,no_planes,N)
    m_set=zeros(N,2);
    %triangular membership function
    if type==1
        m_set(:,1)= linspace(m1,m2,N) ;
        t1= transpose(linspace(0,1,(N+1)/2));
        t2= transpose(linspace(1,0,(N+1)/2));
        t= [t1(1:end-1);t2];
        m_set(:,2)=t;
        plot(m_set(:,1),m_set(:,2));
    end

    if type==2
       m_set(:,1)= linspace(m1,m2,N) ;
       t1= transpose(linspace(0,1,(N+1)/2));
       t2= transpose(linspace(1,0,(N+1)/2));
       t= [t1(1:end-1).^3;t2.^3];
       m_set(:,2)=t;
       % plot(m_set(:,1),m_set(:,2));
    end

    if type==3
       m_set(:,1)= linspace(m1,m2,N) ;
       t1= transpose(linspace(0,1,(N+1)/2));
       t2= transpose(linspace(1,0,(N+1)/2));
       sigma=0.3;
       t1=(1/sigma*sqrt(2*pi))*exp(-0.5*(((t1-1)./sigma).^2));
       t2=(1/sigma*sqrt(2*pi))*exp(-0.5*(((t2-1)./sigma).^2));
       t= [t1(1:end-1);t2];
       t=(t-min(t))/(max(t)-min(t));
       m_set(:,2)=t;
       % plot(m_set(:,1),m_set(:,2));
    end
    %N+1 must be divisible by 9
    if type==4
        m_set(:,1)= linspace(m1,m2,N) ;
        t1= transpose(linspace(0,1,(N+1)/9));
        t2= transpose(linspace(1,0,8*(N+1)/9));
        sigma=0.3;
        t1=(1/sigma*sqrt(2*pi))*exp(-0.5*(((t1-1)./sigma).^2));
        t2=(1/sigma*sqrt(2*pi))*exp(-0.5*(((t2-1)./sigma).^2));
        t= [t1(1:end-1);t2];
        t=(t-min(t))/(max(t)-min(t));
        m_set(:,2)=t;
        % plot(m_set(:,1),m_set(:,2));
    end


    if type==5
        m_set(:,1)= linspace(m1,m2,N) ;
        t1= transpose(linspace(0,1,8*(N+1)/9));
        t2= transpose(linspace(1,0,(N+1)/9));
        sigma=0.3;
        t1=(1/sigma*sqrt(2*pi))*exp(-0.5*(((t1-1)./sigma).^2));
        t2=(1/sigma*sqrt(2*pi))*exp(-0.5*(((t2-1)./sigma).^2));
        t= [t1(1:end-1);t2];
        t=(t-min(t))/(max(t)-min(t));
        m_set(:,2)=t;
        % plot(m_set(:,1),m_set(:,2));
    end
end


function estimate_error(centers,output_labels,labels,no_points,divisions)
    no_clusters=size(divisions,2)-1;
%     fprintf("clusters is ");
%     disp(no_clusters);
    label_map=containers.Map(1,1);
    %disp(label_map);
%     disp (divisions);
%     fprintf("What is size : %d\n",size(divisions,2));
    for j=1:size(divisions,2)-1
%         sprintf("J is %d\n;");
        label_map(mode(output_labels((divisions(j)+1):divisions(j+1))))=j;
    end
%     disp(label_map.keys());
%     disp(label_map.values());

    %calculate errors
    if size(keys(label_map),2)==no_clusters
        error_count=0;
        for i = 1:no_points 
            if label_map(output_labels(i)) ~= labels(i) 
            %             disp([i , label_map(output_labels(i)), labels(i)]);
            %             disp (features(i,:));
                error_count=error_count+1;
            end
        end
        fprintf("Error count = %d & Success Rate = %f \n",error_count,(no_points-error_count)* 100/no_points );
    else
        fprintf("Mode method of mapping fails : Likely too many errors!\n")
    end

end



function centers= calculate_centers(features,labels,feature_dim,no_points,no_clusters)
% Calculates the new centers
    centers=zeros(no_clusters,feature_dim);
    member_counts=zeros(no_clusters,1);
    for i = 1: no_points
        point = features(i,:);
        member_counts(labels(i))=member_counts(labels(i))+1;
        centers(labels(i),:)=centers(labels(i),:)+point;
    end
%     disp(centers);
%     disp(member_counts);
    for i = 1: no_clusters 
        if member_counts(i)~=0
            centers(i,:)= centers(i,:)/member_counts(i);
        else
            disp 'something wrong!';            
        end
    end
%     disp(centers);
    
end


function u_mat = calculate_fcm_memberships(points,new_centers, m)

    no_points=size(points,1);
    no_clusters=size(new_centers,1);
    
    u_mat=zeros(no_points,no_clusters);
    dist=zeros(no_points,no_clusters);
    
    %calculate distances
    for j = 1 : no_clusters 
        t1=points-new_centers(j,:);
        t1=t1.^2;
        t2= sum (t1,2);
        %sum along 2nd index, that is get n x 1 matrix
        t2 = sqrt(t2);
        dist(:,j)= t2;
    end
    
    
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

function eta=calculate_eta(points,u_mat,centers,no_clusters,m)
    u_mat=u_mat.^m;
    eta=zeros(no_clusters,1);
    for j=1:no_clusters
        num=0;
        denom=0;
        for i = 1: size(points,1)
            num=num+u_mat(i,j)*(norm(points(i)-centers(j))^2);
            denom=denom+u_mat(i,j);
        end
        eta(j)=(num/denom);
    end
end
