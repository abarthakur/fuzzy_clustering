# Fuzzy Clustering Algorithms based on K-means

This repo is a collection of fuzzy clustering algorithms, based on (and including) the k-means clustering algorithm. It is implemented in MATLAB.

The algorithms implemented are as follows-

## K-Means 

```[ centers,labels,no_iterations,others] = k_means( points,no_clusters,epsilon,MAX_ITER,init_option,params )```

where 
* **points** = list of feature vectors
* **no_clusters** = K
* **epsilon** = stopping criterion
* **MAX_ITER** = maximum number of iterations
* **init_option** = option for choosing initial centers
      ** 1 (Forgy method : select K random centers)
      ** 2 (Random Partition : randomly partition **points**, and calculate centers)
      ** 3 (Use centers passed in **params**)
* **params** = structure containing additional parameters
      ** init_centers = initial centers

and
* **centers** = list of final K centers of clusters
* **labels** = final partition of points
* **others** = struct containing other useful values
      ** save_centers = list intermediate ( and initial + final) centers

## Fuzzy C-Means       
[Bezdek et al.(1984)](http://www.sciencedirect.com/science/article/pii/0098300484900207)

```[ centers,u_mat,labels,no_iterations,others] = fc_means( points,no_clusters ,m_val,epsilon,MAX_ITER,init_option,params )```

where 
* **m_val** = m
* **init_option** = option for choosing initial centers
      ** 1 (initialize memberships randomly, and calculate centers)
      ** 2 (Use centers passed in **params**)
* **params** = structure containing additional parameters
      ** init_centers = initial centers

and
* **centers** = list of final K centers of clusters
* **u_mat** = final fuzzy partition
* **others** = struct containing other useful values
      ** save_centers = list intermediate ( and initial + final) centers
      

## Possibilistic C-Means
[Krishnapuram,Keller(1996)](http://ieeexplore.ieee.org/document/531779/?denied)

```[ centers,u_mat,labels,no_iterations,others] = pc_means( points,no_clusters ,m_val ,K_val ,alpha, epsilon, MAX_ITER, init_option, params)```

where 
* **init_option** = option for choosing initial centers
      ** 1 (initialize centers and memberships from FCM)
      ** 2 (initialize memberships randomly, and calculate centers(using FCM memberships))
      ** 2 (Use centers and memberships passed in **params**)
* **params** = structure containing additional parameters
      ** init_centers = initial centers
      ** init_u = initial memberships
      ** fcm_m_val = value of m required when initializing with FCM
      ** pcm_variation =
          *** "usual" : usual PCM memberships used
          *** "log" : exponential memberships used (when objective function contains log term)
      ** eta = initial value for eta
      ** eta_option =
          *** 1 (explicitly passed eta is used)
          *** 2 (eta is updated after every iteration (note : this is unstable))
          *** 3 (eta is calculated from inital FCM memberships)
          *** 4 (PCM is run once with inital estimate of eta from FCM memberships, then run again with refined value of eta calculated from memberships obtained after first iteration)

and
* **centers** = list of final K centers of clusters
* **u_mat** = final fuzzy partition
* **others** = struct containing other useful values
      ** save_centers = list intermediate ( and initial + final) centers
      ** fcm_save_centers = list intermediate ( and initial + final) centers of FCM initialization
      ** eta = final value of eta










