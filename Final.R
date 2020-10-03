install.packages("class")
install.packages("MASS")
install.packages("kohonen")

library(class)
library(MASS)
library(kohonen)


Q_100_1=list()
Q_100_2=list()
Q_100_3=list()
Q_100_4=list()

S=sample(1:1000000000,size=100,replace=FALSE)
simulation=matrix(0,200,4)


for(j in 1:100)
{
  s=S[j]
  #Group 1
  set.seed(s)
  Sigma <- diag(c(1,1,1,1))
  simulation[1:50,]=mvrnorm(n = 50, c(-2.5, -(2.5*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 2
  
  simulation[51:100,]=mvrnorm(n = 50, c(2.5, -(2.5*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 3
  
  simulation[101:150,]=mvrnorm(n = 50, c(0, -(2.5*2*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 4
  
  simulation[151:200,]=mvrnorm(n = 50, c(0, 0, (2.5*6*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  
  #모의실험 분석
  
  iris.sc=simulation
  
  ######################################################################################################
  ################################## 1. 사각그리드 / 가상노드 있는 경우 ################################
  ######################################################################################################
  
  Aind=Bind=NULL
  INDEX_rectangular = function(r)
  {
    ind = c(1 : r^2)
    ## Index Set B -> no artificial nodes generated 
    
    Bind = ind[!((ind == 1) | (ind == r) | (ind == r^2-r+1) | (ind == r^2) | (ind > 1 & ind < r ) |
                   (ind > r^2-r+1 & ind < r^2 ) | (ind %in% c(seq(2*r,r^2-r,by=r))) |
                   (ind %in% c(seq(r+1,r^2-2*r+1,by=r))) )]
    
    Aind = ind[-Bind]
    
    ind_mat = matrix(ind, nrow = r, ncol = r, byrow = TRUE)
    index_list = list()
    for(i in 1 : r^2)
      index_list[[i]] = list()
    
    for(j in 1 : r^2)
    {
      index_list[[j]]$adj_index = EVEN(j, r) 
      index_list[[j]]$directions = DIRECTION_TRIANGLE(j, EVEN(j, r)) 
    }
    
    EVEN = function(c, r)
    {
      adjacent = c(c - 1, c + r - 1, c + r , c + r + 1 , c + 1, c - r + 1, c - r, c - r - 1)
      return(matrix(adjacent, nrow = 2, byrow = TRUE))
    }
    
    
    ## Split odd r case and even r case
    
    index_list = EVEN_R(index_list, Aind, r)
    
    ## Renew direction triangle matrix
    for(j in 1 : length(index_list))
      index_list[[j]]$directions = DIRECTION_TRIANGLE(j, index_list[[j]]$adj_index)
    return(index_list)
  }
  
  ## 4개의 direction을 구축하는 rectangular 의 꼭지점들 저장.
  DIRECTION_TRIANGLE = function(j, index_mat) #index_mat=EVEN(c,r)
  {
    direction_mat = matrix(0, nrow = 4, ncol = 4)
    direction_mat[, 1] = j
    ## 4 directions
    direction_mat[1, 2 : 4] = index_mat[1, 1 : 3]
    direction_mat[2, 2 : 4] = as.vector(t(index_mat))[3:5]
    direction_mat[3, 2 : 4] = as.vector(t(index_mat))[c(7:8,1)]
    direction_mat[4, 2 : 4] = index_mat[2, 1 : 3]
    
    return(direction_mat)
  }
  
  
  
  # r이 짝수일 때 artificial node 삽입. 하지만 짝수는 없긴함
  EVEN_R = function(index_list, Aind, r)
  {
    left_bounds = sort(Aind[which(round((Aind - 1) / r) == (Aind - 1) / r)])
    right_bounds = sort(Aind[which((floor(Aind / r) == ceiling(Aind / r)))])
    upper_bounds = sort(Aind[Aind != r^2 & Aind > (r^2 - r + 1)])
    lower_bounds = sort(Aind[Aind != 1 & Aind < r])
    
    
    # 5 artificial node case TYPE1~4
    five = c(1,r,r^2-r+1,r^2)
    index_list[[five[1]]]$adj_index = TYPE1(index_list[[five[1]]]$adj_index)
    index_list[[five[2]]]$adj_index = TYPE2(index_list[[five[2]]]$adj_index)
    index_list[[five[3]]]$adj_index = TYPE3(index_list[[five[3]]]$adj_index)
    index_list[[five[4]]]$adj_index = TYPE4(index_list[[five[4]]]$adj_index)
    
    
    
    ## 3 artifiacial case TYPE5_1~TYPE5_4
    # lower three case TYPE TYPE5_1
    lower_three = lower_bounds
    for (j in 1 : length(lower_three))
    {
      node = lower_three[j]
      index_list[[node]]$adj_index = TYPE5_1(index_list[[node]]$adj_index)
    }
    # upper three case TYPE TYPE5_2
    upper_three = upper_bounds
    for (j in 1 : length(upper_three))
    {
      node = upper_three[j]
      index_list[[node]]$adj_index = TYPE5_2(index_list[[node]]$adj_index)
    }
    # left three case TYPE TYPE5_3
    left_three = left_bounds[-c(1,length(left_bounds))]
    for (j in 1 : length(left_three))
    {
      node = left_three[j]
      index_list[[node]]$adj_index = TYPE5_3(index_list[[node]]$adj_index)
    }
    # right three case TYPE TYPE5_4
    right_three = right_bounds[-c(1,length(right_bounds))]
    for (j in 1 : length(right_three))
    {
      node = right_three[j]
      index_list[[node]]$adj_index = TYPE5_4(index_list[[node]]$adj_index)
    }
    
    return(index_list)
  }
  
  
  # =================================================================
  # Operations - artificial node 생성 유형별 함수.
  # =================================================================
  
  ## 전역변수로 지정되어 지정한 값대로 artificial node 위치에 입력된다.
  
  # 5 artificial node cases. 
  # TYPE1 : 맨 아랫줄 첫번째 노드 ,
  # TYPE2 : 맨 아랫줄 마지막 노드 ,
  # TYPE3 : 맨 윗줄 첫번째 노드 , 
  # TYPE4 : 맨 윗줄 마지막 노드.
  
  TYPE1 = function(adj)
  {
    vadj = as.vector(t(adj))
    vadj[1] = -1
    vadj[2] = -2
    vadj[6] = -6
    vadj[7] = -7
    vadj[8] = -8
    adj = matrix(vadj, nrow = 2, byrow = TRUE)
    adj
  }
  TYPE2 = function(adj)
  {
    vadj = as.vector(t(adj))
    vadj[4] = -4
    vadj[5] = -5
    vadj[6] = -6
    vadj[7] = -7
    vadj[8] = -8
    adj = matrix(vadj, nrow = 2, byrow = TRUE)
    adj
  }
  TYPE3 = function(adj)
  {
    vadj = as.vector(t(adj))
    vadj[1] = -1
    vadj[2] = -2
    vadj[3] = -3
    vadj[4] = -4
    vadj[8] = -8
    adj = matrix(vadj, nrow = 2, byrow = TRUE)
    adj
  }
  
  TYPE4 = function(adj)
  {
    vadj = as.vector(t(adj))
    vadj[2] = -2
    vadj[3] = -3
    vadj[4] = -4
    vadj[5] = -5
    vadj[6] = -6
    adj = matrix(vadj, nrow = 2, byrow = TRUE)
    adj
  }
  
  # 3 artificial node cases TYPE5_1, TYPE5_2, TYPE5_3, TYPE5_4
  # TYPE5_1 : 아랫줄, 
  # TYPE5_2 : 윗줄, 
  # TYPE5_3 : 왼쪽줄, 
  # TYPE5_4 : 오른쪽줄
  TYPE5_1 = function(adj)
  {
    vadj = as.vector(t(adj))
    vadj[6] = -6
    vadj[7] = -7
    vadj[8] = -8
    adj = matrix(vadj, nrow = 2, byrow = TRUE)
    adj
  }
  
  TYPE5_2 = function(adj)
  {
    vadj = as.vector(t(adj))
    vadj[2] = -2
    vadj[3] = -3
    vadj[4] = -4
    adj = matrix(vadj, nrow = 2, byrow = TRUE)
    adj
  }
  
  TYPE5_3 = function(adj)
  {
    vadj = as.vector(t(adj))
    vadj[1] = -1
    vadj[2] = -2
    vadj[8] = -8
    adj = matrix(vadj, nrow = 2, byrow = TRUE)  
    adj
  }
  
  TYPE5_4 = function(adj)
  {
    vadj = as.vector(t(adj))
    vadj[4] = -4
    vadj[5] = -5
    vadj[6] = -6
    adj = matrix(vadj, nrow = 2, byrow = TRUE)  
    adj
  }
  
  EVEN = function(c, r)
  {
    adjacent = c(c - 1, c + r - 1, c + r , c + r + 1 , c + 1, c - r + 1, c - r, c - r - 1)
    return(matrix(adjacent, nrow = 2, byrow = TRUE))
  }
  
  
# =================================================================
  
  Q=matrix(0,200,4)
  Q_final=matrix(0,12,1)
  set.seed(s)
  for(r in 3:14)
  {
    
    aa = INDEX_rectangular(r)
    
    iris.grid = function(r) somgrid(xdim =r, ydim=r, topo="rectangular")
    
    # build model
    iris.som = som(iris.sc, grid=iris.grid(r), rlen=100, alpha=c(0.25,0.001))
    
    #pd 들의 max를 구하기 위해 비교하는 함수 likelihood 반환 x, k, r, l, c : vector (c는 대각선벡터)
    pd=function(x,k,l,r,c){exp(-0.5*euc.dist(x,k))+exp(-0.5*euc.dist(x,l))+exp(-0.5*euc.dist(x,r))+exp(-0.5*euc.dist(x,c))}
    
    # function "euc.dist" : 벡터 x1과 x2 사이의 유클리디안제곱거리
    euc.dist <- function(x1, x2) (sum((x1 - x2) ^ 2)) # 
    
    
    #function "arrow" : rxr 그리드 SOM에서 i 번째 객체의 8개방향에 대한 가상벡터 (승자벡터+방향벡터) 8,7,6,5,4,3,2,1 방향 순
    arrow=function(r,i)rbind(iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][1,]-iris.som$codes[[1]][r+2,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][2,]-iris.som$codes[[1]][r+2,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][2,]-iris.som$codes[[1]][r+1,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][2,]-iris.som$codes[[1]][1,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][r+2,]-iris.som$codes[[1]][1,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][r+2,]-iris.som$codes[[1]][2,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][r+1,]-iris.som$codes[[1]][2,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][1,]-iris.som$codes[[1]][2,])
    
    # function "fincodes" : 원래 생성된 m개의 중량벡터 뒤에 8개의가상벡터연결
    fincodes=function(r,i)rbind(as.matrix(iris.som$codes[[1]]), arrow(r,i))
    
    
    # r노드일때, i 번째 객체의 j번째 방향애 대한 k(node index) 함수 k=1(승자노드), k=2(왼쪽노드),k=3(대각선 노드), k=4(오른쪽노드) 중량벡터
    krl=function(i,j,k,r)fincodes(r,i)[ifelse(aa[[iris.som$unit.classif[i] ]]$directions[j,k]>0,
                                              aa[[iris.som$unit.classif[i] ]]$directions[j,k],
                                              nrow(fincodes(r,i))+aa[[iris.som$unit.classif[i] ]]$directions[j,k]+1),] 
    
    
    #각 개체에 따라 8개의 방향에 대하여 최대인 방향찾고 그 방향 나타내는 index 저장
    
    direction_vector = rep(0, 200)
    
    
    for(i in 1:200)
    {
      vec = function(i,r){c(pd(iris.sc[i,],krl(i,1,1,r),krl(i,1,2,r),krl(i,1,4,r),krl(i,1,3,r)), 
                            pd(iris.sc[i,],krl(i,2,1,r),krl(i,2,2,r),krl(i,2,4,r),krl(i,2,3,r)),
                            pd(iris.sc[i,],krl(i,3,1,r),krl(i,3,2,r),krl(i,3,4,r),krl(i,3,3,r)),
                            pd(iris.sc[i,],krl(i,4,1,r),krl(i,4,2,r),krl(i,4,4,r),krl(i,4,3,r)))}
      maxi_index = which.max(vec(i,r))
      direction_vector[i] = maxi_index
      maxi = max(vec(i,r))
      cat("max direction = ", maxi_index, "value = ", maxi, "\n")
    } 
    
    
    #finwight = 각 개체의 최종 개체가 표현될 승자노드,왼쪽노드,오른쪽 노드 중량벡터의 list
    
    finweight=list()
    
    EACHVECTOR=function(r,i){
      weight_matrix=matrix(0,4,4)
      weight_matrix[1,]=fincodes(r,i)[ifelse(aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1] >0,
                                             aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1] ,
                                             nrow(fincodes(r,i))+aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1]+1),] 
      
      weight_matrix[2,]=fincodes(r,i)[ifelse(aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2] >0,
                                             aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2] ,
                                             nrow(fincodes(r,i))+aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2]+1),] 
      
      weight_matrix[3,]=fincodes(r,i)[ifelse(aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3] >0,
                                             aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3] ,
                                             nrow(fincodes(r,i))+aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3]+1),] 
      
      weight_matrix[4,]=fincodes(r,i)[ifelse(aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],4] >0,
                                             aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],4] ,
                                             nrow(fincodes(r,i))+aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],4]+1),] 
      
      return(weight_matrix)
    }
    
    for(i in 1:200){
      finweight[[i]]= EACHVECTOR(r,i) 
    }
    
    ##  가능도함수 생성 및 p_r, p_l 함수 만들기
    lik=function(i,w,b){
      exp((-0.5/b)*t(i-w)%*%(i-w))
    }
    # function "p_r, p_l" : cxc SOM에서 
    # i 번째 데이터에 대한 b=beta 가능도
    p_h=function(i,b) {
      lik(iris.sc[i,],finweight[[i]][1,],b)/
        (lik(iris.sc[i,],finweight[[i]][4,],b)+lik(iris.sc[i,],finweight[[i]][1,],b)) }
    
    p_v=function(i,b) {
      lik(iris.sc[i,],finweight[[i]][1,],b)/
        (lik(iris.sc[i,],finweight[[i]][2,],b)+lik(iris.sc[i,],finweight[[i]][1,],b)) }
    
    ##########coord 함수사용해서 beta구하기 여기서 a=p_r, b=p_l , k=승자노드, r=오른쪽노드, l=왼쪽노드
    coord=function(i,b){
      p_h(i,b)*p_v(i,b)*finweight[[i]][1,]+
        (1-p_h(i,b))*p_v(i,b)*finweight[[i]][4,]+
        p_h(i,b)*(1-p_v(i,b))*finweight[[i]][2,]+
        (1-p_h(i,b))*(1-p_v(i,b))*finweight[[i]][3,]
    }
    
    #Q_beta 함수
    
    Q_beta=function(b) {
      for(i in 1:200) {
        Q[i,] = (iris.sc[i,]-coord(i,b))
      } 
      return(sum(Q^2))
    }
    Q_final[r-2,1]= optimize(Q_beta,c(0,1000), tol=0.0001)$objective
  }
  
  Q_final=cbind(Q_final,matrix(c(3,4,5,6,7,8,9,10,11,12,13,14),12,1))
  
  Q_100_1[[j]]=Q_final
  
}

#final value of overall 100 simulation mean result may be wrote as below.
Y <- do.call(cbind, Q_100_1)
Y <- array(Y, dim=c(dim(Q_100_1[[1]]), length(Q_100_1)))
Q_1=matrix(0,12,2)
Q_1=apply(Y, c(1, 2), mean)


plot(Q_1[,2],Q_1[,1], xlab="r", ylab="Q_beta",frame = FALSE,axes=FALSE, col = "blue",type="l")
box();axis(1);axis(2)
title("Simulation_1 : Q of rectangular_gridsize ( r ) = 3 : 14")


######################################################################################################
################################## 2. 사각그리드 / 가상노드 없는 경우 ################################
######################################################################################################

Q=matrix(0,200,4)
Q_final=matrix(0,12,1)

for(j in 1:100)
{
  s=S[j]
  #Group 1
  set.seed(s)
  Sigma <- diag(c(1,1,1,1))
  simulation[1:50,]=mvrnorm(n = 50, c(-2.5, -(2.5*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 2
  
  simulation[51:100,]=mvrnorm(n = 50, c(2.5, -(2.5*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 3
  
  simulation[101:150,]=mvrnorm(n = 50, c(0, -(2.5*2*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 4
  
  simulation[151:200,]=mvrnorm(n = 50, c(0, 0, (2.5*6*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  
  #모의실험 분석
  
  iris.sc=simulation
  
  set.seed(s)
  for(r in 3:14)
  {
    
    aa = INDEX_rectangular(r)
    
    iris.grid = function(r) somgrid(xdim =r, ydim=r, topo="rectangular")
    
    # build model
    iris.som = som(iris.sc, grid=iris.grid(r), rlen=100, alpha=c(0.25,0.001))
    
    #pd 들의 max를 구하기 위해 비교하는 함수 likelihood 반환 x, k, r, l, c : vector (c는 대각선벡터)
    pd=function(x,k,l,r,c){exp(-0.5*euc.dist(x,k))+exp(-0.5*euc.dist(x,l))+exp(-0.5*euc.dist(x,r))+exp(-0.5*euc.dist(x,c))}
    
    # function "euc.dist" : 벡터 x1과 x2 사이의 유클리디안제곱거리
    euc.dist <- function(x1, x2) (sum((x1 - x2) ^ 2)) # 
    
    
    # function "fincodes" : 원래 생성된 m개의 중량벡터 뒤에 8개의가상벡터연결
    fincodes=function(r,i)rbind(as.matrix(iris.som$codes[[1]]), matrix(c(100,100,100,100),1,4))
    
    
    # r노드일때, i 번째 객체의 j번째 방향애 대한 k(node index) 함수 k=1(승자노드), k=2(왼쪽노드),k=3(대각선 노드), k=4(오른쪽노드) 중량벡터
    krl=function(i,j,k,r)fincodes(r,i)[ifelse(aa[[iris.som$unit.classif[i] ]]$directions[j,k]>0,
                                              aa[[iris.som$unit.classif[i] ]]$directions[j,k],
                                              nrow(fincodes(r,i))),] 
    
    #각 개체에 따라 8개의 방향에 대하여 최대인 방향찾고 그 방향 나타내는 index 저장
    
    direction_vector = rep(0, 200)
    
    for(i in 1:200)
    {
      vec = function(i,r){c(pd(iris.sc[i,],krl(i,1,1,r),krl(i,1,2,r),krl(i,1,4,r),krl(i,1,3,r)), 
                            pd(iris.sc[i,],krl(i,2,1,r),krl(i,2,2,r),krl(i,2,4,r),krl(i,2,3,r)),
                            pd(iris.sc[i,],krl(i,3,1,r),krl(i,3,2,r),krl(i,3,4,r),krl(i,3,3,r)),
                            pd(iris.sc[i,],krl(i,4,1,r),krl(i,4,2,r),krl(i,4,4,r),krl(i,4,3,r)))}
      maxi_index = which.max(vec(i,r))
      direction_vector[i] = maxi_index
      maxi = max(vec(i,r))
      cat("max direction = ", maxi_index, "value = ", maxi, "\n")
    } 
    
    #finwight = 각 개체의 최종 개체가 표현될 승자노드,왼쪽노드,오른쪽 노드 중량벡터의 list
    
    finweight=list()
    
    EACHVECTOR=function(r,i){
      weight_matrix=matrix(0,4,4)
      weight_matrix[1,]=fincodes(r,i)[ifelse(aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1] >0,
                                             aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1] ,
                                             nrow(fincodes(r,i))+aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1]+1),] 
      
      weight_matrix[2,]=fincodes(r,i)[ifelse(aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2] >0,
                                             aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2] ,
                                             nrow(fincodes(r,i))+aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2]+1),] 
      
      weight_matrix[3,]=fincodes(r,i)[ifelse(aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3] >0,
                                             aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3] ,
                                             nrow(fincodes(r,i))+aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3]+1),] 
      
      weight_matrix[4,]=fincodes(r,i)[ifelse(aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],4] >0,
                                             aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],4] ,
                                             nrow(fincodes(r,i))+aa[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],4]+1),] 
      
      return(weight_matrix)
    }
    
    
    for(i in 1:200){
      finweight[[i]]= EACHVECTOR(r,i) #r 을 매번 바꿔주기 
    }
    
    ##  가능도함수 생성 및 p_r, p_l 함수 만들기
    lik=function(i,w,b){
      exp((-0.5/b)*t(i-w)%*%(i-w))
    }
    # function "p_r, p_l" : cxc SOM에서 
    # i 번째 데이터에 대한 b=beta 가능도
    p_h=function(i,b) {
      lik(iris.sc[i,],finweight[[i]][1,],b)/
        (lik(iris.sc[i,],finweight[[i]][4,],b)+lik(iris.sc[i,],finweight[[i]][1,],b)) }
    
    p_v=function(i,b) {
      lik(iris.sc[i,],finweight[[i]][1,],b)/
        (lik(iris.sc[i,],finweight[[i]][2,],b)+lik(iris.sc[i,],finweight[[i]][1,],b)) }
    
    ##########coord 함수사용해서 beta구하기 여기서 a=p_r, b=p_l , k=승자노드, r=오른쪽노드, l=왼쪽노드
    coord=function(i,b){
      p_h(i,b)*p_v(i,b)*finweight[[i]][1,]+
        (1-p_h(i,b))*p_v(i,b)*finweight[[i]][4,]+
        p_h(i,b)*(1-p_v(i,b))*finweight[[i]][2,]+
        (1-p_h(i,b))*(1-p_v(i,b))*finweight[[i]][3,]
    }
    
    #Q_beta 함수
    
    Q_beta=function(b) {
      for(i in 1:200) {
        Q[i,] = (iris.sc[i,]-coord(i,b))
      } 
      return(sum(Q^2))
    }
    Q_final[r-2,1]= optimize(Q_beta,c(0,1000), tol=0.0001)$objective
  }
  
  Q_final=cbind(Q_final,matrix(c(3,4,5,6,7,8,9,10,11,12,13,14),12,1))
  
  Q_100_2[[j]]=Q_final
  
}
for(i in 1:100){
  Q_100_2[[i]]=Q_100_2[[i]][,1:2]}

#final value of overall 100 simulation mean result may be wrote as below.
Y2 <- do.call(cbind, Q_100_2)
Y2 <- array(Y2, dim=c(dim(Q_100_2[[1]]), length(Q_100_2)))
Q_2=matrix(0,12,2)
Q_2=apply(Y2, c(1, 2), mean)



plot(Q_2[,2],Q_2[,1], xlab="r", ylab="Q_beta",frame = FALSE,axes=FALSE, col = "blue",type="l")
box();axis(1);axis(2)
title("Simulation_2 : Q of rectangular_gridsize ( r ) = 3 : 14")




######################################################################################################
################################## 3. 육각그리드 / 가상노드 있는 경우 ################################
######################################################################################################


INDEX = function(r)
{
  ind = c(1 : r^2)
  ## Index Set B -> no artificial nodes generated
  Bind = ind[(ind > r + 1) & (ind < r * (r - 1)) &
               !( (floor(ind / r) == ind / r) |
                    (floor((ind - 1) / r) == (ind - 1) / r) )]
  Aind = ind[-Bind]
  ind_mat = matrix(ind, nrow = r, ncol = r, byrow = TRUE)
  index_list = list()
  for(i in 1 : r^2)
    index_list[[i]] = list()
  even_row = 2 * c(1 : floor(NROW(ind_mat) / 2))
  even_row_index = as.vector(ind_mat[even_row, ])
  for(i in 1 : length(even_row_index))
  {
    j = even_row_index[i]
    index_list[[j]]$adj_index = EVEN(j, r) #이건무슨함수지
    index_list[[j]]$directions = DIRECTION_TRIANGLE(j, EVEN(j, r)) #이것도무슨함수지
  }
  odd_row_index = as.vector(ind_mat[-even_row, ])
  for(i in 1 : length(odd_row_index))
  {
    j = odd_row_index[i]
    index_list[[j]]$adj_index = ODD(j, r)
    index_list[[j]]$directions = DIRECTION_TRIANGLE(j, ODD(j, r))
  }
  ## Split odd r case and even r case
  if (floor(r / 2) == ceiling(r / 2))
    index_list = EVEN_R(index_list, Aind, r)
  else
    index_list = ODD_R(index_list, Aind, r)
  ## Renew direction triangle matrix
  for(j in 1 : length(index_list))
    index_list[[j]]$directions = DIRECTION_TRIANGLE(j, index_list[[j]]$adj_index)
  return(index_list)
}

# =================================================================
# Operations
# =================================================================

## 홀수번째 행인지, 짝수번째 행인지 따라서 육각형 index를 계산하는 함수
EVEN = function(c, r)
{
  adjacent = c(c - 1, c + r - 1, c + r, c - r - 1, c - r, c + 1)
  return(matrix(adjacent, nrow = 2, byrow = TRUE))
}


ODD = function(c, r)
{
  adjacent = c(c - 1, c + r, c + r + 1, c - r, c - r + 1, c + 1)
  return(matrix(adjacent, nrow = 2, byrow = TRUE))
}

## 6개의 direction을 구축하는 삼각형의 꼭지점들 저장.
DIRECTION_TRIANGLE = function(j, index_mat)
{
  direction_mat = matrix(0, nrow = 6, ncol = 3)
  direction_mat[, 1] = j
  ## 6 directions
  direction_mat[1, 2 : 3] = index_mat[1, 1 : 2]
  direction_mat[2, 2 : 3] = index_mat[1, 2 : 3]
  direction_mat[3, 2 : 3] = index_mat[, 3]
  direction_mat[4, 2 : 3] = rev(index_mat[, 1])
  direction_mat[5, 2 : 3] = index_mat[2, 2 : 1]
  direction_mat[6, 2 : 3] = index_mat[2, 3 : 2]
  return(direction_mat)
}

# r이 짝수일 때 artificial node 삽입
EVEN_R = function(index_list, Aind, r)
{
  left_bounds = sort(Aind[which(round((Aind - 1) / r) == (Aind - 1) / r)])
  right_bounds = sort(Aind[which((floor(Aind / r) == ceiling(Aind / r)))])
  upper_bounds = sort(Aind[which(Aind != r^2 & Aind > (r^2 - r + 1))])
  lower_bounds = sort(Aind[which(Aind != 1 & Aind < r)])
  # 4 artificial node case
  fours = c(min(right_bounds), max(left_bounds))
  index_list[[fours[1]]]$adj_index = TYPE1(index_list[[fours[1]]]$adj_index)
  index_list[[fours[2]]]$adj_index = TYPE2_1(index_list[[fours[2]]]$adj_index)
  # 3 artificial node case
  # Start point, and end point
  startpoint = min(Aind)
  index_list[[startpoint]]$adj_index = TYPE3(index_list[[startpoint]]$adj_index)
  endpoint = max(Aind)
  index_list[[endpoint]]$adj_index = TYPE4_1(index_list[[endpoint]]$adj_index)
  ## Left side three, except 4 node case(제일 윗줄 노드)
  left_three = left_bounds[c(2 * (1 : (r / 2 - 1)))]
  for (j in 1 : length(left_three))
  {
    node = left_three[j]
    index_list[[node]]$adj_index = TYPE5(index_list[[node]]$adj_index)
  }
  ## Right side three, except 4 node case(제일 밑의 노드)
  right_three = right_bounds[c(2 * c(1 : (r / 2)) - 1)][-1]
  for (j in 1 : length(right_three))
  {
    node = right_three[j]
    index_list[[node]]$adj_index = TYPE6(index_list[[node]]$adj_index)
  }
  ## Two artificial node case (윗부분, 아랫부분)
  for(j in 1 : length(upper_bounds))
  {
    node = upper_bounds[j]
    index_list[[node]]$adj_index = TYPE7(index_list[[node]]$adj_index)
  }
  for(j in 1 : length(lower_bounds))
  {
    node = lower_bounds[j]
    index_list[[node]]$adj_index = TYPE8(index_list[[node]]$adj_index)   
  }
  ## One artificial node case
  left_one = left_bounds[c(2 * c(1 : (r / 2)) - 1)][-1]
  for(j in 1 : length(left_one))
  {
    node = left_one[j]
    index_list[[node]]$adj_index = TYPE9(index_list[[node]]$adj_index)   
  }
  right_one = right_bounds[c(2 * (1 : (r / 2 - 1)))]
  for(j in 1 : length(right_one))
  {
    node = right_one[j]
    index_list[[node]]$adj_index = TYPE10(index_list[[node]]$adj_index)   
  }
  return(index_list)
}

# r이 홀수일 때 artificial node 삽입
ODD_R = function(index_list, Aind, r)
{
  left_bounds = sort(Aind[which(round((Aind - 1) / r) == (Aind - 1) / r)])
  right_bounds = sort(Aind[which((floor(Aind / r) == ceiling(Aind / r)))])
  upper_bounds = sort(Aind[which(Aind != r^2 & Aind > (r^2 - r + 1))])
  lower_bounds = sort(Aind[which(Aind != 1 & Aind < r)])
  # 4 artificial node case
  fours = c(r, r^2)
  index_list[[fours[1]]]$adj_index = TYPE1(index_list[[fours[1]]]$adj_index)
  index_list[[fours[2]]]$adj_index = TYPE2_2(index_list[[fours[2]]]$adj_index)
  
  # 3 artificial node case - 좌측 맨위와 맨아래
  threes = c(min(left_bounds), max(left_bounds))
  index_list[[threes[1]]]$adj_index = TYPE3(index_list[[threes[1]]]$adj_index)
  index_list[[threes[2]]]$adj_index = TYPE4_2(index_list[[threes[2]]]$adj_index)
  ## Left side three(특별히 위에서 정의한 노드는 제외하고 짝수 행에서 튀어나온 것들만)
  left_index_three = c(2 * (1 : floor(r / 2)))
  left_three = left_bounds[left_index_three]
  for (j in 1 : length(left_three))
  {
    node = left_three[j]
    index_list[[node]]$adj_index = TYPE5(index_list[[node]]$adj_index)
  }
  ## Right side three, except 4 node case(제일 밑의 노드)
  right_index_three = c(1 : r)[-left_index_three]
  right_index_three = right_index_three[-c(1, length(right_index_three))]
  ## 우측의 경우 첫번째와 마지막 노드는 4개 가지고 있어서 이들을 빼야 함.
  right_three = right_bounds[right_index_three]
  for (j in 1 : length(right_three))
  {
    node = right_three[j]
    index_list[[node]]$adj_index = TYPE6(index_list[[node]]$adj_index)
  }
  ## Two artificial node case (윗부분, 아랫부분)
  for(j in 1 : length(upper_bounds))
  {
    node = upper_bounds[j]
    index_list[[node]]$adj_index = TYPE7(index_list[[node]]$adj_index)
  }
  for(j in 1 : length(lower_bounds))
  {
    node = lower_bounds[j]
    index_list[[node]]$adj_index = TYPE8(index_list[[node]]$adj_index)   
  }
  ## One artificial node case
  left_one = left_bounds[c(2 * c(1 : floor(r / 2)) - 1)][-1]
  for(j in 1 : length(left_one))
  {
    node = left_one[j]
    index_list[[node]]$adj_index = TYPE9(index_list[[node]]$adj_index)   
  }
  right_one = right_bounds[c(2 * (1 : floor(r / 2)))]
  for(j in 1 : length(right_one))
  {
    node = right_one[j]
    index_list[[node]]$adj_index = TYPE10(index_list[[node]]$adj_index)   
  }
  return(index_list)
}

# =================================================================
# Operations - artificial node 생성 유형별 함수.
# =================================================================
## ani = artificial node index.
## 전역변수로 지정되어 지정한 값대로 artificial node 위치에 입력된다.


# 4 artificial node cases. TYPE1 : 맨 아랫줄 마지막 노드,
#TYPE2_1, TYPE2_2 = r이 짝수면 맨 윗줄 첫번째 노드(2_1),
# r이 홀수면 맨 윗줄 마지막 노드(2_2)
TYPE1 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[3] = -3
  vadj[4] = -4
  vadj[5] = -5
  vadj[6] = -6
  adj = matrix(vadj, nrow = 2, byrow = TRUE)
  adj
}
TYPE2_1 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[1] = -1
  vadj[2] = -2
  vadj[3] = -3
  vadj[4] = -4
  adj = matrix(vadj, nrow = 2, byrow = TRUE)
  adj
}
TYPE2_2 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[2] = -2
  vadj[3] = -3
  vadj[5] = -5
  vadj[6] = -6
  adj = matrix(vadj, nrow = 2, byrow = TRUE)
  adj
}

# 3 artificial node cases인데 따로 처리해야 하는 경우. TYPE3 -> first node, 
# TYPE4_1, TYPE4_2 = r이 짝수일 떄 가장 마지막 노드(4_1), r이 홀수일 때 맨 위줄
# 첫번째 노드(4_2)
TYPE3 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[1] = -1
  vadj[4] = -4
  vadj[5] = -5
  adj = matrix(vadj, nrow = 2, byrow = TRUE)
  adj
}

TYPE4_1 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[2] = -2
  vadj[3] = -3
  vadj[6] = -6
  adj = matrix(vadj, nrow = 2, byrow = TRUE)
  adj
}
TYPE4_2 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[1] = -1
  vadj[2] = -2
  vadj[3] = -3
  adj = matrix(vadj, nrow = 2, byrow = TRUE)  
  adj
}

# 좌, 우에서 3개씩 생성될 때(5 -> left, 6 -> right)
TYPE5 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[1] = -1
  vadj[2] = -2
  vadj[4] = -4
  adj = matrix(vadj, nrow = 2, byrow = TRUE)  
  adj
}
TYPE6 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[3] = -3
  vadj[5] = -5
  vadj[6] = -6
  adj = matrix(vadj, nrow = 2, byrow = TRUE)  
  adj
}
# 2개 노드 생성. 7 -> 맨밑줄에서, 8 -> 맨 윗줄에서
TYPE8 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[4] = -4
  vadj[5] = -5
  adj = matrix(vadj, nrow = 2, byrow = TRUE)  
  adj
}
TYPE7 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[2] = -2
  vadj[3] = -3
  adj = matrix(vadj, nrow = 2, byrow = TRUE)  
  adj
}
# 1개 노드 생성. 9 -> 좌, 10 -> 우
TYPE9 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[1] = -1
  adj = matrix(vadj, nrow = 2, byrow = TRUE)  
  adj
}
TYPE10 = function(adj)
{
  vadj = as.vector(t(adj))
  vadj[6] = -6
  adj = matrix(vadj, nrow = 2, byrow = TRUE)  
  adj
}

#################################
##### end of the code ###########
#################################


for(j in 1:100)
{
  s=S[j]
  #Group 1
  set.seed(s)
  Sigma <- diag(c(1,1,1,1))
  simulation[1:50,]=mvrnorm(n = 50, c(-2.5, -(2.5*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 2
  
  simulation[51:100,]=mvrnorm(n = 50, c(2.5, -(2.5*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 3
  
  simulation[101:150,]=mvrnorm(n = 50, c(0, -(2.5*2*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 4
  
  simulation[151:200,]=mvrnorm(n = 50, c(0, 0, (2.5*6*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  
  #모의실험 분석
  
  iris.sc=simulation
  
  set.seed(s)
  Q_2=matrix(0,200,4)
  Q_final2=matrix(0,11,1)
  
  for(r in 4:14)
  {
    bb = INDEX(r)
    
    # build grid
    iris.grid = function(r) somgrid(xdim =  r, ydim = r, topo="hexagonal") #이게잘못된건아님
    
    # build model
    iris.som = som(iris.sc, grid=iris.grid(r), rlen=100, alpha=c(0.25,0.001))
    
    #pd 들의 max를 구하기 위해 비교하는 함수 likelihood 반환 x, k, r, l : vector
    pd=function(x,k,l,r){exp(-0.5*euc.dist(x,k))+exp(-0.5*euc.dist(x,l))+exp(-0.5*euc.dist(x,r))}
    
    # function "euc.dist" : 벡터 x1과 x2 사이의 유클리디안제곱거리
    euc.dist <- function(x1, x2) (sum((x1 - x2) ^ 2)) # 
    
    
    #function "arrow" : rxr 그리드 SOM에서 i 번째 객체의 6개방향에 대한 가상벡터 (승자벡터+방향벡터) 6,5,4,3,2,1 방향 순
    arrow=function(r,i)rbind(iris.som$codes[[1]][iris.som$unit.classif[i],]-iris.som$codes[[1]][1,]+iris.som$codes[[1]][2,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]-iris.som$codes[[1]][r+1,]+iris.som$codes[[1]][1,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]-iris.som$codes[[1]][r+2,]+iris.som$codes[[1]][1,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][r+2,]-iris.som$codes[[1]][1,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][r+1,]-iris.som$codes[[1]][1,],
                             iris.som$codes[[1]][iris.som$unit.classif[i],]+iris.som$codes[[1]][1,]-iris.som$codes[[1]][2,])
    
    # function "fincodes" : 원래 생성된 m개의 중량벡터 뒤에 6개의가상벡터연결
    fincodes=function(r,i)rbind(as.matrix(iris.som$codes[[1]]), arrow(r,i))
    
    
    # r노드일때, i 번째 객체의 j번째 방향애 대한 k(node index) 함수 k=1(승자노드), k=2(왼쪽노드),k=3(오른쪽노드) 중량벡터
    krl=function(i,j,k,r)fincodes(r,i)[ifelse(bb[[iris.som$unit.classif[i] ]]$directions[j,k]>0,
                                              bb[[iris.som$unit.classif[i] ]]$directions[j,k],
                                              nrow(fincodes(r,i))+bb[[iris.som$unit.classif[i] ]]$directions[j,k]+1),] 
    
    
    #각 개체에 따라 6개의 방향에 대하여 최대인 방향찾고 그 방향 나타내는 index 저장
    
    direction_vector = rep(0, 200)
    
    
    for(i in 1:200)
    {
      vec = function(i,r){c(pd(iris.sc[i,],krl(i,1,1,r),krl(i,1,2,r),krl(i,1,3,r)), 
                            pd(iris.sc[i,],krl(i,2,1,r),krl(i,2,2,r),krl(i,2,3,r)),
                            pd(iris.sc[i,],krl(i,3,1,r),krl(i,3,2,r),krl(i,3,3,r)),
                            pd(iris.sc[i,],krl(i,4,1,r),krl(i,4,2,r),krl(i,4,3,r)),
                            pd(iris.sc[i,],krl(i,5,1,r),krl(i,5,2,r),krl(i,5,3,r)),
                            pd(iris.sc[i,],krl(i,6,1,r),krl(i,6,2,r),krl(i,6,3,r)))}
      maxi_index = which.max(vec(i,r))
      direction_vector[i] = maxi_index
      maxi = max(vec(i,r))
      cat("max direction = ", maxi_index, "value = ", maxi, "\n")
    } 
    
    
    
    #finwight = 각 개체의 최종 개체가 표현될 승자노드,왼쪽노드,오른쪽 노드 중량벡터의 list
    
    finweight=list()
    
    EACHVECTOR=function(r,i){
      weight_matrix=matrix(0,3,4)
      weight_matrix[1,]=fincodes(r,i)[ifelse(bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1] >0,
                                             bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1] ,
                                             nrow(fincodes(r,i))+bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1]+1),] 
      
      weight_matrix[2,]=fincodes(r,i)[ifelse(bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2] >0,
                                             bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2] ,
                                             nrow(fincodes(r,i))+bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2]+1),] 
      
      weight_matrix[3,]=fincodes(r,i)[ifelse(bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3] >0,
                                             bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3] ,
                                             nrow(fincodes(r,i))+bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3]+1),] 
      return(weight_matrix)
    }
    
    for(i in 1:200){
      finweight[[i]]= EACHVECTOR(r,i) #r 을 매번 바꿔주기 
    }
    
    ##  가능도함수 생성 및 p_r, p_l 함수 만들기
    lik=function(i,w,b){
      exp((-0.5/b)*t(i-w)%*%(i-w))
    }
    # function "p_r, p_l" : cxc SOM에서 
    # i 번째 데이터에 대한 b=beta 가능도
    p_r=function(i,b) {
      lik(iris.sc[i,],finweight[[i]][1,],b)/
        (lik(iris.sc[i,],finweight[[i]][3,],b)+lik(iris.sc[i,],finweight[[i]][1,],b)) }
    
    p_l=function(i,b) {
      lik(iris.sc[i,],finweight[[i]][1,],b)/
        (lik(iris.sc[i,],finweight[[i]][2,],b)+lik(iris.sc[i,],finweight[[i]][1,],b)) }
    
    
    ##########coord 함수사용해서 beta구하기 여기서 a=p_r, b=p_l , k=승자노드, r=오른쪽노드, l=왼쪽노드
    
    # 원래 coord함수
    coord=function(i,b){
      (((2*p_r(i,b)+2*p_l(i,b)-1)*finweight[[i]][1,])/(4*sqrt(3))+
         ((1-2*p_r(i,b)+p_l(i,b))*finweight[[i]][3,])/(2*sqrt(3))+
         ((1+p_r(i,b)-2*p_l(i,b))*finweight[[i]][2,])/(2*sqrt(3)))/(sqrt(3)/4)
    }
    
    #Q_beta 함수
    
    Q_beta=function(b) {
      for(i in 1:200) {
        Q_2[i,] = (iris.sc[i,]-coord(i,b))
      } 
      return(sum(Q_2^2))
    }
    Q_final2[r-3,1]= optimize(Q_beta,c(0,1000), tol=0.0001)$objective
  }
  Q_final2=cbind(Q_final2,matrix(c(4,5,6,7,8,9,10,11,12,13,14),11,1))
  Q_100_3[[j]]=Q_final2
}

#final value of overall 100 simulation mean result may be wrote as below.
Y3 <- do.call(cbind, Q_100_3)
Y3 <- array(Y3, dim=c(dim(Q_100_3[[1]]), length(Q_100_3)))
Q_3=matrix(0,11,2)
Q_3=apply(Y3, c(1, 2), mean)

plot(Q_3[,2],Q_3[,1], xlab="r", ylab="Q_beta",frame = FALSE,axes=FALSE, col = "blue",type="l")
box();axis(1);axis(2)
title("Simulation_3 : Q_hexagonal of gridsize ( r ) = 4 : 14")

######################################################################################################
################################## 4. 육각그리드 / 가상노드 없는 경우 ################################
######################################################################################################

for(j in 1:100)
{
  s=S[j]
  #Group 1
  set.seed(s)
  Sigma <- diag(c(1,1,1,1))
  simulation[1:50,]=mvrnorm(n = 50, c(-2.5, -(2.5*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 2
  
  simulation[51:100,]=mvrnorm(n = 50, c(2.5, -(2.5*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 3
  
  simulation[101:150,]=mvrnorm(n = 50, c(0, -(2.5*2*sqrt(3))/3, (-2.5*2*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  #Group 4
  
  simulation[151:200,]=mvrnorm(n = 50, c(0, 0, (2.5*6*sqrt(2))/(4*sqrt(3)),0), Sigma)
  
  
  #모의실험 분석
  
  iris.sc=simulation
  set.seed(s)
  Q_2=matrix(0,200,4)
  Q_final2=matrix(0,11,1)
  
  for(r in 4:14)
  {
    bb = INDEX(r)
    
    # build grid
    iris.grid = function(r) somgrid(xdim =  r, ydim=r, topo="hexagonal")
    
    # build model
    iris.som = som(iris.sc, grid=iris.grid(r), rlen=100, alpha=c(0.25,0.001))
    
    #pd 들의 max를 구하기 위해 비교하는 함수 likelihood 반환 x, k, r, l : vector
    pd=function(x,k,l,r){exp(-0.5*euc.dist(x,k))+exp(-0.5*euc.dist(x,l))+exp(-0.5*euc.dist(x,r))}
    
    # function "euc.dist" : 벡터 x1과 x2 사이의 유클리디안제곱거리
    euc.dist <- function(x1, x2) (sum((x1 - x2) ^ 2)) # 
    
    
    # function "fincodes" : 원래 생성된 m개의 중량벡터 뒤에 6개의가상벡터연결
    fincodes=function(r,i) {rbind(as.matrix(iris.som$codes[[1]]),matrix(c(100,100,100,100),1,4))}
    
    # r노드일때, i 번째 객체의 j번째 방향애 대한 k(node index) 함수 k=1(승자노드), k=2(왼쪽노드),k=3(오른쪽노드) 중량벡터
    krl=function(i,j,k,r)fincodes(r,i)[ifelse(bb[[iris.som$unit.classif[i] ]]$directions[j,k]>0,
                                              bb[[iris.som$unit.classif[i] ]]$directions[j,k],
                                              nrow(fincodes(r,i))),]
    
    #각 개체에 따라 6개의 방향에 대하여 최대인 방향찾고 그 방향 나타내는 index 저장
    
    direction_vector = rep(0, 150)
    
    
    for(i in 1:200)
    {
      vec = function(i,r){c(pd(iris.sc[i,],krl(i,1,1,r),krl(i,1,2,r),krl(i,1,3,r)), 
                            pd(iris.sc[i,],krl(i,2,1,r),krl(i,2,2,r),krl(i,2,3,r)),
                            pd(iris.sc[i,],krl(i,3,1,r),krl(i,3,2,r),krl(i,3,3,r)),
                            pd(iris.sc[i,],krl(i,4,1,r),krl(i,4,2,r),krl(i,4,3,r)),
                            pd(iris.sc[i,],krl(i,5,1,r),krl(i,5,2,r),krl(i,5,3,r)),
                            pd(iris.sc[i,],krl(i,6,1,r),krl(i,6,2,r),krl(i,6,3,r)))}
      maxi_index = which.max(vec(i,r))
      direction_vector[i] = maxi_index
      maxi = max(vec(i,r))
      cat("max direction = ", maxi_index, "value = ", maxi, "\n")
    } 
    
    #finwight = 각 개체의 최종 개체가 표현될 승자노드,왼쪽노드,오른쪽 노드 중량벡터의 list
    
    finweight=list()
    
    EACHVECTOR=function(r,i){
      weight_matrix=matrix(0,3,4)
      weight_matrix[1,]=fincodes(r,i)[ifelse(bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1] >0,
                                             bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1] ,
                                             nrow(fincodes(r,i))+bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],1]+1),] 
      
      weight_matrix[2,]=fincodes(r,i)[ifelse(bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2] >0,
                                             bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2] ,
                                             nrow(fincodes(r,i))+bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],2]+1),] 
      
      weight_matrix[3,]=fincodes(r,i)[ifelse(bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3] >0,
                                             bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3] ,
                                             nrow(fincodes(r,i))+bb[[iris.som$unit.classif[i] ]]$directions[direction_vector[i],3]+1),] 
      return(weight_matrix)
    }
    
    for(i in 1:200){
      finweight[[i]]= EACHVECTOR(r,i) #r 을 매번 바꿔주기 
    }
    
    ##  가능도함수 생성 및 p_r, p_l 함수 만들기
    lik=function(i,w,b){
      exp((-0.5/b)*t(i-w)%*%(i-w))
    }
    # function "p_r, p_l" : cxc SOM에서 
    # i 번째 데이터에 대한 b=beta 가능도
    p_r=function(i,b) {
      lik(iris.sc[i,],finweight[[i]][1,],b)/
        (lik(iris.sc[i,],finweight[[i]][3,],b)+lik(iris.sc[i,],finweight[[i]][1,],b)) }
    
    p_l=function(i,b) {
      lik(iris.sc[i,],finweight[[i]][1,],b)/
        (lik(iris.sc[i,],finweight[[i]][2,],b)+lik(iris.sc[i,],finweight[[i]][1,],b)) }
    
    
    ##########coord 함수사용해서 beta구하기 여기서 a=p_r, b=p_l , k=승자노드, r=오른쪽노드, l=왼쪽노드
    
    coord=function(i,b){
      (((2*p_r(i,b)+2*p_l(i,b)-1)*finweight[[i]][1,])/(4*sqrt(3))+
         ((1-2*p_r(i,b)+p_l(i,b))*finweight[[i]][3,])/(2*sqrt(3))+
         ((1+p_r(i,b)-2*p_l(i,b))*finweight[[i]][2,])/(2*sqrt(3)))/(sqrt(3)/4)
    }
    
    
    
    #Q_beta 함수
    
    Q_beta=function(b) {
      for(i in 1:200) {
        Q_2[i,] = (iris.sc[i,]-coord(i,b))
      } 
      return(sum(Q_2^2))
    }
    Q_final2[r-3,1]= optimize(Q_beta,c(0,1000), tol=0.0001)$objective
  }
  Q_final2=cbind(Q_final2,matrix(c(4,5,6,7,8,9,10,11,12,13,14),11,1))
  Q_100_4[[j]]=Q_final2
}



#final value of overall 100 simulation mean result may be wrote as below.
Y4 <- do.call(cbind, Q_100_4)
Y4 <- array(Y4, dim=c(dim(Q_100_4[[1]]), length(Q_100_4)))
Q_4=matrix(0,11,2)
Q_4=apply(Y4, c(1, 2), mean)

plot(Q_4[,2],Q_4[,1], xlab="r", ylab="Q_beta",frame = FALSE,axes=FALSE, col = "blue",type="l")
box();axis(1);axis(2)
title("Simulation_4 : Q_hexagonal of gridsize ( r ) = 4 : 14")

