library(MASS)
library(kohonen)

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
    index_list[[j]]$adj_index = EVEN(j, r) 
    index_list[[j]]$directions = DIRECTION_TRIANGLE(j, EVEN(j, r)) 
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


# =================================================================
# End of the code
# =================================================================

INDEX(4)
r
library(MASS)
data(USAirlines)

data(crabs)
crabs$group=as.matrix(c(rep(1,50),rep(2,50),rep(3,50),rep(4,50)))
crabs=crabs[,-c(1:3)]
pairs(crabs[,1:5])
library(Tsphere)
TransSphere(crabs[,1:5],crabs[,6])


#data("USAirlines")
#iris.sc=as.matrix(USAirlines[,3:6])
#iris.sc = scale(iris.sc)
INDEX(4)
r=5
plot(iris.som$unit.classif)
plot(iris.som$changes)
plot(iris.som,type="mapping")

set.seed(01039114051)
Q_2=matrix(0,150,4)
Q_final2=matrix(0,9,1)

for(r in 4:12)
{
  bb = INDEX(r)
  
  # install the kohonen package
  install.packages("kohonen")
  install.packages("SOMbrero")
  install.packages("distances")
  
  # load the kohonen package
  library("kohonen")
  library("SOMbrero")
  library("distances")
  
  # scale data
  #iris.sc=as.matrix(iris[,1:4])
  #iris.sc = scale(iris[, 1:4])
  
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
  
  
  #하고싶은건, 각 개체에 따라 6개의 방향에 대하여 최대인 방향찾고 그 방향 나타내는 index 저장
  
  direction_vector = rep(0, 150)
  
  
  for(i in 1:150)
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
  
  for(i in 1:150){
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
  
  
  coord=function(i,b){
    (((2*p_r(i,b)+2*p_l(i,b)-1)*finweight[[i]][1,])/(4*sqrt(3))+
       ((1-2*p_r(i,b)+p_l(i,b))*finweight[[i]][3,])/(2*sqrt(3))+
       ((1+p_r(i,b)-2*p_l(i,b))*finweight[[i]][2,])/(2*sqrt(3)))/(sqrt(3)/4)
  }
  
  #Q_beta 함수
  
  Q_beta=function(b) {
    for(i in 1:150) {
      Q_2[i,] = (iris.sc[i,]-coord(i,b))
    } 
    return(sum(Q_2^2))
  }
  Q_final2[r-3,1]= optimize(Q_beta,c(0,1000), tol=0.0001)$objective
}

Q_final2=cbind(Q_final2,matrix(c(4,5,6,7,8,9,10,11,12),9,1))
plot(Q_final2[,2],Q_final2[,1], xlab="r", ylab="Q_beta",frame = FALSE,axes=FALSE, col = "blue",type="l")
box();axis(1);axis(2)
title("Q_hexagoanl of gridsize ( r ) = 4 : 12")



