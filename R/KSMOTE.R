#' @title A dataset used for illustration.
#' @name yeast_tra
#' @description We use this imbalanced dataset to carry out our interpolation algorithm.
NULL

#' @title A dataset used for illustration.
#' @name yeast_test
#' @description We use this imbalanced dataset to carry out our interpolation algorithm.
NULL

#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rbinom rbeta
#' @useDynLib StatComp21092
NULL

#' This is some descriptio of this function.
#' @title kernel function
#' @param x is a vector in feature space
#' @param y is a vector
#' @return kernel inner product
#' @export
#' @examples x=c(1,2);y=c(1,2);print(kernel_poly(x,y))


kernel_poly<-function(x,y)
{
  return((x%*%y+1)^3)
}

#' This is a function define distance in kernel space.
#' @title distance in kernel space
#' @param x is a vertor in feature space
#' @param y is a vertor in feature space
#' @return distance in kernel space
#' @export

distance<-function(x,y)
{
  d_square=kernel_poly(x,x)+kernel_poly(y,y)-2*kernel_poly(x,y)
  return(d_square)
}


#' This a function to generate kernel matrix on train_set and test_set.
#' @title K_SMOTE
#' @param train_set is the train set.
#' @param test_set is the test set.
#' @param P is the number of generating samples.
#' @param k is the k of kNN algorithm.
#' @import stats
#' @return output is a list of two matrix train_K_mat and test_K_mat.
#' @export

K_SMOTE<-function(train_set,test_set,P,k)
{
  #train_set,test_set分别为训练集和测试集，类型为dataframe
  #P为需要生成的少数类样本数
  #k为kNN中近邻的个数


  #记录训练集中少数类和多数类的index
  min_index=(which(train_set['Class']==1))
  max_index=(which(train_set['Class']==0))


  #初始化核矩阵
  N = nrow(train_set)
  K = matrix(0,N+P,N+P)
  D = matrix(0,N+P,N+P)

  #初始化生成对与随机数
  p_index=numeric(length = P)
  q_index=numeric(length = P)
  delta=runif(P,0,1)

  #计算原核矩阵
  for (i in 1:N)
    for (j in 1:N)
      K[i,j]=kernel_poly(data.matrix(train_set)[i,3:ncol(train_set)-1],data.matrix(train_set)[j,3:ncol(train_set)-1])

  #计算距离矩阵(我们这里只关心少数类)
  for (i in min_index)
    for (j in min_index)
      D[i,j]=distance(data.matrix(train_set)[i,3:ncol(train_set)-1],data.matrix(train_set)[j,3:ncol(train_set)-1])



  #用少数类的距离矩阵结合KNN得到生成对
  for (i in 1:P)
  {
    p=sample(min_index,1,replace = TRUE)
    r=min_index[rank(D[p,min_index])]
    q=sample(r[1:k],1,replace = TRUE)
    p_index[i]=p
    q_index[i]=q
  }



  #生成，即填补核矩阵

  #填K2和K2转置：一个为原样本，一个为生成样本
  for(i in 1:P)
    for(j in 1:N)
    {
      p=p_index[i]
      q=q_index[i]
      K[i+N,j]=(1-delta[i])*K[p,j]+delta[i]*K[q,j]
      K[j,i+N]=K[i+N,j]
    }

  #填K3：均为生成样本
  for(i in 1:P)
    for(j in 1:P)
    {
      l=p_index[i]
      m=q_index[i]
      p=p_index[j]
      q=q_index[j]
      K[i+N,j+N]=(1-delta[j])*(1-delta[i])*K[p,l]+(1-delta[j])*delta[i]*K[p,m]+delta[j]*(1-delta[i])*K[q,l]+delta[j]*delta[i]*K[q,m]
    }


  # 下面填测试样本与训练集的核矩阵
  n=nrow(test_set)
  gram_test=matrix(0,n,N+P)

  for(i in 1:n)
    for(j in 1:N)
      gram_test[i,j]=kernel_poly(data.matrix(test_set)[i,3:ncol(train_set)-1],data.matrix(train_set)[j,3:ncol(train_set)-1])

  for(i in 1:n)
    for(j in 1:P)
    {
      p=p_index[j]
      q=q_index[j]
      gram_test[i,j+N]=(1-delta[j])*gram_test[i,p]+delta[j]*gram_test[i,q]
    }

  mat=list("train_km"=K,"test_km"=gram_test)

  return(mat)
}

#' @title modify_label
#' @param train_set is the train set with target.
#' @param test_set is the test set with target.
#' @param P is the number of generating samples.
#' @return output is a list of modified target.
#' @export


modify_label=function(train_set, test_set, P)
{
  train_target=array(1, dim = nrow(train_set)+P)
  test_target=array(1, dim = nrow(test_set))

  for(i in 1:nrow(train_set))
    if(train_set[i,'Class']==0)
    {
      train_target[i]=-1
    }
  if(train_set[i,'Class']==1)
  {
    train_target[i]=1
  }

  for(i in 1:nrow(test_set))
    if(test_set[i,'Class']==0)
    {
      test_target[i]=-1
    }
  if(test_set[i,'Class']==1)
  {
    test_target[i]=1
  }
  target=list("train"=train_target,"test"=test_target)

  return(target)
}

#' @title update alpha
#' @param alpha_2 old alpha.
#' @param L the low.
#' @param H the high.
#' @return output is alpha_new.

update_alpha_2=function(alpha_2,L,H)
{
  if(alpha_2>H)
    return(H)
  else
    if(alpha_2<L)
      return(L)
  else
    return(alpha_2)
}

#' @title SMO_SVC
#' @param train_mat the kernel matrix of train set.
#' @param test_mat the kernel matrix of test set.
#' @param train_target train set target.
#' @param test_target test set target.
#' @param iter_max the maximum of iteration.
#' @param epsilon a next condition.
#' @param C the parameter of wrong.
#' @return the predict result.
#' @export

SMO_SVC=function(train_mat,test_mat,train_target,test_target,iter_max=10,epsilon=0.01,C=1)
{
  #初始化
  N=nrow(train_mat)
  alpha=numeric(length = nrow(train_mat))
  b=1
  k=0#迭代次数
  #定义函数g
  g=function(i)
  {
    return(sum(alpha*train_target*train_mat[i,])+b)
  }

  sample_without_i=function(i)
  {
    j=sample(1:N,1)
    while(j==i)
      j=sample(1:N,1)
    return(j)
  }

  while (k<iter_max)
  {
    change=0
    #取出i，j
    for(i in 1:N)
    {
      a_i=alpha[i]
      y_i=train_target[i]
      x_i=i
      E_i=g(x_i)-y_i

      j=sample_without_i(i)
      a_j=alpha[j]
      y_j=train_target[j]
      x_j=j
      E_j=g(x_j)-y_j

      eta=train_mat[i,i]+train_mat[j,j]-2*train_mat[i,j]
      if(eta<=0)
        next

      #未经剪辑的new_a_j
      new_a_j=a_j+y_j*(E_i-E_j)/eta

      #剪辑a_j_new
      if(y_i!=y_j)
      {
        L=max(0,a_j-a_i)
        H=min(C,C+a_j-a_i)
      }
      else
      {
        L=max(0,a_j+a_i-C)
        H=min(C,a_i+a_j)
      }

      a_j_new=update_alpha_2(new_a_j,L,H)
      a_i_new=a_i+y_i*y_j*(a_j-a_j_new)

      if(abs(a_j_new-a_j)<epsilon)
        next

      #更新alpha
      alpha[i]=a_i_new
      alpha[j]=a_j_new

      #更新b
      b_i=-E_i-y_i*train_mat[i,i]*(a_i_new - a_i) - y_j*train_mat[i,j]*(a_j_new - a_j) + b
      b_j=-E_j-y_i*train_mat[i,j]*(a_i_new - a_i) - y_j*train_mat[j,j]*(a_j_new - a_j) + b

      if(abs(a_i_new-C/2)<C/2)
        b=b_i
      else
        if(abs(a_j_new-C/2)<C/2)
          b=b_j
      else
        b=(b_i+b_j)/2

      change=change+1
    }
    if(change==0)
      k=k+1
    else
      k=0
  }

  #计算测试集结果

  result=numeric(length = nrow(test_mat))
  for (i in 1:nrow(test_mat))
  {
    f=sum(alpha*train_target*test_mat[i,])+b
    if(f>0)
      result[i]=1
    else
      result[i]=-1
  }
  return(result)
}

#' @title Confusion Matrix
#' @param pred the prediction of model.
#' @param target the true label.
#' @return the confusion matrix.
#' @export


ConfusionMatrix=function(pred,target)
{
  N=length(pred)
  TP=FP=FN=TN=0
  for(i in 1:N)
  {
    if(pred[i]==1 &&target[i]==1)
      TP=TP+1
    if(pred[i]==-1&&target[i]==-1)
      TN=TN+1
    if(pred[i]==1&&target[i]==-1)
      FP=FP+1
    if(pred[i]==-1&&target[i]==1)
      FN=FN+1
  }

  mat=matrix(c(TP,FP,FN,TN),ncol=2)
  return(mat)
}
