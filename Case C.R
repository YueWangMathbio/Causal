# This is the code of simulation in "Causal inference in degenerate systems: An impossibility result".
# Four algorithms: Alg.2-AF, Alg.2-KI, Alg.S.1, Alg.S.2 are implemented in four test cases. 
# In Alg.2-AF and Alg.S.1, Alg.1 is used.
# Each algorithm runs for different numbers of observations: 
# n=300, 500, 800, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 10000, 12000, 14000, 17000, 20000, 23000, 26000, 30000.
# The output "ou" is the number of runs that claim there are multiple Markov boundaries.
# In output "ou", the order is Alg.S.2, Alg.2-KI, Alg.2-AF, Alg.S.1.




nn=500;#each algorithm runs 500 times in each case
library("infotheo")#used in calculating CMI
alpha=0.999;#significance level=0.001



####################case C#######################


for(cou in 1:19){
  ou=c(0,0,0,0);#number of runs that produce correct result
  if(cou==1){n=300;}
  if(cou==2){n=500;}
  if(cou==3){n=800;}
  if(cou==4){n=1000;}  
  if(cou==5){n=1500;}
  if(cou==6){n=2000;}
  if(cou==7){n=3000;}
  if(cou==8){n=4000;}
  if(cou==9){n=5000;}
  if(cou==10){n=6000;}
  if(cou==11){n=7000;}
  if(cou==12){n=10000;}
  if(cou==13){n=12000;}
  if(cou==14){n=14000;}
  if(cou==15){n=17000;}
  if(cou==16){n=20000;}
  if(cou==17){n=23000;}
  if(cou==18){n=26000;}
  if(cou==19){n=30000;}
  
  
  
  a1=0;
  a2=0;
  a3=0;
  a4=0;
  
  print(n);

  for(qq in 1:nn){
    
    ################start producing test data##########
    x=matrix(0,nrow=n,ncol=10);
    for(i in 1:10){
      x[,i]=sample(0:1,n,replace=T);
    }
    y=x[,1];
    for(i in 1:n){
      if(x[i,1]==x[i,2]){
        y[i]=0;
      }else{
        y[i]=1;
      }
    }
    x[,10]=y;
    te=runif(n);
    for(i in 1:n){
      if(te[i]<0.05){
        x[i,10]=1-x[i,10];
      }
    }
    x[,9]=x[,10];
    for(i in 1:n){
      k=ceiling(runif(1)*10);
      if(k>2 & k<5){
        y[i]=x[i,k];
      }
    }
    #########finish producing test data#######
    
    ###############start of Alg.S.2 ####################
    
    ch=qchisq(alpha, df=2^(10-1));
    z=x;
    c=runif(10);
    f=0;
    nd1=rep(1,10);
    for(i in 10:1){
      c[i]=mutinformation(y,x[,])-mutinformation(y,x[,-i]);
      c[i]=c[i]*n*2;
      if(c[i]<ch){
        
        if(is.null(ncol(z))==FALSE){
          z=z[,-i];
        }else{z=matrix(0,nrow=n,ncol=1);}
        nd1[i]=0;
        f=f+1;
      }
    }
    
    c2=mutinformation(y,x)-mutinformation(y,z);
    c2=c2*n*2;
    d=(2^f-1)*2^(10-f);
    ch2=qchisq(alpha, df=d);
    
    if(c2>ch2){
      a1=a1+1;
    }else{
      
    }
    #############finish of Alg.S.2 #############
    
    ###############start of Alg.S.1 #################
    
    f2=10;
    z=x;
    s=runif(10);
    p=rep(0,10);
    nd=rep(1,10);
    while(f2>0){
      min=2;
      for(i in 1:10)
      {
        if(nd[i]>0){
          s[i]=mutinformation(y,z)-mutinformation(y,z[,-i]);
          s[i]=s[i]*n*2;
          d=2^(f2-1);
          p[i]=pchisq(s[i],df=d);
          if(p[i]<min)
          {
            min=p[i];
            mn=i;
          }
        }
      }
      if(min>alpha)
      {
        break;
      }else
      {
        z[,mn]=rep(0,n);
        nd[mn]=0;
        f2=f2-1;
      }
    }
    ##########finish Alg.1 with an MB#######
    m=0;
    for(i in 10:1){
      if(nd[i]==1){
        s2=mutinformation(y,x)-mutinformation(y,x[,-i]);
        s2=s2*n*2;
        ch4=qchisq(alpha, df=2^(10-1));
        if(s2<=ch4){
          m=1;
          
          break;
        }
      }
    }
    if(m==0){
      
    }else{
      
      a2=a2+1;
    }
    #################finish of Alg.S.1 #################
    
    ##############start of Alg.2-KI #############
    k0=0.8;
    #KIAMB
    #forward
    m=matrix(0,nrow = n,ncol = 10);
    m2=matrix(0,nrow = n,ncol = 11);
    nd2=rep(0,10);#0 means not M, 1 means M
    pva=rep(0,10);#temp p-value
    b=0;#real col n of M
    while(prod(nd2)==0){#E is not empty
      
      
      r=0;#temp index
      mark=rep(0,10);# 1 means in CanMB, 0 means not
      for(i in 1:10){
        if(nd2[i]==0){
          m2[,11]=x[,i];
          pva[i]=mutinformation(y,m2)-mutinformation(y,m);
          d=qchisq(alpha, df=2^b);
          if(pva[i]*n*2>d){
            mark[i]=1;
          }
        }
      }# finish constructing CanMB
      if(sum(mark)==0){
        break;
      }
      no2=1;
      if(floor(sum(mark)*k0)>no2){
        no2=floor(sum(mark)*k0);
      }
      while(sum(mark)>no2){
        tem8=ceiling(runif(1)*sum(mark));
        for(i in 1:10){
          if(sum(mark[1:i])==tem8){
            mark[i]=0;
            break;
          }
        }
      }
      ma=-1;
      for(i in 1:10){
        if(mark[i]==1){
          if(ma<pva[i]){
            ma=pva[i];
            r=i;
          }
        }
      }
      nd2[r]=1;
      b=b+1;
      m[,r]=x[,r];
      m2[,r]=x[,r];
    }
    #backward
    
    for(i in 1:10){
      if(nd2[i]==1){
        m3=m;
        m3[,i]=rep(0,n);
        te=mutinformation(y,m)-mutinformation(y,m3);
        te=te*n*2;
        d=qchisq(alpha, df=2^(b-1));
        if(te<d){
          b=b-1;
          nd2[i]=0;
          m[,i]=rep(0,n);
        }
      }
    }
    
    #####finish of KIAMB######
    la=0;
    for(l in 1:10){
      if(nd2[l]==1){#start another IAMB
        m4=matrix(0,nrow = n,ncol = 10);
        m5=matrix(0,nrow = n,ncol = 11);
        nd3=rep(0,10);#0 means E, 1 means M
        pva=rep(0,10);#temp p-value
        b=0;#real col n of M
        w=x;
        w[,l]=rep(0,n);
        while(prod(nd3)==0){#E is not empty
          r=0;#temp index
          mark=rep(0,10);# 1 means in CanMB, 0 means not
          for(i in 1:10){
            if(nd3[i]==0){
              m5[,11]=w[,i];
              pva[i]=mutinformation(y,m5)-mutinformation(y,m4);
              d=qchisq(alpha, df=2^b);
              if(pva[i]*n*2>d){
                mark[i]=1;
              }
            }
          }# finish constructing CanMB
          if(sum(mark)==0){
            break;
          }
          no2=1;
          if(floor(sum(mark)*k0)>no2){
            no2=floor(sum(mark)*k0);
          }
          while(sum(mark)>no2){
            tem8=ceiling(runif(1)*sum(mark));
            for(i in 1:10){
              if(sum(mark[1:i])==tem8){
                mark[i]=0;
                break;
              }
            }
          }
          ma=-1;
          for(i in 1:10){
            if(mark[i]==1){
              if(ma<pva[i]){
                ma=pva[i];
                r=i;
              }
            }
          }
          nd3[r]=1;
          b=b+1;
          m4[,r]=w[,r];
          m5[,r]=w[,r];
        }
        #backward
        
        for(i in 1:10){
          if(nd3[i]==1){
            m6=m4;
            m6[,i]=rep(0,n);
            te=mutinformation(y,m4)-mutinformation(y,m6);
            te=te*n*2;
            d=qchisq(alpha, df=2^(b-1));
            if(te<d){
              b=b-1;
              nd3[i]=0;
              m4[,i]=rep(0,n);
            }
          }
        }
        
        m7=matrix(0,nrow=n,ncol=10);
        p1=0;
        
        p2=0;
        for(co in 1:10){
          if(nd2[co]==1||nd3[co]==1){
            p1=p1+1;
            m7[,co]=x[,co];
          }
          if(nd3[co]==1){
            p2=p2+1;
          }
        }
        
        
        te=mutinformation(y,m7)-mutinformation(y,m4);
        te=te*n*2;
        f=(2^(p1-p2)-1)*2^p2;
        d=qchisq(alpha, df=f);
        
        if(te<d){
          
          la=1;
          break
        }
      }#end of another IAMB
      
    }
    if(la==0){
      
    }else{
      
      a3=a3+1;
    }
    
    ##############finish of Alg.2-KI ###############
    
    ############start of Alg.2-AF #########
    f2=10;
    z=x;
    s=runif(10);
    p=rep(0,10);
    nd=rep(1,10);
    while(f2>0){
      min=2;
      for(i in 1:10)
      {
        if(nd[i]>0){
          s[i]=mutinformation(y,z)-mutinformation(y,z[,-i]);
          s[i]=s[i]*n*2;
          d=2^(f2-1);
          p[i]=pchisq(s[i],df=d);
          if(p[i]<min)
          {
            min=p[i];
            mn=i;
          }
        }
      }
      if(min>alpha)
      {
        break;
      }else
      {
        z[,mn]=rep(0,n);
        nd[mn]=0;
        f2=f2-1;
      }
    }
    
    #####finish Alg.1 with an MB######
    
    m=0;
    for(i in 1:10){
      ndt=rep(1,10);
      zn=x;
      if(nd[i]==1){
        ndt[i]=0;
        zn[,i]=rep(0,n);
        f2=9;
        
        s=runif(10);
        p=rep(0,10);
        
        while(f2>0){
          min=2;
          for(i in 1:10)
          {
            if(ndt[i]>0){
              s[i]=mutinformation(y,zn)-mutinformation(y,zn[,-i]);
              s[i]=s[i]*n*2;
              d=2^(f2-1);
              p[i]=pchisq(s[i],df=d);
              if(p[i]<min)
              {
                min=p[i];
                mn=i;
              }
            }
          }
          if(min>alpha)
          {
            break;
          }else
          {
            zn[,mn]=rep(0,n);
            ndt[mn]=0;
            f2=f2-1;
          }
        }
        
        zq=matrix(0,nrow=n,ncol=10);
        for(k in 1:10){
          if(nd[k]==1 | ndt[k]==1){
            zq[,k]=x[,k];
          }
        }
        s2=mutinformation(y,zq)-mutinformation(y,zn);
        s2=s2*n*2;
        nda=nd*ndt;
        q1=sum(nd)-sum(nda);
        q2=sum(ndt);
        d=(2^q1-1)*2^q2;
        ch5=qchisq(alpha, df=d);
        
        if(s2<=ch5){
          m=1;
          break;
        }
      }
    }
    if(m==0){
      
    }else{
      
      a4=a4+1;
    }
    ou=c(a1,a3,a4,a2);
  }
  ###########finish of Alg.2-AF ##########
  
  if(cou==1){write.csv(ou, file = "C300.csv")}
  if(cou==2){write.csv(ou, file = "C500.csv")}
  if(cou==3){write.csv(ou, file = "C800.csv")}
  if(cou==4){write.csv(ou, file = "C1000.csv")}  
  if(cou==5){write.csv(ou, file = "C1500.csv")}
  if(cou==6){write.csv(ou, file = "C2000.csv")}
  if(cou==7){write.csv(ou, file = "C3000.csv")}
  if(cou==8){write.csv(ou, file = "C4000.csv")}
  if(cou==9){write.csv(ou, file = "C5000.csv")}
  if(cou==10){write.csv(ou, file = "C6000.csv")}
  if(cou==11){write.csv(ou, file = "C7000.csv")}
  if(cou==12){write.csv(ou, file = "C10000.csv")}
  if(cou==13){write.csv(ou, file = "C12000.csv")}
  if(cou==14){write.csv(ou, file = "C14000.csv")}
  if(cou==15){write.csv(ou, file = "C17000.csv")}
  if(cou==16){write.csv(ou, file = "C20000.csv")}
  if(cou==17){write.csv(ou, file = "C23000.csv")}
  if(cou==18){write.csv(ou, file = "C26000.csv")}
  if(cou==19){write.csv(ou, file = "C30000.csv")}
  print(ou)
  
}