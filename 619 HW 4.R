# 14.3
y=ts(matrix(0,6,221))
y[,1]=t(c(0,0,0,5,10,20))
for(i in 1:6){
  e=rnorm(220)
  for(j in 2:221){
    y[i,j]=y[i,j-1]+e[j-1]
  }
}
y=ts(t(y[,21:221]))
plot(y,plot.type='s',col=c(1,1,1,2,2,2),main='Random Walks with Different Initial')
y=ts(matrix(0,6,221))
y[,1]=t(c(0,0,0,5,10,20))
gamma=c(0,0,0,0.1,0.1,0.1)
for(i in 1:6){
  e=rnorm(220)
  for(j in 2:221){
    y[i,j]=gamma[i]+y[i,j-1]+e[j-1]
  }
}
y=ts(t(y[,21:221]))
plot(y,plot.type='s',col=c(1,1,1,2,2,2),main='Random Walks with/without drift')
#14.6
tic
rho=c(-0.7,0,0.7,1)
N=1000
r_ratio=matrix(0,16,5)
S=c(50,100,200,400,800)
for(s in 1:5){
  c=0
  for(i in 1:4){
    for(j in 1:4){
      c=c+1
      r=0
      for(n in 1:S[s]){
        e1=rnorm(N)
        e2=rnorm(N)
        y=ts(matrix(0,N,1))
        x=ts(matrix(0,N,1))
        y[1]=e1[1]
        x[1]=e2[1]
        for(t in 2:N){
          y[t]=rho[i]*y[t-1]+e1[t]
          x[t]=rho[j]*x[t-1]+e2[t]
        }
        model=lm(y~x)
        se=sigma(model)/sqrt(sum((x-mean(x))^2))
        if (abs(model$coef[2]/se)>1.96){
          r=r+1
        }
      }
      r_ratio[c,s]=r/S[s]
    }
  }
}
toc
#14.7
rho=c(-0.7,0,0.7,1)
N=10000
r_ratio=matrix(0,16,5)
S=c(50,100,200,400,800)
for(s in 1:5){
  c=0
  for(i in 1:4){
    for(j in 1:4){
      c=c+1
      r=0
      for(n in 1:S[s]){
        e1=rnorm(N)
        e2=rnorm(N)
        y=ts(matrix(0,N,1))
        x=ts(matrix(0,N,1))
        y[1]=e1[1]
        x[1]=e2[1]
        for(t in 2:N){
          y[t]=rho[i]*y[t-1]+e1[t]
          x[t]=rho[j]*x[t-1]+e2[t]
        }
        model=dynlm(y~x+lag(y,-1))
        if (abs(model$coef[2]/sigma(model))>1.96){
          r=r+1
        }
      }
      r_ratio[c,s]=r/S[s]
    }
  }
}
#14.8
rho=0.8
S=c(50,100,200,400,800,1600,3200,6400)
p=round(3*S^0.25)
N=10000
r_ratio=matrix(0,8,2)
for(s in 1:8){
      r=0
        e1=rnorm(N)
        e2=rnorm(N)
        y=ts(matrix(0,N,1))
        x=ts(matrix(0,N,1))
        y[1]=e1[1]
        x[1]=e2[1]
        for(t in 2:N){
          y[t]=rho[i]*y[t-1]+e1[t]
          x[t]=rho[j]*x[t-1]+e2[t]
        }
        model=dynlm(y~x+lag(y,-1))
        u_hat1=unname(model$residuals)
        W1=unname(window(cbind(ita,x,lag(y,-1)),start=2,end=N))
        Ga_hat=array(0,c(3,3,8))
        for(j in 1:8){
          u_hat=u_hat1
          for(i in (j+1):N){
            Ga_hat[,,j]=Ga_hat[,,j]+u_hat[i]*u_hat[i-j]*(t(t(W[i,]))%*%t(W[i-j,]))
          }
          Ga_hat[,,j]=Ga_hat[,,j]/N
        }
        Ga0=matrix(0,3,3)
        for(i in 1:176){
          Ga0=Ga0+u_hat[i]^2*(t(t(W[i,]))%*%t(W[i,]))
        }
        Ga0=Ga0/176
        Sig_NW1=array(0,c(3,3,8))
        Psi=matrix(0,3,3)
        Omega=matrix(0,8,3)
        for(p in 1:8){
          for(j in 1:p){
            Sig_NW1[,,p]=Sig_NW1[,,p]+(1-(j/(p+1)))*(Ga_hat[,,j]+t(Ga_hat[,,j]))
          }
          Sig_NW1[,,p]=Sig_NW1[,,p]+Ga0
          Psi=solve(t(W)%*%W/176)%*%(Sig_NW1[,,p]/176)%*%solve(t(W)%*%W/176)
          Omega[p,]=sqrt(diag(Psi))
        }
        if (abs(model$coef[2]/sigma(model))>1.96){
          r=r+1
        }
      }
      r_ratio[c,s]=r/S[s]
    }
  }
}