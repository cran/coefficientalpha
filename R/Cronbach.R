##############################################
## Function to find missing data patterns   ##
## Function rsem.pattern                    ##
##############################################
rsem.pattern<-function(x,print=FALSE){
  ## This function generate missing data patterns
  ## INPUT
  ## x: data set
  ## OUTPUT
  ## misinfo
  ##   [1,] number of cases in the pattern
  ##   [2,] number of observed variables in the pattern
  ##   [3:(p+2), ] contains a permutation of {1,2,3,...,p} 
  ##                      the first part corresponding to observed variables
  ##                      the remaining corresponding to missing variables
  if (missing(x)) stop("A data set has to be provided!")
  if (!is.matrix(x)) x<-as.matrix(x)
  y<-x
  ## check to see if all data are missing
  M<-is.na(x)
  nM<-max(apply(M,1,sum))
  n<-dim(x)[1]
  p<-dim(x)[2]
  if (nM==p) stop("Some cases have missing data on all variables. Please delete them first.")
  misorder<-rep(0,n)
  for (i in 1:n){
    misorderj<-0
    for (j in 1:p){
      if (is.na(x[i,j])) misorderj<-misorderj+2^(j-1)
    }
    misorder[i]<-misorderj
  }
  ## Combine data with missing pattern indicator
  ## order data according to misorder
  #id<-1:n
  temp<-order(misorder)
  x<-x[temp,]
  misn<-misorder[temp]
  #id<-id[temp]
  
  ##identifying the subscripts of missing variables and put them in misinfo;
  mi<-0; nmi<-0;oi<-0; noi<-0;
  for (j in 1:p){
    if (is.na(x[1,j])){
      mi<-c(mi,j)  ## recording the missing variable subscript in the first case
      nmi<-nmi+1   ## number of missing values in the first case
    }else{
      oi<-c(oi,j)
      noi<-noi+1
    }
  }
  oi<-oi[2:(noi+1)]
  if (nmi==0){
    misinfo_0 = c(noi, oi)
  }else{
    mi<-mi[2:(nmi+1)]
    misinfo_0<-c(noi,oi,mi) ##recording the # observed variables, locations;
  }	
  patnobs <- 0 ## number of observed cases in a pattern
  totpat<-1; ncount<-1;
  t1<-misn[1]
  for (i in 2:n){
    if (misn[i]==t1){
      ncount<-ncount+1
    }else{
      patnobs<-c(patnobs,ncount)
      t1<-misn[i]
      ncount<-1
      totpat<-totpat+1
      mi<-0; nmi<-0;oi<-0; noi<-0;
      for (j in 1:p){
        if (is.na(x[i,j])){
          mi<-c(mi,j)
          nmi<-nmi+1
        }else{
          oi<-c(oi,j)
          noi<-noi+1
        }
      }
      oi<-oi[2:(noi+1)]
      mi<-mi[2:(nmi+1)]
      misinfo_0 <- rbind(misinfo_0, c(noi,oi,mi))
    }
  }
  patnobs<-c(patnobs, ncount)
  patnobs<-patnobs[2:(totpat+1)]
  if (is.vector(misinfo_0)){
    misinfo<-c(patnobs, misinfo_0)
  }else{
    misinfo<-cbind(patnobs, misinfo_0)
  }
  if (!is.matrix(misinfo)){
    misinfo<-matrix(misinfo, nrow=1)
  }
  
  ## Different presentation of missing data patterns
  nr<-nrow(misinfo)
  mispat<-matrix(1, nrow=nr, ncol=p)
  
  for (i in 1:nr){
    if (misinfo[i,2]<p){
      ind<-misinfo[i, (misinfo[i,2]+3):(p+2)]
      mispat[i, ind]<-0
    }
  }
  mispat<-cbind(misinfo[,1:2, drop=FALSE], mispat)
  rownames(mispat)<-paste('Pattern ', 1:nr, sep="")
  colnames(mispat)<-c('n','nvar',colnames(x))
  
  if (print) print(mispat)
  
  invisible(list(misinfo=misinfo, mispat=mispat, x=x, y=y))
}


##############################################
## Weight function for each case            ##
## Function rsem.weight                     ##
##############################################

rsem.weight<-function(x, varphi, mu0, sig0){
	##x<-xpattern$x
	if (!is.matrix(x)) x<-as.matrix(x)
	n<-dim(x)[1]
	p<-dim(x)[2]
	prob<-1-varphi ## chi-square p-value
	wi1all<-wi2all<-NULL
	if (varphi==0){
		wi1all<-wi2all<-rep(1, n)
	}else{
		for (i in 1:n){
			## take out one row of x
			xi<-x[i, ]
			xid<-which(!is.na(xi))
			xic<-xi[xid]
			
			ximu<-mu0[xid]
			xisig<-as.matrix(sig0[xid, xid])
			xidiff<-as.matrix(xic-ximu)
			
			di2<-t(xidiff)%*%solve(xisig)%*%xidiff
			di<-sqrt(di2)
			
			pi<-length(xic)
    		chip<-qchisq(prob, pi)
    		ck<-sqrt(chip)
    		cbeta<-( pi*pchisq(chip,pi+2) + chip*(1-prob) )/pi
    		
			if (di <= ck){
				wi1all<-c(wi1all, 1)
				wi2all<-c(wi2all, 1/cbeta)
			}else{
				wi1<-ck/di
				wi1all<-c(wi1all, wi1)
				wi2all<-c(wi2all, wi1*wi1/cbeta)
			}
		}		
    }       
    return(list(w1=wi1all, w2=wi2all))
}

##############################################
## EM algorithm for unstructured (mu, sigma)##
## Function rsem.emmusig                    ##
##############################################

rsem.emmusig<-function(xpattern, varphi=.1, max.it=1000, st='i'){
  ## x: data set
  ## misinfo: missing data pattern
  ## varphi: 
  if (is.null(xpattern$mispat)) stop("The output from the function rsem.pattern is required")
  x<-xpattern$x
  misinfo<-xpattern$misinfo
  
  ep <- 1e-6  ## precision
  n<-dim(x)[1]
  p<-dim(x)[2]
  
  mu0<-rep(0,p)
  sig0<-diag(p)
  if (st=='mcd'){
    y<-na.omit(x)
    ny<-nrow(y)
    par.st<-cov.rob(y, method='mcd')
    mu0<-par.st$center
    sig0<-par.st$cov
  }
  
  n_it<-0;        
  dt<-1;
  if (varphi==0){
    ck<-10e+10
    cbeta<-1
  }else{
    prob<-1-varphi ## chi-square p-value
    chip<-qchisq(prob, p)
    ck<-sqrt(chip)
    cbeta<-( p*pchisq(chip, p+2) + chip*(1-prob) )/p
  }
  while (dt>ep && n_it <= max.it){
    sumx<-rep(0,p); sumxx<-array(0,dim=c(p,p)); sumw1<-0; sumw2<-0;
    npat<-dim(misinfo)[1]  ## number of missing data patterns
    p1<-misinfo[1,2]       ## number of observed variables in pattern 1
    n1<-misinfo[1,1]       ## number of cases in pattern 1
    if (p1==p){            ## complete data
      sigin <- solve(sig0)  ## matrix inverse
      for (i in 1:n1){
        xi<-x[i,]
        xi0<-xi-mu0
        di2<-xi0%*%sigin%*%xi0
        di<-sqrt(di2)
        ## Huber weight functions
        if (di<=ck){
          wi1<-1
          wi2<-1/cbeta
        }else{
          wi1<-ck/di
          wi2<-wi1*wi1/cbeta
        } ## end Huber weight
        sumw1<-sumw1+wi1;
        xxi0<-xi0%*%t(xi0)
        sumx<-sumx+wi1*xi
        sumxx<-sumxx+c(wi2)*xxi0
        sumw2<-sumw2+wi2
      } ## end for
    }else{ ## end p1==p
      ## with missing data
      if (varphi==0){
        ck1<-1e+10
        cbeta1<-1
      }else{ 
        chip1<-qchisq(prob, p1)
        ck1<-sqrt(chip1)
        cbeta1<-( p1*pchisq(chip1,p1+2) + chip1*(1-prob) )/p1
      }
      o1<-misinfo[1,3:(2+p1)]
      m1<-misinfo[1,(2+p1+1):(p+2)]
      mu_o<-mu0[o1]; mu_m<-mu0[m1]
      sig_oo<-sig0[o1,o1]; sig_om<-sig0[o1,m1];
      if (p1==1) {sig_mo<-sig_om}else{sig_mo<-t(sig_om)} 
      sig_mm<-sig0[m1,m1];
      sigin_oo<-solve(sig_oo)
      beta_mo<-sig_mo%*%sigin_oo
      
      delt <- array(0, dim=c(p,p))
      delt[m1,m1]<-sig_mm - beta_mo%*%sig_om
      for (i in 1:n1){
        xi<-x[i,]
        xi_o<-xi[o1]
        xi0_o<-xi_o-mu_o
        stdxi_o<-sigin_oo%*%xi0_o
        di2<-xi0_o%*%stdxi_o
        di<-sqrt(di2)
        if (di<=ck1){ ##Huber weight
          wi1<-1
          wi2<-1/cbeta1
        }else{
          wi1<-ck1/di
          wi2<-wi1*wi1/cbeta1
        }
        sumw1<-sumw1+wi1
        xm1<-mu_m+sig_mo%*%stdxi_o
        xi[m1]<-xm1
        xi0<-xi-mu0
        xxi0<-xi0%*%t(xi0)
        sumx<-sumx+wi1*xi
        sumxx<-sumxx+c(wi2)*xxi0+delt
        sumw2<-sumw2+wi2
      } ##end for 1:n1  
    }## end of (p1=p)
    ## start from pattern 2	
    if (npat>1){
      snj<-n1	
      for (j in 2:npat){
        nj<-misinfo[j,1]; pj<-misinfo[j,2];
        oj<-misinfo[j, 3:(2+pj)]; mj<-misinfo[j, (2+pj+1):(p+2)];
        mu_o<-mu0[oj]; mu_m<-mu0[mj];
        sig_oo<-sig0[oj,oj]; sig_om<-sig0[oj,mj];
        if (pj==1) {sig_mo<-sig_om}else{sig_mo<-t(sig_om)} 
        sig_mm<-sig0[mj,mj];
        sigin_oo<-solve(sig_oo)
        beta_mo<-sig_mo%*%sigin_oo
        delt <- array(0, dim=c(p,p))
        delt[mj,mj]<-sig_mm - beta_mo%*%sig_om
        if (varphi==0){
          ckj<-10e+10
          cbetaj<-1
        }else{
          chipj<-qchisq(prob,pj)
          ckj<-sqrt(chipj)
          cbetaj<- ( pj*pchisq(chipj, pj+2) + chipj*(1-prob) )/pj
        }
        for (i in ((snj+1):(snj+nj))){
          xi<-x[i,]
          xi_o<-xi[oj]
          xi0_o<-xi_o - mu_o
          stdxi_o<-sigin_oo%*%xi0_o
          di2<-xi0_o%*%stdxi_o
          di<-sqrt(di2)
          if (di<=ckj){ ##Huber weight
            wi1<-1
            wi2<-1/cbetaj
          }else{
            wi1<-ckj/di
            wi2<-wi1*wi1/cbetaj
          }
          sumw1<-sumw1+wi1
          xmj<-mu_m+sig_mo%*%stdxi_o
          xi[mj]<-xmj
          xi0<-xi-mu0
          xxi0<-xi0%*%t(xi0)
          sumx<-sumx+wi1*xi
          sumxx<-sumxx+c(wi2)*xxi0+delt
          sumw2<-sumw2+wi2
        }
        snj<-snj+nj
      } ## for (j in 2:npat)
    }
    mu1<-sumx/sumw1
    sig1<-sumxx/n
    dt<-max(c(max(abs(mu1-mu0)), max(abs(sig1-sig0))));
    mu0<-mu1;
    sig0<-sig1;
    n_it<-n_it+1;
  } ## end while
  if (n_it>=max.it) warning("The maximum number of iteration was exceeded. Please increase max.it in the input.")
  rownames(sig1)<-colnames(sig1)
  names(mu1)<-colnames(sig1)
  weight<-rsem.weight(xpattern$y, varphi, mu1, sig1)
  list(mu=mu1, sigma=sig1, max.it=n_it, weight=weight)
}


##############################################
## Stacking a matrix to a vector            ##
## Function rsem.vec                        ##
##############################################
rsem.vec<-function(x){
  t(t(as.vector(x)))
}

##############################################
## computing the estimator of the asymptotic## 
## covariance of \hat\Omega_{\hat\beta};    ##
## Function rsem.Ascov                      ##
##############################################  
rsem.gname<-function(name){
    temp.name<-NULL
    k<-length(name)
    for (i in 1:k){
      for (j in i:k){
        temp.name<-c(temp.name, paste(name[i], ".", name[j], sep=""))
      }
    }
    temp.name
}

###########
##############################################
## Generate a duplication matrix            ##
## Function rsem.DP                         ##
##############################################
rsem.DP <- function(x){ 	
  mat <- diag(x)
  index <- seq(x*(x+1)/2)
  mat[ lower.tri( mat , TRUE ) ] <- index
  mat[ upper.tri( mat ) ] <- t( mat )[ upper.tri( mat ) ]
  outer(c(mat), index , function( x , y ) ifelse(x==y, 1, 0 ) )
}

rsem.Ascov<-function(xpattern, musig, varphi=.1){
  if (is.null(xpattern$mispat)) stop("The output from the function rsem.pattern is required")
  x<-xpattern$x
  misinfo<-xpattern$misinfo
  
  mu0<-musig$mu
  sig0<-musig$sig
  
  n<-dim(x)[1]; p<-dim(x)[2];
  ps<-p*(p+1)/2; pps<-p+ps;
  dup<-rsem.DP(p) ##duplication matrix
  dupt<-t(dup) 
  i_p<-diag(p)
  B11<-array(0, c(p,p)); B12<-array(0,c(p,ps)); B22<-array(0,c(ps,ps));
  ddl11<-array(0,c(p,p)); ddl12<-array(0,c(p,ps)); 
  ddl21<-array(0,c(ps,p)); ddl22<-array(0,c(ps,ps));
  if (varphi==0){
    ck<-1e+10
    cbeta<-1
  }else{
    prob<-1-varphi
    chip<-qchisq(prob,p);
    ck<-sqrt(chip);
    cbeta<-( p*pchisq(chip,p+2)+ chip*(1-prob) )/p;
  }
  dl<-rep(0,pps)
  npat<-dim(misinfo)[1]
  n1<-misinfo[1,1]; p1<-misinfo[1,2];
  if (p1==p){
    ## complete data
    sigin<-solve(sig0)
    vecsig<-rsem.vec(sig0)
    Wmat<-kronecker(sigin,sigin)/2  ##Kronecker product
    for (i in 1:n1){
      xi<-x[i,]
      xi0<-xi-mu0;
      stdxi<-sigin%*%xi0
      stdxit<-t(stdxi)
      di2<-xi0%*%stdxi
      di<-sqrt(di2)
      ## Huber weight functions;
      if (di<=ck){
        wi1<-1
        wi2<-1/cbeta
        dwi1<-0
        dwi2<-0				
      }else{
        wi1<-ck/di
        wi2<-wi1*wi1/cbeta
        dwi1<-wi1/di2
        dwi2<-wi2/di2
      } ##end Huber weight
      ## computing B_\beta
      dlimu<-c(wi1)*stdxi
      xixi0<-xi0%*%t(xi0)
      vecyi<-rsem.vec(xixi0)
      wvecyi<-c(wi2)*vecyi
      dlisig<-dupt%*%Wmat%*%(wvecyi-vecsig)
      B11<-B11+dlimu%*%t(dlimu)
      B12<-B12+dlimu%*%t(dlisig)
      B22<-B22+dlisig%*%t(dlisig)
      dl_i<-c(dlimu, dlisig)
      dl<-dl+dl_i
      ## computing A_\beta
      Hi<-stdxi%*%stdxit
      tti<-c(wi1)*sigin
      uui<-c(wi2)*sigin
      ddl11<-ddl11 + (-tti+c(dwi1)*Hi)
      ddl22<-ddl22 + dupt%*%( Wmat - kronecker(Hi, (uui-.5*c(dwi2)*Hi) ) )%*%dup
      ddl12<-ddl12 + kronecker( (-tti+.5*c(dwi1)*Hi) ,stdxit)%*%dup
      ddl21<-ddl21 + dupt%*%kronecker( (-uui+c(dwi2)*Hi), stdxi )
    } ## for 1:n1		
  }else{			
    ## missing data
    if (varphi==0){
      ck1<-1e+10
      cbeta1<-1
    }else{
      chip1<-qchisq(prob,p1)
      ck1<-sqrt(chip1)
      cbeta1<-( p1*pchisq(chip1,p1+2) + chip1*(1-prob) )/p1
    }
    o1<-misinfo[1,3:(2+p1)]
    mu_o<-mu0[o1]
    sig_oo<-sig0[o1,o1]
    vecsig_oo<-rsem.vec(sig_oo)
    sigin_oo<-solve(sig_oo)
    E1<-i_p[o1,]; 
    if (o1==1){Et1=E1}else{Et1<-t(E1)}
    F1<-kronecker(E1, E1)%*%dup;
    Ft1<-t(F1)
    Wmat1<-.5*kronecker(sigin_oo, sigin_oo)
    for (i in 1:n1){
      xi<-x[i,]
      xi_o<-xi[o1]
      xi0_o<-xi_o-mu_o
      xi0_ot<-t(xi0_o)
      stdxi_o<-sigin_oo%*%xi0_o
      stdxit_o<-t(stdxi_o)
      di2<-xi0_o%*%stdxi_o
      di<-sqrt(di2)
      ## Huber weight functions;
      if (di<=ck){
        wi1<-1
        wi2<-1/cbeta
        dwi1<-0
        dwi2<-0				
      }else{
        wi1<-ck/di
        wi2<-wi1*wi1/cbeta
        dwi1<-wi1/di2
        dwi2<-wi2/di2
      } ##end Huber weight
      ## computing B_\beta
      dlimu<-c(wi1)*Et1%*%stdxi_o
      xixi0_o<-xi0_o%*%t(xi0_o)
      vecyi<-rsem.vec(xixi0_o)
      wvecyi<-c(wi2)*vecyi
      dlisig<-Ft1%*%Wmat1%*%(wvecyi-vecsig_oo)
      B11<-B11+dlimu%*%t(dlimu)
      B12<-B12+dlimu%*%t(dlisig)
      B22<-B22+dlisig%*%t(dlisig)
      dl_i<-c(dlimu, dlisig)
      dl<-dl+dl_i
      ## computing A_\beta
      Hi<-stdxi_o%*%stdxit_o
      tti<-c(wi1)*sigin_oo
      uui<-c(wi2)*sigin_oo
      ddl11<-ddl11 + Et1%*%(-tti+c(dwi1)*Hi)%*%E1
      ddl22<-ddl22 + Ft1%*%( Wmat1 - kronecker(Hi, (uui-.5*c(dwi2)*Hi) ) )%*%F1
      ddl12<-ddl12 + Et1%*%kronecker( (-tti+.5*c(dwi1)*Hi) ,stdxit_o)%*%F1
      ddl21<-ddl21 + Ft1%*%kronecker( (-uui+c(dwi2)*Hi), stdxi_o )%*%E1
      
    } ## end for 1:n1		
  }##end p1==p
  
  ## for patterns 2 and above
  if (npat>1){
    snj<-n1
    for (j in 2:npat){
      nj<-misinfo[j,1]
      pj<-misinfo[j,2]
      if (varphi==0){
        ckj<-1e+10
        cbetaj<-1
      }else{
        chipj<-qchisq(prob, pj)
        ckj<-sqrt(chipj)
        cbetaj<-( pj*pchisq(chipj,pj+2) + chipj*(1-prob) )/pj
      }
      oj<-misinfo[j, 3:(2+pj)]
      mu_o<-mu0[oj]
      sig_oo<-sig0[oj,oj]
      sigin_oo<-solve(sig_oo)
      vecsig_oo<-rsem.vec(sig_oo)
      Ej<-i_p[oj, ]
      if (pj==1){Etj<-Ej}else{Etj<-t(Ej)}
      Fj<-kronecker(Ej, Ej)%*%dup
      Ftj<-t(Fj)
      Wmati<-0.5*kronecker(sigin_oo, sigin_oo)
      for (i in (snj+1):(snj+nj)){
        xi<-x[i,]
        xi_o<-xi[oj]
        xi0_o<-xi_o-mu_o
        xi0_ot<-t(xi0_o)
        stdxi_o<-sigin_oo%*%xi0_o
        stdxit_o<-t(stdxi_o)
        di2<-xi0_o%*%stdxi_o
        di<-sqrt(di2)
        ## Huber weight functions;
        if (di<=ckj){
          wi1<-1
          wi2<-1/cbetaj
          dwi1<-0
          dwi2<-0				
        }else{
          wi1<-ckj/di
          wi2<-wi1*wi1/cbetaj
          dwi1<-wi1/di2
          dwi2<-wi2/di2
        } ##end Huber weight
        ## computing B_\beta
        dlimu<-c(wi1)*Etj%*%stdxi_o
        xixi0_o<-xi0_o%*%t(xi0_o)
        vecyi<-rsem.vec(xixi0_o)
        wvecyi<-c(wi2)*vecyi
        dlisig<-Ftj%*%Wmati%*%(wvecyi-vecsig_oo)
        B11<-B11+dlimu%*%t(dlimu)
        B12<-B12+dlimu%*%t(dlisig)
        B22<-B22+dlisig%*%t(dlisig)
        dl_i<-c(dlimu, dlisig)
        dl<-dl+dl_i
        ## computing A_\beta
        Hi<-stdxi_o%*%stdxit_o
        tti<-c(wi1)*sigin_oo
        uui<-c(wi2)*sigin_oo
        ddl11<-ddl11 + Etj%*%(-tti+c(dwi1)*Hi)%*%Ej
        ddl22<-ddl22 + Ftj%*%( Wmati - kronecker(Hi, (uui-.5*c(dwi2)*Hi) ) )%*%Fj
        ddl12<-ddl12 + Etj%*%kronecker( (-tti+.5*c(dwi1)*Hi) ,stdxit_o)%*%Fj
        ddl21<-ddl21 + Ftj%*%kronecker( (-uui+c(dwi2)*Hi), stdxi_o )%*%Ej
        
      } ## end for snj+1:snj+nj
      snj<-snj+nj
    } ##end for 2:npat
  }

  ## Constructing B_\bata and A_\beta matrices
  Bbeta<-rbind( cbind(B11, B12), cbind(t(B12), B22) )
  Abeta<-rbind( cbind(ddl11, ddl12), cbind(ddl21, ddl22) )
  Abin<-solve(Abeta)
  Omega<-n*Abin%*%Bbeta%*%t(Abin)
  Gamma<-Omega
  
  xnames<-colnames(x)
  if (is.null(xnames)) xnames<-paste('V', 1:p)
  
  mnames<-rsem.gname(xnames)
  colnames(Abeta)<-colnames(Bbeta)<-colnames(Gamma)<-rownames(Abeta)<-rownames(Bbeta)<-rownames(Gamma)<-c(xnames, mnames)
  
  list(Abeta=Abeta, Bbeta=Bbeta, Gamma=Gamma)
}


## R function for cronbach alpha
cronbach<-function(y, varphi=0.1, se=FALSE, complete=FALSE){
	if (complete) y<-na.omit(y)	
	p<-ncol(y)
	n<-nrow(y)
	rownames(y)<-1:n
	
	pat<-rsem.pattern(y)
	cov1<-rsem.emmusig(pat,varphi=varphi)
	
	weight<-rep(1, n)
	if (varphi>0) weight<-cov1$weight
	
	sigma<-cov1$sigma

	alpha<-p/(p-1)*(1-sum(diag(sigma))/sum(sigma))

	## robust standard error
	se.alpha<-NA
	if (se){
		s<-rsem.Ascov(pat, cov1, varphi=varphi)
		firstd<-NULL
		
		for (i in 1:p){
			for (j in i:p){
				if (i==j){
					temp<--p/(p-1)*(sum(sigma)-sum(diag(sigma)))/((sum(sigma))^2)
				}else{
					temp<- 2*p/(p-1)*sum(diag(sigma))/((sum(sigma))^2)
				}
				firstd<-c(firstd, temp)
			}
		}
		q<-p+p*(p+1)/2
		gamma<-s$Gamma[(p+1):q, (p+1):q]
		se.alpha<-sqrt(firstd%*%gamma%*%t(t(firstd)))/sqrt(n)
	}
	
	if (varphi>0) prop.down<-(n-sum(weight$w1==1))/n*100
	
	cat('The alpha is ', alpha, sep='')
	if (se) cat(' with the standard error ', se.alpha, sep='')
	cat('.\n')
	
	if (varphi>0) cat('About ', round(prop.down,2), '% of cases were downweighted.\n', sep='')
	
	alphaobject<-list(alpha=alpha, se=se.alpha, weight=weight, y=y, varphi=varphi, musig=cov1)
	class(alphaobject)<-'alpha'
	invisible(alphaobject)	
}


plot.alpha<-function(x, type="weight", profile=5, interval=0.01, center=TRUE, scale=FALSE, w1=FALSE, numbered=FALSE, pos='topright', ...){
	## type: weight, profile, diagnosis
	res<-x
	y<-res$y
	outcase.temp<-sum(res$weight$w1<1)
	## find the smallest weight
	if (profile==0){
		outcase<-5
		if (outcase.temp<5) outcase<-outcase.temp
	}
	if (outcase.temp<5) profile<-outcase.temp
		
	## find the smallest cases
	temp<-sort(res$weight$w1)
	
	if (profile==0){
		idx<-which(res$weight$w1 <= temp[outcase])
	}else{
		idx<-which(res$weight$w1 <= temp[profile])
	}
	
	if (substr(type,1,1)=="w"){
		if (w1){
			par(mfrow=c(2,1))
			plot(res$weight$w1, ylim=c(0,1), xlab='Case number', ylab='Weight w1', main='Weights for mean estimates',...)
			text(idx, res$weight$w1[idx], idx, pos=1, col='red', ...)
			plot(res$weight$w2, ylim=c(0,max(res$weight$w2)), xlab='Case number', ylab='Weight w2', main='Weights for covariance estimates',...)		
			text(idx, res$weight$w2[idx], idx, pos=1, col='red', ...)
			par(mfrow=c(1,1))
		}else{
			plot(res$weight$w2, ylim=c(0,max(res$weight$w2)), xlab='Case number', ylab='Weight w2',...)		
			text(idx, res$weight$w2[idx], idx, pos=1, col='red', ...)
		}
	} 
	
	if (substr(type,1,1)=="p"){
		## profile plot		
		if (center){
			if (scale){
				for (i in 1:nrow(y)){
				   y[i, ]<-(y[i, ]-res$musig$mu)/sqrt(diag(res$musig$sigma))
			   }
			}else{
				for (i in 1:nrow(y)){
				   y[i, ]<-y[i, ]-res$musig$mu
				}
			}			
		}
		l.min<-min(y, na.rm=TRUE)
		l.max<-max(y, na.rm=TRUE)
	
		p<-ncol(y)
		
		if (center){
			plot(1:p, rep(0,p), type='l', ylim=c(l.min, l.max), lwd=3, ylab='Score', xlab='', axes=FALSE, col='black', ...)
		}else{
			plot(1:p, res$musig$mu, type='l', ylim=c(l.min, l.max), lwd=3, ylab='Score', xlab='', axes=FALSE, col='black', ...)
		}
		box()
		axis(1, at=1:p, labels=names(y), las=2, ...)
		axis(2, ...)
		
		if (!is.null(pos)){
			ltyno<-1	
			idx.type<-NULL	
			for (i in idx){
				lines(1:p, y[i, ], lty=ltyno, ...)
				if (numbered & center) text(1, y[i,1], i)
				ltyno<-ltyno+1
				## type of outliers
				temp.type<-'(O)'
				if (max(y[i, ], na.rm=TRUE)<0)  temp.type<-'(L-)'
				if (min(y[i, ], na.rm=TRUE)>0)  temp.type<-'(L+)'
				idx.type<- c(idx.type, temp.type)
			}
			legend(pos, legend=paste(idx, idx.type), lty=1:ltyno)
		}
	}
	
	if (substr(type,1,1)=="d"){
		varphi<-res$varphi
		phi<-seq(0, varphi, by=interval)
		alpha.diag<-NULL
		for (i in 1:length(phi)){
			alpha.diag<-c(alpha.diag, cronbach(y, phi[i])$alpha)
		}
		plot(phi, alpha.diag, xlab='varphi', ylab='alpha',...)
	}
}

summary.alpha<-function(object, prob=.95,...){
    if (prob > .5) prob <- 1 - prob
	res<-object
	cat("\nThe estimated alpha is \n")
	
	t0.txt <- sprintf("  %-20s", "alpha")
	t1.txt <- sprintf("  %10.3f", res$alpha)
	cat(t0.txt, t1.txt, "\n", sep="")
  
	if (!is.na(res$se)){
		t0.txt <- sprintf("  %-20s", "se")
		t1.txt <- sprintf("  %10.3f", res$se)
		cat(t0.txt, t1.txt, "\n", sep="")
  
		## output p-value, too
		z<-abs(res$alpha/res$se)
		pvalue<-(1-pnorm(z))
		t0.txt <- sprintf("  %-20s", "p-value (alpha>0)")
		t1.txt <- sprintf("  %10.3f", pvalue)
		cat(t0.txt, t1.txt, "\n", sep="")
  
		## output confidence interval
		ci<-res$alpha+c(-1,1)*abs(qnorm(prob/2))*res$se
		if (ci[2]>1) ci[2]<-1.000
		if (ci[1]<0) ci[1]<-0.000
  
		t0.txt <- sprintf("  %-20s", "Confidence interval")
		t1.txt <- sprintf("  %10.3f", ci)
		cat(t0.txt, t1.txt, "\n\n", sep="")
	}
}

