###########################################################################
##                                                                       ##
## Assess improvements to false-discovery rate control in cFDR method    ##
## Functions for simulation                                              ##
##                                                                       ##
## James Liley, 15/2/18                                                  ##
##                                                                       ##
###########################################################################

###########################################################################
## Packages and scripts ###################################################
###########################################################################

#library(mnormt)
#library(mgcv)
#library(pbivnorm)
#library(MASS)
#library(fields)

###########################################################################
## Functions ##############################################################
###########################################################################

##' Return co-ordinates of L-regions for cFDR method, with or without estimation of Pr(H0|Pj<pj).
##' 
##' Parameter 'mode' defines the way in which L-curves are constructed. L-curves define a mapping 
##' from the unit square [0,1]^2 onto the unit interval [0,1], and we consider the mapped values 
##' of each point (p[i],q[i]). In this method (cFDR1) the mapping depends on the points used to 
##' generate L, and if these points coincide with the points we are using the map for, the 
##' behaviour of the map is very complicated. Parameter 'mode' governs the set of points used in 
##' generating L-curves, and can be used to ensure that the points used to define curves L (and 
##' hence the mapping) are distinct from the points the mapping is used on.
##' 
##' @title vl
##' @param p principal p-values
##' @param q conditional p-values
##' @param adj adjust cFDR values and hence curves L using estimate of Pr(H0|Pj<pj)
##' @param indices instead of at cfdr cutoffs, compute v(L) at indices of points. Overrides parameter at if set.
##' @param at cfdr cutoff/cutoffs. Defaults to null
##' @param mode determines set of variables to use for computing cFDR. Set to 0 to leave all of 'indices' in, 1 to remove each index one-by-one (only when computing the curve L for that value of (p,q)), 2 to remove all indices in variable 'fold' and compute curves L only using points not in 'fold'
##' @param fold If mode=2, only compute L-curves using points not in 'fold'. 
##' @param p_threshold if H0 is to be rejected automatically whenever p<p_threshold, include this in all regions L
##' @param nt number of test points in x-direction, default 5000
##' @param nv resolution for constructing L-curves, default 1000
##' @param scale return curves on the p- or z- plane. Y values are equally spaced on the z-plane.
##' @param closed determines whether curves are closed polygons encircling regions L (closed=T), or lines indicating the rightmost border of regions L
##' @param verbose print progress if mode=1
##' @export
##' @author James Liley
##' @return list containing elements x, y. Assuming n curves are calculated (where n=length(indices) or length(at)) and closed=T, x is a matrix of dimension n x (4+nv), y ix a vector of length (4+nv).
##' @examples 
##' # Generate standardised simulated dataset
##' set.seed(1); n=10000; n1p=100; n1pq=100; n1q=100
##' zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))
##' fold_id=(1:n) %% 3
##' 
##' 
##' # points to generate L-regions for
##' example_indices=c(4262, 268,83,8203)
##' 
##' 
##' 
##' par(mfrow=c(1,3))
##' 
##' v1=vl(p,q,indices=example_indices,mode=0,nv=5000); 
##' plot(p,q,cex=0.5,col="gray",xlim=c(0,0.001),ylim=c(0,1), main="All points in"); 
##' for (i in 1:length(example_indices)) lines(v1$x[i,],v1$y); 
##' for (i in 1:length(example_indices)) points(p[example_indices[i]],q[example_indices[i]],qh=16,col="blue")
##' 
##' # L1_{S-fold} example (fold left out; all points are in fold 2). Spikes disappear.
##' v2=vl(p,q,indices=example_indices,mode=2,fold=which(fold_id==2)); 
##' plot(p,q,cex=0.5,col="gray",xlim=c(0,0.001),ylim=c(0,1),main="Fold removed"); 
##' for (i in 1:length(example_indices)) lines(v2$x[i,],v2$y); 
##' for (i in 1:length(example_indices)) points(p[example_indices[i]],q[example_indices[i]],pch=16,col="blue")
##' 
##' # L1_{S-(p,q)} example (each point left out of points used to generate curve through that point). Spikes disappear.
##' v3=vl(p,q,indices=example_indices,mode=1); 
##' plot(p,q,cex=0.5,col="gray",xlim=c(0,0.001),ylim=c(0,1), main="Single point removed"); 
##' for (i in 1:length(example_indices)) lines(v3$x[i,],v3$y); 
##' for (i in 1:length(example_indices)) points(p[example_indices[i]],q[example_indices [i]],pch=16,col="blue")
##'
vl=function(p,q,adj=T,indices=NULL,at=NULL,mode=0,fold=NULL,nt=5000, nv=2000, p_threshold=0, scale=c("p","z"), closed=T,verbose=F) {
  
  zp=-qnorm(p/2); zq=-qnorm(q/2)
  
  if (any(!is.finite(zp+zq))) stop("P-values p,q must be in [1e-300,1]")
  
  #nt=5000; # number of test point in x direction
  #nv=1000; # resolution
  mx=max(c(10,abs(-qnorm(p/2)))); my=max(c(10,abs(-qnorm(q/2)))); #mx=10; my=10 # maximum limits for integration

  gx=0 #1/1000  # use for 'tiebreaks'- if a point is on a curve with nonzero area, force the L-curve through that point
  
  if (is.null(indices)) {
    if (is.null(at)) stop("One of the parameters 'indices', 'at' must be set")
    ccut=at
    mode=0
  }
    
  if (!is.null(indices)) { # set ccut. NOT equal to cfdr at the points; needs to be adjusted since an additional test point is used in determination of L
    if (mode==1) {
      yval2=seq(0,my,length.out=nv); xval2=outer(rep(1,length(indices)),yval2); pval2=2*pnorm(-yval2)
      xtest=seq(0,mx,length.out=nt); ptest=2*pnorm(-xtest)
      
      ccut=rep(0,length(indices))
      for (i in 1:length(indices)) {
        w=which(zq[-indices[i]] >= zq[indices[i]])
        if (length(w) >= 1) {
          cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[-indices[i]][w])(xtest)); cfsub=cummin(cfsub)
          ccut[i]=approx(xtest,cfsub,zp[indices[i]],rule=2)$y
        } else {
          #if (length(w)==1) 
          #  if (p[-indices[i]][w] < p[indices[i]]) ccut[i]= min(2*p[indices[i]],p[-indices[i]][w]) else  ccut[i]= p[indices[i]] 
            if (length(w)==0) ccut[i]= p[indices[i]] 
        }
      }
    } 
    if (mode==2) {
      
      yval2=seq(0,my,length.out=nv); xval2=outer(rep(1,length(indices)),yval2); pval2=2*pnorm(-yval2)
      xtest=seq(0,my,length.out=nt); ptest=2*pnorm(-xtest)
      
      #if (is.null(indices)) indices=indices
      
      ccut=rep(0,length(indices))
      for (i in 1:length(indices)) {
        w=which(zq[-fold] >= zq[indices[i]])
        if (length(w) >= 1) {
          cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[-fold][w])(xtest)); cfsub=cummin(cfsub)
          ccut[i]=approx(xtest,cfsub-gx*xtest + gx*mx,zp[indices[i]],rule=2)$y
        } else ccut[i]=p[indices[i]]*1/(1+ length(which(zp[-fold][w]>= zp[indices[i]])))
      }
      #    ccut= p[indices]*sapply(indices,function(x) 
      #      (1+length(which(q[-fold] <= q[x])))/(1+length(which(p[-fold] <= p[x] & q[-fold] <= q[x]))))# set ccut
    }
    if (mode==0) {
      yval2=seq(0,my,length.out=nv); xval2=outer(rep(1,length(indices)),yval2); pval2=2*pnorm(-yval2)
      xtest=seq(0,my,length.out=nt); ptest=2*pnorm(-xtest)
      
      ccut=rep(0,length(indices))
      for (i in 1:length(indices)) {
        w=which(zq >= zq[indices[i]])
        if (length(w) >= 2) {
          cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[w])(xtest)); cfsub=cummin(cfsub)
          ccut[i]=approx(xtest,cfsub-gx*xtest + gx*mx,zp[indices[i]],rule=2)$y
        } else ccut[i]=p[indices[i]]*1/(1+ length(which(zp[w]>= zp[indices[i]])))
      }
      #    ccut=p[indices]*sapply(indices,function(x) 
      #      (1+length(which(q <= q[x])))/(1+length(which(p <= p[x] & q <= q[x]))))# set ccut
    }
  }  
  
  ccut=ccut*(1+ 1e-6) # ccut + 1e-8 # prevent floating-point comparison errors

  if (verbose & mode==1) print(paste0(length(ccut)," regions to calculate"))

  out=rep(0,length(ccut))
  
  if (mode==0) {
    yval2=seq(0,my,length.out=nv); xval2=outer(rep(1,length(ccut)),yval2); pval2=2*pnorm(-yval2)
    xtest=seq(0,mx,length.out=nt); ptest=2*pnorm(-xtest)
    zp_ind=ceiling(zp*nt/mx); zp_ind=pmax(1,pmin(zp_ind,nt))
    zq_ind=ceiling(zq*nv/my); zq_ind=pmax(1,pmin(zq_ind,nv))
    
    zq[which(zq>my)]=my; zp[which(zp>mx)]=mx; p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))
    
    if (adj) correct=cummin((1+ecdf(q[which(p>0.5)])(pval2)*length(p))/
                              (1+ ecdf(q)(pval2)*length(p))) else correct=rep(1,length(pval2)) # adjustment factor for q[i]
    if (!is.null(indices)) ccut=ccut*correct[zq_ind[indices]]
    
    
    for (i in 1:length(yval2)) {
      w=which(zq > yval2[i])
      if (length(w)>= 1) {
        cfsub= ptest/(1+(1/length(w))-ecdf(zp[w])(xtest)); cfsub=cummin(cfsub)
        xval2[,i]=approx(cfsub*correct[i]-gx*xtest + gx*mx,xtest,ccut,rule=2,method="const",f=1)$y 
      }  else xval2[,i]=-qnorm(ccut/2)
    }
    
  } 
  if (mode==1) {
    
    yval2=seq(0,my,length.out=nv); xval2=outer(rep(1,length(ccut)),yval2); pval2=2*pnorm(-yval2)
    xtest=seq(0,mx,length.out=nt); ptest=2*pnorm(-xtest)
    
    zp_ind=ceiling(zp*nt/mx); zp_ind=pmax(1,pmin(zp_ind,nt))
    zq_ind=ceiling(zq*nv/my); zq_ind=pmax(1,pmin(zq_ind,nv))
    # populate zp_ind[i],zq_ind[i] with the index of xval2,yval2 closest (least above?) zp[i],zq[i]. Only need to do at 'sub'
    
    xmat=matrix(0,nt,nv) # base
    
    pqmat=xmat
    qvec=rep(0,nv)
  
    for (i in 1:length(p)) {
      pqmat[1:zp_ind[i],1:zq_ind[i]]=1+pqmat[1:zp_ind[i],1:zq_ind[i]]
      qvec[1:zq_ind[i]]=1+qvec[1:zq_ind[i]]
    }
    # matrix [i,j] is 1+ number of SNPs with zp>xval2[i]; zq>yval2[i], only works for i in sub
    pqmat=1+pqmat
    qvec=1+qvec
    
    #qvec=1+  length(zq_ind)*(1-ecdf(zq_ind)(1:nv))
    # qvec[i] is 1+ number of SNPs with zq>yval[i], only for i in sub
    
    cf_mat=outer(ptest,qvec)/pqmat #pmat*qmat/pqmat
    # [i,j] has cFDR value at xtest[i], ytest[i]; only valid for i in sub
    cf_mat=apply(cf_mat,2,cummin) # Smooth shape of L, see paragraph in methods
    
    l_new=rep(0,length(p)) # set to L in new method
    
    for (i in 1:length(indices)) {
      
      if (adj) correctx=cummin((1+ecdf(q[-indices[i]][which(p[-indices[i]]>0.5)])(pval2)*length(p[-indices[i]]))/
                              (1+ ecdf(q[-indices[i]])(pval2)*length(p[-indices[i]]))) else correctx=rep(1,length(pval2)) # adjustment factor for q[i]
      ccut[i]=ccut[i]*correctx[zq_ind[indices[i]]]
      pqnew=pqmat
      pqnew[1:zp_ind[indices[i]],1:zq_ind[indices[i]]]=-1+pqnew[1:zp_ind[indices[i]],1:zq_ind[indices[i]]] #remove contribution of current SNP
      qnew=qvec
      qnew[1:zq_ind[indices[i]]]=-1+qnew[1:zq_ind[indices[i]]]
      cfx=apply(outer(ptest,qnew)/pqnew,2,cummin)
      cfx=t(t(cfx)*correctx)
      #cxi=interp.surface(list(x=xtest,y=yval2,zp=cfx),cbind(zp[i],zq[i])); 
      #if (is.na(cxi)) cxi=cfx[zp_ind[indices[i]],zq_ind[indices[i]]]
      cxi=ccut[i]
     xv=suppressWarnings(apply(cfx,2,function(x) xtest[max(which(x> cxi))]))
      xv[which(xv<0)]=0; 
      xv[which(!is.finite(xv))]=0
      xval2[i,]=xv;
      if (verbose) print(i)
    }
    
  }
  if (mode==2) {
    yval2=seq(0,my,length.out=nv); xval2=outer(rep(1,length(ccut)),yval2); pval2=2*pnorm(-yval2)
    xtest=seq(0,mx,length.out=nt); ptest=2*pnorm(-xtest)
    zp_ind=ceiling(zp[indices]*nt/mx); zp_ind=pmax(1,pmin(zp_ind,nt))
    zq_ind=ceiling(zq[indices]*nv/my); zq_ind=pmax(1,pmin(zq_ind,nv))
 
    if (adj) correct=cummin((1+ecdf(q[which(p[-fold]>0.5)])(pval2)*length(p[-fold]))/
                              (1+ ecdf(q[-fold])(pval2)*length(p[-fold]))) else correct=rep(1,length(pval2)) # adjustment factor for q[i]
    ccut=ccut*correct[zq_ind]
    
    
    for (i in 1:length(yval2)) {
      w=which(zq[-fold] > yval2[i])
      if (length(w) >= 1 ) {
        cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[-fold][w])(xtest)); cfsub=cummin(cfsub)
        xval2[,i]=approx(cfsub*correct[i]-gx*xtest + gx*mx,xtest,ccut,rule=2,f=1)$y 
      }  else xval2[,i]=-qnorm(ccut/2)
    }
    
  } 
  
  xval2[which(xval2> -qnorm(p_threshold/2))]=-qnorm(p_threshold/2)
  if (closed) {
    yval2=c(Inf,0,yval2,Inf,Inf)
    xval2=cbind(Inf,Inf,xval2,xval2[,nv],Inf)
  }
  
  if (scale[1]=="p") {
    X=2*pnorm(-abs(xval2))
    Y=2*pnorm(-abs(yval2))
  } else {
    X=xval2
    Y=yval2
  }
  
  return(list(x=X,y=Y))  
  
}











##' Return co-ordinates of L-regions for cFDR method using four-groups method, with or without estimation
##'  of Pr(H0|Pj<pj).
##'
##' Assumes (P,Q) follows a bivariate mixture-Gaussian distribution with four components, each centred at
##' the origin and with covariance matrices I2, (s1^2,0; 0,1), (1,0; 0,t1^2), (s2^2,0; 0,t2^2). Components
##' have weights pi0,pi1,pi2 and (1-pi0-pi1-pi2) respectively. Function fit.4g fits parameters to data.
##' 
##' This function does not have a 'mode' option (as for function vl) since the mapping from
##' [0,1]^2 -> [0,1] defined by L-curves depends only on 'pars' in this case and not on observed p,q.
##' 
##' @title vlx
##' @param p principal p-values
##' @param q conditional p-values
##' @param pars fitted parameters governing distribution of (P,Q). A seven element-vector: (pi0,pi1,pi2,s1,s2,t1,t2) respectively. See function description.
##' @param adj adjust cFDR values and hence curves L using estimate of Pr(H0|Pj<pj)
##' @param indices compute v(L) at indices of points. Overrides parameter at if set.
##' @param at cfdr cutoff/cutoffs. Defaults to null
##' @param nt number of test points in x-direction, default 5000
##' @param nv resolution for constructing L-curves, default 1000
##' @param p_threshold if H0 is to be rejected automatically whenever p<p_threshold, include this in all regions L
##' @param scale return curves on the p- or z- plane. Y values are equally spaced on the z-plane.
##' @param closed determines whether curves are closed polygons encircling regions L (closed=T), or lines indicating the rightmost border of regions L
##' @author James Liley
##' @export
##' @return list containing elements x, y. Assuming n curves are calculated and closed=T (where n=length(indices) or length(at)), x is a matrix of dimension n x (4+nv), y ix a vector of length (4+nv).
##' @examples 
##' # Generate standardised simulated dataset
##' set.seed(1); n=10000; n1p=100; n1pq=100; n1q=100
##' zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))
##' fold_id=(1:n) %% 3
##' 
##' # estimate parameters of underying dataset
##' fit_pars=fit.4g(cbind(zp,zq))$pars
##' 
##' # estimate parameters of underying dataset, removing fold 1
##' fit_pars_fold23=fit.4g(cbind(zp[which(fold_id!=1)],zq[which(fold_id!=1)]))$pars
##' 
##' 
##' # points to generate L-regions for
##' example_indices=c(4262, 268,83,8203)
##' 
##' par(mfrow=c(1,2))
##' 
##' v=vlx(p,q,pars=fit_pars,indices=example_indices); 
##' plot(p,q,cex=0.5,col="gray",xlim=c(0,0.001),ylim=c(0,1), main="All points in"); 
##' for (i in 1:length(example_indices)) lines(v$x[i,],v$y); 
##' for (i in 1:length(example_indices)) points(p[example_indices[i]],q[example_indices[i]],pch=16,col="blue")
##' 
##' v2=vlx(p,q,pars=fit_pars_fold23,indices=example_indices); 
##' plot(p,q,cex=0.5,col="gray",xlim=c(0,0.001),ylim=c(0,1), main="Fold removed"); 
##' for (i in 1:length(example_indices)) lines(v2$x[i,],v2$y); 
##' for (i in 1:length(example_indices)) points(p[example_indices[i]],q[example_indices[i]],pch=16,col="blue")
##'
##' 
vlx=function(p,q,pars,adj=T,indices=NULL,at=NULL,p_threshold=0,nt=5000, nv=1000,scale=c("p","z"),closed=T) {
  
  zp=-qnorm(p/2); zq=-qnorm(q/2)

  if (any(!is.finite(zp+zq))) stop("P-values p,q must be in [1e-300,1]")
  
  #nt=5000; # number of test point in x direction
  #nv=1000; # resolution
  mx=max(c(10,abs(-qnorm(p/2)))); my=max(c(10,abs(-qnorm(q/2)))); #mx=10; my=10 # maximum limits for integration
  
  
  if (is.null(indices)) {
    if (is.null(at)) stop("One of the parameters 'indices', 'at' must be set")
    ccut=at
    mode=0
  }
    
  # functions for computation
  denom1=function(zp,zq,pars) {
    2*pars[1]*pnorm(-zp)*pnorm(-zq) + 
    2*pars[2]*pnorm(-zp/pars[4])*pnorm(-zq) +
    2*pars[3]*pnorm(-zp)*pnorm(-zq/pars[6]) +
    2*(1-pars[1]-pars[2]-pars[3])*pnorm(-zp/pars[5])*pnorm(-zq/pars[7]); }
  denom2=function(zp,zq,pars) {
    (pars[1]+pars[2])*pnorm(-zq) + pars[3]*pnorm(-zq/pars[6]) +
      (1-pars[1]-pars[2]-pars[3])*pnorm(-zq/pars[7]) ;}
  denom=function(zp,zq,pars) denom1(zp,zq,pars)/denom2(zp,zq,pars)
  raw_cfx=function(p,q) p/denom(-qnorm(p/2),-qnorm(q/2),pars)
  
  # inverse denom
  id=function(ccut,y,pars) 
    uniroot(function(x) raw_cfx(2*pnorm(-x),2*pnorm(-y))-ccut, c(0,20))$root
  
  if (!is.null(indices)) {
    ccut=raw_cfx(p[indices],q[indices])
    #ccut=rep(1,length(indices))
    #for (i in 1:length(indices))
    #  ccut[i]= min(raw_cfx(seq(p[indices[i]],1,length.out=1000),q[indices[i]]))# set ccut
  }
  
  zp=-qnorm(p/2); zq=-qnorm(q/2)
  out=rep(0,length(ccut))
  
  yval2=seq(0,my,length.out=nv); xval2=outer(rep(1,length(ccut)),yval2); pval2=2*pnorm(-yval2)
  xtest=seq(0,mx,length.out=nt); ptest=2*pnorm(-xtest)
  
  if (adj) correct=cummin((1+ecdf(q[which(p>0.5)])(pval2)*length(p))/
                            (1+ ecdf(q)(pval2)*length(p))) else correct=rep(1,length(pval2)) # adjustment factor for q[i]
  
  for (i in 1:length(yval2)) {
    cfx=cummin(raw_cfx(ptest,2*pnorm(-yval2[i])))
    xval2[,i]=approx(cfx,xtest,ccut/correct[i],rule=2,method="const",f=1)$y
  }
  

  xval2[which(xval2> -qnorm(p_threshold/2))]=-qnorm(p_threshold/2)
  if (closed) {
    yval2=c(Inf,0,yval2,Inf,Inf)
    xval2=cbind(Inf,Inf,xval2,xval2[,nv],Inf)
  }
  
  if (scale[1]=="p") {
    X=2*pnorm(-abs(xval2))
    Y=2*pnorm(-abs(yval2))
  } else {
    X=xval2
    Y=yval2
  }
  
  return(list(x=X,y=Y))  
  

}





##' Return co-ordinates of L-regions for cFDR method using kernel density method, with or without estimation
##'  of Pr(H0|Pj<pj).
##'
##' Estimates empirical distribution of (P,Q) by fitting a kernel density estimate to observed values.
##' 
##' 
##' @title vly
##' @param p principal p-values
##' @param q conditional p-values
##' @param adj adjust cFDR values and hence curves L using estimate of Pr(H0|Pj<pj)
##' @param indices compute v(L) at indices of points. Overrides parameter at if set.
##' @param at cfdr cutoff/cutoffs. Defaults to null
##' @param mode set to 0 to leave all of 'indices' in, 1 (NOT CURRENTLY SUPPORTED) to remove each index only when computing L at that point, 2 to remove all of 'indices' from p,q
##' @param fold If mode=2, only compute L-curves using points not in 'fold'. 
##' @param p_threshold if H0 is to be rejected automatically whenever p<p_threshold, include this in all regions L
##' @param nt number of test points in x-direction, default 5000
##' @param nv resolution for constructing L-curves, default 1000
##' @param scale return curves on the p- or z- plane. Y values are equally spaced on the z-plane.
##' @param closed determines whether curves are closed polygons encircling regions L (closed=T), or lines indicating the rightmost border of regions L
##' @param ... other parameters passed to function kde2d. Can be used to set a non-Gaussian kernel.
##' @author James Liley
##' @export
##' @return list containing elements x, y. Assuming n curves are calculated and closed=T (where n=length(indices) or length(at)), x is a matrix of dimension n x (4+nv), y ix a vector of length (4+nv).
##' @examples 
##' # Generate standardised simulated dataset
##' set.seed(1); n=10000; n1p=100; n1pq=100; n1q=100
##' zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))
##' fold_id=(1:n) %% 3
##' 
##' # points to generate L-regions for
##' example_indices=c(4262, 268,83,8203) #c(164,23,74)
##' 
##' par(mfrow=c(1,2))
##' 
##' v1=vly(p,q,indices=example_indices,mode=0); 
##' plot(p,q,cex=0.5,col="gray",xlim=c(0,0.001),ylim=c(0,1), main="All points in"); 
##' for (i in 1:length(example_indices)) lines(v1$x[i,],v1$y); 
##' for (i in 1:length(example_indices)) points(p[example_indices[i]],q[example_indices[i]],pch=16,col="blue")
##' 
##' v2=vly(p,q,indices=example_indices,mode=2,fold=which(fold_id==1)); 
##' plot(p,q,cex=0.5,col="gray",xlim=c(0,0.001),ylim=c(0,1), main="Fold removed"); 
##' for (i in 1:length(example_indices)) lines(v2$x[i,],v2$y); 
##' for (i in 1:length(example_indices)) points(p[example_indices[i]],q[example_indices[i]],pch=16,col="blue")
##' 
vly=function(p,q,adj=T,indices=NULL,at=NULL,mode=0,fold=NULL,p_threshold=0,nt=5000,nv=1000,scale=c("p","z"),closed=T,...) {
  
  res=200
  
  if (mode==1) stop("Leave-one-out method is not directly supported in this function, but can be implemented for individual indices i using mode=2, indices=i, and fold=(1:length(p))[-i]")
  
  mx=max(c(10,abs(-qnorm(p/2)))); my=max(c(10,abs(-qnorm(q/2)))); 
  # mx=10; my=10 # maximum limits for integration
  
  zp=-qnorm(p/2); zq=-qnorm(q/2)
  
  if (is.null(indices)) mode=0
  if (mode==2) sub=setdiff(1:length(zp),fold) else sub=1:length(zp)
  
  zp[which(zp > mx)]=mx
  zq[which(zq > my)]=my
  
  zp[sub][which(zp>mx)]=mx; zq[sub][which(zq>my)]=my; w=which(is.finite(zp[sub]+zq[sub])); rr=max(c(12,mx,my)) 
  kpq=kde2d(zp[sub],zq[sub],n=res,lims=c(0,rr,0,rr),...)
  int2_kpq=t(apply(apply(kpq$zp[res:1,res:1],2,cumsum),1,cumsum))[res:1,res:1]
  int2_kpq[which(int2_kpq==0)]=min(int2_kpq[which(int2_kpq>0)]) # avoid 0/0 errors
  int_kpq=t(outer(int2_kpq[1,],rep(1,res)))
  int_p=outer(2*pnorm(-kpq$x),rep(1,res))
  kgrid=kpq; kgrid$z=int2_kpq/int_kpq; kgrid$z[which(kgrid$z< 1/length(sub))]=1/length(sub)
  cgrid=kgrid; cgrid$z=int_p/kgrid$zp; cgrid$z=apply(cgrid$z,2,cummin)

  
  if (!is.null(indices)) 
    ccut= interp.surface(cgrid,cbind(pmin(zp[indices],mx),pmin(zq[indices],my)))
  out=rep(0,length(ccut))
  
  yval2=seq(0,my,length.out=nv); xval2=outer(rep(1,length(ccut)),yval2); pval2=2*pnorm(-yval2)
  xtest=seq(0,my,length.out=nt); ptest=2*pnorm(-xtest)
  
  if (adj) correct=cummin((1+ecdf(q[which(p>0.5)])(pval2)*length(p))/
                            (1+ ecdf(q)(pval2)*length(p))) else correct=rep(1,length(pval2)) # adjustment factor for q[i]
  
  
  for (i in 1:length(yval2)) {
    w=which(zq > yval2[i])
    xdenom=interp.surface(kgrid,cbind(xtest,rep(yval2[i],length(xtest))))
    if (length(w)>2) {
      cfx=cummin(ptest/xdenom)
      xval2[,i]=approx(cfx,xtest,ccut/correct[i],rule=2,method="const",f=1)$y
    #} else xval2[,i]=-qnorm(ccut/2)    
    } else if (i>1) xval2[,i]=xval2[,i-1] else xval2[,i]=-qnorm(ccut/2)    
  }
  
  xval2[which(xval2> -qnorm(p_threshold/2))]=-qnorm(p_threshold/2)
  if (closed) {
    yval2=c(Inf,0,yval2,Inf,Inf)
    xval2=cbind(Inf,Inf,xval2,xval2[,nv],Inf)
  }
  
  if (scale[1]=="p") {
    X=2*pnorm(-abs(xval2))
    Y=2*pnorm(-abs(yval2))
  } else {
    X=xval2
    Y=yval2
  }
  
  return(list(x=X,y=Y))  
  
}




##' Integrate over L (general). Assumes Zq| H^p=0 has a mixture distribution, being N(0,1) with probability pi0, and taking some other distribution distx with probability (1-pi0)
##' 
##' We generally assume that distx is a Gaussian centred at 0. 
##'
##' @title il
##' @param X either output from functions vl, vlx, or vly, or matrix nk x nc of x-co-ordinates (Z-plane) of regions to integrate over. X[k,] is the set of co-ordinates for the kth region.
##' @param Y vector of length nc of y-co-ordinates to integrate over. Assumed to be constant for all columns of X. Leave as NULL if X is output from vl, vlx, or vly.
##' @param pi0_null parameter for distribution of Q|H^p=0. Can be a vector of parameters of length np.
##' @param sigma_null scale parameter for distribution of Q|H^p=0. Can be a vector of parameters
##' @param rho_null optional parameter governing covariance between Z scores under the null; for instance, from shared controls
##' @param distx distribution type for distribution of Q|H^p=0. Should be a text string which can be appended to 'd' to get PDF and 'p' to get CDF
##' @param ... additional parameters passed to CDF and PDF functions, ie df=3
##' @author James Liley
##' @export
##' @return matrix of dimension nk x np; [k,p]th element is the integral for the kth region using the pth parameter values.
##' @examples 
##' # Generate standardised simulated dataset
##' set.seed(1); n=10000; n1p=100; n1pq=100; n1q=100
##' zp=c(rnorm(n1p,sd=2), rnorm(n1q,sd=1),rnorm(n1pq,sd=2), rnorm(n-n1p-n1q-n1pq,sd=1))
##' zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=2),rnorm(n1pq,sd=2), rnorm(n-n1p-n1q-n1pq,sd=1))
##' p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))
##'
##' # Generate some L regions
##' example_indices=c(70,67,226)
##' v1=vl(p,q,indices=example_indices,mode=0); 
##' 
##' 
##' # The true distribution of Zq|H^p=0 is N(0,1) with weight n1q/(n-n1p-n1pq), and N(0,2^2) with weight 1- (n1q/(n-n1p-n1pq))
##' true_q0_pars=c(n1q/(n-n1p-n1pq),2)
##' 
##' # Estimate parameters:
##' est_q0_pars=fit.2g(q[which(p>0.5)])$pars
##' 
##' ############### Integrals ################
##' int_true=il(v1,pi0_null=true_q0_pars[1],sigma_null=true_q0_pars[2],distx="norm")
##' int_est=il(v1,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
##' ##########################################
##'
##' ############# Check integral #############
##' # Sample values zp0,zq0 and p0,q0 following distribution of Zp,Zq|H^p=0
##' nsc=1000000 # generate nsc sample values
##' n0=round(nsc*true_q0_pars[1]); n1=nsc-n0
##' zp0=rnorm(nsc,sd=1)
##' zq0=c(rnorm(n0,sd=1), rnorm(n1,sd=true_q0_pars[2]))
##' p0=2*pnorm(-abs(zp0)); q0=2*pnorm(-abs(zq0))
##' 
##' 
##' # Proportion of values p0,q0 falling in region with co-ordinates v1$x[2,],v1$y
##' length(which(in.out(cbind(v1$x[2,],v1$y),cbind(p0,q0))))/nsc
##' 
##' # Value of integral over the region
##' int_true[2] 
##' ##########################################
il=function(X,Y=NULL,pi0_null=NULL,sigma_null=rep(1,length(pi0_null)),rho_null=0,distx=c("norm","t","cauchy"),...) {
  
  if (is.null(Y) & length(X)!=2) stop("If Y is null, X must be a two-element list with components x and y") 
  if (is.null(Y)) {
    Y=X$y; X=X$x 
  }
  
  if (is.null(pi0_null | is.null(sigma_null))) stop("Parameters pi0_null and sigma_null must be set")
  
  # if X,Y are on the p-value scale
  if (all(Y<= 1)) {
    X=-qnorm(X/2); Y=-qnorm(Y/2)
  }
  
  # if X,Y define closed polygons
  if (any(!is.finite(Y))) {
    ntemp=length(Y)
    X=matrix(X[,3:(ntemp-3)],nrow=nrow(X)); Y=Y[3:(ntemp-3)]
  }
  
  pdf=get(paste0("d",distx[1]))
  cdf=get(paste0("p",distx[1]))
  
  nv=length(Y)
  
  # integrate
  nn=length(pi0_null); out=c()
  
  ysc=(range(Y)[2]-range(Y)[1])/length(Y)
  yw=which.max(Y)
  if (max(abs(rho_null))<0.001) {
    for (ii in 1:nn) {
      ypart= ysc*colSums(4*pnorm(-t(X))*(pi0_null[ii]*dnorm(Y) + 
                                           (1-pi0_null[ii])*pdf(Y/sigma_null[ii],...)/sigma_null[ii])) 
      infpart=  4*pnorm(-X[,yw])*(pi0_null[ii]*pnorm(-Y[yw]) + 
                                    (1-pi0_null[ii])*cdf(-Y[yw]/sigma_null[ii],...))
      out=c(out,ypart+infpart)
    }
  } else {
    if (length(rho_null)==1) rho_null=rep(rho_null,length(pi0_null))
    for (ii in 1:nn) {
      m1a=rho_null[ii]*Y; m1b=rho_null[ii]*Y/(sigma_null[ii]^2) # means of conditional distributions given y=Y
      s1a=sqrt(1-(rho_null[ii]^2)); s1b=sqrt(1-(rho_null[ii]/sigma_null[ii])^2)
      ypart1=ysc* colSums(2*(pnorm(-t(X),mean=m1a,sd=s1a)*pi0_null[ii]*dnorm(Y) +
                               pnorm(-t(X),mean=m1b,sd=s1b)*(1-pi0_null[ii])*
                               pdf(Y/sigma_null[ii],...)/sigma_null[ii])) # positive rho part
      ypart2=ysc* colSums(2*(pnorm(-t(X),mean=-m1a,sd=s1a)*pi0_null[ii]*dnorm(Y) +
                               pnorm(-t(X),mean=-m1b,sd=s1b)*(1-pi0_null[ii])*
                               pdf(Y/sigma_null[ii],...)/sigma_null[ii])) # negative rho part
      infpart1=2*pi0_null[ii]*pbivnorm(-X[,yw],-Y[yw],rho=rho_null[ii]) +    
        2*(1-pi0_null[ii])*pbivnorm(-X[,yw],-Y[yw]/sigma_null[ii],rho=rho_null[ii]/sigma_null[ii])
      infpart2=2*pi0_null[ii]*pbivnorm(-X[,yw],-Y[yw],rho=-rho_null[ii]) +    
        2*(1-pi0_null[ii])*pbivnorm(-X[,yw],-Y[yw]/sigma_null[ii],rho=-rho_null[ii]/sigma_null[ii])
      out=c(out,ypart1+ypart2 + infpart1+infpart2)
    }
  }
  out
}






##' Fit a specific two Guassian mixture distribution to a set of P or Z values.
##'
##' Assumes 
##' Z ~ N(0,1) with probability pi0, Z ~ N(0,sigma^2) with probability 1-pi0
##'
##' Returns MLE for pi0 and sigma. Uses R's optim function. Can weight observations.
##' 
##' @title fit.2g
##' @param P numeric vector of observed data, either p-values or z-scores. If rho=0, should be one-dimensional vector; if rho is set, should be bivariate observations (P,Q)
##' @param pars initial values for parameters
##' @param weights optional weights for parameters
##' @param sigma_range range of possible values for sigma (closed interval). Default [1,100]
##' @return a list containing parameters pars, likelihoods under h1 (Z distributed as above), likelihood under h0 (Z~N(0,1)) and likelihood ratio lr.
##' @export
##' @author James Liley
##' @examples
##' sigma=2; pi0 <- 0.8
##' 
##' n=10000; n0=round(pi0*n); n1=n-n0
##' Z = c(rnorm(n0,0,1),rnorm(n1,0,sqrt(1+ (sigma^2))))
##' fit=fit.2g(Z)
##' fit$pars
fit.2g=function(P, pars = c(0.5, 1.5), weights = rep(1, min(length(Z),dim(Z)[1])), 
                sigma_range = c(1,100),rho=0,...) {
  if (all(P<=1) & all(P>= 0)) Z=-qnorm(P/2) else Z=P # P-values or Z scores
  
  pars = as.numeric(pars)
  #Z = as.numeric(Z)
  if (length(pars) != 2) 
    stop("Parameter 'pars' must be a two-element vector containing values pi0, s")
  if (pars[1] >= 1 | pars[1] <= 0 | pars[2] <= sigma_range[1]) 
    stop("Initial value of pi0 must all be between 0 and 1 and initial value of s must be greater than 0 (or sigma_range[1] if set)")
  if (abs(rho)< 0.001) {
    if (length(dim(Z))>1) Z=Z[,2]
    w = which(!is.na(Z * weights))
    Z = Z[w]
    weights = weights[w]
    w1=which(abs(Z)<30); w2=setdiff(1:length(Z),w1)
    l2 = function(pars = c(0.5, 1.5)) -(sum(weights[w1] * log(pars[1] * 
                                                                dnorm(Z[w1], sd = 1) + (1 - pars[1]) * dnorm(Z[w1], sd = pars[2])))) +
      -(sum(weights[w2] * (log(1 - pars[1]) + dnorm(Z[w2], sd = pars[2],log=T)))) 
  } else {
    w1=which(abs(Z[,2])<30); w2=setdiff(1:dim(Z)[1],w1)
    w = which(!is.na(Z[,1]*Z[,2] * weights))
    Z = Z[w,]
    weights = weights[w]
    l2=function(pars = c(0.5, 1.5)) {
      v1=rbind(c(1,rho),c(rho,1))
      v1r=rbind(c(1,-rho),c(-rho,1))
      v2=rbind(c(1,rho),c(rho,pars[2]^2))
      v2r=rbind(c(1,-rho),c(-rho,pars[2]^2))
      pi0=pars[1]
      fw1= pi0*dmnorm(Z[w1,],varcov=v1) + (1-pi0)*dmnorm(Z[w1,],varcov=v2) 
      fw2= sum(weights[w2]*(log(1-pi0) + dmnorm(Z[w2,],varcov=v2,log=T) ))
      fw1r= pi0*dmnorm(Z[w1,],varcov=v1r) + (1-pi0)*dmnorm(Z[w1,],varcov=v2r)  
      fw2r= sum(weights[w2]*(log(1-pi0) + dmnorm(Z[w2,],varcov=v2r,log=T) ))
      -sum(weights[w1]*log(fw1+fw1r)) - (fw2+ fw2r)
    }
  }
  zx = optim(pars, function(p) l2(p), lower = c(1e-05, sigma_range[1]), 
             upper = c(1 - (1e-05), sigma_range[2]), method = "L-BFGS-B", control = list(factr = 10), 
             ...)
  h1 = -zx$value
  h0 = -l2()
  yy = list(pars = zx$par, h1value = h1, h0value = h0, lr = h1 - 
              h0)
  return(yy)
}







##' Fit a four-part mixture normal model to bivariate data. Assumes that data are distributed as one of
##'  N(0,1) x N(0,1)        with probability pi0
##'  N(0,s1^2) x N(0,1)     with probability pi1
##'  N(0,1) x N(0,t1^2)     with probability pi2
##'  N(0,s2^2) x N(0,t2^2)  with probability 1-pi0-pi1-pi2
##' Fits the set of parameters (pi0,pi1,pi2,s1,s2,t1,t2) using an E-M algorithm
##' 
##' @title fit.4g
##' @param P matrix N x 2 of data points (Z scores or P-values)
##' @param pars initial parameter values
##' @param weights weights for points
##' @param C include additive term C*log(pi0*pi1*pi2*(1-pi0-pi1-pi2)) in objective function to ensure identifiability of model
##' @param maxit maximum number of iterations (supersedes tol)
##' @param tol stop after increment in log-likelihood is smaller than this
##' @param sgm force s1,s2,t1,t2 to be at least this value
##' @export
##' @author James Liley
##' @return list with elements pars (fitted parameters), lhood (log likelihood) and hist (fitted parameters during algorithm
##' @examples 
##' pi0=0.5; pi1=0.15; pi2=0.25
##' s1=3; s2=2; t1=2; t2=3
##' true_pars=c(pi0,pi1,pi2,s1,s2,t1,t2)
##' 
##' 
##' n=100000; n0=round(pi0*n); n1=round(pi1*n); n2=round(pi2*n); n3=n-n0-n1-n2
##' 
##' zs=c(rnorm(n0,sd=1),rnorm(n1,sd=s1),rnorm(n2,sd=1),rnorm(n3,sd=s2))
##' zt=c(rnorm(n0,sd=1),rnorm(n1,sd=1),rnorm(n2,sd=t1),rnorm(n3,sd=t2))
##' 
##' Z=cbind(zs,zt)
##' 
##' f=fit.4g(Z)
##' f$pars
fit.4g= function (P, pars = c(0.7, 0.1,0.1, 2, 2, 2, 2), weights = rep(1, 
    dim(Z)[1]), C = 1, maxit = 10000, tol = 1e-04, 
    sgm = c(0.8,100)) {

  if (all(P<=1) & all(P>= 0)) Z=-qnorm(P/2) else Z=P # P-values or Z scores
  
  zs=Z[,1]; zt=Z[,2]; wsum=sum(weights)
  
  pars[which(pars< 1e-8)]=1e-8
  if (sum(pars[1:3]) >= 1) pars[1:3]=pars[1:3]/(sum(pars[1:3]) + 1e-5)
  
  lhood1=function(pars) {
    pars[1]*dnorm(zs)*dnorm(zt)
  }
  
  lhood2=function(pars) {
    pars[2]*dnorm(zs,sd=pars[4])*dnorm(zt)
  }

  lhood3=function(pars) {
    pars[3]*dnorm(zs)*dnorm(zt,sd=pars[6])
  }

  lhood4=function(pars) {
    (1-pars[1]-pars[2]-pars[3])*dnorm(zs,sd=pars[5])*dnorm(zt,sd=pars[7])
  }

  lhood=function(pars) {
    #lh=matrix(0,length(zs),4)
    lh[,1]=lhood1(pars); lh[,2]=lhood2(pars); lh[,3]=lhood3(pars); lh[,4]=lhood4(pars)
    lh
    #cbind(lhood1(pars),lhood2(pars),lhood3(pars),lhood4(pars))
  }
  
  lh=matrix(0,length(zs),4)
  hist=matrix(0,maxit,8)
  diff=Inf
  
  i=1; lh0=-Inf
  while(i < maxit & diff>tol) {
  
  lp=lhood(pars)
  lp[which(lp==0)]=min(lp[which(lp>0)])
  lx=lp/rowSums(lp)
  if (any(!is.finite(lx))) lx[which(!is.finite(lx))] = 0
  p=(colSums(lx * weights) + C)/(wsum + 3 * C)

  # E step
  pars[1] = p[1]; pars[2] = p[2]; pars[3]=p[3]
  
  # M step
  #pars[4]=sqrt(sum(weights*lx[,2]* zs^2)/sum(weights*lx[,2]))
  #pars[6]=sqrt(sum(weights*lx[,3]* zt^2)/sum(weights*lx[,3]))
  #pars[5]=sqrt(sum(weights*lx[,4]* zs^2)/sum(weights*lx[,4]))
  #pars[7]=sqrt(sum(weights*lx[,4]* zt^2)/sum(weights*lx[,4]))

  lw=lx*weights; cc=colSums(lw)
  pars[4]=min(sgm[2],max(sgm[1],sqrt(sum(lw[,2]* zs^2)/cc[2])))
  pars[6]=min(sgm[2],max(sgm[1],sqrt(sum(lw[,3]* zt^2)/cc[3])))
  pars[5]=min(sgm[2],max(sgm[1],sqrt(sum(lw[,4]* zs^2)/cc[4])))
  pars[7]=min(sgm[2],max(sgm[1],sqrt(sum(lw[,4]* zt^2)/cc[4])))

    
  logl=sum(weights*log(rowSums(lp)))
  hist[i,]=c(pars,logl)
  
  i=i+1
  diff=logl-lh0
  lh0=logl
}
  hist=hist[1:(i-1),]  
  return(list(pars=pars,lhood=hist[i-1,8],hist=hist))
}



##' Run the Benjamini-Hochberg procedure
##' 
##' @title bh
##' @param P list of p-values
##' @param alpha FDR control level
##' @export
##' @return indices of p-values to reject
##' @examples 
##' # no
bh=function(P,alpha) {
  ox=rank(P); n=length(P)
  w=which(P/ox <= alpha/n)
  if (length(w)>0) return(which(ox<= max(ox[w]))) else return(c())
}





##' Estimate cFDR at a set of points using counting-points method (cFDR1 or cFDR1s)
##' 
##' @title cfdr
##' @param p vector of p-values for dependent variable of interest
##' @param q vector of p-values from other dependent variable#
##' @param sub list of indices at which to compute cFDR estimates
##' @param exclude list of indices to exclude (each point (p[i],q[i]) is still automatically included in the computation of its own cFDR value)
##' @param adj include estimate of Pr(H^p=0|Q<q) in estimate
##' @return vector of cFDR values; set to 1 if index is not in 'sub'
##' @export
##' @author James Liley
##' @examples 
##' # Generate standardised simulated dataset
##' set.seed(1); n=10000; n1p=100; n1pq=100; n1q=100
##' zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))
##'
##' cx=cfdr(p,q)
##' 
##' plot(p,q,cex=0.5,xlim=c(0,0.05)); points(p[which(cx<0.5)],q[which(cx<0.5)],col="red") 
cfdr=function(p,q,sub=which(qnorm(p/2)^2 + qnorm(q/2)^2 > 4),exclude=NULL,adj=F) {
  cf=rep(1,length(p))
  for (i in 1:length(sub)) {
    ww=which(q[-exclude]<q[sub[i]]); 
    qq=(1+length(which(p[-exclude][ww]<p[sub[i]])))/(1+length(ww))
    cf[sub[i]]=p[sub[i]]/qq
  }
  
  if (adj) {
    correction_num=1+ (ecdf(q[which(p>0.5)])(q)*length(q))
    correction_denom=(1+ rank(q))
    correct=cummin((correction_num/correction_denom)[order(-q)])[order(order(-q))]
    
    cf[which(cf<1)]=(cf*correct)[which(cf<1)]
  }
  
  cf[which(cf>1)]=1
  cf
}










##' Estimate cFDR at a set of points using parametrisation (cFDR2 or cFDR2s)
##' 
##' @title cfdrx
##' @param p vector of p-values for dependent variable of interest
##' @param q vector of p-values from other dependent variable
##' @param pars parameters governing fitted distribution of P,Q; get from function fit.4g
##' @param sub list of indices at which to compute cFDR estimates
##' @param adj include estimate of Pr(H^p=0|Q<q) in estimate
##' @return vector of cFDR values; set to 1 if index is not in 'sub'
##' @export
##' @author James Liley
##' @examples 
##' # Generate standardised simulated dataset
##' set.seed(1); n=10000; n1p=100; n1pq=100; n1q=100
##' zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))
##'
##' fit_pars=fit.4g(cbind(zp,zq))$pars
##'
##' cx=cfdrx(p,q,pars=fit_pars)
##' 
##' plot(p,q,cex=0.5,xlim=c(0,0.05)); points(p[which(cx<0.5)],q[which(cx<0.5)],col="red") 
cfdrx=function(p,q,pars,sub=1:length(p),adj=F) {
  denom1=function(zp,zq,pars) {
    pars[1]*pbivnorm(-zp,-zq) + 
      pars[2]*pbivnorm(-zp,-zq/pars[4]) +
      (1-pars[1]-pars[2])*pbivnorm(-zp/pars[3],-zq/pars[5],rho=pars[6]/(pars[3]*pars[5])) ;}
  denom2=function(zp,zq,pars) {
    pars[1]*pnorm(-zq) + pars[2]*pnorm(-zq/pars[4]) + (1-pars[1]-pars[2])*pnorm(-zq/pars[5]) ;}
  denom=function(zp,zq,pars) 2*denom1(zp,zq,pars)/denom2(zp,zq,pars)
  if (adj) {
    correction_num=1+ (ecdf(q[which(p>0.5)])(q)*length(q))
    correction_denom=(1+ rank(q))
    correct=cummin((correction_num/correction_denom)[order(-q)])[order(order(-q))]
  } else correct=rep(1,length(p))
  out=rep(1,length(p))
  out[sub]=(p*correct/denom(-qnorm(p/2),-qnorm(q/2),pars))[sub]
}



##' Estimate cFDR at a set of points using kernel density estimate (cFDR3 or cFDR3s)
##' 
##' @title cfdry
##' @param p vector of p-values for dependent variable of interest
##' @param q vector of p-values from other dependent variable
##' @param sub list of indices at which to compute cFDR estimates
##' @param exclude list of indices to exclude (each point (p[i],q[i]) is still automatically included in the computation of its own cFDR value)
##' @param adj include estimate of Pr(H^p=0|Q<q) in estimate
##' @param ... other parameters passed to kde2d
##' @return vector of cFDR values; set to 1 if index is not in 'sub'
##' @export
##' @author James Liley
##' @examples 
##' # Generate standardised simulated dataset
##' set.seed(1); n=10000; n1p=100; n1pq=100; n1q=100
##' zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
##' p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))
##'
##'
##' cx=cfdry(p,q)
##' 
##' plot(p,q,cex=0.5,xlim=c(0,0.05)); points(p[which(cx<0.5)],q[which(cx<0.5)],col="red")
cfy=function(p,q,sub=1:length(p),exclude=NULL,adj=F,...) {

  res=200
  mx=max(c(10,abs(-qnorm(p/2)))); my=max(c(10,abs(-qnorm(q/2)))); 

  
  zp=-qnorm(p/2); zq=-qnorm(q/2)
  zp[which(zp>mx)]=mx; zq[which(zq>my)]=my; w=which(is.finite(zp+zq)); 
  kpq=kde2d(zp[-exclude],zq[-exclude],n=res,lims=c(0,mx,0,my),...)
  int2_kpq=t(apply(apply(kpq$z[res:1,res:1],2,cumsum),1,cumsum))[res:1,res:1]
  int_kpq=t(outer(int2_kpq[1,],rep(1,res)))
  kgrid=kpq; kgrid$z=int2_kpq/int_kpq
  kdenom=interp.surface(kgrid,cbind(zp[sub],zq[sub]))
  
  if (adj) {
    correction_num=1+ (ecdf(q[which(p>0.5)])(q)*length(q))
    correction_denom=(1+ rank(q))
    correct=cummin((correction_num/correction_denom)[order(-q)])[order(order(-q))]
  } else correct=rep(1,length(p))
  
  out=rep(1,length(p))
  out[sub]= (p*correct)[sub]/kdenom
}





###########################################################################
## Incidental functions ###################################################
###########################################################################

px=function(x,add=F,...) if (!add) plot(-log10((1:length(x))/(1+length(x))),-log10(sort(x)),...) else points(-log10((1:length(x))/(1+length(x))),-log10(sort(x)),...) 
ab=function() abline(0,1,col="red")


