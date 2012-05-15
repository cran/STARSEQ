#include <algorithm> 
#include <iostream>
#include <ctime>
#include <vector>
#include <cmath>
#include <R.h>
#include <Rmath.h>
extern "C"{
  double kbacFunc(std::vector<std::vector<double> >& marker, std::vector<double>& y, std::vector<double>& yBinaryUb, std::vector<double>& yBinaryLb, std::vector<int> & ixPerm, double markerMax, double no1, double no0);
  double vtFunc(std::vector<std::vector<double> >& marker, std::vector<double>& y, std::vector<double>& yBinaryUb, std::vector<double>& yBinaryLb, std::vector<int> & ixPerm,double markerMax, double no1, double no0);
  double wssFunc(std::vector<std::vector<double> >& marker, std::vector<double>& y, std::vector<double>& yBinaryUb, std::vector<double>& yBinaryLb, std::vector<int> & ixPerm,double markerMax, double no1, double no0);
  double skatFunc(std::vector<std::vector<double> >& marker, std::vector<double>& y, std::vector<double>& yBinaryUb, std::vector<double>& yBinaryLb, std::vector<int> & ixPerm,double markerMax, double no1, double no0);
  double rarecoverFunc(std::vector<std::vector<double> >& marker, std::vector<double>& y, std::vector<double>& yBinaryUb, std::vector<double>& yBinaryLb, std::vector<int> & ixPerm,double markerMax, double no1, double no0);
  double skatoSumstatFunc(std::vector<double >& X_T_times_Y, std::vector<std::vector<double> >& X_T_times_X,double varY, int alternative, std::vector<double>& weight,std::vector<std::vector<double> >& lambdaMat);
  double matrixByCol(double* mat,int m,int n,int nrow,int ncol);
  double pMixChisq(double,std::vector<double> lambda);
  double covXY(double* x, double*y, std::vector<int>& ixPerm);
  double corXY(double* x, double*y, std::vector<int>& ixPerm, double& denominator)
  {
    double xMean=0.0,yMean=0.0,xyMean=0.0,noPerson=((double) ixPerm.size());
    for(int ii=0;ii<ixPerm.size();ii++)
      {
	xMean+=x[ii];
	yMean+=y[ii];
	xyMean+=(x[ii]*y[ixPerm[ii]]);
      }
    xMean/=((double) noPerson);
    yMean/=((double) noPerson);
    xyMean/=((double) noPerson);
    return (xyMean-xMean*yMean)/denominator;
  }
  
  double covXY(double* x, double*y, std::vector<int>& ixPerm)
  {
    double xMean=0.0,yMean=0.0,xyMean=0.0,noPerson=((double) ixPerm.size());
    for(int ii=0;ii<ixPerm.size();ii++)
      {
	xMean+=x[ii];
	yMean+=y[ii];
	xyMean+=(x[ii]*y[ixPerm[ii]]);
      }
    xMean/=((double) noPerson);
    yMean/=((double) noPerson);
    xyMean/=((double) noPerson);
    return (xyMean-xMean*yMean);
  }
  double vectorCor(std::vector<double>& x, std::vector<double>& y,std::vector<int>& ixPerm)
  {
    double xMean=0.0,yMean=0.0,xyMean=0.0,xxMean=0.0,yyMean=0.0,noPerson=((double) ixPerm.size());
    for(int ii=0;ii<ixPerm.size();ii++)
      {
	xMean+=x[ii];
	xxMean+=(x[ii])*(x[ii]);
	yyMean+=(y[ii])*(y[ii]);
	yMean+=y[ii];
	xyMean+=(x[ii]*y[ixPerm[ii]]);
      }
    xMean/=noPerson;
    yMean/=noPerson;
    xxMean/=noPerson;
    yyMean/=noPerson;
    xyMean/=noPerson;    
    return (xyMean-xMean*yMean)/sqrt((xxMean-xMean*xMean)*(yyMean-yMean*yMean));
  }
  
  double vectorMax(std::vector<double>& x)
  {
    return(*max_element(x.begin(),x.end()));
  }

  void genericPerm1(double* xPt, double* yPt,int* noPersonPt,int* noPermPt,double* pValuePt,double* statisticPt)
  {  
    std::vector<double> statPerm(noPermPt[0]);  
    std::vector<int> ixPerm(noPersonPt[0]);
    double denominator=0.0,xMean=0.0,yMean=0.0,xSumSq=0.0,ySumSq=0.0,pValue=0.0;
    int ii=0;
    for(ii=0;ii<noPersonPt[0];ii++) {
      ixPerm[ii]=ii;
    }
    for(ii=0;ii<noPersonPt[0];ii++)
      {
	xMean+=xPt[ii];
	yMean+=yPt[ii];
	xSumSq+=(xPt[ii])*(xPt[ii]);
	ySumSq+=(yPt[ii])*(yPt[ii]);  
      }
    xMean/=((double) noPersonPt[0]);
    yMean/=((double) noPersonPt[0]);
    xSumSq/=((double) noPersonPt[0]);
    ySumSq/=((double) noPersonPt[0]);

    denominator=sqrt((xSumSq-xMean*xMean)*(ySumSq-yMean*yMean));
    double statData=corXY(xPt,yPt,ixPerm,denominator);

    statData*= ((statData)*noPersonPt[0]);

    for(ii=0;ii<noPermPt[0];ii++)
      {
	random_shuffle(ixPerm.begin(),ixPerm.end());
	statPerm[ii]=corXY(xPt,yPt,ixPerm,denominator);
	statPerm[ii] *= (statPerm[ii]*(noPersonPt[0]));
	pValue += (statData<(statPerm[ii]) ? 1.0 : 0.0);  
      }

    pValuePt[0]=pValue/((double) noPermPt[0]);  
    statisticPt[0]=statData;
  }
  
  void genericPerm2(double* markerPt,double*yPt,double* yUbPt, double* yLbPt, int* noPersonPt, int* noMarkerPt, int* noPermPt, double* pValuePt, int* testPt, double* markerMaxPt,double* statisticPt,double* extraPar)
  {
    double alpha=extraPar[0];
    int ii=0,jj=0,kk=0,ll=0;
    double pValue=0.0,pValueTmp=0.0;
    std::vector<std::vector<double> > marker(noMarkerPt[0],std::vector<double>(noPersonPt[0]));
    std::vector<double> y(noPersonPt[0]);
    std::vector<int> ixPerm(noPersonPt[0]);
    //noPerm must be a multiple of 1000; need to check in R code;
    int jjMax=((int) (noPermPt[0])/1000);
    if(jjMax==0) jjMax=1;
    double pValueMin=0.0,pValueMax=0.0;
    std::vector<double> statPerm(jjMax*1000,0.0);
    double yUb=yUbPt[0],yLb=yLbPt[0];
    std::vector<double> yBinaryUb(noPersonPt[0]), yBinaryLb(noPersonPt[0]);
    double markerMax=markerMaxPt[0];
    for(jj=0;jj< noPersonPt[0];jj++)
      {
	y[jj]=yPt[jj];
	yBinaryUb[jj] = ((yPt[jj] > yUb) ? 1.0 : 0.0);
	yBinaryLb[jj] = ((yPt[jj] < yLb) ? 1.0 : 0.0);
	for(ii=0;ii< noMarkerPt[0];ii++)
	  {
	    marker[ii][jj]=markerPt[ii*noPersonPt[0]+jj];
	  }
      }
    double noUb1=0.0,noLb1=0.0;
    for(ii=0;ii<noPersonPt[0];ii++) {
      ixPerm[ii]=ii;
      noUb1+=yBinaryUb[ii];
      noLb1+=yBinaryLb[ii];
    }
    if(testPt[0]!=6)
      {
	double (*statFunc)(std::vector<std::vector<double> >& , std::vector<double>& , std::vector<double>& , std::vector<double>&, std::vector<int>&, double , double , double );
	
	if(testPt[0]==1)
	  statFunc=&kbacFunc;
	if(testPt[0]==2)
	  statFunc=&vtFunc;
	if(testPt[0]==3)
	  statFunc=&wssFunc;
	if(testPt[0]==4)
	  statFunc=&skatFunc;
	if(testPt[0]==5)
	  statFunc=&rarecoverFunc;
	double statData=(*statFunc)(marker,y,yBinaryUb,yBinaryLb, ixPerm, markerMax,noUb1,noLb1);
	statisticPt[0]=statData;
	double statPermTmp=0.0;
	for(jj=0;jj<jjMax;jj++)
	  {
	    for(ii=0;ii<1000;ii++)
	      {
		random_shuffle(ixPerm.begin(),ixPerm.end());
		statPermTmp=(*statFunc) (marker,y,yBinaryUb,yBinaryLb, ixPerm,markerMax,noUb1,noLb1);
		statisticPt[ii+jj*1000+1]=statPermTmp;
		//pValue += (statData<(statPerm[ii]) ? 1.0 : 0.0);
		pValue += (statData<=(statPermTmp) ? 1.0 : 0.0);
	      }	
	    pValueTmp=pValue/((double) (jj+1)*1000);
	    pValueMin=pValueTmp-1.96*sqrt(pValueTmp*(1-pValueTmp)/((double) (jj+1)*1000));
	    pValueMax=pValueTmp+1.96*sqrt(pValueTmp*(1-pValueTmp)/((double) (jj+1)*1000));
	    //it could be that the p-value is 0 and the permutation should be continued;
	    if(pValueMin>alpha) jj=jjMax;//if the p-value is too large terminate the permutation procedure;
	  }
	pValuePt[0]=pValueTmp;	
      }
    if(testPt[0]==6)
      {
	std::vector<std::vector<double> > X_T_times_X(noMarkerPt[0],std::vector<double> (noMarkerPt[0],0.0));
	std::vector<double>  X_T_times_Y(noMarkerPt[0],0.0);
	double yMean=0.0;
	for(ii=0;ii<y.size();ii++) yMean+=y[ii];
	yMean=yMean/((double) y.size());
	for(ii=0;ii<y.size();ii++) y[ii]-=yMean;
	//centralize marker;
	std::vector<double> markerMean(noMarkerPt[0],0.0);
	
	for(ii=0;ii<noMarkerPt[0];ii++)
	  {  
	    for(jj=0;jj<y.size();jj++)
	      {
		markerMean[ii]+=(marker[ii][jj]);
	      }
	    markerMean[ii]/=((double) y.size());
	  }

	for(ii=0;ii<noMarkerPt[0];ii++)
	  {  
	    for(jj=0;jj<y.size();jj++)
	      {
		marker[ii][jj]-=markerMean[ii];
	      }
	  }
	for(ii=0;ii<noMarkerPt[0];ii++)
	  {  
	    for(jj=0;jj<noMarkerPt[0];jj++)	      
	      {
		for(kk=0;kk<noPersonPt[0];kk++)
		  {
		    X_T_times_X[ii][jj]+=((marker[ii][kk])*(marker[jj][kk]));
		  }
	      }
	  }
	for(jj=0;jj<noPersonPt[0];jj++)
	  {  
	    for(ii=0;ii<noMarkerPt[0];ii++)
	      {
		X_T_times_Y[ii]+=(marker[ii][jj]*y[jj]);
	      }
	  }
	std::vector<double> weight(noMarkerPt[0],0.0);
	std::vector<std::vector<double> > lambdaMat(11,std::vector<double>(noMarkerPt[0],0.0));	
	for(ii=0;ii<noMarkerPt[0];ii++)
	  {
	    weight[ii]=extraPar[ii+1];
	  }
	for(ii=0;ii<11;ii++)
	  {
	    for(jj=0;jj<noMarkerPt[0];jj++)
	      {
		lambdaMat[ii][jj]=matrixByCol(&(extraPar[noMarkerPt[0]+1]),ii,jj,11,noMarkerPt[0]);
	      }
	  }
	for(ii=0;ii<noMarkerPt[0];ii++)
	  {
	    extraPar[noMarkerPt[0]+1+11*(noMarkerPt[0])+ii]=lambdaMat[1][ii];
	  }
	double varY=covXY(&(y[0]),&(y[0]),ixPerm);
	double statData=skatoSumstatFunc(X_T_times_Y,X_T_times_X,varY,0,weight,lambdaMat);
	statisticPt[0]=statData;	
	for(jj=0;jj<jjMax;jj++)
	  {
	    for(ii=0;ii<1000;ii++)
	      {
		random_shuffle(ixPerm.begin(),ixPerm.end());
		std::vector<double> X_T_times_YPerm(noMarkerPt[0],0.0);
		for(kk=0;kk<noPersonPt[0];kk++)
		  {  
		    for(ll=0;ll<noMarkerPt[0];ll++)
		      {
			X_T_times_YPerm[ll]+=(marker[ll][kk]*y[ixPerm[kk]]);
		      }
		  }
		statPerm[ii+jj*1000]=skatoSumstatFunc(X_T_times_YPerm,X_T_times_X,varY,0,weight,lambdaMat);
		statisticPt[ii+jj*1000+1]=statPerm[ii+jj*1000];
		pValue += (statData>=(statPerm[ii+jj*1000]) ? 1.0 : 0.0);
	      }	
	    pValueTmp=pValue/((double) (jj+1)*1000);
	    pValueMin=pValueTmp-1.96*sqrt(pValueTmp*(1-pValueTmp)/((double) (jj+1)*1000));
	    pValueMax=pValueTmp+1.96*sqrt(pValueTmp*(1-pValueTmp)/((double) (jj+1)*1000));
	    if(pValueMin>alpha) jj=jjMax;//if the p-value is too large terminate the permutation procedure;
	  }
	pValuePt[0]=pValueTmp;		
      }
  }
  double matrixByCol(double* mat,int m,int n,int nrow,int ncol)
  {
    if(m <nrow & n<ncol)
      {
	return mat[m+n*nrow];
      }
    return(0.0);
  }

  double kbacFunc(std::vector<std::vector<double> >& marker, std::vector<double>& y, std::vector<double>& yBinaryUb, std::vector<double>& yBinaryLb, std::vector<int> & ixPerm, double markerMax, double noUb1, double noLb1)
  {
    std::vector<int> markerCountUb1(markerMax);
    std::vector<int> markerCountLb1(markerMax);
    std::vector<int> markerCountAll(markerMax);
    std::vector<double> weight1(markerMax,0), weight0(markerMax,0);
    int noPerson=y.size();
    std::vector<double> x1(noPerson),x0(noPerson);
    for(int ii=0;ii<markerMax;ii++)
      {
	for(int jj=0;jj<noPerson;jj++)      
	  {
	    if((marker[0][jj])==ii && yBinaryUb[ixPerm[jj] ]==1.0) markerCountUb1[ii]+=1;
	    if((marker[0][jj])==ii && yBinaryLb[ixPerm[jj] ]==1.0) markerCountLb1[ii]+=1;
	    if((marker[0][jj])==ii) markerCountAll[ii]+=1;
	  }
      }  
    double noUb0,noLb0;
    noUb0=noPerson-noUb1;
    noLb0=noPerson-noLb1;

    for(int ii=0;ii<markerMax;ii++)
      {
	weight1[ii]=phyper((double) markerCountUb1[ii],(double) noUb1, (double) noUb0,(double) markerCountAll[ii],1,0);
	weight0[ii]=phyper((double) markerCountLb1[ii],(double) noLb1,(double) noLb0,(double) markerCountAll[ii],1,0);
      }

    weight1[0]=0;weight0[0]=0;
    //assign genotype coding;
    for(int ii=0;ii<noPerson;ii++)
      {
	x1[ii]=weight1[(int) marker[0][ii]];
	x0[ii]=weight0[(int) marker[0][ii]];
      }
    double stat1,stat0;
    stat1=vectorCor(x1,y,ixPerm);   
    stat0=vectorCor(x0,y,ixPerm);
    stat1*=(stat1*noPerson);
    stat0*=(stat0*noPerson);
    return stat1>stat0 ? stat1 : stat0;
  }
  
  double vtFunc(std::vector<std::vector<double> >& marker, std::vector<double>& y, std::vector<double>& yBinaryUb, std::vector<double>& yBinaryLb, std::vector<int> & ixPerm,double markerMax, double noUb1, double noLb1)
  {
    int noMarker=marker.size();
    std::vector<double> stat_vec(noMarker);
    double tmp;
    double noPerson=y.size();
    //the calculation of the likelihood has to be 
    for(int ii=0;ii<noMarker;ii++)
      {
	tmp=vectorCor(marker[ii],y,ixPerm);		
	stat_vec[ii]=tmp*tmp;
      }
    double vt_stat=vectorMax(stat_vec);
    return vt_stat*noPerson;
  }


  double skatFunc(std::vector<std::vector<double> >& marker, std::vector<double>& y, std::vector<double>& yBinaryUb, std::vector<double>& yBinaryLb, std::vector<int> & ixPerm,double markerMax, double noUb1, double noLb1)
  {
    //marker here is a weighted marker which is given by GWG^T
    double skatStat=0.0;
    for(int ii=0;ii<marker.size();ii++)
      {
	for(int jj=0;jj<marker.size();jj++)
	  {
	    skatStat+=(y[ixPerm[ii]])*(marker[ii][jj])*(y[ixPerm[jj]]);
	  }
      }
    return skatStat;
  }

  double wssFunc(std::vector<std::vector<double> >& marker, std::vector<double>& y, std::vector<double>& yBinaryUb, std::vector<double>& yBinaryLb, std::vector<int> & ixPerm,double markerMax, double noUb1, double noLb1)
  {
    int noMarker=marker.size();
    std::vector<int> markerCountUb1(noMarker);
    std::vector<int> markerCountLb1(noMarker);
    std::vector<int> markerCountAll(noMarker);
    std::vector<double> weight1(noMarker,0), weight0(noMarker,0);
    int noPerson=y.size();
    std::vector<double> x1(noPerson),x0(noPerson);
    for(int ii=0;ii<noMarker;ii++)
      {
	for(int jj=0;jj<noPerson;jj++)      
	  {
	    if(yBinaryUb[ixPerm[jj] ]==0.0) markerCountUb1[ii]+=(marker[ii][jj]);
	    if(yBinaryLb[ixPerm[jj] ]==0.0) markerCountLb1[ii]+=(marker[ii][jj]);
	  }
      }  
    double noUb0,noLb0;
    noUb0=noPerson-noUb1;
    noLb0=noPerson-noLb1;
    
    for(int ii=0;ii<noMarker;ii++)
      {
	weight1[ii]=1.0/sqrt(noPerson*(markerCountUb1[ii]+1)/(2*noUb0+2)*(1-(markerCountUb1[ii]+1)/(2*noUb0+2)));
	weight0[ii]=1.0/sqrt(noPerson*(markerCountLb1[ii]+1)/(2*noLb0+2)*(1-(markerCountLb1[ii]+1)/(2*noLb0+2)));
      }
    
    // weight1[0]=0;weight0[0]=0;
    //assign genotype coding;
    for(int ii=0;ii<noMarker;ii++)
      {
	for(int jj=0;jj<noPerson;jj++)
	  {
	    x1[jj]+=marker[ii][jj]*weight1[ii];
	    x0[jj]+=marker[ii][jj]*weight0[ii];
	  }
      }
    double stat1,stat0;
    stat1=vectorCor(x1,y,ixPerm);   
    stat0=vectorCor(x0,y,ixPerm);
    stat1*=(stat1*noPerson);
    stat0*=(stat0*noPerson);
    return stat1>stat0 ? stat1 : stat0;
  } 

  double rarecoverFunc(std::vector<std::vector<double> >& marker, std::vector<double>& y, std::vector<double>& yBinaryUb, std::vector<double>& yBinaryLb, std::vector<int> & ixPerm,double markerMax, double noUb1, double noLb1)
  {
    int cont=1;
    double Q=0.0;
    std::vector<int> ixIn(marker.size(),0);
    std::vector<int> ixOut(marker.size(),1);
    double statOld=0.0,statNewTemp=0.0;
    std::vector<double> xii(y.size(),0.0),xiiTemp(y.size(),0.0);
    int noOut=marker.size(), noIn=0,ixMax=0;
    std::vector<double> statVec(marker.size(),0.0);
    //for(int ii=0;ii<ixPerm.size();ii++) ixPerm[ii]=ii;
    while(noOut>0 & cont==1)
      {
	for(int ii=0;ii<ixOut.size();ii++)
	  {
	    if(ixOut[ii]==1)
	      {
		for(int jj=0;jj<y.size();jj++)
		  xiiTemp[jj]=xii[jj]+(marker[ii][jj]);
		statVec[ii]=vectorCor(xiiTemp,y,ixPerm);
		statVec[ii]*=(y.size()*(statVec[ii]));
	      }	    
	  }
	ixMax=max_element(statVec.begin(),statVec.end())-statVec.begin();
	statNewTemp=*max_element(statVec.begin(),statVec.end());
	if(statNewTemp-statOld>=Q)
	  {
	    ixIn[ixMax]=1;
	    ixOut[ixMax]=0;
	    statOld=statNewTemp;
	    for(int ii=0;ii<xii.size();ii++)
	      xii[ii]+=marker[ixMax][ii];
	    cont=1;
	    noOut--;
	    for(int ii=0;ii<statVec.size();ii++) 
	      {
		statVec[ii]=0.0;
	      }
	  }
	if(statNewTemp-statOld<Q)
	  cont=0;
      }
    return statOld;
  }
  double skatoSumstatFunc(std::vector<double >& X_T_times_Y, std::vector<std::vector<double> >& X_T_times_X,double varY, int alternative, std::vector<double>& weight,std::vector<std::vector<double> >& lambdaMat)
  {
    int ii=0,jj=0,kk=0,noMarker=X_T_times_Y.size();
    std::vector<double> rho(11,0.0),stat(11,0.0),pValue(11,0.0);    
    for(ii=0;ii<11;ii++) rho[ii]=0.1*ii;
    for(ii=0;ii<11;ii++)
      {
	for(jj=0;jj<noMarker;jj++)
	  {
	    for(kk=0;kk<noMarker;kk++)
	      {
		stat[ii]+=(1-rho[ii])*(X_T_times_Y[jj])*(X_T_times_Y[kk]);
		if(jj==kk) stat[ii]+=(rho[ii])*(X_T_times_Y[jj])*weight[kk]*(X_T_times_Y[kk]);
	      }
	  }
	pValue[ii]=pMixChisq(stat[ii],lambdaMat[ii]);
      }
    return (*min_element(pValue.begin(),pValue.end()));
  }
  
}

