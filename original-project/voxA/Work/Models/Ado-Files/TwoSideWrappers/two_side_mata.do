

mata:

mata clear
mata set matastrict on

real scalar bo_k(real scalar L, real scalar U, real scalar y, real scalar d)
{
	real scalar k
	k=0
	if (d<=(y-U)){	
		k=U
	}
	else if (((y-U)<d)&(d<(y-L))){	
		k=(y-d)
	}
	else {
		k=L
	}
	return(k)
}

real scalar bo_dk(real scalar L, real scalar U, real scalar y, real scalar d)
{
	real scalar k
	k=0
	if (d<=(y-U)){	
		k=0
	}
	else if (((y-U)<d)&(d<(y-L))){	
		k=-1
	}
	else {
		k=0
	}
	return(k)
}

real scalar bo_K(real scalar L, real scalar U, real scalar y, real scalar d)
{
	real scalar k
	
	k=0
	if (d<=(y-U)){
		k=2*y*U-2*d*U-U*U
	}
	else if (((y-U)<d)&(d<(y-L))){
		k=(y-d)^2
	}
	else {
		k=2*y*L-2*d*L-L^2
	}
	return(k)
}

real scalar bo_S(real scalar L1, real scalar L2, real scalar U1, real scalar U2, real scalar y1, real scalar y2, real scalar d)
{

	return(bo_K(L2,U2,y1,d)+bo_K(L1,U1,y2,-d)-d^2)
}

real scalar bo_V(real scalar L1, real scalar L2, real scalar U1, real scalar U2, real scalar y1, real scalar y2, real scalar d)
{
	real scalar v
	v=0
	
	 
		if (d<=(L1-U2)){
			v=bo_S(L1,L2,U1,U2,y1,y2,L1-U2)
		}
		else if (((L1-U2)<d)&(d<(U1-L2))){
			v=bo_S(L1,L2,U1,U2,y1,y2,d)
		}
		else {
			v=bo_S(L1,L2,U1,U2,y1,y2,U1-L2)
		}
	 
	return(v)
	
}

real scalar bo_u(real scalar L1, real scalar L2, real scalar U1, real scalar U2, real scalar y1, real scalar y2, real scalar d)
{
	real scalar u, a1,a2,a3
	u=0
	if (((L1-U2)<d)&(d<(U1-L2))){
		u=bo_k(L2,U2,y1,d)-bo_k(L1,U1,y2,-d)+d;
	}
	return(u)
}

real scalar bo_du(real scalar L1, real scalar L2, real scalar U1, real scalar U2, real scalar y1, real scalar y2, real scalar d)
{
	real scalar u, a1,a2,a3
	u=0
	if (((L1-U2)<d)&(d<(U1-L2))){
		u=bo_dk(L2,U2,y1,d)+bo_dk(L1,U1,y2,-d)+1;
	}
	return(u)
}

void bo_fct2side_new(real vector todo,
			real vector bet,
			real vector y,
			real matrix x,
			real vector Li,
			real vector Ui,
			real matrix nn,
			real scalar npers,
			real vector wei,
			real matrix scores,
			real vector nused,
			real scalar f,
			real vector g,
			real matrix H
)
{	
	real vector fi, xb, dx
	real scalar i, j1, j2, de, k
	real matrix di
	

	fi=J(npers,1,0)
	k=cols(bet)
	H=J(k,k,0)
	di=J(npers,k,0)
	nused=0\0

	xb=x*bet'
	for (i=1; i<=npers; i++) {
	if (wei[i]>0) {
		nused[1]=nused[1]+1
		for (j1=nn[i,1]; j1<nn[i,2]; j1++) {
		for (j2=j1+1; j2<=nn[i,2];  j2++) {
		de=xb[j1]-xb[j2]

		
		fi[i]=fi[i]+bo_V(Li[j1], Li[j2], Ui[j1], Ui[j2], y[j1], y[j2], de)

		if (todo >=1){
			dx=x[j1,.]-x[j2,.]
			di[i,.]=di[i,.]+bo_u(Li[j1], Li[j2], Ui[j1], Ui[j2], y[j1], y[j2], de)*dx

		}
		
		if (todo >=2){
			H=H+bo_du(Li[j1], Li[j2], Ui[j1], Ui[j2], y[j1], y[j2], de)*dx'*dx/wei[i]
		}
				}
		}
	fi[i]=-0.5*fi[i]/wei[i]
	di[i,.]=di[i,.]/wei[i]
	}
	}
	f=colsum(fi)
	g=colsum(di)
	scores=di
}

void bo_fct2side(real vector todo,
			real vector bet,
			real vector y,
			real matrix x,
			real matrix nn,
			real scalar npers,
			real vector wei,
			real matrix scores,
			real vector nused,
			real scalar f,
			real vector g,
			real matrix H
)
{	
	real vector fi, xb, dx
	real scalar i, j1, j2, de, k
	real matrix di
	

	fi=J(npers,1,0)
	k=cols(bet)
	H=J(k,k,0)
	di=J(npers,k,0)
	nused=0\0

	xb=x*bet'
	for (i=1; i<=npers; i++) {
	if (wei[i]>0) {
		nused[1]=nused[1]+1
		for (j1=nn[i,1]; j1<nn[i,2]; j1++) {
		for (j2=j1+1; j2<=nn[i,2];  j2++) {
		de=xb[j1]-xb[j2]
		
		fi[i]=fi[i]+bo_rrr(y[j1],de)+bo_rrr(y[j2],-de)
			
		if (todo >=1){
			dx=x[j1,.]-x[j2,.]
			di[i,.]=di[i,.]+(bo_uuu(y[j1],de)-bo_uuu(y[j2],-de))*dx

		}
		if (todo >=2){
			H=H+(bo_ddd(y[j1],de)+bo_ddd(y[j2],-de))*dx'*dx/wei[i]
		}
				}
		}
	fi[i]=fi[i]/wei[i]
	di[i,.]=di[i,.]/wei[i]
	}
	}
	f=colsum(fi)
	g=colsum(di)
	scores=di
}

void two_sidemata(string scalar varname, string scalar touse, real scalar estmet, real scalar pd, real scalar bs, real scalar vc)
{
	real matrix V, bodat, x, xx, nn, bboot, H, scores, gamma 
	real vector b, y, ident, idsort, wei,h,p,pp,g1,g2, g, Lit,Uit
	real scalar n1, n2, n3, n4, n5, k, i, j, fval, chi2stat, chi2prob  
	real scalar nused, fcen1,fcen2,ymin,ymax, simdel,f, f1,f2, errorm, converged, nboor, improved
	real scalar todo, delt, L1, U1, d1, y1, L2, U2,y2,df
	transmorphic S

	

	
	
	
	
	ymin=0
	ymax=1
		
	st_view(bodat=.,.,tokens(varname),touse)
	bodat=sort(bodat,cols(bodat))

	y=bodat[.,1]
	if (vc==1){
		Lit=bodat[.,(cols(bodat)-2)]
		Uit=bodat[.,(cols(bodat)-1)]
		x=bodat[.,2::(cols(bodat)-3)]
	}
	else {
		Lit=J(rows(bodat),1,ymin)
		Uit=J(rows(bodat),1,ymax)
		x=bodat[.,2::(cols(bodat)-1)]
	}
	ident=bodat[.,cols(bodat)]
	idsort=uniqrows(ident)

	k=cols(x)
	n1=rows(x)
	n2=rows(idsort)
	fcen1=mean(y:<=J(n1,1,ymin))
	fcen2=mean(y:>=J(n1,1,ymax))
	nn=1,0
	i=2
	j=1
	for (i=2; i<=n1; i++) {
		if (ident[i] != ident[i-1]) {
			nn[j,2]=i-1
			nn = nn \(0,0)
			j=j+1
			nn[j,1]=i
		}
	} 
	nn[j,2]=n1
	wei=nn[.,2]-nn[.,1]

	if (rows(nn)!=rows(idsort)) {
	"something went wrong with sorting the data"
	return
	}
	

	xx=x,J(rows(x),1,1)
	b=invsym(xx'*xx)*xx'y
	b=b[1::k,.]'
	


	fval=0
	H=0
	g=0
	todo=0
	scores=0
	nused=0\0

	S=optimize_init()

	if (pd==0) {
		optimize_init_verbose(S,0)
		optimize_init_tracelevel(S, "none")
	}

	converged=0
	while (converged==0) {

		optimize_init_evaluator(S, &bo_fct2side_new())
		optimize_init_technique(S, "nr")
		optimize_init_evaluatortype(S, "d2")
		optimize_init_params(S,b)
		optimize_init_argument(S, 1, y)
		optimize_init_argument(S, 2, x)
		optimize_init_argument(S, 3, Lit)
		optimize_init_argument(S, 4, Uit)
		optimize_init_argument(S, 5, nn)
		optimize_init_argument(S, 6, n2)
		optimize_init_argument(S, 7, wei)
		optimize_init_argument(S, 8, scores)
		optimize_init_argument(S, 9, nused)
		optimize_init_which(S,"max")
		p=optimize(S)

		converged=1
		todo=0
		bo_fct2side_new(todo,p,y,x, Lit,Uit,nn,n2,wei,scores,nused,fval,g,H)
		b=p
		h=J(k,1,0)
		h=h:+0.0001
		g1=0
		f1=0
		
		for  (i=1; i<=k; i++) {
			improved=1
			while (improved==1) {
				improved=0
				pp=b
				pp[i]=pp[i]+h[i]
				bo_fct2side_new(todo,pp,y,x, Lit,Uit,nn,n2,wei,scores,nused,f1,g1,H)
				if (f1>fval) {
"Function improved" 
f1
					improved=1
					fval=f1
					b=pp
					converged=0
				}
			}
			improved=1
			while (improved==1) {
				improved=0
				pp=b
				pp[i]=pp[i]-h[i]
				bo_fct2side_new(todo,pp,y,x,Lit,Uit,nn,n2,wei,scores,nused,f1,g1,H)
				if (f1>fval) {
"Function improved" 
f1
					improved=1
					fval=f1
					b=pp
					converged=0
				}
			}
		}

	}


	
/* calculating Variance */

	n3=n1/n2
	n4=colmin(wei)+1
	n5=colmax(wei)+1
	todo=2
	bo_fct2side_new(todo,p,y,x, Lit,Uit,nn,n2,wei,scores,nused,fval,g,H)
	gamma=scores'*scores


	H=invsym(-H)
	V=H*gamma*H'



	chi2stat=p*invsym(V)*p'
	chi2prob=chi2tail(cols(p),chi2stat)

	st_matrix("r(beta)", p)
	st_matrix("r(V)", V)
	st_matrix("r(N)", (n1\n2\n3\n4\n5\nused[1]\fcen1\fcen2))
	st_matrix("r(chi2)", (chi2stat\chi2prob))

}

real scalar function bo_ddd(real scalar y, real scalar d)
/* 
Derivative of bo_uuu
*/
{
	real scalar dd

	dd=0
	if (d>y) dd=0
	else if (d>0) dd=-1
	else if (d>y-1) dd=0
	else if (d>-1) dd=1
	else dd=0
	
	return(dd)
}

real scalar function bo_uuu(real scalar y, real scalar d)
/* 
Note that 
		u1(y,d)=u(y,d). 
		u2(y,d)=u(y,-d)
*/
{
	real scalar uu

	if ((d>1)|(d<-1)){
		uu=0
		}
	else{
		if (d>0) uu=max(((y-d)\0))
		else uu=min((y\(1+d)))
		}
	return(uu)
}

real scalar function bo_rrr(real scalar y, real scalar dd)
/* 
Note that 
		r1(y,d)=r(y,d). 
		r2(y,d)=-r(y,-d)
*/
{
	real scalar rr, d

	d=dd
	if (d>1) d=1
	if (d<-1) d=-1

	if (d<y-1) rr=d+(d^2+(1-y)^2)/2
	else if (d<0) rr=d*y
	else if (d<y) rr=d*y-(d^2)/2
	else rr=(y^2)/2
	
	return(rr)
}

end
