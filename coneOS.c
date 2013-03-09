#include "coneOS.h"

Sol * coneOS(Data * d, Cone * k)
{
	int i;
	double err = -1, EPS_PRI = -1;  
	Work * w = initWork(d);
  for (i=0;i<d->MAX_ITERS;i++){             
		projectLinSys(d,w);
		relax(d,w);
		projectCones(d,w,k);
		updateDualVars(w);

    // err = calcPriResid(d,w);
    // EPS_PRI = sqrt(w->l)*d->EPS_ABS + 
    //       d->EPS_REL*fmax(calcNorm(w->uv,w->l),calcNorm(w->uv_t,w->l));
    // if (err < EPS_PRI) break;
		if (i % 10 == 0) printSummary(d,w,i,err,EPS_PRI);
	}
	Sol * sol = malloc(sizeof(Sol));
  setx(d,w,sol);
  sety(d,w,sol);
	//getSolution(d,w,sol);
	printSummary(d,w,i,err,EPS_PRI);
	printSol(d,sol);
	freeWork(w);
	return sol;
}

void freeWork(Work * w){
  freePriv(w);
  free(w->z_half);
  free(w->z);
  free(w->u);
  free(w->ztmp);
  free(w);
}

void printSol(Data * d, Sol * sol){
	int i;
	printf("%s\n",sol->status); 
	if (0 && sol->x != NULL){
		for ( i=0;i<d->n; i++){
			printf("x[%i] = %4f\n",i, sol->x[i]);
		}
	}
	if (0 && sol->y != NULL){
		for ( i=0;i<d->m; i++){
			printf("y[%i] = %4f\n",i, sol->y[i]);
		}
	}
}

void updateDualVars(Work * w){
	addScaledArray(w->u,w->ztmp,w->l,1);
	addScaledArray(w->u,w->z,w->l,-1);
}

void projectCones(Data *d,Work * w,Cone * k){
	memcpy(w->z,w->ztmp,w->l*sizeof(double));
	addScaledArray(w->z,w->u,w->l,1);
	/* s onto K */
	projCone(w->z + w->si,k);

	/* r onto 0 */
	memset(w->z + w->ri,0,(sizeof(double)*d->n));
  
  /* y onto K^* */
	projDualCone(w->z + w->yi,k);
}

// void getSolution(Data * d,Work * w,Sol * sol){
//   double tau = (w->uv[d->n+d->m]+w->uv_t[d->n+d->m])/2;
//   double kap = (w->uv[w->l-1]+w->uv_t[w->l-1])/2;
//   setx(d,w,sol);
//   setz(d,w,sol);
//   if (tau > d->UNDET_TOL && tau > kap){
//     sol->status = strdup("Solved");
//     scaleArray(sol->x,1/tau,d->n);
//     scaleArray(sol->z,1/tau,d->m);
//   }
//   else{ 
//     if (calcNorm(w->uv,w->l)<d->UNDET_TOL*sqrt(w->l)){
//       sol->status = strdup("Undetermined");
//       scaleArray(sol->x,NAN,d->n);
//       scaleArray(sol->z,NAN,d->m);
//     }
//     else {
//       double ip_z = innerProd(d->b,sol->z,d->m);  
//       double ip_x = innerProd(d->c,sol->x,d->n);  
//       if (ip_z < ip_x){
//         sol->status = strdup("Infeasible");
//         scaleArray(sol->z,-1/ip_z,d->m);
//         scaleArray(sol->x,NAN,d->n);
//       }
//       else{
//         sol->status = strdup("Unbounded");
//         scaleArray(sol->x,-1/ip_x,d->n);
//         scaleArray(sol->z,NAN,d->m);
//       }
//     }
//   }
// }

void sety(Data * d,Work * w, Sol * sol){
	sol->y = malloc(sizeof(double)*d->m);
	memcpy(sol->y, w->z + 2*(d->n) + d->m, d->m*sizeof(double));
}

void setx(Data * d,Work * w, Sol * sol){
	sol->x = malloc(sizeof(double)*d->n);
	memcpy(sol->x, w->z, d->n*sizeof(double));
}

void relax(Data * d,Work * w){   
	int j;
	for( j=0;j<w->l;j++){
		w->ztmp[j] = d->ALPH*w->z_half[j]+(1-d->ALPH)*w->z[j];
	}
}

void printSummary(Data * d,Work * w,int i, double err, double EPS_PRI){
	printf("Iteration %i, primal residual %4f, primal tolerance %4f\n",i,err,EPS_PRI);
}
