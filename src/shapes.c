#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SDF.h"


double circleSDF(struct vec3 *v, struct vec3 *p){
        return genSDFsphere(newVec(0.0, 0.0, -1.0), 1.2, v, p);
}


struct vec3 *circleSDFGrad(struct vec3 *p){
        return normalize(vecsubtract(p, newVec(0.0, 1.2, -1.0)));
}

double stemSDF(struct vec3 *v, struct vec3 *p){
	struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(-1.0, 0.0, 0.25);
        double dims[3] = {0.07, 0.07, 0.25};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        double result = genSDFbox(norms, corn, dims, v, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

struct vec3 *stemgrad(struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(-1.0, 0.0, 0.25);
        double dims[3] = {0.07, 0.7, 0.25};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
	struct vec3 *result = boxGrad(norms, corn, dims, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;	
}


double planeSDF(struct vec3 *v, struct vec3 *p){
        struct vec3 *pointon = newVec(0, 0, -6);
        struct vec3 *n = newVec(0.0, 0.0, 1.0);
        struct vec3 *to = vecsubtract(p, pointon);
        if (dot3(v,n) == 0) return 1111.1;
        double result = -1*dot3(to, n)/dot3(v, n);
        if (result >= 0) return result;
        else return 1111.1;
}

struct vec3 *planeSDFGrad(struct vec3 *p){
        return newVec(0,0,1);
}


double leg1SDF(struct vec3 *v, struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(0.0, 0.0, -2.0);
        double dims[3] = {0.3, 0.3, 4};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        double result = genSDFbox(norms, corn, dims, v, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

struct vec3 *leg1SDFgrad(struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(0.0, 0.0, -2.0);
        double dims[3] = {0.3, 0.3, 4};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        struct vec3 *result = boxGrad(norms, corn, dims, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}


double leg2SDF(struct vec3 *v, struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(-2.5284, -2.5284, -2.0);
        double dims[3] = {0.3, 0.3, 4};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        double result = genSDFbox(norms, corn, dims, v, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

struct vec3 *leg2SDFgrad(struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
	struct vec3 *n3 = newVec(0.0,0.0,1.0);
	struct vec3 *corn = newVec(-2.5284, -2.5284, -2.0);
	double dims[3] = {0.3, 0.3, 4};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
	norms[2] = n3;
        struct vec3 *result = boxGrad(norms, corn, dims, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

double leg3SDF(struct vec3 *v, struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(-2.5284, 2.5284, -2.0);
        double dims[3] = {0.3, 0.3, 4};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        double result = genSDFbox(norms, corn, dims, v, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

struct vec3 *leg3SDFgrad(struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(-2.5284, 2.5284, -2.0);
        double dims[3] = {0.3, 0.3, 4};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        struct vec3 *result = boxGrad(norms, corn, dims, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

double leg4SDF(struct vec3 *v, struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(-5.057, 0.0, -2.0);
        double dims[3] = {0.3, 0.3, 4};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        double result = genSDFbox(norms, corn, dims, v, p);
        free(n1);
        free(n2);
	free(n3);
        //free(dims);
        free(corn);
        return result;
}

struct vec3 *leg4SDFgrad(struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(-5.057, 0.0, -2.0);
        double dims[3] = {0.3, 0.3, 4};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        struct vec3 *result = boxGrad(norms, corn, dims, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

double boxSDF(struct vec3 *v, struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(0.0, 0.0, -0.5);
        double dims[3] = {2, 2, 1.5};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        double result = genSDFbox(norms, corn, dims, v, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

struct vec3 *boxSDFgrad(struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(0.0, 0.0, -0.5);
        double dims[3] = {2, 2, 1.5};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        struct vec3 *result = boxGrad(norms, corn, dims, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

double roundBoxSDF(struct vec3 *v, struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(0.0, -2.0, -1.0);
        double dims[3] = {2, 2, 1.5};
	struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
	struct vec3 *dif = vecsubtract(p,corn);
	double result = genSDFbox(norms, corn, dims, v, p);
	struct vec3 *onbox = vecadd(p, vecmult(result-0.001, v));
	if (result != 1111.1 && vecEq(boxGrad(norms, corn, dims, onbox), newVec(0,0,0))!=1){
		//printf("%f, ", result);
		result = result - (0.5/dot3(boxGrad(norms, corn, dims, onbox), v));
		//printf("%f\n", -1/dot3(boxGrad(norms, corn, dims, onbox), v));
	}
	free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

struct vec3 *roundboxSDFgrad(struct vec3 *p){
        struct vec3 *n1 = normalize(newVec(1.0,-1.0,0.0));
        struct vec3 *n2 = normalize(newVec(1.0,1.0,0.0));
        struct vec3 *n3 = newVec(0.0,0.0,1.0);
        struct vec3 *corn = newVec(0.0, -2.0, -1.0);
        double dims[3] = {2, 2, 1.5};
        struct vec3 *norms[3];
        norms[0] = n1;
        norms[1] = n2;
        norms[2] = n3;
        struct vec3 *result = RoundboxGrad(0.5,norms, corn, dims, p);
        free(n1);
        free(n2);
        free(n3);
        //free(dims);
        free(corn);
        return result;
}

double stdBoxSDF(struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *p){
	double result = 0;
        struct vec3 *dif = vecsubtract(p,corner);
	for (int i = 0; i<3; i++){
                if (dot3(normals[i],dif) > 0) result = result + pow(dot3(normals[i],dif),2);
	}
	return result;
	
}

double cylinder1(struct vec3 *v, struct vec3 *p){
        double h = 2.7;
        double r = 1;
        struct vec3 *b = newVec(0.0,-1.0,-1.0);
        struct vec3 *a = normalize(newVec(-0.6,1.0,0.0));
        struct vec3 *c = vecadd(vecmult(dot3(vecsubtract(p,b),a), a), b);
        double d1 = (cylSDF(v,p,c,r,a));
        double d2 = (genSDFplane(a, vecadd(b, vecmult(h, a)), v, p));
        double d3 = -1.0*(genSDFplane(a, b, v, p));
        double max = -1111.1;
        if (d1 == d2 && d2 == d3) max = d1;
	if (d1 > max && d1 != 1111.1) max = d1;
        if (d2 > max && d2 != 1111.1) max = d2;
        if (d3 > max && d3 != 1111.1) max = d3;
        if (max == 1111.1) return 1111.1;
        struct vec3 *direct = vecmult(max, v);
        struct vec3 *point = vecadd(p, direct);
        int zerocount = 0;
        struct vec3 *axdist = vecsubtract(vecsubtract(point,b), vecmult(dot3(vecsubtract(point, b), a), a));
        if (max != d1 && dot3(axdist,axdist) > r*r) max = 1111.1;
        free(axdist);
        free(direct);
        free(point);
        return max;
}

double cylinder2(struct vec3 *v, struct vec3 *p){
        double h = 0.2;
        double r = 4.5*tan(3.14/18.0);
        struct vec3 *a = normalize(newVec(0.5,0.0,1.0));
        struct vec3 *b = vecadd(newVec(0.0,0.0, -3.5), vecmult(4.5, a));
	struct vec3 *c = vecadd(vecmult(dot3(vecsubtract(p,b),a), a), b);
        double d1 = (cylSDF(v,p,c,r,a));
        double d2 = (genSDFplane(a, vecadd(b, vecmult(h, a)), v, p));
        double d3 = -1.0*(genSDFplane(a, b, v, p));
	double height = dot3(vecsubtract( vecadd(vecmult(d1, v), p), b), a); 
	double max;
	if (dot3(vecsubtract( vecadd(vecmult(d1, v), p), b), a) > 0 && dot3(vecsubtract( vecadd(vecmult(d1, v), p), b), a) < h) max = d1;
	else max = d2;
	struct vec3 *dif = vecsubtract(p,b);
        if (d1 == d2 && d2 == d3) max = d1;
        if (d2 > d1 && d2 != 1111.1) max = d2;
        if (d3 > max && d3 != 1111.1) max= d3;
        if (max == 1111.1) return 1111.1;
        struct vec3 *direct = vecmult(max, v);
        struct vec3 *point = vecadd(p, direct);
        int truthCount = 0;
        int zerocount = 0;
        struct vec3 *axdist = vecsubtract(vecsubtract(point,b), vecmult(dot3(vecsubtract(point, b), a), a));
        if (max != d1 && dot3(axdist,axdist) > r*r) max = 1111.1;
        free(axdist);
        free(direct);
        free(point);
        return max;
}

struct vec3 *cylinderGrad1(struct vec3 *p){
	double h = 2.7;
	struct vec3 *b = newVec(0.0,-1.0,-1.0);
	struct vec3 *a = normalize(newVec(-0.6,1.0,0.0));
	if (0.000001 < dot3(vecsubtract(p,b), a) && dot3(vecsubtract(p,b), a) < h) return cylgrad(p, b, a);
	else  if ((h-0.000001) <= dot3(vecsubtract(p,b), a)) return a;
	else return vecmult(-1,a);	
}

struct vec3 *cylinderGrad2(struct vec3 *p){
        double h = 0.2;
        struct vec3 *a = normalize(newVec(0.5,0.0,1.0));
        struct vec3 *b = vecadd(newVec(0.0,0.0, -3.5), vecmult(4.5, a));
	if (0 < dot3(vecsubtract(p,b), a) && dot3(vecsubtract(p,b), a) < (h-0.0001)) return cylgrad(p, b, a);
        else  if ((h-0.0001) <= dot3(vecsubtract(p,b), a)) return a;
        else return vecmult(-1,a);      
}

double infCyl(struct vec3 *v, struct vec3 *p){
	double r = 1;
	struct vec3 *a = normalize(newVec(1.0,1.0,0.0));
	struct vec3 *b = newVec(0.0,0.0,-1.0);
	struct vec3 *c = vecadd(vecmult(dot3(vecsubtract(p,b),a), a), b);
	return cylSDF(v,p,c,r,a);
}

struct vec3 *infCylGrad(struct vec3 *p){
        struct vec3 *a = normalize(newVec(1.0,1.0,0.0));
	struct vec3 *b = newVec(0.0,0.0,-1.0);
	return cylgrad(p,b,a);
}

double prismgenSDF(struct vec3 *v, struct vec3 *p){
	struct vec3 *base = newVec(0,0,-2); 
	struct vec3 *a = newVec(0,0,1);
	struct vec3 *n0 = newVec(-1,0,0);
	int faces = 10;
	double height = 1;
	double facewidth = 1; 
	return prismSDF(base, a, n0, faces, height, facewidth, v, p);
}

struct vec3 *prismSDFgrad(struct vec3 *p){
	struct vec3 *base = newVec(0,0,-2);
        struct vec3 *a = newVec(0,0,1);
        struct vec3 *n0 = newVec(-1,0,0);
        int faces = 10;
        double height = 1;
        double facewidth = 1;
	return prismGrad(base, a, n0, faces, height, facewidth, p);
}

double coneGenSDF(struct vec3 *v, struct vec3 *p){
	struct vec3 *a = newVec(0.0,0.0,1.0);
	double h  = 2.0;
	double t = 0.52359877559;
	struct vec3 *b = newVec(0,0,-1);
	return coneSDF(v, p, b, a, t);
}

struct vec3 *coneSDFGrad(struct vec3 *p){
	struct vec3 *a = newVec(0.0,0.0,1.0);
        double h  = 2.0;
        double t = 0.52359877559;
	struct vec3 *b = newVec(0,0,-1);
	return coneGrad(p, b, a, t);
}

struct contact *cylEval(struct vec3 *v, struct vec3 *p){
	struct vec3 *a = newVec(0.0,0.0,1.0);
	struct vec3 *b = newVec(0.0,0.0,-1.0);	
	double h = 2.0;
	double r = 0.5;
	struct vec3 *c = vecadd(component(vecsubtract(p,b),a),b);
	return evalCylinder(v, p, c, r, a, h, b);
}

double infcone(struct vec3 *v, struct vec3 *p){
	struct vec3 *b = newVec(0.0,0.0, -1.0);
	struct vec3 *a = normalize(newVec(-0.75 ,0.0, -1.0));
	double t = 3.14/9.0;
	return coneSDF(v, p, b, a, t);
//	if (bottom < top) return bottom;
//	else return top;
}

struct vec3 *infconeGrad(struct vec3 *p){
	struct vec3 *b = newVec(0.0,0.0, -1.0);
        struct vec3 *a = normalize(newVec(-0.75 ,0.0, -1.0));
        double t = 3.14/9.0;
	struct vec3 *dif = vecsubtract(p,b);
	return coneGrad(p, b, a, t);
}

double iceCreamCone(struct vec3 *v, struct vec3 *p){
	struct vec3 *b = newVec(0.0,0.0, -3.5);
        struct vec3 *a = normalize(newVec(0 ,0.5, 1.2));
        double t = 3.14/18.0;
        struct vec3 *dif = vecsubtract(p,b);
	double h = 4.5;
	struct vec3 *c = vecadd(b, vecmult(h, a));
	double d1 = genSDFplane(a, c, v, p);
	double d2 = coneSDF(v, p, b, a, t);
	double cone;
	double d3 = genSDFsphere(c, h*tan(t), v, p);
	struct vec3 *direc = vecadd(p, vecmult(d1, v));
	double ret;
	struct vec3 *difDirec = vecsubtract(vecadd(p, vecmult(d2, v)), b);
	if (d1>d2 && dot3(vecsubtract(direc, c), vecsubtract(direc, c)) <= h*h*tan(t)*tan(t)) return d1;
	else if (d2>d1 && h > dot3(difDirec, a) && dot3(difDirec, a)>0) return d2;
	else return 1111.1;
//	return d2;
}

struct vec3 *iceCreamConeGrad(struct vec3 *p){
	struct vec3 *b = newVec(0.0,0.0, -3.5);
        struct vec3 *a = normalize(newVec(0 ,0.5, 1.2));
        double t = 3.14/18.0;
        struct vec3 *dif = vecsubtract(p,b);
        double h = 4.5;
//	printf("difdot: %f, %f \n",dot3(normalize(dif), a),cos(t));
	if (valCompare(dot3(normalize(dif), a), cos(t)) == 1) {
//		printf("CONEGRAD\n");
		return coneGrad(p, b, a, t);
	}
	else return a;

//	return coneGrad(p, b, a, t);
}

double iceCream(struct vec3 *v, struct vec3 *p){
	struct vec3 *b = newVec(0.0,0.0, -3.5);
        struct vec3 *a = normalize(newVec(0.5 ,0.0, 1.0));
	double h = 4.5;
	struct vec3 *c = vecadd(b, vecmult(h, a));
	double t = 3.14/18.0;
	return genSDFsphere(c, h*tan(t), v, p);
}

struct vec3 *iceCreamGrad(struct vec3 *p){
	struct vec3 *b = newVec(0.0,0.0, -3.5);
        struct vec3 *a = normalize(newVec(0.5 ,0.0, 1.0));
        double h = 4.5;
        struct vec3 *c = vecadd(b, vecmult(h, a));
	return normalize(vecsubtract(p,c));
}
