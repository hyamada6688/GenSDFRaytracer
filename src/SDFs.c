#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SDF.h"

/*int Render(struct Node *head, struct Camera *camera, double hsteps, double vsteps, FILE *fp);

double genSDFbox(struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *v, struct vec3 *p);

struct vec3 *boxGrad(struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *p);

struct vec3 *RoundboxGrad(double r, struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *p);

double cylSDF(struct vec3 *v, struct vec3 *p, struct vec3 *c, double r, struct vec3 *a, double h);

struct vec3 *cylgrad(struct vec3 *p, struct vec3 *b, struct vec3 *a);
  
*/
double circleSDF(struct vec3 *v, struct vec3 *p){
        return genSDFsphere(newVec(-1.0, 0.0, -0.5), 0.5, v, p);
}


struct vec3 *circleSDFGrad(struct vec3 *p){
        return normalize(vecsubtract(p, newVec(-1.0, 0.0, -0.5)));
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
        struct vec3 *corn = newVec(0.0, 2.0, -1.0);
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
        struct vec3 *corn = newVec(0.0, 2.0, -1.0);
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


double cylinder(struct vec3 *v, struct vec3 *p){
	double h = 2.0;
	double r = 0.5;
        struct vec3 *b = newVec(0.0,0.0,-1.0);
        struct vec3 *a = newVec(0.0,0.0,1.0);
	struct vec3 *c = vecadd(vecmult(dot3(vecsubtract(p,b),a), a), b);
	double d1 = (cylSDF(v,p,c,r,a,h));
        double d2 = (genSDFplane(a, vecadd(b, vecmult(h, a)), v, p));
        double d3 = (genSDFplane(a, b, v, p));
        double max = d1;
	if ( 0 < dot3(vecsubtract(p,b), a) && dot3(vecsubtract(p,b), a) < (h-0.01)){
		d3 = -d3;
		d2 = -d2;
	} 
	else if (dot3(vecsubtract(p,b), a) >= (h-0.01)) d3 = -d3;
	else d2 = -d2;
        if (d1 == d2 && d2 == d3) return d1;
        if (d2 > d1) max = d2;
        if (d3 > max) max= d3;
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

struct vec3 *cylinderGrad(struct vec3 *p){
	double h = 2.0;
	struct vec3 *b = newVec(0.0,0.0,-1.0);
	struct vec3 *a = newVec(0.0,0.0,1.0);
	if (0 < dot3(vecsubtract(p,b), a) && dot3(vecsubtract(p,b), a) < (h-0.01)) return cylgrad(p, b, a);
	else  if ((h-0.1) <= dot3(vecsubtract(p,b), a)) return a;
	else return vecmult(-1,a);	
}

double prismgenSDF(struct vec3 *v, struct vec3 *p){
	struct vec3 *base = newVec(0,0,-2); 
	struct vec3 *a = newVec(0,0,1);
	struct vec3 *n0 = newVec(-1,0,0);
	int faces = 100;
	double height = 1;
	double facewidth = 0.1; 
	return prismSDF(base, a, n0, faces, height, facewidth, v, p);
}

struct vec3  *prismSDFgrad(struct vec3 *p){
	struct vec3 *base = newVec(0,0,-2);
        struct vec3 *a = newVec(0,0,1);
        struct vec3 *n0 = newVec(-1,0,0);
        int faces = 100;
        double height = 1;
        double facewidth = 0.1;
	return prismGrad(base, a, n0, faces, height, facewidth, p);

}
