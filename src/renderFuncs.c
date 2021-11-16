#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <png.h>
#include "SDF.h"

struct Camera *newCam(struct vec3 *focus, struct vec3 *corner[4], struct vec3 *light){
        struct Camera *newGuy = malloc(sizeof(struct Camera));
        newGuy->focus = focus;
        newGuy->corner[0] = corner[0];
        newGuy->corner[1] = corner[1];
        newGuy->corner[2] = corner[2];
        newGuy->corner[3] = corner[3];
	newGuy->light = light;
        return newGuy;
}

struct vec3 *newVec(double x, double y, double z){
        struct vec3 *vec = malloc(sizeof(struct vec3));
        //struct vec3 *vec;
	vec->x = x;
        vec->y = y;
        vec->z = z;
        return vec;
}

struct Sphere *newSphere(double radius, struct vec3 *center){
	struct Sphere *new = malloc(sizeof(struct Sphere));
	new->radius = radius;
	new->center = center;
	return new;
}

struct Plane *newPlane(struct vec3 *normal, struct vec3 *pointon){
	struct Plane *new = malloc(sizeof(struct Plane));
	new->normal = normal;
	new->pointOn = pointon;
	return new;
}

double dot3(struct vec3 *v1, struct vec3 *v2){
        return  (v1->x * v2->x) + (v1->y * v2->y) + (v1->z * v2->z);
}

struct vec3 *cross3(struct vec3 *v1, struct vec3 *v2){
	double e1 = v1->y*v2->z - v1->z*v2->y;
	double e2 = v1->z*v2->x - v1->x*v2->z;
	double e3 = v1->x*v2->y - v1->y*v2->x;
	struct vec3 *result = newVec(e1, e2, e3);
	return result;
}

struct vec3 *vecsubtract(struct vec3 *v1, struct vec3 *v2){
        struct vec3 *diff = malloc(sizeof(struct vec3));
        //struct vec3 *diff;
	diff->x = v1->x - v2->x;
        diff->y = v1->y - v2->y;
        diff->z = v1->z - v2->z;
        return diff;
}

struct vec3 *vecadd(struct vec3 *v1, struct vec3 *v2){
        struct vec3 *sum = malloc(sizeof(struct vec3));
        //struct vec3 *sum;
	sum->x = v1->x + v2->x;
        sum->y = v1->y + v2->y;
        sum->z = v1->z + v2->z;
        return sum;
}

struct vec3 *vecmult(double scalar, struct vec3 *vector){
        struct vec3 *prod = malloc(sizeof(struct vec3));
        //struct vec3 *prod;
	prod->x = scalar*(vector->x);
        prod->y = scalar*(vector->y);
        prod->z = scalar*(vector->z);
        return prod;
}

struct vec3 *normalize(struct vec3 *OG){
        double norm = pow((double) pow(OG->x, 2.0) + pow(OG->y, 2.0) + pow(OG->z, 2.0), 0.5);
        OG->x = (OG->x)/norm;
        OG->y = (OG->y)/norm;
        OG->z = (OG->z)/norm;
        return OG;
}

double vecnorm(struct vec3 *OG){
	return pow((double) pow(OG->x, 2.0) + pow(OG->y, 2.0) + pow(OG->z, 2.0), 0.5);
}

int vecEq(struct vec3 *v1, struct vec3 *v2){
	if (pow(v1->x - v2->x, 2) + pow(v1->y - v2->y, 2) + pow(v1->z - v2->z, 2) < 0.0001*0.0001) return 1;
	else return 0;
}

struct vec3 *antiComponent(struct vec3 *pb, struct vec3 *a){
	double len = dot3(pb, a);
	struct vec3 *component = vecmult(len, a);
	struct vec3 *antico = vecsubtract(pb, component);
	free(component);
	return antico;

}

struct vec3 *component(struct vec3 *pb, struct vec3 *a){
	return vecmult(dot3(pb, a), a);
}

double genSDFsphere(struct vec3 *center, double r, struct vec3 *v, struct vec3 *p){
        struct vec3 *a = vecsubtract(p, center);
        double result;
	if (pow(dot3(a, v), 2.0) + pow(r, 2.0) - dot3(a, a) < 0.0) return 1111.1;
        else result = -dot3(v, a) - pow((double) pow(dot3(a, v), 2.0) + pow(r, 2.0) - dot3(a, a), 0.5);
	//printf("%f\n", result);
        free(a);
	if (result > 0.0) return result;
	else return 1111.1;
}

double genSDFplane(struct vec3 *normal, struct vec3 *pointOn, struct vec3 *v, struct vec3 *p){
        struct vec3 *a = vecsubtract(p,pointOn);
	if (dot3(a, normal) == 0 || dot3(a, normal) == -0) return 0;
	if (dot3(v, normal) == 0) return 1111.1;
	double result = -1.0*fabs(dot3(a, normal))/dot3(v, normal);
	free(a);
	//if (result < 0) return 1111.1;
	return result;
}

double genSDFbox(struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *v, struct vec3 *p){
	double d1 = (genSDFplane(normals[0], corner, v, p));
	double d2 = (genSDFplane(normals[1], corner, v, p));
	double d3 = (genSDFplane(normals[2], corner, v, p));
	double max;
	if (d1 != 1111.1) max = d1;
	else max = 0;
	//printf("dists: %f, %f, %f\n", d1, d2, d3);
	if (d1 == d2 && d2 == d3) return d1;
	if (d2 > max && d2 != 1111.1) max = d2;
	if (d3 > max && d3 != 1111.1) max= d3;
	struct vec3 *direct = vecmult(max, v);
	struct vec3 *point = vecadd(p, direct);
	struct vec3 *onBox = vecsubtract(point, corner);
	int truthCount = 0;
	int zerocount = 0;
	for (int i = 0; i < 3; i++){
		if (dot3(onBox, normals[i]) >=  -1*dims[i] && dot3(onBox, normals[i]) <= 0.0001) truthCount++;
		printf("%f\n", (dot3(onBox, normals[i])));
	}
	//printf("\n");
	free(direct);
	free(point);
	free(onBox);
	//if (truthCount == 3 && zerocount > 1) printf("|");
	if (truthCount == 3) return max;
	else return 1111.1;
	
}

struct vec3 *boxGrad(struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *p){
	struct vec3 *a = vecsubtract(p, corner);
	struct vec3 *grad = newVec(0,0,0);
	int count = 0;
	double dists[3] = {dot3(normals[0],a), dot3(normals[1],a), dot3(normals[2],a)};
	for (int i = 0; i < 3; i++){
		if ((dists[i] < 0.0001 && dists[i] > -0.0001) || ((dists[i]-dims[i]) < 0.0001 && (dists[i]-dims[i]) > -0.0001))  {
			grad = vecadd(grad, normals[i]);
			count++;
			
		}
		else if (dists[i] > 0.0001 || dists[i] < -1*dims[i] - 0.0001) {
			count++;
			grad  = vecadd(grad, component(a,normals[i]));
		}
	}
	if (count ==  0){
		if (dists[0] == dists[1] &&dists[1] ==  dists[2]) grad = a;
	        if (dists[1] < dists[0] && dists[1]<dists[2]) grad = component(a, normals[1]);
        	if (dists[2] < dists[0] && dists[2]<dists[1]) grad = component(a, normals[2]);
		if (dists[0] < dists[1] && dists[0]<dists[2]) grad = component(a, normals[0]);
	}
	return normalize(grad);
}

struct vec3 *RoundboxGrad(double r, struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *p){
        struct vec3 *a = vecsubtract(p, corner);
        struct vec3 *grad = newVec(0,0,0);
        struct vec3 *newGrad = grad;
        struct vec3 *proj;
        int count = 0;
        double dists[3] = {-1*dot3(normals[0],a),  -1*dot3(normals[1],a), -1*dot3(normals[2],a)};
        for (int i = 0; i < 3; i++){
                if (dot3(normals[i],a) > r || dot3(normals[i],a) < r-1*dims[i]+r) {
                        count++;
                        grad  = vecadd(grad, component(a,normals[i]));
                }
		else if (dot3(normals[i],a) == r || dot3(normals[i],a) == r-1*dims[i]){
                        grad = vecadd(grad, normals[i]);
		}
        }
        if (count == 0){
                if (dists[0] == dists[1] &&dists[1] ==  dists[2]) grad = a;
                if (dists[1] > dists[0] && dists[1]>dists[2]) grad = vecmult(-1,normals[1]);
                if (dists[2] > dists[0] && dists[2]>dists[1]) grad = vecadd(grad, component(a, normals[2]));
                if (dists[0] > dists[1] && dists[0]>dists[2]) grad = vecadd(grad, component(a, normals[0]));
        }
         if (count == 0) return newGrad;
        return normalize(grad);
}

double cylSDF(struct vec3 *v, struct vec3 *p, struct vec3 *c, double r, struct vec3 *a, double h){
        struct vec3 *vperp = normalize(vecsubtract(v, vecmult(dot3(a,v),a)));
        return genSDFsphere(c, r, vperp, p)/(dot3(v,vperp));
}

struct vec3 *cylgrad(struct vec3 *p, struct vec3 *b, struct vec3 *a){
	struct vec3 *pb = vecsubtract(p,b);
        return normalize(vecsubtract(pb, vecmult(dot3(pb,a), a)));
}

double coneSDF(struct vec3 *v, struct vec3 *p, struct vec3 *b, struct vec3 *a, double h, double theta){
	struct vec3 *dif = vecsubtract(p,b);
	double al = dot3(dif, v);
	double bet = dot3(dif,a);
	double gam = dot3(v,a);
	double d =vecnorm(dif);
	double k = cos(theta);
	free(dif);
	return (k*k*al-bet*gam+2*k*sqrt(k*k*(al*al-d*d) + bet*bet + d*d*gam*gam - 2*al*bet*gam))/(gam*gam - k*k);
}

/*struct vec3 *coneGrad(struct vec3 *p, struct vec3 *b, struct vec3 *a, double h, double theta){
	struct vec3 *dif = vecsubtract(p,b);
	struct vec3 *toax = antiComponent(dif, a);
	double difa = dot3(dif, a);
	struct vec3 *sincos = malloc(sizeof(struct vec3));
	sincos->x = sin(theta);
	sincos->y = cos(theta);
	if (difa < 0) struct vec3 *grad = normalize(dif);
	else if (0 < difa < h) {
		struct vec3 *grad = malloc(sizeof(struct vec3));
		grad->x = 2*costheta*dot3(p,
	}
}*/

double prismSDF(struct vec3 *base, struct vec3 *a, struct vec3 *n0, int faces, double height, double facewidth, struct vec3 *v, struct vec3 *p){
	n0 = normalize(n0);
	struct vec3 *norm = malloc(sizeof(struct vec3));
	struct vec3 *n1 = cross3(a, n0);
	double max = 0;
	double curdist;
	struct vec3 *pointon;
	struct vec3 *directed;
	int count = 0;
	struct vec3 *b;
	double thickness = facewidth/(2*tan(3.141592/faces));
	for (int i=0; i < faces; i++){
		norm = vecadd(vecmult(cos(i*6.283/faces), n0), vecmult(sin(i*6.283/faces),n1));
		pointon = vecadd(vecmult(thickness, norm), base);
		curdist = genSDFplane(norm, pointon,v,p);
		directed = vecadd(p, vecmult(curdist, v));
		b = vecsubtract(directed, pointon);
		if (curdist == 1111.1) break;
		if (curdist > max && curdist != 1111.1 && 0 <= dot3(b, a) && dot3(b, a) <= height && vecnorm(antiComponent(b,a)) <= facewidth/2 && vecnorm(antiComponent(b,a)) >= -1*facewidth/2) {
			max = curdist; 
		}
	}
	pointon = vecadd(base, vecmult(height, a));
	curdist = genSDFplane(a, pointon, v, p);
	if (max == 0 && curdist > max && curdist != 1111.1) {
        	directed = vecadd(p, vecmult(curdist, v));
		for (int i=0; i < faces; i++) {
			norm = vecadd(vecmult(cos(i*6.283/faces), n0), vecmult(sin(i*6.283/faces),n1)); 
			if (dot3(vecsubtract(directed, pointon), norm) <= thickness) count ++;
			else break;
			}
		if (count == faces) max = curdist;
	}
	count = 0;
	pointon = base;
        curdist = genSDFplane(vecmult(-1.0,a), pointon, v, p);
        if (max == 0 && curdist > max && curdist != 1111.1) {
                directed = vecadd(p, vecmult(curdist, v));
                for (int i=0; i < faces; i++) {
                        norm = vecadd(vecmult(cos(i*6.283/faces), n0), vecmult(sin(i*6.283/faces),n1));
                        if (dot3(vecsubtract(directed, pointon), norm) <= thickness) count ++;
                        else break;
                        }
                if (count == faces) max = curdist;
        }
	
	if (max == 0) max = 1111.1;
	
	
	free(pointon);
	free(norm);
	free(n1);
	return max;

}

struct vec3 *prismGrad(struct vec3 *base, struct vec3 *a, struct vec3 *n0, int faces, double height, double facewidth, struct vec3 *p){
	n0 = normalize(n0);
	struct vec3 *grad = newVec(0,0,0);
	struct vec3 *norm;
	struct vec3 * b;
	struct vec3 *n1 = cross3(a, n0);
	double thickness = facewidth/(2*tan(3.141592/faces));
	double max = -1*(height + 2*thickness);
	for (int i=0; i < faces; i++){
                norm = vecadd(vecmult(cos(i*6.283/faces), n0), vecmult(sin((i*6.283/faces)),n1));
		b = vecsubtract(p,vecadd(vecmult(thickness, norm), base));
		if (dot3(b,norm) > max){
			max = dot3(b,norm);
			grad = norm;
		}
	}
	b = vecsubtract(p,vecadd(vecmult(height, a), base));
	if (dot3(b, a) > max) {
		max = dot3(b, a);
		grad = a;
	}
	b = vecsubtract(p, base);
	
	if (-1.0*dot3(a,b) > max) grad = vecmult(-1,a);
	
	return grad;
}
	
struct vec3 *shade(struct Node *head, struct Camera *camera, struct SD *sd){
	if (sd->dist == 1111.1) return newVec(0, 0, 0);
	struct vec3 *grad = sd->gradient(sd->pos);
	struct vec3 *lightnorm = normalize(vecsubtract(camera->light, sd->pos));
	struct vec3 *lightv = vecmult(-1, lightnorm);
	if (vecEq(sd->pos, calcSDF(head, lightv, camera->light)->pos)==0) return newVec(0, 0, 0);
	double scale = dot3(grad, lightnorm);
	if (scale < 0) scale = 0;
	struct vec3 *shade = vecmult(scale, sd->color);
	free(grad);
	free(lightnorm);
	free(sd);
	return shade;
}

struct SD *calcSDF(struct Node *head, struct vec3 *v, struct vec3 *p){
	struct Node *cur = head;
	struct Node *minNode = head;
	double min = 1111.1;
	while (cur->next != 0){
		if (cur->shape->SDF(v, p) < min) {
			min = cur->shape->SDF(v, p);
			minNode = cur;
		}
		cur = cur->next;
		
	}
	if (cur->shape->SDF(v, p) < min) {
		min = cur->shape->SDF(v, p);
		minNode = cur;
		
	}
	struct SD *sd = malloc(sizeof(struct SD));
	sd->dist = min;
	sd->gradient = minNode->shape->gradSDF;
	sd->color = minNode->shape->color;
	struct vec3 *direct = vecmult(min, v);
	sd->pos = vecadd(p, direct);
	free(direct);	
	return sd;
}

void newNode(struct Node *head, struct Shape *shape){
	struct Node *cur = head;
	while (cur->next !=0){
		cur = cur->next;
		printf("step\n");
	}
	cur->next = malloc(sizeof(struct Node));
	cur->next->shape = shape;
	cur->next->next = 0;
}

struct Shape *newShape(struct vec3 *color){
        struct Shape *newShape = malloc(sizeof(struct Shape));
        newShape->color = color;
        return newShape;
}


int Render(struct Node *head, struct Camera *camera, double hsteps, double vsteps, FILE *fp){
	//fprintf(fp, "a");
	printf("MadeItHere\n");
	struct vec3 *sheetpoint;
	struct vec3 *k = vecsubtract(camera->corner[1], camera->corner[0]);
	struct vec3 *h = vecsubtract(camera->corner[3], camera->corner[0]);
	struct vec3 *hnow;
	struct vec3 *know;
	struct vec3 *v;
	struct SD *sdcur;
	struct vec3 *curcolor;
	struct vec3 *hraw;
	struct vec3 *pointon;
	struct vec3 *lightvec;
	struct vec3 * grad;
	//fprintf(fp, "[");
	for (int i = 0; i < vsteps; i ++){
		hraw = vecmult(((double) i)/((double) vsteps), h);
		hnow = vecadd(hraw, camera->corner[0]);
		//fprintf(fp, "[");
		free(hraw);
		for (int j = 0; j < hsteps; j++){
			know = vecmult(((double) j)/((double) hsteps), k);
			sheetpoint = vecadd(hnow, know);
			v = normalize(vecsubtract(sheetpoint, camera->focus));
			sdcur = calcSDF(head, v, camera->focus);
			pointon = vecadd(camera->focus, vecmult(sdcur->dist, v));
			grad = sdcur->gradient(pointon);
			double dis = sdcur->dist;
			lightvec = normalize(vecsubtract(camera->light, sdcur->pos));
			curcolor = shade(head, camera, sdcur);
			//if (dis != 1111.1 && curcolor->x == 0) printf("%f, %f, %f at %f, %f, %f\n", grad->x, grad->y, grad->z, pointon->x, pointon->y, pointon->z);
//			if (dis != 1111.1) printf("%f\n", pointon->z);
			if (j == 0) fprintf(fp, "%f,%f,%f", curcolor->x, curcolor->y, curcolor->z);
			else fprintf(fp, ";%f,%f,%f", curcolor->x, curcolor->y, curcolor->z);
			//if (curcolor->y != 0) printf("o");
			//else printf(" ");
			//printf("%f\n", curcolor->y);
			//if (curcolor->y != 0 && !isnan(curcolor->y)) printf("o %f,%f,%f\n", pointon->x, pointon->y, pointon->z);
			//if ((curcolor->y)==0) printf("* %f,%f,%f\n", pointon->x, pointon->y, pointon->z);
                        //if (curcolor->y != 0 && !isnan(curcolor->y)) printf("o ");
			//else printf("* ");
			if (dis != 1111.1) printf("* ");
                        else printf("  ");
			//else printf("\n. point on: %f, %f, %f   gradient: %f, %f, %f   light: %f, %f, %f  dot: %f\n", pointon->x, pointon->y, pointon->z, grad->x, grad->y, grad->z, lightvec->x, lightvec->y, lightvec->z, dot3(grad, lightvec));
			//else printf(". ");
			//free(sheetpoint);
			//free(know);
			//free(v);
			//free(curcolor);
		}
		printf("\n");
		fprintf(fp, ";\n");
		fprintf(fp, "\n");
		//free(hnow);
	}	
	fprintf(fp, "\n");
	printf("\n");
	free(h);
	free(k);
	return 0;
}


