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

double maxim(double arr[], int size){
	double maxi=0;
	for (int i = 0; i < size; i++) if (arr[i] > maxi) maxi = arr[i];
	return maxi;
	
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

int valCompare(double a, double b){
//	printf("%f\n", (a-b)*(a-b));
	if ((a-b)*(a-b) < 0.000001 ) return 1;
	else return 0;
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

double cylSDF(struct vec3 *v, struct vec3 *p, struct vec3 *c, double r, struct vec3 *a){
        a = normalize(a);
	struct vec3 *vperp = normalize(vecsubtract(v, vecmult(dot3(a,v),a)));
	if (dot3(v, vperp)==0) return 1111.1;
        return genSDFsphere(c, r, vperp, p)/(dot3(v,vperp));
}

struct vec3 *cylgrad(struct vec3 *p, struct vec3 *b, struct vec3 *a){
	a = normalize(a);
	struct vec3 *pb = vecsubtract(p,b);
        return normalize(vecsubtract(pb, vecmult(dot3(pb,a), a)));
}

struct contact *evalCylinder(struct vec3 *v, struct vec3 *p, struct vec3 *c, double r, struct vec3 *a, double h, struct vec3 *b){
	struct contact *point = malloc(sizeof(struct contact));
	double d3;
	a = normalize(a);
        struct vec3 *vperp = normalize(vecsubtract(v, vecmult(dot3(a,v),a)));
        if (dot3(v, vperp)==0) d3 = 1111.1;
        else d3 = genSDFsphere(c, r, vperp, p)/(dot3(v,vperp));
	double d1 = genSDFplane(a, vecadd(b, vecmult(h,a)), v, p);
	double d2 = genSDFplane(vecmult(-1,a), b, v, p);
	double SDFs[3] = {d1, d2, d3};
	double maxi = maxim(SDFs, 3);
	if (maxi == 1111.1) {
		point->dist = 1111.1;
		return point;
	}
	struct vec3 *pointOn = vecadd(vecmult(maxi, v), p);
	if (maxi == d3 && dot3(vecsubtract(p, b),a)<=h && dot3(vecsubtract(p, b),a)>=0){
		point->grad = normalize(antiComponent(vecsubtract(pointOn,b),a));
		printf("Top cap : %f, Bottom cap: %f, cylinder : %f, maxi = %f\n", d1, d2, d3, maxi);
		printf("%f %f %f\n", point->grad->x, point->grad->y, point->grad->z);
		point->dist = d3;
		printf("cylinderHIt");
	}
	else if (maxi == d3 || vecnorm(antiComponent(vecsubtract(pointOn, b), a)) > r){
		point->dist = 1111.1;
                return point;
	}
	else if (maxi ==d1 ){
		point->grad = a;
		point->dist = d1;
	}
	else{
		point->grad = vecmult(-1,a);
                point->dist = d2;
	}
	return point;
}

double coneSDF(struct vec3 *v, struct vec3 *p, struct vec3 *b, struct vec3 *a, double theta){
	struct vec3 *dif = vecsubtract(p,b);
	double al = dot3(dif, v);
	double bet = dot3(dif,a);
	double gam = dot3(v,a);
	double d =vecnorm(dif);
	double k = cos(theta);
	free(dif);
	if (gam == k || (bet*bet - 2*al*bet*gam + al*al*k*k + d*d*(gam*gam-k*k)) < 0) return 1111.1;
	//else return (k*k*al-bet*gam - 2*k*sqrt(k*k*(al*al-d*d) + bet*bet + d*d*gam*gam - 2*al*bet*gam))/(gam*gam - k*k);
	else return (-bet*gam + al*k*k - fabs(k)*sqrt(bet*bet - 2*al*bet*gam + al*al*k*k + d*d*(gam*gam-k*k)))/(gam*gam - k*k);
}

struct vec3 *coneGrad(struct vec3 *p, struct vec3 *b, struct vec3 *a, double theta){
	struct vec3 *dif = vecsubtract(p,b);
	struct vec3 *toax = antiComponent(dif, a);
//	struct vec3 *grad = normalize(vecadd(vecmult(cos(theta),a), vecmult(sin(theta),toax)));
//	if (dot3(dif, a) < h) return grad;
//	else return a;
	struct vec3 *tangent = cross3(a, toax);
	struct vec3 *grad = normalize(cross3(dif, tangent));
//	printf("%f, %f, %f\n", grad->x, grad->y, grad->z);
//	printf("%f\n", dot3(grad, tangent));
//	printf(" toax^2 = %f    tangent^2 = %f   Grad^2 = %f", dot3(toax,toax), dot3(tangent,tangent), dot3(grad, grad));
	return grad;
}

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

/*struct contact *evalPrism(struct vec3 *base, struct vec3 *a, struct vec3 *n0, int faces, double height, double facewidth, struct vec3 *v, struct vec3 *p){
	n0 = normalize(n0);
	struct contact *result = malloc(sizeof(struct *contact));
        struct vec3 *norm = malloc(sizeof(struct vec3));
        struct vec3 *n1 = cross3(a, n0);
        double max = -1111.1;
        double curdist;
        struct vec3 *pointon;
        struct vec3 *directed;
        int count = 0;
        struct vec3 *b;
        double thickness = facewidth/(2*tan(3.141592/faces));
	double maxOff = -1*(height + 2*thickness);
	struct vec3 *bestNorm;
	struct vec3 *bestPointon;
	struct vec3 *bestDirect;
        for (int i=0; i < faces; i++){
                norm = vecadd(vecmult(cos(i*6.283/faces), n0), vecmult(sin(i*6.283/faces),n1));
                pointon = vecadd(vecmult(thickness, norm), base);
                curdist = genSDFplane(norm, pointon,v,p);
                directed = vecadd(p, vecmult(curdist, v));
                b = vecsubtract(directed, pointon);
                if (curdist == 1111.1) break;
                if (curdist > max) {
                        max = curdist;
			bestNorm = norm;
			bestB = b;
			bestDirect = directed;
                }
        }
	if (0 <= dot3(bestB, a) && dot3(bestB, a) <= height && abs(vecnorm(antiComponent(bestB,a))) <= facewidth/2) {
		result->dist = max;
		result->grad = bestNorm;
		result->p = bestDirect;
		return result;
		
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
}*/

struct vec3 *shade(struct Node *head, struct Camera *camera, struct SD *sd){
	if (sd->dist == 1111.1) return newVec(0, 0, 0);
	struct vec3 *grad = sd->gradient(sd->pos);
	struct vec3 *lightnorm = normalize(vecsubtract(camera->light, sd->pos));
	struct vec3 *lightv = vecmult(-1, lightnorm);
	struct SD *sdlight = calcSDF(head, lightv, camera->light);
	if (vecEq(sd->pos, sdlight->pos)==0) return newVec(0, 0, 0);
//	printf("matches");
	double scale = dot3(grad, lightnorm);
	if (scale < 0) scale = 0;
	struct vec3 *shade = vecmult(scale, sd->color);
	shade->x = round(shade->x);
	shade->y = round(shade->y);
	shade->z = round(shade->z);
//	printf("scale = %f    shade = %f %f %f\n", scale, shade->x, shade->y, shade->z);
	free(grad);
	free(lightnorm);
	free(sd);
	return shade;
}


struct SD *calcSDF(struct Node *head, struct vec3 *v, struct vec3 *p){
	struct Node *cur = head;
	struct Node *minNode = head;
	double min = head->shape->SDF(v,p);
	while (cur->next != 0){
	 	cur = cur->next;
		if (cur->shape->SDF(v, p) < min) {
			min = cur->shape->SDF(v, p);
			minNode = cur;
//			if (cur != head) printf("Cyl2 selected\n");
		}
		
	}
	
//	if (head->shape->SDF(v,p) != 1111.1 && head->next->shape->SDF(v,p) != 1111.1 && head->next->shape->SDF(v,p) > head->shape->SDF(v,p)) printf("cone: %f and sph: %f\n", head->shape->SDF(v,p), head->next->shape->SDF(v,p));
	struct SD *sd = malloc(sizeof(struct SD));
	sd->dist = min;
	sd->gradient = minNode->shape->gradSDF;
	sd->color = minNode->shape->color;
	struct vec3 *direct = vecmult(min, v);
	sd->pos = vecadd(p, direct);
	free(direct);	
	return sd;
}

/*struct vec3 *shade(struct Node *head, struct Camera *camera, struct SD *sd){
        if (sd->dist == 1111.1) return newVec(0, 0, 0);
        struct vec3 *lightnorm = normalize(vecsubtract(camera->light, sd->pos));
        struct vec3 *lightv = vecmult(-1, lightnorm);
        if (vecEq(sd->pos, calcSDF(head, lightv, camera->light)->pos)==0) return newVec(0, 0, 0);
        double scale = dot3(sd->grad, lightnorm);
        if (scale < 0) scale = 0;
        struct vec3 *shade = vecmult(scale, sd->color);
        shade->x = round(shade->x);
        shade->y = round(shade->y);
        shade->z = round(shade->z);
        free(lightnorm);
        free(sd);
        return shade;
}

struct SD *calcSDF(struct Node *head, struct vec3 *v, struct vec3 *p){
        struct Node *cur = head;
        struct Node *minNode = head;
	struct contact *minCont = head->shape->eval(v,p);
	printf("dis = %f\n", minCont->dist);
        double min = minCont->dist;
	struct vec3 *minGrad = minCont->grad;
        while (cur->next != 0){
                cur = cur->next;
                if (cur->shape->eval(v, p)->dist < min) {
                        minCont = cur->shape->eval(v, p);
			min = minCont->dist;
			minNode = cur;
//                      if (cur != head) printf("Cyl2 selected\n");
                }
                
        }
        
//      printf("cyl1: %f and cyl2: %f\n", head->shape->SDF(v,p), head->next->shape->SDF(v,p));
        struct SD *sd = malloc(sizeof(struct SD));
        sd->dist = min;
        sd->grad = minCont->grad;
        sd->color = minNode->shape->color;
        struct vec3 *direct = vecmult(min, v);
        sd->pos = vecadd(p, direct);
        free(direct);   
        return sd;
}*/

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
	for (int i = 0; i < vsteps; i ++){
		hraw = vecmult(((double) i)/((double) vsteps), h);
		hnow = vecadd(hraw, camera->corner[0]);
		free(hraw);
		for (int j = 0; j < hsteps; j++){
			know = vecmult(((double) j)/((double) hsteps), k);
			sheetpoint = vecadd(hnow, know);
			v = normalize(vecsubtract(sheetpoint, camera->focus));
			sdcur = calcSDF(head, v, camera->focus);
			pointon = vecadd(camera->focus, vecmult(sdcur->dist, v));
			grad = sdcur->grad;
			double dis = sdcur->dist;
			lightvec = normalize(vecsubtract(camera->light, sdcur->pos));
			curcolor = shade(head, camera, sdcur);
//			if (dis != 1111.1 && curcolor->x == 0) printf("%f, %f, %f at %f, %f, %f\n", grad->x, grad->y, grad->z, pointon->x, pointon->y, pointon->z);
//			if (dis != 1111.1) printf("%f\n", pointon->z);
			if (j == 0) fprintf(fp, "%f,%f,%f", curcolor->x, curcolor->y, curcolor->z);
			else fprintf(fp, ";%f,%f,%f", curcolor->x, curcolor->y, curcolor->z);
			if (dis != 1111.1) printf("* ");
                        else printf("  ");
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

