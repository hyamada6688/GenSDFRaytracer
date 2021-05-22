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

double genSDFsphere(struct vec3 *center, double r, struct vec3 *v, struct vec3 *p){
        struct vec3 *a = vecsubtract(p, center);
        double result;
	if (pow(dot3(a, v), 2.0) + pow(r, 2.0) - dot3(a, a) < 0.0) return 1111.1;
        else result = -dot3(v, a) + pow((double) pow((dot3(a, v)), 2.0) + pow(r, 2.0) - dot3(a, a), 0.5);
	//printf("%f\n", result);
        free(a);
	if (result > 0.0) return result;
	else return 1111.1;
}

double genSDFplane(struct vec3 *normal, struct vec3 *pointOn, struct vec3 *v, struct vec3 *p){
        struct vec3 *a = vecsubtract(pointOn, p);
	if (dot3(a, normal) == 0 || dot3(a, normal) == -0) return 0;
	if (dot3(v, normal) == 0) return 1111.1;
	double result = dot3(a, normal)/dot3(v, normal);
	free(a);
	if (result < 0) return 1111.1;
	else return result;
}

double genSDFbox(struct vec3 *normals[3], struct vec3 *corner, double dims[3], struct vec3 *v, struct vec3 *p){
	double d1 = (genSDFplane(normals[0], corner, v, p));
	double d2 = (genSDFplane(normals[1], corner, v, p));
	double d3 = (genSDFplane(normals[2], corner, v, p));
	double max = d1;
	//printf("dists: %f, %f, %f\n", d1, d2, d3);
	if (d1 == d2 && d2 == d3) return d1;
	if (d2 > d1) max = d2;
	if (d3 > max) max= d3;
	if (max == 1111.1) return 1111.1;
	struct vec3 *direct = vecmult(max, v);
	struct vec3 *point = vecadd(p, direct);
	struct vec3 *onBox = vecsubtract(point, corner);
	int truthCount = 0;
	int zerocount = 0;
	for (int i = 0; i < 3; i++){
		if (-1*dot3(onBox, normals[i]) <= dims[i] && 0.00001 >= dot3(onBox, normals[i])) truthCount++;
		//if (-1*dot3(onBox, normals[i]) <= dims[i] && 0 < dot3(onBox, normals[i]) && pow(dot3(onBox, normals[i]), 2) < 0.001) printf("'");
		//printf("TRUTH: %d (dot : %f, dim : %f), ", truthCount, -1*dot3(onBox, normals[i]), dims[i]);
		//if (pow(-1*dot3(onBox, normals[i]), 2) < 0.001) zerocount++;
	}
	//printf("\n");
	free(direct);
	free(point);
	free(onBox);
	//if (truthCount == 3 && zerocount > 1) printf("|");
	if (truthCount == 3) return max;
	else return 1111.1;
	
}

struct vec3 *shade(struct Camera *camera, struct vec3 *v, struct vec3 *pos, struct SD *sd){
	if (sd->dist == 1111.1) return newVec(0, 0, 0);
	struct vec3 *contact = vecadd(pos, vecmult(sd->dist, v));
	struct vec3 *lightnorm =  normalize(vecsubtract(camera->light, contact));
	struct vec3 *grad = sd->gradient(pos);
	double scale = dot3(grad, lightnorm);
	if (scale < 0) scale = 0;
	struct vec3 *shade = vecmult(scale, sd->color);
	free(contact);
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
	return sd;
}

void newNode(struct Node *head, struct Shape *shape){
	struct Node *cur = head;
	while (cur->next !=0){
		cur = cur->next;
	}
	cur->next = malloc(sizeof(struct Node));
	cur->next->shape = shape;
}

struct Shape *newShape(struct vec3 *color){
        struct Shape *newShape = malloc(sizeof(struct Shape));
        newShape->color = color;
        return newShape;
}


int Render(struct Node *head, struct Camera *camera, double hsteps, double vsteps, FILE *fp){
//int Render(double radius, struct vec3 *center, struct Camera *camera, double hsteps, double vsteps){
	/*struct vec3 *center = malloc(sizeof(struct vec3));
	center->x = 0.0;
	center->y = 0.0;
	center->z = 2.0;
	double radius = 2;*/

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
			//sdcur = calcSDF(head, v, camera->focus);
			//curcolor = shade(camera, v, camera->focus, sdcur);
			//if (j == 0) fprintf(fp, "%f, %f, %f", curcolor->x, curcolor->y, curcolor->z);
			//else fprintf(fp, ";%f, %f, %f", curcolor->x, curcolor->y, curcolor->z);
			//if (curcolor->y != 0) printf("o");
			//else printf(" ");
			if (head->shape->SDF(v,camera->focus) != 1111.1) printf("o");
                        else printf(".");
			//free(sheetpoint);
			//free(know);
			//free(v);
			//free(curcolor);
		}
		printf("\n");
		//fprintf(fp, ";\n");
		//fprintf(fp, "\n");
		//free(hnow);
	}	
	//fprintf(fp, "\n");
	printf("\n");
	free(h);
	free(k);
	return 0;
}


