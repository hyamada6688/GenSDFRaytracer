#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct Camera{
	struct vec3 *focus;
	struct vec3 *corner[4];
	struct vec3 *light;
};

struct Node{
	struct Node *next;
	struct Shape *shape;
};


//typedef struct vec3 *(*r3r3)(struct Sphere *shape, struct vec3 *v, struct vec3 *pos);

//typedef double (*r3r1)(struct Sphere *shape, struct vec3 *pos);

typedef struct vec3 *(*r3r3)(struct vec3 *);

typedef double (*r3r1)(struct vec3 *, struct vec3 *);


struct Shape {
        struct vec3 *color;
	r3r3 gradSDF;
	r3r1 SDF;
};


struct Sphere {
	struct Shape base;
	double radius;
	struct vec3 *center;
	r3r1 SDF;
	struct vec3 *color;
};

struct vec3{
	double x;
	double y;
	double z;
};

struct Plane {
	struct Shape base;
	struct vec3 *normal;
	struct vec3 *pointOn;
};


struct SD {
	double dist;
	r3r3 gradient;
	struct vec3 *color;
};

struct Camera *newCam(struct vec3 *focus, struct vec3 *corner[4], struct vec3 *light);

struct vec3 *newVec(double x, double y, double z);

double dot3(struct vec3 *v1, struct vec3 *v2);

struct vec3 *vecsubtract(struct vec3 *v1, struct vec3 *v2);

struct vec3 *vecmult(double scalar, struct vec3 *vector);

struct vec3 *normalize(struct vec3 *OG);

double genSDFsphere(struct vec3 *center, double r, struct vec3 *v, struct vec3 *p);

double genSDFplane(struct vec3 *normal, struct vec3 *pointOn, struct vec3 *v, struct vec3 *p);

struct vec3 *shade(struct Camera *camera, struct vec3 *v, struct vec3 *pos, struct SD *sd);

struct SD *calcSDF(struct Node *head, struct vec3 *v, struct vec3 *p);

void newNode(struct Node *head, struct Shape *shape);

struct Plane *newPlane(struct vec3 *normal, struct vec3 *pointon);

struct Sphere *newSphere(double radius, struct vec3 *center);

struct Shape *newShape(struct vec3 *color);
