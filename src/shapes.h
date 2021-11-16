#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double circleSDF(struct vec3 *v, struct vec3 *p);


struct vec3 *circleSDFGrad(struct vec3 *p);

double stemSDF(struct vec3 *v, struct vec3 *p);

struct vec3 *stemgrad(struct vec3 *p);

double boxSDF(struct vec3 *v, struct vec3 *p);

struct vec3 *boxSDFgrad(struct vec3 *p);


double planeSDF(struct vec3 *v, struct vec3 *p);

struct vec3 *planeSDFGrad(struct vec3 *p);

double leg1SDF(struct vec3 *v, struct vec3 *p);

struct vec3 *leg1SDFgrad(struct vec3 *p);

double leg2SDF(struct vec3 *v, struct vec3 *p);

struct vec3 *leg2SDFgrad(struct vec3 *p);

double leg3SDF(struct vec3 *v, struct vec3 *p);

struct vec3 *leg3SDFgrad(struct vec3 *p);

double leg4SDF(struct vec3 *v, struct vec3 *p);

struct vec3 *leg4SDFgrad(struct vec3 *p);

double roundBoxSDF(struct vec3 *v, struct vec3 *p);

struct vec3 *roundboxSDFgrad(struct vec3 *p);

double cylinder(struct vec3 *v, struct vec3 *p);

struct vec3 *cylinderGrad(struct vec3 *p);

struct vec3 *prismSDFgrad(struct vec3 *p);

double prismgenSDF(struct vec3 *v, struct vec3 *p);


