#include "triangle.h"

Hit Triangle::intersect(Ray const &ray) {
  double denominator = N.dot(ray.D);
  if (denominator == 0) {
    return Hit::NO_HIT();
  }

  double t = N.dot(v0 - ray.O) / denominator;
  if (t < 0) {
    return Hit::NO_HIT();
  }

  Point intersection = ray.at(t);

  double side0 = (v1 - v0).cross(intersection - v0).dot(N);
  double side1 = (v2 - v1).cross(intersection - v1).dot(N);
  double side2 = (v0 - v2).cross(intersection - v2).dot(N);

  if (side0 < 0 || side1 < 0 || side2 < 0) {
    return Hit::NO_HIT();
  }

  Vector hitNormal = N;
  if (hitNormal.dot(ray.D) > 0) {
    hitNormal = -hitNormal;
  }
  return Hit(t, hitNormal);
}

Triangle::Triangle(Point const &v0, Point const &v1, Point const &v2) : v0(v0), v1(v1), v2(v2), N() {
  Vector edge1 = v1 - v0;
  Vector edge2 = v2 - v0;
  N = edge1.cross(edge2).normalized();
}
