#ifndef MESH_H_
#define MESH_H_

#include <string>
#include <vector>

#include "../object.h"

class Mesh : public Object {
  struct TriangleUV {
    Point v0;
    Point v1;
    Point v2;
    Vector uv0;
    Vector uv1;
    Vector uv2;
  };

  std::vector<ObjectPtr> d_tris;
  std::vector<TriangleUV> d_uvTriangles;

 public:
  Mesh(std::string const &filename, Point const &position, Vector const &rotation, Vector const &scale);

  Hit intersect(Ray const &ray) override;
  Vector toUV(Point const &hit) override;
  void tangentFrame(Point const &hit, Vector const &normal, Vector &tangentU, Vector &tangentV) const override;
};

#endif
