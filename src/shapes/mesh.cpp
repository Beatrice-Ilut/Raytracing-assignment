#include "mesh.h"

#include <cmath>
#include <iostream>
#include <limits>

#include "../objloader.h"
#include "../vertex.h"
#include "triangle.h"

using namespace std;

Hit Mesh::intersect(Ray const &ray) {
  Hit closestHit(numeric_limits<double>::infinity(), Vector());

  // for all triangles in the mesh check if the ray intersect a tringle 
  // the first intersections is the closest
  for (ObjectPtr const &triangle : d_tris) {
    Hit hit = triangle->intersect(ray);
    if (hit.t >= 0 && hit.t < closestHit.t) {
      closestHit = hit;
    }
  }

   if (std::isinf(closestHit.t)) {
    return Hit::NO_HIT();
  }

  return closestHit;
}

Vector Mesh::toUV(Point const &hit) {
  for (TriangleUV const &triangle : d_uvTriangles) {
    Vector const edge1 = triangle.v1 - triangle.v0;
    Vector const edge2 = triangle.v2 - triangle.v0;
    Vector const offset = hit - triangle.v0;
    Vector const normal = edge1.cross(edge2).normalized();
    if (abs(offset.dot(normal)) > 1e-3) continue;

    double const d00 = edge1.dot(edge1);
    double const d01 = edge1.dot(edge2);
    double const d11 = edge2.dot(edge2);
    double const d20 = offset.dot(edge1);
    double const d21 = offset.dot(edge2);
    double const denominator = d00 * d11 - d01 * d01;
    if (abs(denominator) < numeric_limits<double>::epsilon()) continue;

    double const u = (d11 * d20 - d01 * d21) / denominator;
    double const v = (d00 * d21 - d01 * d20) / denominator;
    if (u < -1e-8 || v < -1e-8 || u + v > 1.0 + 1e-8) continue;

    return (1.0 - u - v) * triangle.uv0 + u * triangle.uv1 + v * triangle.uv2;
  }

  return Vector();
}

void Mesh::tangentFrame(Point const &hit, Vector const &normal, Vector &tangentU, Vector &tangentV) const {
  for (TriangleUV const &triangle : d_uvTriangles) {
    Vector const edge1 = triangle.v1 - triangle.v0;
    Vector const edge2 = triangle.v2 - triangle.v0;
    Vector const offset = hit - triangle.v0;
    Vector const triangleNormal = edge1.cross(edge2).normalized();
    if (abs(offset.dot(triangleNormal)) > 1e-3) continue;

    double const d00 = edge1.dot(edge1);
    double const d01 = edge1.dot(edge2);
    double const d11 = edge2.dot(edge2);
    double const d20 = offset.dot(edge1);
    double const d21 = offset.dot(edge2);
    double const denominator = d00 * d11 - d01 * d01;
    if (abs(denominator) < numeric_limits<double>::epsilon()) continue;

    double const baryU = (d11 * d20 - d01 * d21) / denominator;
    double const baryV = (d00 * d21 - d01 * d20) / denominator;
    if (baryU < -1e-8 || baryV < -1e-8 || baryU + baryV > 1.0 + 1e-8) continue;

    double const uv1U = triangle.uv1.x - triangle.uv0.x;
    double const uv1V = triangle.uv1.y - triangle.uv0.y;
    double const uv2U = triangle.uv2.x - triangle.uv0.x;
    double const uv2V = triangle.uv2.y - triangle.uv0.y;
    double const uvDeterminant = uv1U * uv2V - uv1V * uv2U;
    if (abs(uvDeterminant) < numeric_limits<double>::epsilon()) break;

    tangentU = (uv2V * edge1 - uv1V * edge2) / uvDeterminant;
    tangentV = (-uv2U * edge1 + uv1U * edge2) / uvDeterminant;
    if (normal.dot(triangleNormal) < 0.0) tangentV = -tangentV;
    return;
  }

  Object::tangentFrame(hit, normal, tangentU, tangentV);
}

Mesh::Mesh(string const &filename, Point const &position, Vector const &rotation, Vector const &scale) {
  OBJLoader model(filename);
  d_tris.reserve(model.numTriangles());
  vector<Vertex> vertices = model.vertex_data();

  double cosX = cos(rotation.x);
  double sinX = sin(rotation.x);
  double cosY = cos(rotation.y);
  double sinY = sin(rotation.y);
  double cosZ = cos(rotation.z);
  double sinZ = sin(rotation.z);

  auto transform = [&](Point const &vertex) {
    //scaling
    Point point(vertex.x * scale.x, vertex.y * scale.y, vertex.z * scale.z);

    //rotation
    double rotatedY = point.y * cosX - point.z * sinX;
    double rotatedZ = point.y * sinX + point.z * cosX;
    point.y = rotatedY;
    point.z = rotatedZ;

    double rotatedX = point.x * cosY + point.z * sinY;
    rotatedZ = -point.x * sinY + point.z * cosY;
    point.x = rotatedX;
    point.z = rotatedZ;

    rotatedX = point.x * cosZ - point.y * sinZ;
    rotatedY = point.x * sinZ + point.y * cosZ;
    point.x = rotatedX;
    point.y = rotatedY;

    //translation
    return point + position;
  };
  // Create triangles from the vertices and apply transformations
  for (size_t tri = 0; tri != model.numTriangles(); ++tri) {
    Vertex one = vertices[tri * 3];
    Point v0(one.x, one.y, one.z);

    Vertex two = vertices[tri * 3 + 1];
    Point v1(two.x, two.y, two.z);

    Vertex three = vertices[tri * 3 + 2];
    Point v2(three.x, three.y, three.z);

    v0 = transform(v0);
    v1 = transform(v1);
    v2 = transform(v2);

    d_tris.push_back(ObjectPtr(new Triangle(v0, v1, v2)));
    d_uvTriangles.push_back(TriangleUV{v0, v1, v2, Vector(one.u, one.v, 0), Vector(two.u, two.v, 0),
                                       Vector(three.u, three.v, 0)});
  }

  cout << "Loaded model: " << filename << " with " << model.numTriangles() << " triangles.\n";
}
