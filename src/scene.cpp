#include "scene.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "hit.h"
#include "image.h"
#include "material.h"
#include "ray.h"

using namespace std;

pair<ObjectPtr, Hit> Scene::castRay(Ray const &ray) const {
  // Find hit object and distance
  Hit min_hit(numeric_limits<double>::infinity(), Vector());
  ObjectPtr obj = nullptr;
  for (unsigned idx = 0; idx != objects.size(); ++idx) {
    Hit hit(objects[idx]->intersect(ray));
    if (hit.t < min_hit.t) {
      min_hit = hit;
      obj = objects[idx];
    }
  }

  return pair<ObjectPtr, Hit>(obj, min_hit);
}

Color Scene::trace(Ray const &ray, unsigned depth) {
  pair<ObjectPtr, Hit> mainhit = castRay(ray);
  ObjectPtr obj = mainhit.first;
  Hit min_hit = mainhit.second;

  // No hit? Return background color.
  if (!obj) return Color(0.0, 0.0, 0.0);

  Material const &material = obj->material;
  Point hit = ray.at(min_hit.t);
  Vector V = -ray.D;

  // Pre-condition: For closed objects, N points outwards.
  Vector N = min_hit.N;

  // The shading normal always points in the direction of the view,
  // as required by the Phong illumination model.
  Vector shadingN;
  if (N.dot(V) >= 0.0)
    shadingN = N;
  else
    shadingN = -N;

  Color matColor = material.color;

  // Add ambient once, regardless of the number of lights.
  Color color = material.ka * matColor;

  // Add diffuse and specular components.
  for (auto const &light : lights) {
    Vector L = (light->position - hit).normalized();

    // render shadows if enabled, offset to avaoid shadow acne
    if (renderShadows) {
      Vector offsetNormal = N.dot(L) >= 0.0 ? N : -N;
      Point shadowOrigin = hit + epsilon * offsetNormal;
      Vector shadowDirection = (light->position - shadowOrigin).normalized();
      pair<ObjectPtr, Hit> shadowhit = castRay(Ray(shadowOrigin, shadowDirection));
      double lightDistance = (light->position - shadowOrigin).length();

      if (shadowhit.first && shadowhit.second.t < lightDistance) continue;
    }

    // Add diffuse.
    double dotNormal = shadingN.dot(L);
    double diffuse = std::max(dotNormal, 0.0);
    color += diffuse * material.kd * light->color * matColor;

    // Add specular.
    if (dotNormal > 0) {
      Vector reflectDir = reflect(-L, shadingN);
      double specAngle = std::max(reflectDir.dot(V), 0.0);
      double specular = std::pow(specAngle, material.n);

      color += specular * material.ks * light->color;
    }
  }

  if (depth > 0 && material.isTransparent) {
    // The object is transparent, and thus refracts and reflects light.

    // use Snell's law to compute the refraction direction and trace it recursively.
    bool entering = N.dot(ray.D) < 0.0;
    Vector refractionNormal = entering ? N : -N;
    double n1 = entering ? 1.0 : material.nt;
    double n2 = entering ? material.nt : 1.0;
    double cosTheta = max(0.0, min(1.0, -ray.D.dot(refractionNormal)));
    double eta = n1 / n2;
    double sin2ThetaT = eta * eta * (1.0 - cosTheta * cosTheta);

    Vector reflectionDirection = reflect(ray.D, refractionNormal).normalized();
    Vector reflectionOffsetNormal = N.dot(reflectionDirection) >= 0.0 ? N : -N;
    Point reflectionOrigin = hit + epsilon * reflectionOffsetNormal;
    Color reflectedColor = trace(Ray(reflectionOrigin, reflectionDirection), depth - 1);

    double reflectance = 1.0;
    if (sin2ThetaT < 1.0) {
      double cosThetaT = sqrt(1.0 - sin2ThetaT);
      double r0 = (n1 - n2) / (n1 + n2);
      r0 *= r0;
      // schlick's approximation for reflectance
      reflectance = r0 + (1.0 - r0) * pow(1.0 - cosTheta, 5.0);
      
      Vector refractionDirection = (eta * ray.D + (eta * cosTheta - cosThetaT) * refractionNormal).normalized();
      Vector refractionOffsetNormal = N.dot(refractionDirection) >= 0.0 ? N : -N;
      Point refractionOrigin = hit + epsilon * refractionOffsetNormal;
      color += (1.0 - reflectance) * trace(Ray(refractionOrigin, refractionDirection), depth - 1);
    }

    color += reflectance * reflectedColor;
  } else if (depth > 0 && material.ks > 0.0) {
    // The object is not transparent, but opaque.
    // Reflect the ray and trace it recursively.
    Vector reflectionDirection = reflect(ray.D, shadingN);
    Vector offsetNormal = N.dot(reflectionDirection) >= 0.0 ? N : -N;
    Point reflectionOrigin = hit + epsilon * offsetNormal;
    color += material.ks * trace(Ray(reflectionOrigin, reflectionDirection), depth - 1);
  }

  return color;
}

void Scene::render(Image &img) {
  unsigned w = img.width();
  unsigned h = img.height();

  for (unsigned y = 0; y < h; ++y)
    for (unsigned x = 0; x < w; ++x) {
      Point pixel(x + 0.5, h - 1 - y + 0.5, 0);
      Ray ray(eye, (pixel - eye).normalized());
      Color col = trace(ray, recursionDepth);
      col.clamp();
      img(x, y) = col;
    }
}

// --- Misc functions ----------------------------------------------------------

// Defaults
Scene::Scene() : objects(), lights(), eye(), renderShadows(false), recursionDepth(0), supersamplingFactor(1) {}

void Scene::addObject(ObjectPtr obj) { objects.push_back(obj); }

void Scene::addLight(Light const &light) { lights.push_back(LightPtr(new Light(light))); }

void Scene::setEye(Triple const &position) { eye = position; }

unsigned Scene::getNumObject() { return objects.size(); }

unsigned Scene::getNumLights() { return lights.size(); }

void Scene::setRenderShadows(bool shadows) { renderShadows = shadows; }

void Scene::setRecursionDepth(unsigned depth) { recursionDepth = depth; }

void Scene::setSuperSample(unsigned factor) { supersamplingFactor = factor; }
