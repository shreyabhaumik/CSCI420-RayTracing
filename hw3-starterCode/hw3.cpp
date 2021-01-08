/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Shreya Bhaumik
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <iostream>
#include <cmath>
#include "myStructures.h"

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

#define PI 3.141592
#define ε 1e-13
#define MAX_Z -1e13
// for softshadow
#define MAX_POINTS 20


unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
bool isAntialias = false; // for Supersampling Anti-aliasing
bool isSoftshadow = false; // for Soft Shadows
// for recursive reflection
int totReflect = 0; 
double toRad = ((double)PI / 180.0);

/*---------------------------------------------------------Ray start---------------------------------------------------------*/
struct Ray
{
  Vec3 origin;
  Vec3 direction;

  Ray() {}
  Ray(const Vec3& o, const Vec3& d): origin(o), direction(d) {}  
  /*Overloaded function - for a sphere and a triangle
  to see if ray intersects any sphere/triangle - this function is to be executed per sphere/triangle. 
  If yes, returns true, puts point of intersection in hitpoint, 
  and distance of origin of ray till point of intersection in dist.
  If no, returns false.
  */
  bool checkSingleIntersect(const Sphere& s, Vec3& hitpoint) // for sphere
  {
    Vec3 do2c = (origin - Vec3(s.position[0], s.position[1], s.position[2])); // origin of ray - center of sphere
    double r = (double)s.radius;
    double a = direction.dot(direction);  //needing a because I wanted to compute roots differently
    double b = 2.0 * direction.dot(do2c);
    double c = do2c.dot(do2c) - r*r;
    double discriminant = b*b - 4.0*a*c;
    double t0, t1, t;
    t0 = t1 = t = MAX_Z;
    // if discriminant is 0 or very close to 0 (here abs(discriminant) <= 1e-13) we treat it like it has touched the sphere
    // so both roots will be equal. This is to take care of spheres with rough surfaces.
    if(discriminant < -ε)  return false; // (imaginary roots)
    double absD = std::abs(discriminant);
    if(absD < ε) 
    { // since ε is so close to 0, these two roots will actually be almost the same - 1 real root
      t0 = (- b + std::sqrt(absD)) * 0.5 / a; // (-b (+/-) sqrt(absD))/2a, but a = 1
      t1 = (- b - std::sqrt(absD)) * 0.5 / a;
      t = t0;
    }
    else // if discriminant > ε - there are 2 real roots now
    { // to avoid huge problem in computation if b and sqrt(discriminant) are very close in value, I won't use the regular formula
      // https://en.wikipedia.org/wiki/Loss_of_significance - see 'Instability of the quadratic equation'
      if(b < 0) // Can safely assume b!=0 because if it was then determinant would be -ve and we have already handled that case
        t0 = (- b + std::sqrt(discriminant)) * (double)0.5; // Trying to ensure that only numbers of the same sign are added
      else
        t0 = (- b - std::sqrt(discriminant)) * (double)0.5;
      t1 = c/t0;
      t0 = t0/a;
      t = t0;
    }
    if (t0 < 0 && t1 < 0) return false; // backside of ray - part beyond origin of ray is intersecting sphere
    if (t0 > t1 && t1 > 0) t = t1;
    else t = t0;
    hitpoint = origin + direction * t;
    return true;
  }

  bool checkSingleIntersect(const Triangle& tri, Vec3& hitpoint)  // for triangle
  { // the 3 vertices of triangle being tested
    Vec3 C0(tri.v[0].position[0], tri.v[0].position[1], tri.v[0].position[2]);
    Vec3 C1(tri.v[1].position[0], tri.v[1].position[1], tri.v[1].position[2]);
    Vec3 C2(tri.v[2].position[0], tri.v[2].position[1], tri.v[2].position[2]);
    Vec3 triN; // Normal of the triangle being tested
    triN.cross((C1 - C0),(C2 - C0)).Normalize();
    // If the ray is parallel or almost totally parallel (ε) to the plane of triangle, there's no intersection
    double triNdotDir = triN.dot(direction);
    if(std::abs(triNdotDir) < ε) return false; // [-ε, ε] then no intersection
    // calculating t - point of intersection
    double d = - (triN.dot(C0)); // finding value of d by putting value of a vertex(point) on the equation of the triangle(plane)
    double t = - (triN.dot(origin) + d) / triNdotDir; // Already ruled out the case where triNdotDir could be in [-ε, ε]
    // if t <= 0, the intersection is behind ray origin (or at ray origin i.e. on the camera when t=0 - which is too near)
    if(t <= 0) return false;
    hitpoint = origin + direction * t;
    /* Point-in-triangle testing(using concept of Barycentric coordinates)
    now that the intersection is found, need to check if it is inside or outside the triangle
    Forming 3 new triangles with 2 vertices of triangle being testing and the intersection point.
    For each of the new triangles, will check if it's normal is in same direction as the normal of the triangle being tested 
    - that we can see by checking sign of dot product of the two normals.
    If any one of the new triangles fail this test, the intersection point lies outside the triangle being tested.
    */
    Vec3 newTriN;
    double newTri0 = newTriN.cross((C1 - C0),(hitpoint - C0)).dot(triN);
    double newTri1 = newTriN.cross((C2 - C1),(hitpoint - C1)).dot(triN);
    double newTri2 = newTriN.cross((C0 - C2),(hitpoint - C2)).dot(triN);
    if(newTri0 < 0 || newTri1 < 0 || newTri2 < 0) return false;
    return true;
  }
};
/*----------------------------------------------------------Ray end----------------------------------------------------------*/

Ray * primaryRay(double pixelx, double pixely)
{
  double imageAspectRatio = (double)WIDTH/(double)HEIGHT;
  double tanhalfFOV = std::tan((fov * 0.5) * toRad);
  if(!isAntialias)
  {
    Ray* primary = new Ray[1];
    // normalized coordinates of the pixels defined in NDC (Normalized Device Coordinates) space
    double pixelNDCx = (pixelx+0.5)/(double)WIDTH; // +0.5 makes it ray pass through the center of the pixel
    double pixelNDCy = (pixely+0.5)/(double)HEIGHT;
    // Pixel coordinates expressed in NDC space are in the range [0,1]. But, film or image plane is centred around the world's origin.
    // So we need to map it to the range [-1,1].
    double pixelScreenx = 2.0*pixelNDCx - 1;
    // double pixelScreeny = 1 - 2.0*pixelNDCy; // because the [0,1] in y runs from up to down which should now be [1,-1] from up to down
    double pixelScreeny = 2.0*pixelNDCy - 1;
    // The upper calculations were done assuming that the image is square - so we have to account for aspect ratio
    // We are assuming that the image plane lies at z=-1. So from trig, the amount by which x and y is up or down from center is tan(FOV/2).
    double pixelCamerax = pixelScreenx * tanhalfFOV * imageAspectRatio; // x would change
    double pixelCameray = pixelScreeny * tanhalfFOV; // the y would remain unchanged
    // Now finally the point's x and y are in camera space. see slide 9(Ray Tracing) the final point is (pixelCamerax, pixelCameray, -1). z=-1 or -f doesn't matter.
    Vec3 direc(pixelCamerax, pixelCameray, -1.0);
    primary[0] = Ray(Vec3(0.0,0.0,0.0),direc.Normalize()); // normalizing direction as it is a direction
    return primary;
  }
  else  // Supersampling Anti-aliasing - 4 rays per pixel
  {
    Ray* primaries = new Ray[4];
    // 1st Ray
    double pixelNDCx = (pixelx+0.25)/(double)WIDTH;
    double pixelNDCy = (pixely+0.25)/(double)HEIGHT;
    double pixelScreenx = 2.0*pixelNDCx - 1;
    double pixelScreeny = 2.0*pixelNDCy - 1;
    double pixelCamerax = pixelScreenx * tanhalfFOV * imageAspectRatio;
    double pixelCameray = pixelScreeny * tanhalfFOV;
    Vec3 direc(pixelCamerax, pixelCameray, -1.0);
    primaries[0] = Ray(Vec3(0.0,0.0,0.0),direc.Normalize());
    // 2nd Ray
    pixelNDCx = (pixelx+0.75)/(double)WIDTH;
    pixelNDCy = (pixely+0.25)/(double)HEIGHT;
    pixelScreenx = 2.0*pixelNDCx - 1;
    pixelScreeny = 2.0*pixelNDCy - 1;
    pixelCamerax = pixelScreenx * tanhalfFOV * imageAspectRatio;
    pixelCameray = pixelScreeny * tanhalfFOV;
    direc = Vec3(pixelCamerax, pixelCameray, -1.0);
    primaries[1] = Ray(Vec3(0.0,0.0,0.0),direc.Normalize());
    // 3rd Ray
    pixelNDCx = (pixelx+0.75)/(double)WIDTH;
    pixelNDCy = (pixely+0.75)/(double)HEIGHT;
    pixelScreenx = 2.0*pixelNDCx - 1;
    pixelScreeny = 2.0*pixelNDCy - 1;
    pixelCamerax = pixelScreenx * tanhalfFOV * imageAspectRatio;
    pixelCameray = pixelScreeny * tanhalfFOV;
    direc = Vec3(pixelCamerax, pixelCameray, -1.0);
    primaries[2] = Ray(Vec3(0.0,0.0,0.0),direc.Normalize());
    // 4th Ray
    pixelNDCx = (pixelx+0.25)/(double)WIDTH;
    pixelNDCy = (pixely+0.75)/(double)HEIGHT;
    pixelScreenx = 2.0*pixelNDCx - 1;
    pixelScreeny = 2.0*pixelNDCy - 1;
    pixelCamerax = pixelScreenx * tanhalfFOV * imageAspectRatio;
    pixelCameray = pixelScreeny * tanhalfFOV;
    direc = Vec3(pixelCamerax, pixelCameray, -1.0);
    primaries[3] = Ray(Vec3(0.0,0.0,0.0),direc.Normalize());

    return primaries;
  }  
}

//Phong shading - overloaded functions for spheres and triangles
void CLAMP(double& d)
{
  if(d > 1.0) d = 1.0;
  else if(d < 0.0) d = 0.0;
}
Color PhongShading(const Sphere& s, const Light& l, Vec3& intersection) // for spheres
{
  // Diffuse part - L.N
  Vec3 l_Pos(l.position[0], l.position[1], l.position[2]);
  Vec3 L = (l_Pos - intersection).Normalize();  // Light direction vector - L
  Vec3 N = (intersection - Vec3(s.position[0], s.position[1], s.position[2])).Normalize();  // Normal - N
  double LdotN = L.dot(N);  CLAMP(LdotN);
  // Specular part - R.V
  Vec3 R = (N * (2.0 * (LdotN)) - L).Normalize(); // Reflection direction - R = 2(l.n)n - l
  Vec3 V = (-intersection).Normalize();  // Eye direction - V (0 - intersection, since eye/camera is at 0)
  double RdotV = R.dot(V);  CLAMP(RdotV);
  // Phong Equation Calculation
  Color kd(s.color_diffuse[0], s.color_diffuse[1], s.color_diffuse[2]);
  Color ks(s.color_specular[0], s.color_specular[1], s.color_specular[2]);
  double spec = std::pow(RdotV, s.shininess);
  Color c;
  c.r = l.color[0] * ((kd.r * LdotN) + (ks.r * spec));
  c.g = l.color[1] * ((kd.g * LdotN) + (ks.g * spec));
  c.b = l.color[2] * ((kd.b * LdotN) + (ks.b * spec));
  return c;
}

Color PhongShading(const Triangle& t, const Light& l, Vec3& intersection) // for triangles
{
  Vec3 C0(t.v[0].position[0], t.v[0].position[1], t.v[0].position[2]);
  Vec3 C1(t.v[1].position[0], t.v[1].position[1], t.v[1].position[2]);
  Vec3 C2(t.v[2].position[0], t.v[2].position[1], t.v[2].position[2]);
  Vec3 C = intersection;
  // to find Barycentric Coordinates of triangle
  Vec3 crossProd;
  double AreaC0C1C2 = crossProd.cross((C1 - C0),(C2 - C0)).len(); // dividing area by 2 is unnecessary as we will ultimately divide 2 areas to get alpha, beta or gamma
  double AreaCC1C2 = crossProd.cross((C1 - C),(C2 - C)).len();
  double alpha = AreaCC1C2 / AreaC0C1C2; // alpha
  double AreaC0CC2 = crossProd.cross((C - C0),(C2 - C0)).len();
  double beta = AreaC0CC2 / AreaC0C1C2; // beta
  double gamma = 1.0 - alpha - beta; // gamma
  // interpolating normal with Barycentric Coordinates of triangle
  Vec3 N = Vec3((alpha * t.v[0].normal[0]) + (beta * t.v[1].normal[0]) + (gamma * t.v[2].normal[0]),
                (alpha * t.v[0].normal[1]) + (beta * t.v[1].normal[1]) + (gamma * t.v[2].normal[1]),
                (alpha * t.v[0].normal[2]) + (beta * t.v[1].normal[2]) + (gamma * t.v[2].normal[2])).Normalize();
  // Diffuse part - L.N
  Vec3 l_Pos(l.position[0], l.position[1], l.position[2]);
  Vec3 L = (l_Pos - intersection).Normalize();  // Light direction vector - L
  double LdotN = L.dot(N);  CLAMP(LdotN);
  // Specular part - R.V
  Vec3 R = (N * (2.0 * (LdotN)) - L).Normalize(); // Reflection direction - R = 2(l.n)n - l
  Vec3 V = (-intersection).Normalize();  // Eye direction - V (0 - intersection, since eye/camera is at 0)
  double RdotV = R.dot(V);  CLAMP(RdotV);
  // Phong Equation Calculation
  Color kd((alpha * t.v[0].color_diffuse[0]) + (beta * t.v[1].color_diffuse[0]) + (gamma * t.v[2].color_diffuse[0]),
            (alpha * t.v[0].color_diffuse[1]) + (beta * t.v[1].color_diffuse[1]) + (gamma * t.v[2].color_diffuse[1]),
            (alpha * t.v[0].color_diffuse[2]) + (beta * t.v[1].color_diffuse[2]) + (gamma * t.v[2].color_diffuse[2]));
  Color ks((alpha * t.v[0].color_specular[0]) + (beta * t.v[1].color_specular[0]) + (gamma * t.v[2].color_specular[0]),
            (alpha * t.v[0].color_specular[1]) + (beta * t.v[1].color_specular[1]) + (gamma * t.v[2].color_specular[1]),
            (alpha * t.v[0].color_specular[2]) + (beta * t.v[1].color_specular[2]) + (gamma * t.v[2].color_specular[2]));
  double shininess = (alpha * t.v[0].shininess) + (beta * t.v[1].shininess) + (gamma * t.v[2].shininess);
  double spec = std::pow(RdotV, shininess);
  Color c;
  c.r = l.color[0] * ((kd.r * LdotN) + (ks.r * spec));
  c.g = l.color[1] * ((kd.g * LdotN) + (ks.g * spec));
  c.b = l.color[2] * ((kd.b * LdotN) + (ks.b * spec));
  return c;
}

Color checkAllSphereIntersect(Ray& r, const Color& c, Vec3& intersect, int& s_ID, double& nearest) // spheres
{
  Color pixel_c = c;
  s_ID = -1;
  for(int i = 0; i < num_spheres; i++) // checking all spheres
  {
    Vec3 intersection = Vec3(0.0, 0.0, (double)MAX_Z);
    if(r.checkSingleIntersect(spheres[i], intersection) && (intersection.z > nearest))
    {
      pixel_c = Color();
      s_ID = i;
      for(int j = 0; j < num_lights; j++)
      {
        // creating shadow ray
        Vec3 l_Pos(lights[j].position[0], lights[j].position[1], lights[j].position[2]);
        Vec3 shadowRayDir = (l_Pos - intersection).Normalize(); // 'intersection' is the origin for the shadow ray
        Ray shadow(intersection, shadowRayDir); // shadow ray

        bool isInShadow = false; // to see if intersection point is in shadow of other object(s)
        // checking shadow ray's intersection against all other objects
        for(int k = 0; k < num_spheres; k++) // checking all spheres
        {
          Vec3 hitpoint;
          if(shadow.checkSingleIntersect(spheres[k], hitpoint) && (k != i)) // if shadow ray hits a sphere and it is not the sphere that made the initial intersection
          {
            Vec3 shadowOriginToHitpoint = hitpoint - intersection;  // 'intersection' is the origin for the shadow ray
            Vec3 shadowOriginToLight = l_Pos - intersection;
            if(shadowOriginToLight.len() > shadowOriginToHitpoint.len())
            {
              isInShadow = true;
              break;
            }
          }
        }
        for(int k = 0; k < num_triangles; k++) // checking all triangles
        {
          Vec3 hitpoint;
          if(shadow.checkSingleIntersect(triangles[k], hitpoint)) // if shadow ray hits a triangle
          {
            Vec3 shadowOriginToHitpoint = hitpoint - intersection;  // 'intersection' is the origin for the shadow ray
            Vec3 shadowOriginToLight = l_Pos - intersection;
            if(shadowOriginToLight.len() > shadowOriginToHitpoint.len())
            {
              isInShadow = true;
              break;
            }
          }
        }
        if(!isInShadow) // if point is not in shadow, it will show it's color
          pixel_c += PhongShading(spheres[i], lights[j], intersection);
      }
      nearest = intersection.z;
      intersect = intersection;
    }
  }
  return pixel_c;
}

Color checkAllTriangleIntersect(Ray& r, const Color& c, Vec3& intersect, int& t_ID, double& nearest) // triangles
{
  Color pixel_c = c;
  t_ID = -1;
  for(int i = 0; i < num_triangles; i++) // checking all spheres
  {
    Vec3 intersection = Vec3(0.0, 0.0, (double)MAX_Z);
    double dist = (double)MAX_Z;
    if(r.checkSingleIntersect(triangles[i], intersection) && (intersection.z > nearest))
    {
      pixel_c = Color();
      t_ID = i;
      for(int j = 0; j < num_lights; j++)
      {
        // creating shadow ray
        Vec3 l_Pos(lights[j].position[0], lights[j].position[1], lights[j].position[2]);
        Vec3 shadowRayDir = (l_Pos - intersection).Normalize(); // 'intersection' is the origin for the shadow ray
        Ray shadow(intersection, shadowRayDir); // shadow ray

        bool isInShadow = false; // to see if intersection point is in shadow of other object(s)
        // checking shadow ray's intersection against all other objects
        for(int k = 0; k < num_spheres; k++) // checking all spheres
        {
          Vec3 hitpoint;
          if(shadow.checkSingleIntersect(spheres[k], hitpoint)) // if shadow ray hits a sphere
          {
            Vec3 shadowOriginToHitpoint = hitpoint - intersection;  // 'intersection' is the origin for the shadow ray
            Vec3 shadowOriginToLight = l_Pos - intersection;
            if(shadowOriginToLight.len() > shadowOriginToHitpoint.len())
            {
              isInShadow = true;
              break;
            }
          }
        }
        for(int k = 0; k < num_triangles; k++) // checking all triangles
        {
          Vec3 hitpoint;
          if(shadow.checkSingleIntersect(triangles[k], hitpoint) && (k != i)) // if shadow ray hits a triangle and it is not the triangle that made the initial intersection
          {
            Vec3 shadowOriginToHitpoint = hitpoint - intersection;  // 'intersection' is the origin for the shadow ray
            Vec3 shadowOriginToLight = l_Pos - intersection;
            if(shadowOriginToLight.len() > shadowOriginToHitpoint.len())
            {
              isInShadow = true;
              break;
            }
          }
        }
        if(!isInShadow) // if point is not in shadow, it will show it's color
          pixel_c += PhongShading(triangles[i], lights[j], intersection);
      }
      nearest = intersection.z;
      intersect = intersection;
    }
  }
  return pixel_c;
}

double MAX_POINTSinv = 1.0 / (double)MAX_POINTS;
Light * pointToSphereLights(Light L)  // to convert a point light into a spherical light source with multiple point lights
{
  Light * spherePoints = new Light[MAX_POINTS];
  for(int i=0; i<MAX_POINTS; i++)
  {
    // to find a point in a sphere
    // reference used: https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
    double distFromCenter = ((double)rand()/(double)RAND_MAX);  // distFromCenter is clamped to between [0,1]
    double phi = (double)(rand() % 360) * toRad;
    double theta = (double)(rand() % 360) * toRad;
    double sinPhi = std::sin(phi);  double cosPhi = std::sqrt(1.0 - sinPhi*sinPhi);
    double sinTheta = std::sin(theta);  double cosTheta = std::sqrt(1.0 - sinTheta*sinTheta);
    double xDisp = distFromCenter * sinPhi * cosTheta;
    double yDisp = distFromCenter * sinPhi * sinTheta;
    double zDisp = distFromCenter * cosPhi;
    // to calculate the position on the point wrt the world
    spherePoints[i].position[0] = L.position[0] + xDisp;
    spherePoints[i].position[1] = L.position[1] + yDisp;
    spherePoints[i].position[2] = L.position[2] + zDisp;
    // to reduce contribution by all new lights
    spherePoints[i].color[0] = L.color[0] * MAX_POINTSinv;
    spherePoints[i].color[1] = L.color[1] * MAX_POINTSinv;
    spherePoints[i].color[2] = L.color[2] * MAX_POINTSinv;
  }
  return spherePoints;
}

Color checkAllSphereIntersectSoftShadow(Ray& r, const Color& c, Vec3& intersect, int& s_ID, double& nearest) // spheres
{
  Color pixel_c = c;
  s_ID = -1;
  for(int i = 0; i < num_spheres; i++) // checking all spheres
  {
    Vec3 intersection = Vec3(0.0, 0.0, (double)MAX_Z);
    if(r.checkSingleIntersect(spheres[i], intersection) && (intersection.z > nearest))
    {
      pixel_c = Color();
      s_ID = i;
      for(int j = 0; j < num_lights; j++)
      {
        Light * subLights = pointToSphereLights(lights[j]);
        for(int a = 0; a < MAX_POINTS; a++) // all sublights from one light
        {
          // creating shadow ray
          Vec3 l_Pos(subLights[a].position[0], subLights[a].position[1], subLights[a].position[2]);
          Vec3 shadowRayDir = (l_Pos - intersection).Normalize(); // 'intersection' is the origin for the shadow ray
          Ray shadow(intersection, shadowRayDir); // shadow ray

          bool isInShadow = false; // to see if intersection point is in shadow of other object(s)
          // checking shadow ray's intersection against all other objects
          for(int k = 0; k < num_spheres; k++) // checking all spheres
          {
            Vec3 hitpoint;
            if(shadow.checkSingleIntersect(spheres[k], hitpoint) && (k != i)) // if shadow ray hits a sphere and it is not the sphere that made the initial intersection
            {
              Vec3 shadowOriginToHitpoint = hitpoint - intersection;  // 'intersection' is the origin for the shadow ray
              Vec3 shadowOriginToLight = l_Pos - intersection;
              if(shadowOriginToLight.len() > shadowOriginToHitpoint.len())
              {
                isInShadow = true;
                break;
              }
            }
          }
          for(int k = 0; k < num_triangles; k++) // checking all triangles
          {
            Vec3 hitpoint;
            if(shadow.checkSingleIntersect(triangles[k], hitpoint)) // if shadow ray hits a triangle
            {
              Vec3 shadowOriginToHitpoint = hitpoint - intersection;  // 'intersection' is the origin for the shadow ray
              Vec3 shadowOriginToLight = l_Pos - intersection;
              if(shadowOriginToLight.len() > shadowOriginToHitpoint.len())
              {
                isInShadow = true;
                break;
              }
            }
          }
          if(!isInShadow) // if point is not in shadow, it will show it's color
            pixel_c += PhongShading(spheres[i], subLights[a], intersection);
        }
      }
      nearest = intersection.z;
      intersect = intersection;
    }
  }
  return pixel_c;
}

Color checkAllTriangleIntersectSoftShadow(Ray& r, const Color& c, Vec3& intersect, int& t_ID, double& nearest) // triangles
{
  Color pixel_c = c;
  t_ID = -1;
  for(int i = 0; i < num_triangles; i++) // checking all spheres
  {
    Vec3 intersection = Vec3(0.0, 0.0, (double)MAX_Z);
    double dist = (double)MAX_Z;
    if(r.checkSingleIntersect(triangles[i], intersection) && (intersection.z > nearest))
    {
      pixel_c = Color();
      t_ID = i;
      for(int j = 0; j < num_lights; j++)
      {
        Light * subLights = pointToSphereLights(lights[j]);
        for(int a = 0; a < MAX_POINTS; a++) // all sublights from one light
        {
          // creating shadow ray
          Vec3 l_Pos(subLights[a].position[0], subLights[a].position[1], subLights[a].position[2]);
          Vec3 shadowRayDir = (l_Pos - intersection).Normalize(); // 'intersection' is the origin for the shadow ray
          Ray shadow(intersection, shadowRayDir); // shadow ray

          bool isInShadow = false; // to see if intersection point is in shadow of other object(s)
          // checking shadow ray's intersection against all other objects
          for(int k = 0; k < num_spheres; k++) // checking all spheres
          {
            Vec3 hitpoint;
            if(shadow.checkSingleIntersect(spheres[k], hitpoint)) // if shadow ray hits a sphere
            {
              Vec3 shadowOriginToHitpoint = hitpoint - intersection;  // 'intersection' is the origin for the shadow ray
              Vec3 shadowOriginToLight = l_Pos - intersection;
              if(shadowOriginToLight.len() > shadowOriginToHitpoint.len())
              {
                isInShadow = true;
                break;
              }
            }
          }
          for(int k = 0; k < num_triangles; k++) // checking all triangles
          {
            Vec3 hitpoint;
            if(shadow.checkSingleIntersect(triangles[k], hitpoint) && (k != i)) // if shadow ray hits a triangle and it is not the triangle that made the initial intersection
            {
              Vec3 shadowOriginToHitpoint = hitpoint - intersection;  // 'intersection' is the origin for the shadow ray
              Vec3 shadowOriginToLight = l_Pos - intersection;
              if(shadowOriginToLight.len() > shadowOriginToHitpoint.len())
              {
                isInShadow = true;
                break;
              }
            }
          }
          if(!isInShadow) // if point is not in shadow, it will show it's color
            pixel_c += PhongShading(triangles[i], subLights[a], intersection);
        }
      }
      nearest = intersection.z;
      intersect = intersection;
    }
  }
  return pixel_c;
}

Color recursiveComputeColor(Ray& r, int countReflect) // this handles the primary ray as well as recursing for reflection
{
  if(countReflect > totReflect) return Color(0.0, 0.0, 0.0);
  double R = 0.0, G = 0.0, B = 0.0;
  double nearest = (double)MAX_Z;
  Vec3 intersect;
  Color fromLightOrShadow = Color(1.0, 1.0, 1.0);
  Color fromReflect = Color();
  int s_ID, t_ID;
  if(isSoftshadow)
  {
    fromLightOrShadow = checkAllSphereIntersectSoftShadow(r, fromLightOrShadow, intersect, s_ID, nearest);
    fromLightOrShadow = checkAllTriangleIntersectSoftShadow(r, fromLightOrShadow, intersect, t_ID, nearest);
  }
  else
  {
    fromLightOrShadow = checkAllSphereIntersect(r, fromLightOrShadow, intersect, s_ID, nearest);
    fromLightOrShadow = checkAllTriangleIntersect(r, fromLightOrShadow, intersect, t_ID, nearest);
  }
  Color ks;
  Ray reflect;
  // reflection
  if(s_ID != -1) // means we have an intersection with a sphere
  {
    Sphere s = spheres[s_ID];
    Vec3 N = (intersect - Vec3(s.position[0], s.position[1], s.position[2])).Normalize();  // Normal - N
    Vec3 L = -r.direction;
    double LdotN = L.dot(N);  CLAMP(LdotN);
    Vec3 R = (N * (2.0 * (LdotN)) - L).Normalize(); // Reflection direction - R = 2(l.n)n - l 
    Vec3 origin = intersect + (R * ε); // (R * ε) so that we don't repeat same
    reflect = Ray(origin, R.Normalize());
    ks = Color(s.color_specular[0], s.color_specular[1], s.color_specular[2]);
  }
  else if(t_ID != -1) // means we have an intersection with a triangle
  {
    Triangle t = triangles[t_ID];
    Vec3 C0(t.v[0].position[0], t.v[0].position[1], t.v[0].position[2]);
    Vec3 C1(t.v[1].position[0], t.v[1].position[1], t.v[1].position[2]);
    Vec3 C2(t.v[2].position[0], t.v[2].position[1], t.v[2].position[2]);
    Vec3 C = intersect;
    // to find Barycentric Coordinates of triangle
    Vec3 crossProd;
    double AreaC0C1C2 = crossProd.cross((C1 - C0),(C2 - C0)).len(); // dividing area by 2 is unnecessary as we will ultimately divide 2 areas to get alpha, beta or gamma
    double AreaCC1C2 = crossProd.cross((C1 - C),(C2 - C)).len();
    double alpha = AreaCC1C2 / AreaC0C1C2; // alpha
    double AreaC0CC2 = crossProd.cross((C - C0),(C2 - C0)).len();
    double beta = AreaC0CC2 / AreaC0C1C2; // beta
    double gamma = 1.0 - alpha - beta; // gamma
    // interpolating normal with Barycentric Coordinates of triangle
    Vec3 N = Vec3((alpha * t.v[0].normal[0]) + (beta * t.v[1].normal[0]) + (gamma * t.v[2].normal[0]),
                    (alpha * t.v[0].normal[1]) + (beta * t.v[1].normal[1]) + (gamma * t.v[2].normal[1]),
                    (alpha * t.v[0].normal[2]) + (beta * t.v[1].normal[2]) + (gamma * t.v[2].normal[2])).Normalize();
    Vec3 L = -r.direction;
    double LdotN = L.dot(N);  CLAMP(LdotN);
    Vec3 R = (N * (2.0 * (LdotN)) - L).Normalize(); // Reflection direction - R = 2(l.n)n - l 
    Vec3 origin = intersect + (R * ε); // (R * ε) so that we don't repeat same
    reflect = Ray(origin, R.Normalize());
    ks = Color((alpha * t.v[0].color_specular[0]) + (beta * t.v[1].color_specular[0]) + (gamma * t.v[2].color_specular[0]),
            (alpha * t.v[0].color_specular[1]) + (beta * t.v[1].color_specular[1]) + (gamma * t.v[2].color_specular[1]),
            (alpha * t.v[0].color_specular[2]) + (beta * t.v[1].color_specular[2]) + (gamma * t.v[2].color_specular[2]));
  }
  double fact = 0.05;
  if(totReflect == 0) // only the primary ray - no reflections
  {
    if(s_ID == -1 && t_ID == -1)  // misses all objects
      return Color(1.0, 1.0, 1.0);
    else  // hits an object
    {
      fromReflect = Color();
      R = (1.0 - fact) * fromLightOrShadow.r + fact * fromReflect.r;
      G = (1.0 - fact) * fromLightOrShadow.g + fact * fromReflect.g;
      B = (1.0 - fact) * fromLightOrShadow.b + fact * fromReflect.b;
      return Color(R,G,B) + recursiveComputeColor(reflect, countReflect + 1)*(fact*2.0);
    }
  }
  else  // when there's the primary ray as well as reflections
  {
    if(countReflect == 0) // the primary ray
    {
      if(s_ID == -1 && t_ID == -1)  // misses all objects
        return Color(1.0, 1.0, 1.0);
      else  // hits an object
      {
        fromReflect = recursiveComputeColor(reflect, countReflect);
        R = (1.0 - fact) * fromLightOrShadow.r + fact * fromReflect.r;
        G = (1.0 - fact) * fromLightOrShadow.g + fact * fromReflect.g;
        B = (1.0 - fact) * fromLightOrShadow.b + fact * fromReflect.b;
        return Color(R,G,B) + recursiveComputeColor(reflect, countReflect + 1)*(fact*2.0);
      }
    }
    else // all the reflection rays
    {
      if(s_ID == -1 && t_ID == -1)  // misses all objects
        return Color(0.0, 0.0, 0.0);
      else  // hits an object
      {
        fromReflect = recursiveComputeColor(reflect, countReflect + 1);
        R = (1.0 - ks.r) * fromLightOrShadow.r + ks.r * fromReflect.r;
        G = (1.0 - ks.g) * fromLightOrShadow.g + ks.g * fromReflect.g;
        B = (1.0 - ks.b) * fromLightOrShadow.b + ks.b * fromReflect.b;
        return Color(R,G,B)*(1.0 - fact) + recursiveComputeColor(reflect, countReflect + 2)*(fact);
      }
    }
  }
}

Color handlePrimaryRay_s(int x, int y) // to handle antialiasing per pixel enabled or disabled for primary rays
{
  Color finCol = Color();
  double r = 0.0, g = 0.0, b = 0.0;
  if(isAntialias)
  {
    Ray * primaries = primaryRay(x,y);
    for(int i=0; i<4; i++)
    {
      Color c = recursiveComputeColor(primaries[i], 0);
      r += c.r;
      g += c.g;
      b += c.b;
    }
    r /= 4.0; // averaging over the color computed by the 4 rays for a pixel
    g /= 4.0;
    b /= 4.0;
    finCol = Color(r,g,b);
  }
  else
  {
    Ray * primary = primaryRay(x,y);
    finCol = recursiveComputeColor(primary[0], 0);
  }
  return finCol;
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//MODIFY THIS FUNCTION
void draw_scene()
{
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      // plot_pixel(x, y, x % 256, y % 256, (x+y) % 256);
      Color finCol = handlePrimaryRay_s(x,y);
      finCol += Color(ambient_light[0], ambient_light[1], ambient_light[2]); // finally adding the ambient color
      finCol.colClamp();
      plot_pixel(x, y, finCol.r * 255, finCol.g * 255, finCol.b * 255);
    }
    glEnd();
    glFlush();
  }
  printf("Ray Tracing Completed!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

void keyboardFunc(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 27: // ESC key
      exit(0); // exit the program
      break;

    default:
      break;
  }
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 6))
  {  
    printf ("Usage: %s <input scenefile> [y/n - to enable/disable SSAA Anti-aliasing] [y/n - to enable/disable Soft Shadows] [numberOfReflectionsInInteger] [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 6)
  {
    mode = MODE_JPEG;
    filename = argv[5];

    totReflect = atoi(argv[4]);

    char op = (argv[3])[0];
    if(op == 121 || op == 89) // ASCII for y - 121, Y - 89
      isSoftshadow = true;

    op = (argv[2])[0];
    if(op == 121 || op == 89) // ASCII for y - 121, Y - 89
      isAntialias = true;
  }
  else if(argc == 5) // mode = MODE_DISPLAY; by default
  {
    totReflect = atoi(argv[4]);

    char op = (argv[3])[0];
    if(op == 121 || op == 89) // ASCII for y - 121, Y - 89
      isSoftshadow = true;

    op = (argv[2])[0];
    if(op == 121 || op == 89) // ASCII for y - 121, Y - 89
      isAntialias = true;
  }
  else if(argc == 4) // mode = MODE_DISPLAY; by default
  {
    char op = (argv[3])[0];
    if(op == 121 || op == 89) // ASCII for y - 121, Y - 89
      isSoftshadow = true;

    op = (argv[2])[0];
    if(op == 121 || op == 89) // ASCII for y - 121, Y - 89
      isAntialias = true;
  }
  else if(argc == 3) // mode = MODE_DISPLAY; by default
  {
    char op = (argv[2])[0];
    if(op == 121 || op == 89) // ASCII for y - 121, Y - 89
      isAntialias = true;
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);
  if(isAntialias)
    std::cout << "Supersampling Anti-aliasing enabled." << std::endl;
  else
    std::cout << "Supersampling Anti-aliasing disabled." << std::endl;
  if(isSoftshadow)
    std::cout << "Soft Shadows enabled." << std::endl;
  else
    std::cout << "Soft Shadows disabled." << std::endl;
  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  glutKeyboardFunc(keyboardFunc); // callback for pressing the keys on the keyboard
  init();
  glutMainLoop();
}
