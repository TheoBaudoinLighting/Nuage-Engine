   #pragma once
   
   #include <glm/glm.hpp>
   #include <memory>
   
   #include "nei_aabb.h"
   #include "nei_material.h"
   
   struct Ray;
   
   struct Hittable;
   typedef std::shared_ptr<Hittable> HittablePtr;
   
   struct HitRecord {
       float m_T; // Intersection time
       glm::vec3 m_Point; // Intersection point
       glm::vec3 m_Normal; // Normal at the intersection
       MaterialPtr m_Material; // Material of the intersected object
       HittablePtr m_Object; // Object of the intersection
   };
   
   struct Hittable {
       // Check if the ray hits the object
       virtual bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& record) const = 0;
       
       // Get the bounding box of the object
       virtual bool BoundingBox(AABB& output_box) const = 0;
       
       // Get the material of the object
       virtual MaterialPtr GetMaterial() const { return nullptr; } 
       
       // Get the area of the object
       virtual float Area() const {
           return 0.0f;
       }
   };