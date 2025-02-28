use super::traits::{Bounded, Surface};
use super::{Point3D, Ray, Vec3D};

#[derive(Debug)]
pub struct Sphere {
    pub center: Point3D,
    pub radius: f32,
}

impl Sphere {
    pub fn new(center: Point3D, radius: f32) -> Self {
        Self { center, radius }
    }

    #[inline]
    fn normal(&self, point: Point3D) -> Vec3D {
        return (point - self.center) / self.radius;
    }
}

impl Surface for Sphere {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)> {
        let r_oc: Vec3D = ray.origin - self.center;
        let a = ray.direction.norm_squared();
        let half_b = Vec3D::dot(ray.direction, r_oc);
        let c = r_oc.norm_squared() - self.radius * self.radius;
        let d = half_b * half_b - a * c;

        if d < 0.0 {
            return None;
        }
        let sqrt_d = d.sqrt();

        let mut t = -(half_b + sqrt_d) / a;
        if t < t_max && t > t_min {
            return Some((t, self.normal(ray.at(t))));
        }

        t = (-half_b + sqrt_d) / a;
        if t < t_max && t > t_min {
            return Some((t, self.normal(ray.at(t))));
        }

        return None;
    }
}

impl Bounded for Sphere {
    fn bounding_box(&self) -> AxisAlignedBoundingBox {
        let v = Vec3D::new(self.radius, self.radius, self.radius);
        let upper = self.center + v;
        let lower = self.center - v;
        return AxisAlignedBoundingBox { upper, lower };
    }
}

#[derive(Debug)]
pub struct Plane {
    d: f32, // distance from origo to plane along the normal direction
    normal: Vec3D,
}

impl Plane {
    pub fn new(p: Point3D, normal: Vec3D) -> Self {
        let normal = normal.normalise();
        let d = -Vec3D::dot(p, normal);
        Self { d, normal }
    }
}

impl Surface for Plane {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)> {
        let denom = Vec3D::dot(self.normal, ray.direction);
        if denom.abs() < f32::EPSILON {
            return None;
        }

        let t = -(Vec3D::dot(self.normal, ray.origin) + self.d) / denom;
        if t < t_max && t > t_min {
            return Some((t, self.normal));
        } else {
            return None;
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct AxisAlignedBoundingBox {
    pub upper: Point3D,
    pub lower: Point3D,
}

impl AxisAlignedBoundingBox {
    #[inline]
    pub fn new(a: Point3D, b: Point3D) -> Self {
        let upper = Point3D::max(a, b);
        let lower = Point3D::min(a, b);
        return Self { upper, lower };
    }

    #[inline]
    pub fn empty() -> Self {
        Self {
            upper: Vec3D::fill(f32::NEG_INFINITY),
            lower: Vec3D::fill(f32::INFINITY),
        }
    }

    #[inline]
    pub fn centroid(&self) -> Point3D {
        (self.upper + self.lower) / 2.0
    }

    #[inline]
    pub fn area(&self) -> f32 {
        self.half_area() * 2.0
    }

    #[inline]
    pub fn half_area(&self) -> f32 {
        let w = self.upper - self.lower; // widths
        return w.x * w.y + w.y * w.z + w.z * w.x;
    }

    #[inline]
    pub fn combine(a: &Self, b: &Self) -> Self {
        let upper = Point3D::max(a.upper, b.upper);
        let lower = Point3D::min(a.lower, b.lower);
        return Self { upper, lower };
    }

    #[inline]
    pub fn grow(&self, p: Point3D) -> Self {
        Self {
            upper: Point3D::max(self.upper, p),
            lower: Point3D::min(self.lower, p),
        }
    }
}

impl Bounded for AxisAlignedBoundingBox {
    fn bounding_box(&self) -> AxisAlignedBoundingBox {
        return self.clone();
    }
}

impl Surface for AxisAlignedBoundingBox {
    fn hit(&self, ray: &Ray, mut t_min: f32, mut t_max: f32) -> Option<(f32, Vec3D)> {
        for axis in 0..3 {
            let d_inv = ray.recip_direction[axis];
            let mut t0: f32 = (self.lower[axis] - ray.origin[axis]) * d_inv;
            let mut t1: f32 = (self.upper[axis] - ray.origin[axis]) * d_inv;
            if d_inv < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }

            t_min = t_min.max(t0);
            t_max = t_max.min(t1);

            if t_max <= t_min {
                return None;
            }
        }
        return Some((t_min, Vec3D::ZERO)); // Should never need normal vector of AABB
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type AABB = AxisAlignedBoundingBox;

    #[test]
    fn hit_sphere() {
        let sphere = Sphere::new(Point3D::new(1.0, 1.0, 1.0), 1.0);
        let ray = Ray::new(Point3D::ZERO, Vec3D::new(1.0, 1.0, 1.0) / f32::sqrt(3.0));
        let (t, normal) = sphere.hit(&ray, 0.0, f32::INFINITY).unwrap();
        assert!((f32::sqrt(3.0) - 1.0 - t).abs() < 1e-6);
        assert!((normal.norm() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn hit_sphere_tangent() {
        let sphere = Sphere::new(Point3D::new(1.0, 1.0, 1.0), 1.0);
        let ray = Ray::new(Point3D::ZERO, Vec3D::new(1.0, 1.0, 0.0) / f32::sqrt(2.0));
        let (t, normal) = sphere.hit(&ray, 0.0, f32::INFINITY).unwrap();
        assert!((f32::sqrt(2.0) - t).abs() < 1e-6);
        assert!((normal.norm() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn miss_sphere() {
        let sphere = Sphere::new(Point3D::new(1.0, 1.0, 1.0), 1.0);
        let ray = Ray::new(Point3D::ZERO, Vec3D::Z);
        let result = sphere.hit(&ray, 0.0, f32::INFINITY);
        assert_eq!(result, None);
    }

    #[test]
    fn hit_bounding_box() {
        let bounding_box = AABB::new(Vec3D::X, Vec3D::new(2.0, 2.0, 2.0));

        let ray = Ray::new(Vec3D::ZERO, Vec3D::ONES);
        assert_eq!(bounding_box.hit(&ray, 0.0, 10.0), Some((1.0, Vec3D::ZERO)));

        let ray = Ray::new(Vec3D::ZERO, -Vec3D::ONES);
        assert_eq!(bounding_box.hit(&ray, 0.0, 10.0), None);
    }

    #[test]
    fn combine_bounding_boxes() {
        let box1 = AABB::new(Vec3D::X, Vec3D::new(2.0, 2.0, 2.0));
        let box2 = AABB::new(Vec3D::ZERO, Vec3D::new(3.0, 0.5, 0.5));
        let result = AABB::combine(&box1, &box2);
        let expected = AABB::new(Vec3D::ZERO, Vec3D::new(3.0, 2.0, 2.0));

        assert_eq!(result, expected);
    }
}
