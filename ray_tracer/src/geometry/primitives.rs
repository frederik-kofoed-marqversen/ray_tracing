use std::f32::consts::PI;

use crate::traits::*;
use crate::Ray;
use crate::Vec3D;

#[derive(Debug)]
pub struct Point {
    pub position: Vec3D,
}

impl Point {
    #[inline]
    pub fn new(position: Vec3D) -> Self {
        Self { position }
    }
}

impl Surface for Point {
    fn hit(&self, _ray: &Ray, _t_min: f32, _t_max: f32) -> Option<(f32, Vec3D)> {
        None
    }
}

#[derive(Debug)]
pub struct Sphere {
    pub center: Vec3D,
    pub radius: f32,
}

impl Sphere {
    #[inline]
    pub fn new(center: Vec3D, radius: f32) -> Self {
        Self { center, radius }
    }

    #[inline]
    pub fn area(&self) -> f32 {
        4.0 * PI * self.radius * self.radius
    }

    #[inline]
    fn normal(&self, point: Vec3D) -> Vec3D {
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

    fn centroid(&self) -> Vec3D {
        self.center
    }
}

#[derive(Debug, Clone)]
pub struct Triangle {
    pub v0: Vec3D,
    pub v1: Vec3D,
    pub v2: Vec3D,
}

impl Triangle {
    #[inline]
    pub fn new(v0: Vec3D, v1: Vec3D, v2: Vec3D) -> Self {
        Self { v0, v1, v2 }
    }

    #[inline]
    pub fn area_weighted_normal(&self) -> Vec3D {
        let e1 = self.v1 - self.v0;
        let e2 = self.v2 - self.v0;
        0.5 * Vec3D::cross(e1, e2)
    }

    #[inline]
    pub fn area(&self) -> f32 {
        self.area_weighted_normal().norm()
    }
}

#[inline]
pub fn triangle_intersection(
    v0: Vec3D,
    v1: Vec3D,
    v2: Vec3D,
    ray: &Ray,
    t_min: f32,
    t_max: f32,
) -> Option<(f32, Vec3D)> {
    // Möller-Trumbore intersection algorithm. More or less copied from Wikipedia
    let e1 = v1 - v0;
    let e2 = v2 - v0;

    let ray_cross_e2 = Vec3D::cross(ray.direction, e2);
    let det = Vec3D::dot(e1, ray_cross_e2);

    if det.abs() < f32::EPSILON {
        return None; // Ray is parallel to triangle.
    }

    let inv_det = 1.0 / det;
    let s = ray.origin - v0;
    let u = inv_det * Vec3D::dot(s, ray_cross_e2);
    if u < 0.0 || u > 1.0 {
        // Ray misses triangle in u-dir
        return None;
    }

    let s_cross_e1 = Vec3D::cross(s, e1);
    let v = inv_det * Vec3D::dot(ray.direction, s_cross_e1);
    if v < 0.0 || u + v > 1.0 {
        // Ray misses triangle in v-dir
        return None;
    }

    let t = inv_det * Vec3D::dot(e2, s_cross_e1);
    if t > t_min && t < t_max {
        // Hit!
        let normal = Vec3D::cross(e1, e2).normalise();
        return Some((t, normal));
    } else {
        // Ray hits outside of given interval
        return None;
    }
}

impl Surface for Triangle {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)> {
        triangle_intersection(self.v0, self.v1, self.v2, ray, t_min, t_max)
    }
}

impl Bounded for Triangle {
    #[inline]
    fn bounding_box(&self) -> AxisAlignedBoundingBox {
        AxisAlignedBoundingBox::empty()
            .grow(self.v0)
            .grow(self.v1)
            .grow(self.v2)
    }

    #[inline]
    fn centroid_axis(&self, axis: usize) -> f32 {
        (self.v0[axis] + self.v1[axis] + self.v2[axis]) / 3.0
    }

    #[inline]
    fn centroid(&self) -> Vec3D {
        (self.v0 + self.v1 + self.v2) / 3.0
    }
}

#[derive(Debug)]
pub struct Plane {
    d: f32, // distance from origo to plane along the normal direction
    normal: Vec3D,
}

impl Plane {
    #[inline]
    pub fn new(p: Vec3D, normal: Vec3D) -> Self {
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
    pub upper: Vec3D,
    pub lower: Vec3D,
}

impl AxisAlignedBoundingBox {
    #[inline]
    pub fn new(a: Vec3D, b: Vec3D) -> Self {
        let upper = Vec3D::max(a, b);
        let lower = Vec3D::min(a, b);
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
    pub fn centroid(&self) -> Vec3D {
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
        let upper = Vec3D::max(a.upper, b.upper);
        let lower = Vec3D::min(a.lower, b.lower);
        return Self { upper, lower };
    }

    #[inline]
    pub fn grow(&self, p: Vec3D) -> Self {
        Self {
            upper: Vec3D::max(self.upper, p),
            lower: Vec3D::min(self.lower, p),
        }
    }
}

impl AxisAlignedBoundingBox {
    pub fn hit(&self, ray: &Ray, mut t_min: f32, mut t_max: f32) -> Option<f32> {
        for axis in 0..3 {
            let d_inv = ray.direction[axis].recip();
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
        return Some(t_min);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type AABB = AxisAlignedBoundingBox;

    #[test]
    fn hit_sphere() {
        let sphere = Sphere::new(Vec3D::new(1.0, 1.0, 1.0), 1.0);
        let ray = Ray::new(Vec3D::ZERO, Vec3D::new(1.0, 1.0, 1.0) / f32::sqrt(3.0));
        let (t, normal) = sphere.hit(&ray, 0.0, f32::INFINITY).unwrap();
        assert!((f32::sqrt(3.0) - 1.0 - t).abs() < 1e-6);
        assert!((normal.norm() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn hit_sphere_tangent() {
        let sphere = Sphere::new(Vec3D::new(1.0, 1.0, 1.0), 1.0);
        let ray = Ray::new(Vec3D::ZERO, Vec3D::new(1.0, 1.0, 0.0) / f32::sqrt(2.0));
        let (t, normal) = sphere.hit(&ray, 0.0, f32::INFINITY).unwrap();
        assert!((f32::sqrt(2.0) - t).abs() < 1e-6);
        assert!((normal.norm() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn miss_sphere() {
        let sphere = Sphere::new(Vec3D::new(1.0, 1.0, 1.0), 1.0);
        let ray = Ray::new(Vec3D::ZERO, Vec3D::Z);
        let result = sphere.hit(&ray, 0.0, f32::INFINITY);
        assert_eq!(result, None);
    }

    #[test]
    fn hit_bounding_box() {
        let bounding_box = AABB::new(Vec3D::X, Vec3D::new(2.0, 2.0, 2.0));

        let ray = Ray::new(Vec3D::ZERO, Vec3D::ONES);
        assert_eq!(bounding_box.hit(&ray, 0.0, 10.0), Some(1.0));

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
