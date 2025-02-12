use super::traits::{Bounded, Surface};
use super::{Point3D, Ray, Vec3D};

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
        let upper = self.center + &v;
        let lower = self.center - v;
        return AxisAlignedBoundingBox { upper, lower };
    }
}

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
        if denom == 0.0 {
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

#[derive(Debug)]
pub struct Triangle {
    pub p1: Point3D,
    pub p2: Point3D,
    pub p3: Point3D,
    d: f32, // distance from origo to plane along the normal direction
    normal: Vec3D,
}

impl Triangle {
    pub fn new(p1: Point3D, p2: Point3D, p3: Point3D) -> Self {
        let v12 = &p2 - &p1;
        let v23 = &p3 - &p2;
        let normal = Vec3D::cross(v12, v23).normalise();
        let d = -Vec3D::dot(p1, normal);
        Self {
            p1,
            p2,
            p3,
            d,
            normal,
        }
    }
}

impl Surface for Triangle {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)> {
        let denom = Vec3D::dot(self.normal, ray.direction);
        if denom == 0.0 {
            return None;
        }

        let t = -(Vec3D::dot(self.normal, ray.origin) + self.d) / denom;
        if t > t_max || t < t_min {
            return None;
        }

        let e1 = self.p2 - self.p1;
        let e2 = self.p3 - self.p2;
        let e3 = self.p1 - self.p3;

        let p = ray.at(t);
        let c1 = p - self.p1;
        let c2 = p - self.p2;
        let c3 = p - self.p3;

        if Vec3D::dot(self.normal, Vec3D::cross(e1, c1)) < 0.0
            || Vec3D::dot(self.normal, Vec3D::cross(e2, c2)) < 0.0
            || Vec3D::dot(self.normal, Vec3D::cross(e3, c3)) < 0.0
        {
            return None;
        } else {
            return Some((t, self.normal));
        }
    }
}

impl Bounded for Triangle {
    fn bounding_box(&self) -> AxisAlignedBoundingBox {
        let mut upper = Vec3D::ZERO;
        let mut lower = Vec3D::ZERO;
        for i in 0..3 {
            upper[i] = self.p1[i].max(self.p2[i]).max(self.p3[i]);
            lower[i] = self.p1[i].min(self.p2[i]).min(self.p3[i]);
        }
        return AxisAlignedBoundingBox { upper, lower };
    }
}

pub type AABB = AxisAlignedBoundingBox;

#[derive(Debug, PartialEq, Clone)]
pub struct AxisAlignedBoundingBox {
    pub upper: Point3D,
    pub lower: Point3D,
}

impl AxisAlignedBoundingBox {
    pub fn new(a: Point3D, b: Point3D) -> Self {
        let lower = Point3D::new(f32::min(a.x, b.x), f32::min(a.y, b.y), f32::min(a.z, b.z));
        let upper = Point3D::new(f32::max(a.x, b.x), f32::max(a.y, b.y), f32::max(a.z, b.z));
        return Self { lower, upper };
    }

    pub fn centroid(&self) -> Point3D {
        (self.upper + self.lower) / 2.0
    }

    pub fn widths(&self) -> [f32; 3] {
        let mut widths = [0.0; 3];
        for i in 0..3 {
            widths[i] = self.upper[i] - self.lower[i];
        }
        return widths;
    }

    pub fn combine(a: &Self, b: &Self) -> Self {
        let lower = Point3D::new(
            f32::min(a.lower.x, b.lower.x),
            f32::min(a.lower.y, b.lower.y),
            f32::min(a.lower.z, b.lower.z),
        );
        let upper = Point3D::new(
            f32::max(a.upper.x, b.upper.x),
            f32::max(a.upper.y, b.upper.y),
            f32::max(a.upper.z, b.upper.z),
        );
        return Self { lower, upper };
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
            let d_inv = 1.0 / ray.direction[axis];
            let mut t0: f32 = (self.lower[axis] - ray.origin[axis]) * d_inv;
            let mut t1: f32 = (self.upper[axis] - ray.origin[axis]) * d_inv;
            if d_inv < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }

            t_min = f32::max(t0, t_min);
            t_max = f32::min(t1, t_max);

            if t_max <= t_min {
                return None;
            }
        }
        return Some((t_min, Vec3D::ZERO)); // think about normal vector?
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(
            bounding_box.hit(&ray, 0.0, 10.0),
            Some((1.0, Vec3D::ZERO))
        );

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
