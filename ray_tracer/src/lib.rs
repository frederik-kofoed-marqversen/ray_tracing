#[macro_use]
extern crate impl_ops;

mod vec3d;
pub use vec3d::Vec3D;
pub type Point3D = Vec3D;

pub mod traits;
pub mod primitives;
pub mod materials;
pub mod stl;
pub mod accelerators;
pub mod engine;

#[derive(Debug, Copy, Clone)]
pub struct Ray {
    pub origin: Point3D,
    pub direction: Vec3D,
}

impl Ray {
    pub fn new(origin: Point3D, direction: Vec3D) -> Self {
        Self{origin, direction}
    }
    
    pub fn at(&self, t: f32) -> Point3D {
        self.origin + self.direction * t
    }
}

pub struct Object {
    pub surface: Box<dyn traits::Surface>,
    pub material: materials::Material
}

impl Object {
    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D, &Object)> {
        return match self.surface.hit(ray, t_min, t_max) {
            Some((t, normal)) => {
                Some((t, normal, &self))
            },
            None => None,
        }
    }
}

pub struct Camera {
    origin: Point3D,
    horizontal: Vec3D,
    vertical: Vec3D,
    lower_left_corner: Point3D,
}

impl Camera {
    pub fn new(origin: Point3D, mut direction: Vec3D, aspect_ratio: f32) -> Self {
        let window_height: f32 = 2.0;
        let window_width: f32 = window_height * aspect_ratio;
        let focal_length: f32 = 2.0;  // Horizontal FOV of 46 degrees

        direction.normalise();
        let horizontal = Vec3D::cross(
            &direction, 
            &Vec3D::e3()
        ).normalise() * window_width;
        let vertical = Vec3D::cross(
            &horizontal,
            &direction
        ).normalise() * window_height;
        let lower_left_corner: Point3D = &origin + direction*focal_length - &horizontal / 2.0 - &vertical / 2.0;

        Self {origin, horizontal, vertical, lower_left_corner}
    }

    pub fn ray(&self, u: f32, v: f32) -> Ray {
        Ray::new(
            self.origin.clone(), 
            self.lower_left_corner + u*self.horizontal + v*self.vertical - self.origin
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ray_at() {
        let point = Point3D::new(1.0, 1.0, 1.0);
        let vec = Vec3D::new(0.0, 1.0, 2.0);
        let ray = Ray::new(point, vec);
        let point = ray.at(2.0);
        assert_eq!(point.x, 1.0);
        assert_eq!(point.y, 3.0);
        assert_eq!(point.z, 5.0);
    }
}