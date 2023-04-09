use super::{Point3D, Vec3D, Ray};

pub struct Camera {
    origin: Point3D,
    horizontal: Vec3D,
    vertical: Vec3D,
    lower_left_corner: Point3D,
}

impl Camera {
    pub fn new(origin: Point3D, mut direction: Vec3D, aspect_ratio: f64) -> Self {
        let window_height: f64 = 2.0;
        let window_width: f64 = window_height * aspect_ratio;
        let focal_length: f64 = 2.0;  // Horizontal FOV of 46 degrees

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

    pub fn ray(&self, u: f64, v: f64) -> Ray {
        Ray::new(
            self.origin.clone(), 
            self.lower_left_corner + u*self.horizontal + v*self.vertical - self.origin
        )
    }
}