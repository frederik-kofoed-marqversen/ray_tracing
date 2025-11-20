use crate::Vec3D;
use crate::Ray;

pub struct Camera {
    origin: Vec3D,
    direction: Vec3D,
    aspect_ratio: f32,
    horizontal: Vec3D,
    vertical: Vec3D,
    lower_left_corner: Vec3D,
}

impl Camera {
    pub fn new(origin: Vec3D, direction: Vec3D, aspect_ratio: f32) -> Self {
        let mut result = Self {
            origin,
            direction,
            aspect_ratio,
            horizontal: Vec3D::default(),
            vertical: Vec3D::default(),
            lower_left_corner: Vec3D::default(),
        };
        result.set_direction(direction);
        result
    }

    pub fn get_direction(&self) -> Vec3D {
        self.direction
    }

    pub fn get_pos(&self) -> Vec3D {
        self.origin
    }

    pub fn set_pos(&mut self, point: Vec3D) {
        let translation = point - self.origin;
        self.origin = point;
        self.lower_left_corner += translation;
    }

    pub fn set_direction(&mut self, direction: Vec3D) {
        let window_height: f32 = 2.0;
        let window_width: f32 = window_height * self.aspect_ratio;
        let focal_length: f32 = 2.0; // Horizontal FOV of 46 degrees

        let direction = direction.normalise();

        self.direction = direction;
        self.horizontal = Vec3D::cross(direction, Vec3D::Z).normalise() * window_width;
        self.vertical = Vec3D::cross(self.horizontal, direction).normalise() * window_height;
        self.lower_left_corner =
            self.origin + direction * focal_length - self.horizontal / 2.0 - self.vertical / 2.0;
    }

    pub fn look_at(&mut self, point: Vec3D) {
        self.set_direction(point - self.origin);
    }

    #[inline]
    pub fn ray(&self, u: f32, v: f32) -> Ray {
        let direction =
            self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin;
        Ray::new(self.origin, direction.normalise())
    }
}