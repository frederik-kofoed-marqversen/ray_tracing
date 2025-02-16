mod vec3d;
pub use vec3d::Vec3D;
pub type Point3D = Vec3D;

pub mod accelerators;
pub mod engine;
pub mod materials;
pub mod primitives;
pub mod stl;
pub mod traits;
pub mod subdivision_surface;
pub mod triangle_mesh;

#[derive(Debug, Copy, Clone)]
pub struct Ray {
    pub origin: Point3D,
    pub direction: Vec3D,
}

impl Ray {
    pub fn new(origin: Point3D, direction: Vec3D) -> Self {
        Self { origin, direction }
    }

    pub fn at(&self, t: f32) -> Point3D {
        self.origin + self.direction * t
    }
}

pub struct Object {
    pub surface: Box<dyn traits::Surface>,
    pub material: materials::Material,
}

impl Object {
    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D, &Object)> {
        return match self.surface.hit(ray, t_min, t_max) {
            Some((t, normal)) => Some((t, normal, &self)),
            None => None,
        };
    }
}

pub struct Camera {
    pub origin: Point3D,
    pub direction: Vec3D,
    aspect_ratio: f32,
    horizontal: Vec3D,
    vertical: Vec3D,
    lower_left_corner: Point3D,
}

impl Camera {
    pub fn new(origin: Point3D, direction: Vec3D, aspect_ratio: f32) -> Self {
        let mut result = Self {
            origin,
            direction,
            aspect_ratio,
            horizontal: Vec3D::default(),
            vertical: Vec3D::default(),
            lower_left_corner: Point3D::default(),
        };
        result.set_direction(direction);
        result
    }

    pub fn set_direction(&mut self, direction: Point3D) {
        let window_height: f32 = 2.0;
        let window_width: f32 = window_height * self.aspect_ratio;
        let focal_length: f32 = 2.0; // Horizontal FOV of 46 degrees

        let direction = direction.normalise();

        self.direction = direction.clone();
        self.horizontal = Vec3D::cross(direction, Vec3D::Z).normalise() * window_width;
        self.vertical = Vec3D::cross(self.horizontal, direction).normalise() * window_height;
        self.lower_left_corner =
            self.origin + direction * focal_length - self.horizontal / 2.0 - self.vertical / 2.0;
    }

    pub fn look_at(&mut self, point: Point3D) {
        self.set_direction(point - self.origin);
    }

    pub fn ray(&self, u: f32, v: f32) -> Ray {
        let direction = self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin;
        Ray::new(
            self.origin.clone(),
            direction / direction.norm(),
        )
    }
}
