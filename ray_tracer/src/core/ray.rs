use crate::Vec3D;

#[derive(Debug, Copy, Clone)]
pub struct Ray {
    pub origin: Vec3D,
    pub direction: Vec3D,
}

impl Ray {
    #[inline]
    pub fn new(origin: Vec3D, direction: Vec3D) -> Self {
        Self { origin, direction }
    }

    #[inline]
    pub fn at(&self, t: f32) -> Vec3D {
        self.origin + self.direction * t
    }
}