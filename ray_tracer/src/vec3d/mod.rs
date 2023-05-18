extern crate fastrand;
use fastrand::Rng;

mod ops;

#[derive(Debug, Copy, Clone)]
pub struct Vec3D {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vec3D {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self {x, y, z}
    }

    #[inline]
    pub fn norm(&self) -> f32 {
        self.norm_squared().sqrt()
    }

    #[inline]
    pub fn norm_squared(&self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    #[inline]
    pub fn normalise(&mut self) -> &Self {
        *self /= self.norm();
        self
    }

    #[inline]
    pub fn almost_zero(&self, abs: f32) -> bool {
        return self.norm_squared() < abs;
    }

    #[inline]
    pub fn almost_equal(vec1: &Self, vec2: &Self, abs: f32) -> bool {
        (vec1 - vec2).almost_zero(abs)
    }

    #[inline]
    pub fn dot(vec1: &Self, vec2: &Self) -> f32 {
        vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z
    }

    #[inline]
    pub fn cross(vec1: &Self, vec2: &Self) -> Self {
        Self::new(
            vec1.y*vec2.z - vec1.z*vec2.y,
            vec1.z*vec2.x - vec1.x*vec2.z,
            vec1.x*vec2.y - vec1.y*vec2.x,
        )
    }

    #[inline]
    pub fn interpolate(vec1: &Self, vec2: &Self, t: f32) -> Self {
        (1.0 - t) * vec1 + t * vec2
    }
}

impl Vec3D {
    #[inline]
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    pub fn ones() -> Self {
        Self::new(1.0, 1.0, 1.0)
    }

    #[inline]
    pub fn e1() -> Self {
        Self::new(1.0, 0.0, 0.0)
    }

    #[inline]
    pub fn e2() -> Self {
        Self::new(0.0, 1.0, 0.0)
    }

    #[inline]
    pub fn e3() -> Self {
        Self::new(0.0, 0.0, 1.0)
    }

    #[inline]
    pub fn unit_vector(direction: &Self) -> Self {
        direction.clone() / direction.norm()
    }

    #[inline]
    pub fn random_rejection(rng: &mut Rng) -> Self {
        loop {
            let (x, y, z) = (rng.f32()-0.5, rng.f32()-0.5, rng.f32()-0.5);
            let norm_sq = x*x + y*y + z*z;
            if norm_sq < 0.25 {
                let norm = norm_sq.sqrt();
                return Self::new(x / norm, y / norm, z / norm)
            }
        }
    }

    #[inline]
    pub fn random_direct(rng: &mut Rng) -> Self {
        let phi = 2.0 * std::f32::consts::PI * rng.f32();
        let cos = 1.0 - 2.0 * rng.f32();
        let sin = f32::sqrt(1.0 - cos*cos);
        let result = Self::new(phi.cos() * sin, phi.sin() * sin, cos);
        return result
    }
}

impl std::fmt::Display for Vec3D {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[{:?}, {:?}, {:?}]", self.x, self.y, self.z)
    }
}

impl Default for Vec3D {
    fn default() -> Self {
        Self::zero()
    }
}

impl std::cmp::PartialEq for Vec3D {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}
impl Eq for Vec3D {}

impl std::ops::Index<usize> for Vec3D {
    type Output = f32;
    fn index(&self, index: usize) -> &f32 {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Index {} out of range", index),
        }
    }
}

impl std::ops::IndexMut<usize> for Vec3D {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Index {} out of range", index),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construct_vector() {
        let vec = Vec3D::new(1.0, 2.0, 3.0);
        assert_eq!(vec.x, 1.0);
        assert_eq!(vec.y, 2.0);
        assert_eq!(vec.z, 3.0);
    }

    #[test]
    fn normalise() {
        let mut vec = Vec3D::new(1.0, 2.0, 3.0);
        let vec2: Vec3D = vec.normalise() + Vec3D::e1();
        assert!((vec.x - 0.2672612419124243).abs() < 1e-6);
        assert!((vec.y - 0.5345224838248487).abs() < 1e-6);
        assert!((vec.z - 0.8017837257372731).abs() < 1e-6);
        assert!((vec2.x - 1.2672612419124243).abs() < 1e-6);
    }

    #[test]
    fn unit_vector() {
        let vec = Vec3D::new(1.0, 2.0, 3.0);
        let u = Vec3D::unit_vector(&vec);
        assert!((u.x - 0.2672612419124243).abs() < 1e-6);
        assert!((u.y - 0.5345224838248487).abs() < 1e-6);
        assert!((u.z - 0.8017837257372731).abs() < 1e-6);
    }

    #[test]
    fn random_rejection() {
        let mut rng = fastrand::Rng::new();
        let vec = Vec3D::random_rejection(&mut rng);
        assert!((vec.norm() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn random_direct() {
        let mut rng = fastrand::Rng::new();
        let vec = Vec3D::random_direct(&mut rng);
        assert!((vec.norm() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn dot() {
        let vec1 = Vec3D::new(1.0, 2.0, 3.0);
        let vec2 = Vec3D::new(1.0, -1.0, 1.0);
        assert_eq!(Vec3D::dot(&vec1, &vec2), 2.0);
    }

    #[test]
    fn cross() {
        let vec1 = Vec3D::new(1.0, 2.0, 3.0);
        let vec2 = Vec3D::new(1.0, -1.0, 1.0);
        assert_eq!(Vec3D::cross(&vec1, &vec2), Vec3D::new(5.0, 2.0, -3.0));
        assert_eq!(Vec3D::cross(&vec2, &vec1), Vec3D::new(-5.0, -2.0, 3.0));
    }
}