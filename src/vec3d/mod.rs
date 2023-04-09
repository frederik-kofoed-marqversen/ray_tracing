use super::{Rng, StandardNormal};
mod ops;

#[derive(Debug, Copy, Clone)]
pub struct Vec3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3D {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self {x, y, z}
    }

    #[inline]
    pub fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }

    #[inline]
    pub fn norm_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    #[inline]
    pub fn normalise(&mut self) -> &Self {
        *self /= self.norm();
        self
    }

    #[inline]
    pub fn almost_zero(&self, abs: f64) -> bool {
        return self.norm_squared() < abs;
    }

    #[inline]
    pub fn dot(vec1: &Self, vec2: &Self) -> f64 {
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
    pub fn interpolate(vec1: &Self, vec2: &Self, t: f64) -> Self {
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
    pub fn random(length: f64, rng: &mut impl Rng) -> Self {
        let v = Self::new(
            rng.sample(StandardNormal), 
            rng.sample(StandardNormal), 
            rng.sample(StandardNormal)
        );
        v * (length / v.norm())
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

#[cfg(test)]
mod tests {
    extern crate rand;
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
        assert!((vec.x - 0.2672612419124243).abs() < 1e-15);
        assert!((vec.y - 0.5345224838248487).abs() < 1e-15);
        assert!((vec.z - 0.8017837257372731).abs() < 1e-15);
        assert!((vec2.x - 1.2672612419124243).abs() < 1e-15);
    }

    #[test]
    fn unit_vector() {
        let vec = Vec3D::new(1.0, 2.0, 3.0);
        let u = Vec3D::unit_vector(&vec);
        assert!((u.x - 0.2672612419124243).abs() < 1e-15);
        assert!((u.y - 0.5345224838248487).abs() < 1e-15);
        assert!((u.z - 0.8017837257372731).abs() < 1e-15);
    }

    #[test]
    fn random_unit_vector() {
        let mut rng = rand::thread_rng();
        let vec = Vec3D::random(2.0, &mut rng);
        assert!((vec.norm() - 2.0).abs() < 1e-15);
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