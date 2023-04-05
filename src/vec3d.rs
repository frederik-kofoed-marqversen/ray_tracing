use std::ops;

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

    pub fn normalise(&mut self) {
        *self /= self.norm();
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
}

impl std::fmt::Display for Vec3D {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[{:?}, {:?}, {:?}]", self.x, self.y, self.z)
    }
}

impl std::cmp::PartialEq for Vec3D {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}
impl Eq for Vec3D {}

impl ops::Add<Vec3D> for Vec3D {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self::new(
            self.x + other.x, 
            self.y + other.y, 
            self.z + other.z
        )
    }
}

impl ops::AddAssign<Vec3D> for Vec3D {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl ops::Sub<Vec3D> for Vec3D {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self {
        Self::new(
            self.x - other.x, 
            self.y - other.y, 
            self.z - other.z
        )
    }
}

impl ops::SubAssign<Vec3D> for Vec3D {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl ops::Neg for Vec3D {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl ops::Mul<f64> for Vec3D {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: f64) -> Self {
        Self::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl ops::MulAssign<f64> for Vec3D {
    #[inline]
    fn mul_assign(&mut self, scalar: f64) {
        self.x *= scalar;
        self.y *= scalar;
        self.z *= scalar;
    }
}

impl ops::Div<f64> for Vec3D {
    type Output = Self;

    #[inline]
    fn div(self, scalar: f64) -> Self {
        Self::new(self.x / scalar, self.y / scalar, self.z / scalar)
    }
}

impl ops::DivAssign<f64> for Vec3D {
    #[inline]
    fn div_assign(&mut self, scalar: f64) {
        self.x /= scalar;
        self.y /= scalar;
        self.z /= scalar;
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
        vec.normalise();
        assert!((vec.x - 0.2672612419124243).abs() < 1e-15);
        assert!((vec.y - 0.5345224838248487).abs() < 1e-15);
        assert!((vec.z - 0.8017837257372731).abs() < 1e-15);
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

    #[test]
    fn add() {
        let vec1 = Vec3D::new(1.0, 2.0, 3.0);
        let vec2 = Vec3D::new(1.0, -1.0, 1.0);
        assert_eq!(vec1 + vec2, Vec3D::new(2.0, 1.0, 4.0));
    }

    #[test]
    fn add_assign() {
        let mut vec1 = Vec3D::new(1.0, 2.0, 3.0);
        let vec2 = Vec3D::new(1.0, -1.0, 1.0);
        vec1 += vec2;
        assert_eq!(vec1, Vec3D::new(2.0, 1.0, 4.0));
    }

    #[test]
    fn sub() {
        let vec1 = Vec3D::new(1.0, 2.0, 3.0);
        let vec2 = Vec3D::new(1.0, -1.0, 1.0);
        assert_eq!(vec1 - vec2, Vec3D::new(0.0, 3.0, 2.0));
    }

    #[test]
    fn sub_assign() {
        let mut vec1 = Vec3D::new(1.0, 2.0, 3.0);
        let vec2 = Vec3D::new(1.0, -1.0, 1.0);
        vec1 -= vec2;
        assert_eq!(vec1, Vec3D::new(0.0, 3.0, 2.0));
    }

    #[test]
    fn neg() {
        let vec = Vec3D::new(1.0, 2.0, 3.0);
        assert_eq!(-vec, Vec3D::new(-1.0, -2.0, -3.0));
    }

    #[test]
    fn mul() {
        let vec = Vec3D::new(1.0, 2.0, 3.0);
        assert_eq!(vec * 2.0, Vec3D::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn mul_assign() {
        let mut vec = Vec3D::new(1.0, 2.0, 3.0);
        vec *= 2.0;
        assert_eq!(vec, Vec3D::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn div() {
        let vec = Vec3D::new(2.0, 4.0, 6.0);
        assert_eq!(vec / 2.0, Vec3D::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn div_assign() {
        let mut vec = Vec3D::new(2.0, 4.0, 6.0);
        vec /= 2.0;
        assert_eq!(vec, Vec3D::new(1.0, 2.0, 3.0));
    }
}