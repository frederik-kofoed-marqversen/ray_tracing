#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vec3D {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vec3D {
    pub const ZERO: Self = Self {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    pub const ONES: Self = Self {
        x: 1.0,
        y: 1.0,
        z: 1.0,
    };
    pub const X: Self = Self {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    pub const Y: Self = Self {
        x: 0.0,
        y: 1.0,
        z: 0.0,
    };
    pub const Z: Self = Self {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    };

    #[inline]
    #[must_use]
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }

    #[inline]
    #[must_use]
    pub fn fill(val: f32) -> Self {
        Self {
            x: val,
            y: val,
            z: val,
        }
    }

    #[inline]
    #[must_use]
    pub fn new_polar(r: f32, theta: f32, phi: f32) -> Self {
        Self {
            x: r * f32::sin(theta) * f32::cos(phi),
            y: r * f32::sin(theta) * f32::sin(phi),
            z: r * f32::cos(theta),
        }
    }

    #[inline]
    #[must_use]
    pub fn norm_squared(self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    #[inline]
    #[must_use]
    pub fn norm(self) -> f32 {
        self.norm_squared().sqrt()
    }

    #[inline]
    #[must_use]
    pub fn norm_recip(self) -> f32 {
        self.norm_squared().sqrt().recip()
    }

    #[inline]
    #[must_use]
    pub fn normalise(self) -> Self {
        let norm_recip = self.norm_recip();
        assert!(norm_recip.is_finite());
        return self * norm_recip;
    }

    #[inline]
    #[must_use]
    pub fn try_normalise(self) -> Option<Self> {
        let norm_recip = self.norm_recip();
        if norm_recip.is_finite() {
            return Some(self * norm_recip);
        } else {
            return None;
        }
    }

    #[inline]
    #[must_use]
    pub fn almost_zero(self, abs: f32) -> bool {
        return self.norm_squared() < abs;
    }

    #[inline]
    #[must_use]
    pub fn almost_equal(vec1: Self, vec2: Self, abs: f32) -> bool {
        (vec1 - vec2).almost_zero(abs)
    }

    #[inline]
    #[must_use]
    pub fn is_zero(self) -> bool {
        return self.x == 0.0 && self.y == 0.0 && self.z == 0.0;
    }

    #[inline]
    #[must_use]
    pub fn dot(vec1: Self, vec2: Self) -> f32 {
        vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z
    }

    #[inline]
    #[must_use]
    pub fn cross(vec1: Self, vec2: Self) -> Self {
        Self::new(
            vec1.y * vec2.z - vec1.z * vec2.y,
            vec1.z * vec2.x - vec1.x * vec2.z,
            vec1.x * vec2.y - vec1.y * vec2.x,
        )
    }

    #[inline]
    #[must_use]
    pub fn min(vec1: Self, vec2: Self) -> Self {
        Self {
            x: vec1.x.min(vec2.x),
            y: vec1.y.min(vec2.y),
            z: vec1.z.min(vec2.z),
        }
    }

    #[inline]
    #[must_use]
    pub fn max(vec1: Self, vec2: Self) -> Self {
        Self {
            x: vec1.x.max(vec2.x),
            y: vec1.y.max(vec2.y),
            z: vec1.z.max(vec2.z),
        }
    }

    #[inline]
    #[must_use]
    pub fn mul_elemwise(vec1: Self, vec2: Self) -> Self {
        Self {
            x: vec1.x * vec2.x,
            y: vec1.y * vec2.y,
            z: vec1.z * vec2.z,
        }
    }

    #[inline]
    #[must_use]
    pub fn min_elem(self) -> f32 {
        self.x.min(self.y).min(self.z)
    }

    #[inline]
    #[must_use]
    pub fn max_elem(self) -> f32 {
        self.x.max(self.y).max(self.z)
    }

    #[inline]
    #[must_use]
    pub fn map(self, mut fun: impl FnMut(f32) -> f32) -> Self {
        Self {
            x: fun(self.x),
            y: fun(self.y),
            z: fun(self.z),
        }
    }

    #[inline]
    #[must_use]
    pub fn interpolate(vec1: Self, vec2: Self, t: f32) -> Self {
        (1.0 - t) * vec1 + t * vec2
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
        let vec = Vec3D::new(1.0, 2.0, 3.0).normalise();
        assert!((vec.x - 0.2672612419124243).abs() < 1e-6);
        assert!((vec.y - 0.5345224838248487).abs() < 1e-6);
        assert!((vec.z - 0.8017837257372731).abs() < 1e-6);
    }

    #[test]
    fn unit_vector() {
        let vec = Vec3D::new(1.0, 2.0, 3.0);
        let u = vec.normalise();
        assert!((u.x - 0.2672612419124243).abs() < 1e-6);
        assert!((u.y - 0.5345224838248487).abs() < 1e-6);
        assert!((u.z - 0.8017837257372731).abs() < 1e-6);
    }

    #[test]
    fn dot() {
        let vec1 = Vec3D::new(1.0, 2.0, 3.0);
        let vec2 = Vec3D::new(1.0, -1.0, 1.0);
        assert_eq!(Vec3D::dot(vec1, vec2), 2.0);
    }

    #[test]
    fn cross() {
        let vec1 = Vec3D::new(1.0, 2.0, 3.0);
        let vec2 = Vec3D::new(1.0, -1.0, 1.0);
        assert_eq!(Vec3D::cross(vec1, vec2), Vec3D::new(5.0, 2.0, -3.0));
        assert_eq!(Vec3D::cross(vec2, vec1), Vec3D::new(-5.0, -2.0, 3.0));
    }
}
