use std::ops;
use super::Vec3D;

impl_op_ex!(+ |lhs: &Vec3D, rhs: &Vec3D| -> Vec3D {
    Vec3D::new(
        lhs.x + rhs.x, 
        lhs.y + rhs.y, 
        lhs.z + rhs.z
    )
});

impl_op_ex!(+= |lhs: &mut Vec3D, rhs: &Vec3D| {
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;
});

impl_op_ex!(- |lhs: &Vec3D, rhs: &Vec3D| -> Vec3D {
    Vec3D::new(
        lhs.x - rhs.x, 
        lhs.y - rhs.y, 
        lhs.z - rhs.z
    )
});

impl_op_ex!(-= |lhs: &mut Vec3D, rhs: &Vec3D| {
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    lhs.z -= rhs.z;
});

impl_op_ex!(- |vec: &Vec3D| -> Vec3D {
    Vec3D::new(-vec.x, -vec.y, -vec.z)
});

impl_op_ex!(* |lhs: &Vec3D, rhs: &Vec3D| -> Vec3D {
    Vec3D::new(
        lhs.x * rhs.x, 
        lhs.y * rhs.y, 
        lhs.z * rhs.z
    )
});

impl_op_ex!(*= |lhs: &mut Vec3D, rhs: &Vec3D| {
    lhs.x *= rhs.x;
    lhs.y *= rhs.y;
    lhs.z *= rhs.z;
});

impl_op_ex_commutative!(* |vec: &Vec3D, scalar: f64| -> Vec3D {
    Vec3D::new(vec.x * scalar, vec.y * scalar, vec.z * scalar)
});

impl_op_ex!(*= |lhs: &mut Vec3D, scalar: f64| {
    lhs.x *= scalar;
    lhs.y *= scalar;
    lhs.z *= scalar;
});

impl_op_ex!(/ |vec: &Vec3D, scalar: f64| -> Vec3D {
    Vec3D::new(vec.x / scalar, vec.y / scalar, vec.z / scalar)
});

impl_op_ex!(/= |lhs: &mut Vec3D, scalar: f64| {
    lhs.x /= scalar;
    lhs.y /= scalar;
    lhs.z /= scalar;
});

/* POSSIBLY WRITE OWN THAT ARE OPTIMISED WRT. MEMORY ALLOCATION WHEN TAKING OWNERSHIP
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

impl ops::Mul<Vec3D> for f64 {
    type Output = Vec3D;

    #[inline]
    fn mul(self, vec: Vec3D) -> Vec3D {
        vec * self
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
} */

#[cfg(test)]
mod tests {
    use super::*;

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
    fn elementwise_mul() {
        let vec1 = Vec3D::new(0.0, 1.0, 2.0);
        let vec2 = Vec3D::new(1.0, -1.0, 1.0);
        assert_eq!(vec1 * vec2, Vec3D::new(0.0, -1.0, 2.0));
    }

    #[test]
    fn elementwise_mul_assign() {
        let mut vec1 = Vec3D::new(0.0, 1.0, 2.0);
        let vec2 = Vec3D::new(1.0, -1.0, 1.0);
        vec1 *= vec2;
        assert_eq!(vec1, Vec3D::new(0.0, -1.0, 2.0));
    }

    #[test]
    fn scalar_mul() {
        let vec = Vec3D::new(1.0, 2.0, 3.0);
        assert_eq!(vec * 2.0, Vec3D::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn scalar_left_mul() {
        let vec = Vec3D::new(1.0, 2.0, 3.0);
        assert_eq!(2.0 * vec, Vec3D::new(2.0, 4.0, 6.0));
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