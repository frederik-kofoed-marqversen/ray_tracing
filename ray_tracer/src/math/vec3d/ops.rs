use super::Vec3D;
use core::iter::Sum;
use core::ops::*;

impl Div<f32> for Vec3D {
    type Output = Self;
    #[inline]
    fn div(self, rhs: f32) -> Self {
        self * rhs.recip()
    }
}

impl Div<&f32> for Vec3D {
    type Output = Vec3D;
    #[inline]
    fn div(self, rhs: &f32) -> Vec3D {
        self.div(*rhs)
    }
}

impl Div<&f32> for &Vec3D {
    type Output = Vec3D;
    #[inline]
    fn div(self, rhs: &f32) -> Vec3D {
        (*self).div(*rhs)
    }
}

impl Div<f32> for &Vec3D {
    type Output = Vec3D;
    #[inline]
    fn div(self, rhs: f32) -> Vec3D {
        (*self).div(rhs)
    }
}

impl DivAssign<f32> for Vec3D {
    #[inline]
    fn div_assign(&mut self, rhs: f32) {
        self.x.div_assign(rhs);
        self.y.div_assign(rhs);
        self.z.div_assign(rhs);
    }
}

impl DivAssign<&f32> for Vec3D {
    #[inline]
    fn div_assign(&mut self, rhs: &f32) {
        self.div_assign(*rhs)
    }
}

impl Div<Vec3D> for f32 {
    type Output = Vec3D;
    #[inline]
    fn div(self, rhs: Vec3D) -> Vec3D {
        Vec3D {
            x: self.div(rhs.x),
            y: self.div(rhs.y),
            z: self.div(rhs.z),
        }
    }
}

impl Div<&Vec3D> for f32 {
    type Output = Vec3D;
    #[inline]
    fn div(self, rhs: &Vec3D) -> Vec3D {
        self.div(*rhs)
    }
}

impl Div<&Vec3D> for &f32 {
    type Output = Vec3D;
    #[inline]
    fn div(self, rhs: &Vec3D) -> Vec3D {
        (*self).div(*rhs)
    }
}

impl Div<Vec3D> for &f32 {
    type Output = Vec3D;
    #[inline]
    fn div(self, rhs: Vec3D) -> Vec3D {
        (*self).div(rhs)
    }
}

impl Mul<f32> for Vec3D {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f32) -> Self {
        Self {
            x: self.x.mul(rhs),
            y: self.y.mul(rhs),
            z: self.z.mul(rhs),
        }
    }
}

impl Mul<&f32> for Vec3D {
    type Output = Vec3D;
    #[inline]
    fn mul(self, rhs: &f32) -> Vec3D {
        self.mul(*rhs)
    }
}

impl Mul<&f32> for &Vec3D {
    type Output = Vec3D;
    #[inline]
    fn mul(self, rhs: &f32) -> Vec3D {
        (*self).mul(*rhs)
    }
}

impl Mul<f32> for &Vec3D {
    type Output = Vec3D;
    #[inline]
    fn mul(self, rhs: f32) -> Vec3D {
        (*self).mul(rhs)
    }
}

impl MulAssign<f32> for Vec3D {
    #[inline]
    fn mul_assign(&mut self, rhs: f32) {
        self.x.mul_assign(rhs);
        self.y.mul_assign(rhs);
        self.z.mul_assign(rhs);
    }
}

impl MulAssign<&f32> for Vec3D {
    #[inline]
    fn mul_assign(&mut self, rhs: &f32) {
        self.mul_assign(*rhs)
    }
}

impl Mul<Vec3D> for f32 {
    type Output = Vec3D;
    #[inline]
    fn mul(self, rhs: Vec3D) -> Vec3D {
        Vec3D {
            x: self.mul(rhs.x),
            y: self.mul(rhs.y),
            z: self.mul(rhs.z),
        }
    }
}

impl Mul<&Vec3D> for f32 {
    type Output = Vec3D;
    #[inline]
    fn mul(self, rhs: &Vec3D) -> Vec3D {
        self.mul(*rhs)
    }
}

impl Mul<&Vec3D> for &f32 {
    type Output = Vec3D;
    #[inline]
    fn mul(self, rhs: &Vec3D) -> Vec3D {
        (*self).mul(*rhs)
    }
}

impl Mul<Vec3D> for &f32 {
    type Output = Vec3D;
    #[inline]
    fn mul(self, rhs: Vec3D) -> Vec3D {
        (*self).mul(rhs)
    }
}

impl Add<Vec3D> for Vec3D {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x.add(rhs.x),
            y: self.y.add(rhs.y),
            z: self.z.add(rhs.z),
        }
    }
}

impl Add<&Vec3D> for Vec3D {
    type Output = Vec3D;
    #[inline]
    fn add(self, rhs: &Vec3D) -> Vec3D {
        self.add(*rhs)
    }
}

impl Add<&Vec3D> for &Vec3D {
    type Output = Vec3D;
    #[inline]
    fn add(self, rhs: &Vec3D) -> Vec3D {
        (*self).add(*rhs)
    }
}

impl Add<Vec3D> for &Vec3D {
    type Output = Vec3D;
    #[inline]
    fn add(self, rhs: Vec3D) -> Vec3D {
        (*self).add(rhs)
    }
}

impl AddAssign<Vec3D> for Vec3D {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x.add_assign(rhs.x);
        self.y.add_assign(rhs.y);
        self.z.add_assign(rhs.z);
    }
}

impl AddAssign<&Self> for Vec3D {
    #[inline]
    fn add_assign(&mut self, rhs: &Self) {
        self.add_assign(*rhs)
    }
}

impl Sub<Vec3D> for Vec3D {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x.sub(rhs.x),
            y: self.y.sub(rhs.y),
            z: self.z.sub(rhs.z),
        }
    }
}

impl Sub<&Vec3D> for Vec3D {
    type Output = Vec3D;
    #[inline]
    fn sub(self, rhs: &Vec3D) -> Vec3D {
        self.sub(*rhs)
    }
}

impl Sub<&Vec3D> for &Vec3D {
    type Output = Vec3D;
    #[inline]
    fn sub(self, rhs: &Vec3D) -> Vec3D {
        (*self).sub(*rhs)
    }
}

impl Sub<Vec3D> for &Vec3D {
    type Output = Vec3D;
    #[inline]
    fn sub(self, rhs: Vec3D) -> Vec3D {
        (*self).sub(rhs)
    }
}

impl SubAssign<Vec3D> for Vec3D {
    #[inline]
    fn sub_assign(&mut self, rhs: Vec3D) {
        self.x.sub_assign(rhs.x);
        self.y.sub_assign(rhs.y);
        self.z.sub_assign(rhs.z);
    }
}

impl SubAssign<&Self> for Vec3D {
    #[inline]
    fn sub_assign(&mut self, rhs: &Self) {
        self.sub_assign(*rhs)
    }
}

impl Neg for Vec3D {
    type Output = Vec3D;
    #[inline]
    fn neg(self) -> Self::Output {
        Self {
            x: self.x.neg(),
            y: self.y.neg(),
            z: self.z.neg(),
        }
    }
}

impl Neg for &Vec3D {
    type Output = Vec3D;
    #[inline]
    fn neg(self) -> Self::Output {
        Vec3D {
            x: self.x.neg(),
            y: self.y.neg(),
            z: self.z.neg(),
        }
    }
}

impl Sum<Self> for Vec3D {
    #[inline]
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Self::ZERO, |a, b| a + b)
    }
}

impl<'a> Sum<&'a Self> for Vec3D {
    #[inline]
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
    {
        iter.fold(Self::ZERO, |a, b| a + b)
    }
}

impl From<[f32; 3]> for Vec3D {
    #[inline]
    fn from(a: [f32; 3]) -> Self {
        Self::new(a[0], a[1], a[2])
    }
}

impl From<Vec3D> for [f32; 3] {
    #[inline]
    fn from(v: Vec3D) -> Self {
        [v.x, v.y, v.z]
    }
}

impl From<(f32, f32, f32)> for Vec3D {
    #[inline]
    fn from(t: (f32, f32, f32)) -> Self {
        Self::new(t.0, t.1, t.2)
    }
}

impl std::fmt::Display for Vec3D {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[{:?}, {:?}, {:?}]", self.x, self.y, self.z)
    }
}

impl Default for Vec3D {
    fn default() -> Self {
        Self::ZERO
    }
}

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
