use super::utils::safe_sqrt;
use core::ops::*;

#[derive(Debug, Clone, Copy)]
pub struct Complex {
    pub real: f32,
    pub imag: f32,
}

impl Complex {
    pub fn conj(self) -> Self {
        Self {
            real: self.real,
            imag: -self.imag,
        }
    }

    pub fn recip(self) -> Self {
        self.conj() / self.abs_square()
    }

    pub fn abs_square(self) -> f32 {
        self.real * self.real + self.imag * self.imag
    }

    pub fn norm(self) -> f32 {
        f32::sqrt(self.abs_square())
    }

    pub fn sqrt(z: Self) -> Self {
        let r = z.norm();
        let a = safe_sqrt((r + z.real) / 2.0);
        let b = safe_sqrt((r - z.real) / 2.0) * z.imag.signum();
        return Complex { real: a, imag: b };
    }

    pub fn sqrt_from_pbrt(z: Self) -> Self {
        let r = z.norm();

        if r == 0.0 {
            return Complex {
                real: 0.0,
                imag: 0.0,
            };
        }

        let t1 = f32::sqrt((r + z.real.abs()) / 2.0);
        let t2 = z.imag / t1 / 2.0;

        if z.real >= 0.0 {
            return Complex { real: t1, imag: t2 };
        } else {
            return Complex {
                real: t2.abs(),
                imag: t1.copysign(z.imag),
            };
        }
    }
}

impl Mul<Complex> for Complex {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Complex) -> Self {
        Self {
            real: self.real.mul(rhs.real) - self.imag.mul(rhs.imag),
            imag: self.real.mul(rhs.imag) + self.imag.mul(rhs.real),
        }
    }
}

impl Mul<f32> for Complex {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f32) -> Self {
        Self {
            real: self.real.mul(rhs),
            imag: self.imag.mul(rhs),
        }
    }
}

impl Mul<Complex> for f32 {
    type Output = Complex;
    #[inline]
    fn mul(self, rhs: Complex) -> Complex {
        Complex {
            real: self.mul(rhs.real),
            imag: self.mul(rhs.imag),
        }
    }
}

impl Div<Complex> for Complex {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Complex) -> Self {
        self.mul(rhs.recip())
    }
}

impl Div<f32> for Complex {
    type Output = Self;
    #[inline]
    fn div(self, rhs: f32) -> Self {
        self * rhs.recip()
    }
}

impl Div<Complex> for f32 {
    type Output = Complex;
    #[inline]
    fn div(self, rhs: Complex) -> Complex {
        self * rhs.recip()
    }
}

impl Add<Complex> for Complex {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            real: self.real.add(rhs.real),
            imag: self.imag.add(rhs.imag),
        }
    }
}

impl Add<f32> for Complex {
    type Output = Self;
    #[inline]
    fn add(self, rhs: f32) -> Self {
        Self {
            real: self.real.add(rhs),
            imag: self.imag,
        }
    }
}

impl Add<Complex> for f32 {
    type Output = Complex;
    #[inline]
    fn add(self, rhs: Complex) -> Complex {
        Complex {
            real: rhs.real.add(self),
            imag: rhs.imag,
        }
    }
}

impl Sub<Complex> for Complex {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            real: self.real.sub(rhs.real),
            imag: self.imag.sub(rhs.imag),
        }
    }
}

impl Sub<f32> for Complex {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: f32) -> Self {
        Self {
            real: self.real.sub(rhs),
            imag: self.imag,
        }
    }
}

impl Sub<Complex> for f32 {
    type Output = Complex;
    #[inline]
    fn sub(self, rhs: Complex) -> Complex {
        Complex {
            real: rhs.real.sub(self),
            imag: rhs.imag,
        }
    }
}

impl Neg for Complex {
    type Output = Complex;
    #[inline]
    fn neg(self) -> Self::Output {
        Self {
            real: self.real.neg(),
            imag: self.imag.neg(),
        }
    }
}
