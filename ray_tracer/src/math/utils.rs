use std::f32::consts::PI;

use super::Vec3D;
use fastrand::Rng;

#[inline]
pub fn safe_sqrt(v: f32) -> f32 {
    if v <= 0.0 {
        0.0
    } else {
        v.sqrt()
    }
}

/// Given unit vector, returns two unit vectors such that all three are pairwise orthogonal.
#[inline]
pub fn tangent_space(vec: Vec3D) -> (Vec3D, Vec3D) {
    let s = f32::signum(vec.z);
    let a = -(s + vec.z).recip();
    let b = vec.x * vec.y * a;

    let t1 = Vec3D::new(1.0 + s * vec.x * vec.x * a, s * b, -s * vec.x);
    let t2 = Vec3D::new(b, s + vec.y * vec.y * a, -vec.y);

    return (t1, t2);
}

/// Uniform area sampling of the unit sphere
#[inline]
pub fn sample_unit_sphere(rng: &mut Rng) -> Vec3D {
    let z = 1.0 - 2.0 * rng.f32();
    let rho = f32::sqrt(1.0 - z * z);
    let phi = 2.0 * PI * rng.f32();
    return Vec3D::new(rho * phi.cos(), rho * phi.sin(), z);
}

/// Samples a distance from an exponential distribution with rate a (pdf: a * exp(-a * x))
#[inline]
pub fn sample_exponential(a: f32, rng: &mut Rng) -> f32 {
    let u = rng.f32();
    return -f32::ln(1.0 - u) / a;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tangent_space() {
        let z = Vec3D::new(1.0, 2.0, 3.0).normalise();
        let (x, y) = tangent_space(z);

        let vecs = [x, y, z];
        for i in 0..3 {
            let j = (i + 1) % 3;
            assert!((vecs[i].norm_squared() - 1.0).abs() < 1e-6);
            assert!(Vec3D::dot(vecs[i], vecs[j]).abs() < 1e-6)
        }
    }
}
