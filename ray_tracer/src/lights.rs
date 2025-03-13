use std::f32::consts::PI;

use fastrand::Rng;

use super::primitives::Sphere;
use super::traits::Light;
use super::vec3d::utils::{sample_unit_sphere, tangent_space};
use super::Vec3D;

#[derive(Debug)]
pub struct PointLight {
    pos: Vec3D,
}

impl Light for PointLight {
    fn sample_light(&self, ref_point: Vec3D, _: &mut Rng) -> Option<(Vec3D, Vec3D, f32)> {
        Some((self.pos, (self.pos - ref_point).normalise(), 1.0))
    }
}

#[derive(Debug)]
pub struct Spotlight {
    pos: Vec3D,
    dir: Vec3D,
    cos_angle: f32,
}

impl Light for Spotlight {
    fn sample_light(&self, ref_point: Vec3D, _: &mut Rng) -> Option<(Vec3D, Vec3D, f32)> {
        let direction = (self.pos - ref_point).normalise();
        if -Vec3D::dot(direction, self.dir) > self.cos_angle {
            return Some((self.pos, direction, 1.0));
        } else {
            return None;
        }
    }
}

impl Light for Sphere {
    fn sample_light(&self, ref_point: Vec3D, rng: &mut Rng) -> Option<(Vec3D, Vec3D, f32)> {
        let z = ref_point - self.center;
        let d2 = z.norm_squared();

        if d2 < self.radius * self.radius {
            // Reference point is inside the sphere => Use uniform area sampling
            let normal = sample_unit_sphere(rng);
            let point = self.radius * normal + self.center;
            let direction = point - ref_point;
            let norm_2 = direction.norm_squared();
            let direction = direction / norm_2.sqrt();
            let pdf = norm_2 / (self.area() * Vec3D::dot(direction, normal).abs());
            return Some((point, direction, pdf));
        }

        let sin_theta_max = self.radius * self.radius / d2;
        let sin_theta_max_2 = sin_theta_max * sin_theta_max;
        let cos_theta_max = f32::sqrt(1.0 - sin_theta_max_2);

        let cos_theta = (cos_theta_max - 1.0) * rng.f32() + 1.0;
        let sin_theta_2 = 1.0 - cos_theta * cos_theta;

        let cos_alpha = sin_theta_2 / sin_theta_max
            + cos_theta * f32::sqrt(1.0 - sin_theta_2 / sin_theta_max_2);
        let sin_alpha = f32::sqrt(1.0 - cos_alpha * cos_alpha);
        let phi = 2.0 * PI * rng.f32();

        let z = z / f32::sqrt(d2);
        let (x, y) = tangent_space(z);

        let point =
            self.radius * (cos_alpha * phi.cos() * x + cos_alpha * phi.sin() * y + sin_alpha * z);
        let direction = (point - ref_point).normalise();
        let pdf = 1.0 / (2.0 * PI * (1.0 - cos_theta_max));
        return Some((point, direction, pdf));
    }
}
