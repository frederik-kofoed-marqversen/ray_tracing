use std::f32::consts::PI;

use fastrand::Rng;

use super::primitives::Sphere;
use super::traits::Surface;
use super::vec3d::utils::{sample_unit_sphere, tangent_space};
use super::{Ray, Vec3D};

pub trait Light {
    fn sample(&self, ref_point: Vec3D, rng: &mut Rng) -> Option<LightSample>;
    fn pdf(&self, ref_point: Vec3D, dir: Vec3D) -> f32;
}

#[derive(Debug)]
pub struct LightSample {
    pub spectrum: Vec3D,
    pub point: Vec3D,
    pub direction: Vec3D,
    pub pdf: f32,
}

#[derive(Debug)]
pub struct PointLight {
    pos: Vec3D,
    intensity: Vec3D,
}

impl Light for PointLight {
    fn sample(&self, ref_point: Vec3D, _: &mut Rng) -> Option<LightSample> {
        let direction = self.pos - ref_point;
        let d2 = direction.norm_squared();

        Some(LightSample {
            spectrum: self.intensity / d2,
            point: self.pos,
            direction: direction / d2.sqrt(),
            pdf: 1.0,
        })
    }

    fn pdf(&self, _: Vec3D, _: Vec3D) -> f32 {
        0.0
    }
}

#[derive(Debug)]
pub struct Spotlight {
    pos: Vec3D,
    dir: Vec3D,
    cos_angle: f32,
    intensity: Vec3D,
}

impl Light for Spotlight {
    fn sample(&self, ref_point: Vec3D, _: &mut Rng) -> Option<LightSample> {
        let direction = (self.pos - ref_point).normalise();
        let d2 = direction.norm_squared();
        let direction = direction / d2.sqrt();

        if -Vec3D::dot(direction, self.dir) > self.cos_angle {
            return Some(LightSample {
                spectrum: self.intensity / d2,
                point: self.pos,
                direction,
                pdf: 1.0,
            });
        } else {
            return None;
        }
    }

    fn pdf(&self, _: Vec3D, _: Vec3D) -> f32 {
        0.0
    }
}

pub struct AreaLight {
    surface: Sphere,  // Currently only SPHERE area lights are supported
    emmited_radiance: Vec3D,
}

impl Light for AreaLight {
    fn sample(&self, ref_point: Vec3D, rng: &mut Rng) -> Option<LightSample> {
        let z = ref_point - self.surface.center;
        let d2 = z.norm_squared();
        let sin_theta_max_2 = self.surface.radius * self.surface.radius / d2;

        if sin_theta_max_2 > 1.0 {
            // d2 < r2
            // Reference point is inside the sphere => Use uniform area sampling
            let normal = sample_unit_sphere(rng);
            let point = self.surface.radius * normal + self.surface.center;
            let direction = point - ref_point;
            let norm_2 = direction.norm_squared();
            let direction = direction / norm_2.sqrt();
            let pdf = norm_2 / (self.surface.area() * Vec3D::dot(direction, normal).abs());
            return Some(LightSample {
                spectrum: self.emmited_radiance,
                point,
                direction,
                pdf,
            });
        }

        // Reference point is outside the sphere => Use uniform solid angle sampling
        let sin_theta_max = sin_theta_max_2.sqrt();
        let cos_theta_max = f32::sqrt(1.0 - sin_theta_max_2);

        let cos_theta = (cos_theta_max - 1.0) * rng.f32() + 1.0;
        let sin_theta_2 = 1.0 - cos_theta * cos_theta;

        let cos_alpha = sin_theta_2 / sin_theta_max
            + cos_theta * f32::sqrt(1.0 - sin_theta_2 / sin_theta_max_2);
        let sin_alpha = f32::sqrt(1.0 - cos_alpha * cos_alpha);
        let phi = 2.0 * PI * rng.f32();

        let z = z / d2.sqrt();
        let (x, y) = tangent_space(z);

        let point = self.surface.radius
            * (cos_alpha * phi.cos() * x + cos_alpha * phi.sin() * y + sin_alpha * z);
        let direction = (point - ref_point).normalise();
        let pdf = 1.0 / (2.0 * PI * (1.0 - cos_theta_max));
        return Some(LightSample {
            spectrum: self.emmited_radiance,
            point,
            direction,
            pdf,
        });
    }

    fn pdf(&self, ref_point: Vec3D, dir: Vec3D) -> f32 {
        let ray = Ray {
            origin: ref_point,
            direction: dir,
        };

        if let Some((t, normal)) = self.surface.hit(&ray, 0.0, f32::INFINITY) {
            let d2 = (ref_point - self.surface.center).norm_squared();
            let sin_theta_max_2 = self.surface.radius * self.surface.radius / d2;

            if sin_theta_max_2 > 1.0 {
                // Uniform area sampling
                return t * t / (self.surface.area() * Vec3D::dot(dir, normal).abs());
            } else {
                // Uniform solid angle sampling
                let cos_theta_max = f32::sqrt(1.0 - sin_theta_max_2);
                return 1.0 / (2.0 * PI * (1.0 - cos_theta_max));
            }
        } else {
            // Light not sampled by given ray
            return 0.0;
        }
    }
}
