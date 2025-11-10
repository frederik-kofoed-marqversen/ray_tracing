use std::f32::consts::PI;

use fastrand::Rng;

use super::primitives::Sphere;
use super::traits::Surface;
use super::vec3d::utils::{sample_unit_sphere, tangent_space};
use super::{Ray, Vec3D};

#[inline]
fn safe_sqrt(v: f32) -> f32 {
    if v <= 0.0 {
        0.0
    } else {
        v.sqrt()
    }
}

pub trait Light {
    fn sample(&self, ref_point: Vec3D, rng: &mut Rng) -> LightSample;
    fn pdf(&self, _ray: &Ray) -> f32 {
        0.0
    }
    fn hit(&self, _ray: &Ray, _t_min: f32, _t_max: f32) -> Option<(f32, Vec3D, Vec3D)> {
        None
    }
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
    pub pos: Vec3D,
    pub intensity: Vec3D,
}

impl Light for PointLight {
    fn sample(&self, ref_point: Vec3D, _: &mut Rng) -> LightSample {
        let direction = self.pos - ref_point;
        let d2 = direction.norm_squared();

        LightSample {
            spectrum: self.intensity / d2,
            point: self.pos,
            direction: direction / d2.sqrt(),
            pdf: 1.0,
        }
    }
}

#[derive(Debug)]
pub struct Spotlight {
    pub pos: Vec3D,
    pub dir: Vec3D,
    pub cos_angle: f32,
    pub intensity: Vec3D,
}

impl Light for Spotlight {
    fn sample(&self, ref_point: Vec3D, _: &mut Rng) -> LightSample {
        let direction = (self.pos - ref_point).normalise();
        let d2 = direction.norm_squared();
        let direction = direction / d2.sqrt();
        let mut sample = LightSample {
            spectrum: Vec3D::ZERO,
            point: self.pos,
            direction,
            pdf: 1.0,
        };

        if -Vec3D::dot(direction, self.dir) > self.cos_angle {
            sample.spectrum = self.intensity / d2;
        }

        return sample;
    }
}

pub struct AreaLight {
    pub surface: Sphere, // Currently only SPHERE area lights are supported
    pub emmited_radiance: Vec3D,
}

impl Light for AreaLight {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D, Vec3D)> {
        self.surface
            .hit(ray, t_min, t_max)
            .map(|(t, normal)| (t, normal, self.emmited_radiance))
    }

    fn sample(&self, ref_point: Vec3D, rng: &mut Rng) -> LightSample {
        let z = self.surface.center - ref_point;
        let d2 = z.norm_squared();
        let sin_theta_max_2 = self.surface.radius * self.surface.radius / d2;

        if sin_theta_max_2 > 1.0 {
            // d2 < r2
            // Reference point is inside the sphere => Use uniform area sampling
            let normal = sample_unit_sphere(rng);
            let point = self.surface.center + self.surface.radius * normal;
            let direction = point - ref_point;
            let norm_2 = direction.norm_squared();
            let direction = direction / norm_2.sqrt();
            let pdf = norm_2 / (self.surface.area() * Vec3D::dot(direction, normal).abs());

            return LightSample {
                spectrum: self.emmited_radiance,
                point,
                direction,
                pdf,
            };
        }

        // Reference point is outside the sphere => Use uniform solid angle sampling
        let sin_theta_max = sin_theta_max_2.sqrt();
        let cos_theta_max = safe_sqrt(1.0 - sin_theta_max_2);

        let cos_theta = (cos_theta_max - 1.0) * rng.f32() + 1.0;
        let sin_theta_2 = 1.0 - cos_theta * cos_theta;

        let cos_alpha = sin_theta_2 / sin_theta_max
            + cos_theta * safe_sqrt(1.0 - sin_theta_2 / sin_theta_max_2);
        let sin_alpha = safe_sqrt(1.0 - cos_alpha * cos_alpha);
        let phi = 2.0 * PI * rng.f32();

        let z = z / d2.sqrt();
        let (x, y) = tangent_space(z);

        let point = self.surface.center
            + self.surface.radius
                * (cos_alpha * phi.cos() * x + cos_alpha * phi.sin() * y + sin_alpha * z);
        let direction = (point - ref_point).normalise();
        let pdf = 1.0 / (2.0 * PI * (1.0 - cos_theta_max));
        return LightSample {
            spectrum: self.emmited_radiance,
            point,
            direction,
            pdf,
        };
    }

    fn pdf(&self, ray: &Ray) -> f32 {
        if let Some((t, normal)) = self.surface.hit(&ray, 0.0, f32::INFINITY) {
            let d2 = (ray.origin - self.surface.center).norm_squared();
            let sin_theta_max_2 = self.surface.radius * self.surface.radius / d2;

            if sin_theta_max_2 > 1.0 {
                // Uniform area sampling
                return t * t / (self.surface.area() * Vec3D::dot(ray.direction, normal).abs());
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
