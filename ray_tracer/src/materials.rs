use fastrand::Rng;

use super::vec3d::utils::sample_unit_sphere;
use super::Vec3D;

#[derive(Default, Clone, Copy)]
pub struct Material {
    pub albedo: Vec3D,
    pub roughness: f32,
    pub emission: Option<(Vec3D, f32)>,
    pub refractive_index: Option<f32>,
}

impl Material {
    pub fn sample_bsdf(&self, dir_out: Vec3D, mut normal: Vec3D, rng: &mut Rng) -> Vec3D {
        // Compute angle of incidence and normal vector pointing up from surface in direction of where
        // the ray came from. Also compute eta = n1/n2 where n1 is the refractive index of material
        // where the ray comes from and n2 that of the material that ray goes into.
        let mut cos_theta = Vec3D::dot(dir_out, normal); // cosine of angle of incidence
        let mut eta = self.refractive_index;
        if cos_theta < 0.0 {
            cos_theta = -cos_theta;
            normal = -normal;
        } else {
            eta = eta.and_then(|n| Some(1.0 / n));
        }

        // Compute reflection coefficient and compute emitted ray by reflection or refraction
        let (eta, reflection_coefficient) = if let Some(eta) = eta {
            (eta, reflection_coefficient(eta, cos_theta))
        } else {
            (1.0, 1.0)
        };
        let mut dir_sampled;
        if rng.f32() < reflection_coefficient {
            // Reflect
            dir_sampled = reflect(dir_out, normal, cos_theta);
        } else {
            // Refract
            dir_sampled = refract(dir_out, normal, eta, cos_theta);
            // Refracted rays scatter along the negative of the normal vector
            normal = -normal;
        };

        // Sample Lambertian scattering direction
        let dir_scattered = scatter(normal, rng);
        // Add scattering direction to ray direction according to roughness
        // roughness=0 => reflection, roughness=1 => diffuse
        dir_sampled = Vec3D::interpolate(dir_sampled, dir_scattered, self.roughness).normalise();

        return dir_sampled;
    }
}

#[inline]
fn reflection_coefficient(eta: f32, cos_theta: f32) -> f32 {
    let sin_theta = f32::sqrt(1.0 - cos_theta * cos_theta);
    if eta * sin_theta > 1.0 {
        return 1.0;
    }
    // Schlick's approximation
    let r0 = f32::powi((eta - 1.0) / (eta + 1.0), 2);
    return r0 + (1.0 - r0) * f32::powi(1.0 - cos_theta, 5);
}

#[inline]
fn refract(dir_in: Vec3D, normal: Vec3D, eta: f32, cos_theta_in: f32) -> Vec3D {
    let sin_theta_2 = (1.0 - cos_theta_in * cos_theta_in).max(0.0);
    let sin_theta_out_2 = eta * eta * sin_theta_2;
    let cos_theta_out = f32::sqrt(1.0 - sin_theta_out_2);
    return -eta * dir_in + (eta * cos_theta_in - cos_theta_out) * normal;
}

#[inline]
fn reflect(dir_out: Vec3D, normal: Vec3D, cos_theta: f32) -> Vec3D {
    -dir_out + 2.0 * cos_theta * normal
}

#[inline]
fn scatter(normal: Vec3D, rng: &mut Rng) -> Vec3D {
    (normal + sample_unit_sphere(rng))
        .try_normalise()
        .unwrap_or(normal)
}

// Constructors
impl Material {
    #[inline]
    pub fn diffuse(albedo: Vec3D) -> Self {
        Self {
            albedo,
            roughness: 1.0,
            emission: None,
            refractive_index: None,
        }
    }

    #[inline]
    pub fn light_source(emission: Vec3D, intensity: f32) -> Self {
        Self {
            albedo: Vec3D::ZERO,
            roughness: 1.0,
            emission: Some((emission, intensity)),
            refractive_index: None,
        }
    }

    #[inline]
    pub fn mirror() -> Self {
        Self {
            albedo: Vec3D::ONES,
            roughness: 0.0,
            emission: None,
            refractive_index: None,
        }
    }

    #[inline]
    pub fn metal(albedo: Vec3D, roughness: f32) -> Self {
        Self {
            albedo,
            roughness,
            emission: None,
            refractive_index: None,
        }
    }

    #[inline]
    pub fn dielectric(albedo: Vec3D, refractive_index: f32) -> Self {
        Self {
            albedo,
            roughness: 0.0,
            emission: None,
            refractive_index: Some(refractive_index),
        }
    }
}
