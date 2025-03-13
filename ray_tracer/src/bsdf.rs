use super::complex::Complex;
use super::vec3d::utils::sample_unit_sphere;
use super::Vec3D;
use fastrand::Rng;
use std::f32::consts::PI;

const PI_RECIP: f32 = 1.0 / PI;

pub struct MicroFacetModel {}

pub struct BSDFSample {
    pub spectrum: Vec3D,
    pub pdf: f32,
    pub dir_in: Vec3D,
}

pub trait BSDF {
    // Assume shading coordinates where normal vector is +ẑ
    fn eval_bsdf(&self, dir_in: Vec3D, dir_out: Vec3D) -> Vec3D;
    fn pdf(&self, dir_in: Vec3D, dir_out: Vec3D) -> f32;
    fn sample_bsdf(&self, dir_out: Vec3D, rng: &mut Rng) -> BSDFSample;
}

pub struct Diffuse {
    pub reflectance: Vec3D,
}

impl BSDF for Diffuse {
    fn eval_bsdf(&self, dir_in: Vec3D, dir_out: Vec3D) -> Vec3D {
        if dir_in.z * dir_out.z < 0.0 {
            Vec3D::ZERO
        } else {
            self.reflectance * PI_RECIP
        }
    }

    fn pdf(&self, dir_in: Vec3D, dir_out: Vec3D) -> f32 {
        if dir_in.z * dir_out.z < 0.0 {
            0.0
        } else {
            dir_in.z.abs() * PI_RECIP
        }
    }

    fn sample_bsdf(&self, dir_out: Vec3D, rng: &mut Rng) -> BSDFSample {
        let mut dir_in = sample_unit_sphere(rng).try_normalise().unwrap_or(Vec3D::Z);

        if dir_out.z < 0.0 {
            dir_in.z = -dir_in.z;
        }

        // Small optimisation to not call functions on self since only valid samples are generated
        let spectrum = self.reflectance * PI_RECIP;
        let pdf = dir_in.z.abs() * PI_RECIP;

        return BSDFSample {
            spectrum,
            pdf,
            dir_in,
        };
    }
}

pub struct Dielectric {
    refractive_index: f32,
    roughness: Option<MicroFacetModel>,
}

impl BSDF for Dielectric {
    fn eval_bsdf(&self, _dir_in: Vec3D, _dir_out: Vec3D) -> Vec3D {
        assert!(self.roughness.is_none()); // currently only smooth surfaces
        Vec3D::ZERO
    }

    fn sample_bsdf(&self, dir_out: Vec3D, rng: &mut Rng) -> BSDFSample {
        assert!(self.roughness.is_none()); // currently only smooth surfaces

        let cos_out = dir_out.z.abs();
        let reflectance = fresnell_reflectance(cos_out, self.refractive_index);

        if rng.f32() < reflectance {
            // Specular reflection
            let dir_in = Vec3D::new(-dir_out.x, -dir_out.y, dir_out.z);
            let cos_in = dir_in.z.abs();
            return BSDFSample {
                spectrum: Vec3D::fill(reflectance / cos_in),
                pdf: reflectance,
                dir_in,
            };
        } else {
            // Specular transmission
            let (dir_in, _eta) = refract(dir_out, self.refractive_index);
            let cos_in = dir_in.z.abs();
            let transmission = 1.0 - reflectance;
            let spectrum = Vec3D::fill(transmission / cos_in);
            // if transport_mode is radiance { // when path tracing starting from the camera
            //     spectrum *= eta * eta;
            // }
            return BSDFSample {
                spectrum,
                pdf: transmission,
                dir_in,
            };
        }
    }

    fn pdf(&self, _dir_in: Vec3D, _dir_out: Vec3D) -> f32 {
        assert!(self.roughness.is_none()); // currently only smooth surfaces
        0.0
    }
}

pub struct Conductor {
    refractive_index: Vec3D,       // Frequency dependent
    extinction_coefficient: Vec3D, // Frequency dependent
    roughness: Option<MicroFacetModel>,
}

impl BSDF for Conductor {
    fn eval_bsdf(&self, _dir_in: Vec3D, _dir_out: Vec3D) -> Vec3D {
        assert!(self.roughness.is_none()); // currently only smooth surfaces
        Vec3D::ZERO
    }

    fn pdf(&self, _dir_in: Vec3D, _dir_out: Vec3D) -> f32 {
        assert!(self.roughness.is_none()); // currently only smooth surfaces
        0.0
    }

    fn sample_bsdf(&self, dir_out: Vec3D, _rng: &mut Rng) -> BSDFSample {
        assert!(self.roughness.is_none()); // currently only smooth surfaces

        // Perfect specular reflection!
        let dir_in = Vec3D::new(-dir_out.x, -dir_out.y, dir_out.z);
        let cos_in = dir_in.z.abs();
        let spectrum = complex_fresnell_reflectance_spectrum(
            cos_in,
            self.refractive_index,
            self.extinction_coefficient,
        ) / cos_in;

        return BSDFSample {
            spectrum,
            pdf: 1.0,
            dir_in,
        };
    }
}

#[inline]
fn refract(dir_in: Vec3D, mut eta: f32) -> (Vec3D, f32) {
    // eta is the relative index of refraction: eta_- / eta_+
    // where eta_- is the refractive index of the dielectric which exist
    // on the negative normal direction, and eta_+ is that of the outside
    // in the positive normal direction.

    // The following calculations assumes that eta is: eta_t / eta_in
    // where eta_in is the index of refraction on the incident side
    // and eta_t is that on the transmitted side.

    // If cos_in is negative, the incident ray is currently on the inside of
    // the dielectric medium, so we must flip the given eta.
    let mut cos_in = dir_in.z; // Signed cosine
    let mut normal = 1.0;
    if cos_in < 0.0 {
        eta = eta.recip();
        cos_in = -cos_in;
        normal = -1.0;
    }

    let sin_in_2 = (1.0 - cos_in * cos_in).max(0.0);
    let sin_out_2 = sin_in_2 / (eta * eta);
    let cos_out = f32::sqrt(1.0 - sin_out_2);

    let mut dir_out = -dir_in / eta;
    dir_out.z += (cos_in / eta - cos_out) * normal;

    return (dir_out, eta);
}

#[inline]
fn fresnell_reflectance(cos_in: f32, mut eta: f32) -> f32 {
    // eta is the relative index of refraction: eta_- / eta_+
    // where eta_- is the refractive index of the dielectric which exist
    // on the negative normal direction, and eta_+ is that of the outside
    // in the positive normal direction.

    // The following calculations assumes that eta is: eta_t / eta_in
    // where eta_in is the index of refraction on the incident side
    // and eta_t is that on the transmitted side.

    // If cos_in is negative, the incident ray is currently on the inside of
    // the dielectric medium, so we must flip the given eta.
    let mut cos_in = cos_in.clamp(-1.0, 1.0);
    if cos_in < 0.0 {
        eta = eta.recip();
        cos_in = -cos_in;
    }

    // Compute "angles" using Snells law
    let sin_in_2 = 1.0 - cos_in * cos_in;
    let sin_t_2 = sin_in_2 / (eta * eta);

    // Handle total internal reflection
    if sin_t_2 > 1.0 {
        return 1.0;
    }

    let cos_t = f32::sqrt(1.0 - sin_t_2);

    // Compute reflection coefficients
    let r_parallel = (eta * cos_in - cos_t) / (eta * cos_in + cos_t);
    let r_perpendicular = (cos_in - eta * cos_t) / (cos_in + eta * cos_t);

    // Assuming unpolarised light
    let r = (r_parallel * r_parallel + r_perpendicular * r_perpendicular) / 2.0;
    return r;
}

#[inline]
fn complex_fresnell_reflectance_spectrum(
    cos_in: f32,
    refractive_index: Vec3D,
    extinction_coefficient: Vec3D,
) -> Vec3D {
    [0, 1, 2]
        .map(|i| Complex {
            real: refractive_index[i],
            imag: extinction_coefficient[i],
        })
        .map(|eta| complex_fresnell_reflectance(cos_in, eta))
        .into()
}

#[inline]
fn complex_fresnell_reflectance(cos_in: f32, eta: Complex) -> f32 {
    // eta is the relative index of refraction: eta_t / eta_in
    // where eta_in is the index of refraction on the incident side
    // and eta_t is that on the transmitted side.
    // Light does not propagate on the inside of a conductor, so
    // we never have to flip eta.

    let cos_in = cos_in.clamp(0.0, 1.0);

    // Compute "angles" using Snells law
    let sin_in_2 = 1.0 - cos_in * cos_in;
    let sin_t_2 = sin_in_2 / (eta * eta);
    let cos_t = Complex::sqrt(1.0 - sin_t_2);

    // Compute reflection coefficients
    let r_parallel = (eta * cos_in - cos_t) / (eta * cos_in + cos_t);
    let r_perpendicular = (cos_in - eta * cos_t) / (cos_in + eta * cos_t);

    // Assuming unpolarised light
    let r = (r_parallel.abs_square() + r_perpendicular.abs_square()) / 2.0;
    return r;
}
