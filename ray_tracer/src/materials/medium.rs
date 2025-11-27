use crate::Ray;
use fastrand::Rng;
use std::f32::consts::PI;

use crate::materials::bsdf::to_world_coords;
use crate::math::utils::{safe_sqrt, sample_exponential};
use crate::Vec3D;
use std::rc::Rc;

pub struct MediumInterface {
    pub outside: Option<Rc<dyn Medium>>,
    pub inside: Option<Rc<dyn Medium>>,
}

impl MediumInterface {
    pub fn new(inside: Option<Rc<dyn Medium>>, outside: Option<Rc<dyn Medium>>) -> Self {
        Self { outside, inside }
    }
}

pub struct ConstantMajorantSegment {
    pub t_start: f32,
    pub t_end: f32,
    pub sigma_majorant: f32,
}

pub struct MediumProperties {
    pub absorbtion_coefficient: f32,
    pub scattering_coefficient: f32,
    pub emission: Vec3D,
}

pub enum MediumInteraction {
    Null { transmittance: f32 },
    Absorption { emission: Vec3D },
    Scatter { ray: Ray, pdf: f32, f: f32 },
}

pub trait Medium {
    fn evaluate(&self, point: Vec3D, wavelength: usize) -> MediumProperties;

    fn majorant_iterator(
        &self,
        ray: &Ray,
        t_max: f32,
    ) -> Box<dyn Iterator<Item = ConstantMajorantSegment> + '_>;

    // Sample a scattering direction given outgoing direction wo at point in the medium.
    // Returns (wi, pdf, f) where f is the value of the phase function.
    fn sample_scattering_direction(
        &self,
        point: Vec3D,
        wo: Vec3D,
        rng: &mut Rng,
    ) -> (Vec3D, f32, f32);

    #[deprecated(note = "currently assume non-spectral media using only wavelength index 0")]
    fn sample_segment(
        &self,
        ray: &Ray,
        segment: &ConstantMajorantSegment,
        rng: &mut Rng,
    ) -> MediumInteraction {
        let majorant = segment.sigma_majorant;
        if majorant == 0.0 {
            return MediumInteraction::Null { transmittance: 1.0 };
        }

        let mut t = segment.t_start;
        // Propagate through segment in discrete steps until interaction or exit
        loop {
            let dt = sample_exponential(majorant, rng);
            t += dt;

            // Check if ray exits segment
            if t >= segment.t_end {
                let dt = segment.t_end - segment.t_start;
                let transmittance = (-majorant * dt).exp();
                return MediumInteraction::Null { transmittance };
            }

            // Ray has advanced by dt onto t which is within the segment
            // Sample interaction at point ray(t)
            // Should be one of absorption, scattering, or null collision
            let p = ray.at(t);
            // CURRENTLY ONLY SUPPORTING NON-SPECTRAL MEDIA
            let mp = self.evaluate(p, 0);

            let u = rng.f32();
            let p_absorption = mp.absorbtion_coefficient / majorant;
            let p_scattering = mp.scattering_coefficient / majorant;
            // let _p_null = 1.0 - p_absorption - p_scattering;
            if u < p_absorption {
                // Absorbtion
                return MediumInteraction::Absorption {
                    emission: mp.emission,
                };
            } else if u < p_absorption + p_scattering {
                // Scattering
                // In principle these computations could be omitted if depth limit is reached
                // with this scattering event, slightly improving performance
                let (wi, pdf, f) = self.sample_scattering_direction(p, ray.direction, rng);
                let ray = Ray {
                    origin: p,
                    direction: wi,
                };
                return MediumInteraction::Scatter { ray, pdf, f };
            } else {
                // Null collision
                // Continue propagation
                continue;
            }
        }
    }

    fn sample_interaction(&self, ray: &Ray, t_max: f32, rng: &mut Rng) -> (MediumInteraction, f32) {
        // Maybe normalise ray direction here although it already should be?

        let mut transmittance = 1.0;

        for segment in self.majorant_iterator(ray, t_max) {
            let interaction = self.sample_segment(ray, &segment, rng);
            match interaction {
                MediumInteraction::Null {
                    transmittance: segment_transmittance,
                } => {
                    transmittance *= segment_transmittance;
                    continue;
                }
                MediumInteraction::Absorption { .. } => return (interaction, transmittance),
                MediumInteraction::Scatter { .. } => return (interaction, transmittance),
            }
        }
        // Null scattering through the entire medium
        return (MediumInteraction::Null { transmittance }, transmittance);
    }
}

pub struct HomogeneousMedium {
    absorbtion_coefficient: Vec3D,
    scattering_coefficient: Vec3D,
    emission: Vec3D,
    // g parameter for Henyey-Greenstein phase function
    g: f32,
}

impl HomogeneousMedium {
    pub fn new(
        absorbtion_coefficient: Vec3D,
        scattering_coefficient: Vec3D,
        g: Option<f32>,
        emission: Vec3D,
    ) -> Self {
        Self {
            absorbtion_coefficient,
            scattering_coefficient,
            g: g.unwrap_or(0.0),
            emission,
        }
    }
}

impl Medium for HomogeneousMedium {
    fn evaluate(&self, _point: Vec3D, wavelength: usize) -> MediumProperties {
        MediumProperties {
            absorbtion_coefficient: self.absorbtion_coefficient[wavelength],
            scattering_coefficient: self.scattering_coefficient[wavelength],
            emission: self.emission,
        }
    }

    fn majorant_iterator(
        &self,
        _ray: &Ray,
        t_max: f32,
    ) -> Box<dyn Iterator<Item = ConstantMajorantSegment>> {
        let sigma_t = self.absorbtion_coefficient + self.scattering_coefficient;
        let sigma_maj = sigma_t.max_elem();

        let segment = ConstantMajorantSegment {
            t_start: 0.0,
            t_end: t_max,
            sigma_majorant: sigma_maj,
        };

        Box::new(std::iter::once(segment))
    }

    fn sample_scattering_direction(&self, _: Vec3D, wo: Vec3D, rng: &mut Rng) -> (Vec3D, f32, f32) {
        sample_henyey_greenstein(self.g, wo, rng)
    }
}

fn henyey_greenstein_phase_function(cos_theta: f32, g: f32) -> f32 {
    let denom = 1.0 + g * g - 2.0 * g * cos_theta;
    (1.0 - g * g) / (4.0 * PI * denom * safe_sqrt(denom))
}

fn sample_henyey_greenstein_direction(wo: Vec3D, g: f32, rng: &Rng) -> Vec3D {
    let u: f32 = rng.f32();
    let v: f32 = rng.f32();

    let cos_theta: f32;
    if g.abs() < 1e-3 {
        cos_theta = 1.0 - 2.0 * u;
    } else {
        let term: f32 = (1.0 - g * g) / (1.0 + g - 2.0 * g * u);
        cos_theta = -(1.0 + g * g - term * term) / (2.0 * g);
    }
    let sin_theta = safe_sqrt(1.0 - cos_theta * cos_theta);
    let phi = 2.0 * PI * v;

    // Sampled direction in coordinates where wo = (0,0,1)
    let wi = Vec3D::new(sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta);
    // Transform to world coordinates
    let wi = to_world_coords(wi, wo);
    return wi;
}

pub fn sample_henyey_greenstein(g: f32, wo: Vec3D, rng: &mut Rng) -> (Vec3D, f32, f32) {
    let wi = sample_henyey_greenstein_direction(wo, g, rng);
    let cos_theta = Vec3D::dot(wi, wo);
    let f = henyey_greenstein_phase_function(cos_theta, g);
    let pdf = f;
    return (wi, pdf, f);
}
