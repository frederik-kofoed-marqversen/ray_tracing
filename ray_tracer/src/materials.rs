use super::Vec3D;

#[derive(Default, Clone, Copy)]
pub struct Material {
    pub albedo: Vec3D,
    pub roughness: f32,
    pub emission: Option<(Vec3D, f32)>,
    pub refractive_index: Option<f32>,
}

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