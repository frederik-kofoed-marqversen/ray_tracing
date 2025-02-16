pub type Colour = super::Vec3D;

#[derive(Default, Clone)]
pub struct Material {
    pub albedo: Colour,
    pub roughness: f32,
    pub emission: Option<(Colour, f32)>,
    pub refractive_index: Option<f32>,
}

impl Material {
    #[inline]
    pub fn diffuse(albedo: Colour) -> Self {
        Self {
            albedo,
            roughness: 1.0,
            emission: None,
            refractive_index: None,
        }
    }

    #[inline]
    pub fn light_source(emission: Colour, intensity: f32) -> Self {
        Self {
            albedo: Colour::ZERO,
            roughness: 1.0,
            emission: Some((emission, intensity)),
            refractive_index: None,
        }
    }

    #[inline]
    pub fn mirror() -> Self {
        Self {
            albedo: Colour::ONES,
            roughness: 0.0,
            emission: None,
            refractive_index: None,
        }
    }

    #[inline]
    pub fn metal(albedo: Colour, roughness: f32) -> Self {
        Self {
            albedo,
            roughness,
            emission: None,
            refractive_index: None,
        }
    }

    #[inline]
    pub fn dielectric(albedo: Colour, refractive_index: f32) -> Self {
        Self {
            albedo,
            roughness: 0.0,
            emission: None,
            refractive_index: Some(refractive_index),
        }
    }
}