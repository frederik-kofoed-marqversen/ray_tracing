pub type Colour = super::Vec3D;

#[derive(Default)]
pub struct Material {
    pub albedo: Colour,
    pub roughness: f32,
    pub emission: Colour,
    pub intensity: f32,
}

impl Material {
    #[inline]
    pub fn diffuse(albedo: Colour) -> Self {
        Self {
            albedo,
            roughness: 1.0,
            emission: Colour::zero(),
            intensity: 0.0,
        }
    }

    #[inline]
    pub fn light_source(emission: Colour, intensity: f32) -> Self {
        Self {
            albedo: Colour::zero(),
            roughness: 1.0,
            emission,
            intensity,
        }
    }

    #[inline]
    pub fn mirror() -> Self {
        Self {
            albedo: Colour::ones(),
            roughness: 0.0,
            emission: Colour::zero(),
            intensity: 0.0,
        }
    }

    #[inline]
    pub fn metal(albedo: Colour, roughness: f32) -> Self {
        Self {
            albedo,
            roughness,
            emission: Colour::zero(),
            intensity: 0.0,
        }
    }
}