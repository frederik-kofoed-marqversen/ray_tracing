use fastrand::Rng;
use crate::Object;
use crate::materials::BSDF;
use crate::Vec3D;
use crate::Ray;

pub struct Material {
    pub emission: Option<Box<dyn Emission>>,
    pub bsdf: Option<Box<dyn BSDF>>,
}

pub trait Emission {
    fn sample(&self, obj: &Object, ref_point: Vec3D, rng: &mut Rng) -> LightSample;
    fn pdf(&self, obj: &Object, ray: &Ray) -> f32;
}

#[derive(Debug)]
pub struct LightSample {
    pub spectrum: Vec3D,
    pub point: Vec3D,
    pub direction: Vec3D,
    pub pdf: f32,
}

#[derive(Debug)]
pub struct Spotlight {
    pub pos: Vec3D,
    pub dir: Vec3D,
    pub cos_angle: f32,
    pub intensity: Vec3D,
}

impl Emission for Spotlight {
    fn sample(&self, _obj: &Object, ref_point: Vec3D, _rng: &mut Rng) -> LightSample {
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

    fn pdf(&self, _obj: &Object, _ray: &Ray) -> f32 {
        0.0
    }
}

// impl Light for Spotlight {
//     fn sample(&self, ref_point: Vec3D, _: &mut Rng) -> LightSample {
//         let direction = (self.pos - ref_point).normalise();
//         let d2 = direction.norm_squared();
//         let direction = direction / d2.sqrt();
//         let mut sample = LightSample {
//             spectrum: Vec3D::ZERO,
//             point: self.pos,
//             direction,
//             pdf: 1.0,
//         };

//         if -Vec3D::dot(direction, self.dir) > self.cos_angle {
//             sample.spectrum = self.intensity / d2;
//         }

//         return sample;
//     }
// }