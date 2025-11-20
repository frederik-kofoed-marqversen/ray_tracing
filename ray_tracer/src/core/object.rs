use crate::materials::bsdf;
use crate::traits;
use crate::AffineTransform;
use crate::Ray;
use crate::Vec3D;
use std::rc::Rc;

pub struct Object {
    pub surface: Rc<dyn traits::Surface>,
    pub material: Rc<dyn bsdf::BSDF>,
    pub affine_transform: Option<AffineTransform>,
}

impl Object {
    #[inline]
    pub fn new(
        surface: Rc<dyn traits::Surface>,
        material: Rc<dyn bsdf::BSDF>,
        affine_transform: Option<AffineTransform>,
    ) -> Self {
        Self {
            surface,
            affine_transform: affine_transform.map(|at| at.inverse()),
            material,
        }
    }

    #[inline]
    pub fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)> {
        if self.affine_transform.is_none() {
            return self.surface.hit(&ray, t_min, t_max);
        }

        let transform = &self.affine_transform.unwrap();

        let origin = transform.matrix.mul(ray.origin) + transform.translation;
        let direction = transform.matrix.mul(ray.direction);
        let norm_recip = direction.norm_recip();
        let obj_ray = Ray {
            origin,
            direction: direction * norm_recip,
        };

        match self.surface.hit(&obj_ray, t_min, t_max) {
            None => return None,
            Some((obj_t, obj_normal)) => {
                // A normal vector n' in obj. space is defined as a vector orthogonal
                // to the tangent space. Let t' be a tangent vector in obj. space. Then
                // 0 = t'^T * n' = (M^{-1} * t)^T * n' = t^T * (M^{-T} * n')
                // So normal vectors in world space are given by (transpose transform):
                // n = M^{-T} * n'
                // The direction of a normal vector is flipped if det(M) < 0.
                let world_t = obj_t * norm_recip;
                let world_normal = transform.matrix.transpose_mul(obj_normal).normalise();
                return Some((world_t, world_normal));
            }
        }
    }

    #[inline]
    pub fn hit_bool(&self, ray: &Ray, t_min: f32, t_max: f32) -> bool {
        if self.affine_transform.is_none() {
            return self.surface.hit_bool(&ray, t_min, t_max);
        }

        let transform = &self.affine_transform.unwrap();

        let origin = transform.matrix.mul(ray.origin) + transform.translation;
        let direction = transform.matrix.mul(ray.direction);
        let norm_recip = direction.norm_recip();
        let obj_ray = Ray {
            origin,
            direction: direction * norm_recip,
        };

        return self.surface.hit_bool(&obj_ray, t_min, t_max);
    }

    pub fn apply_transform(&mut self, transform: AffineTransform) {
        // Invert transform since we store world-to-object transform. For affine transformations
        // one can show (B o A)^{-1} = A^{-1} o B^{-1}, so when adding a transform B to an object
        // we must add the inverse transform on the right.
        let transform = transform.inverse();
        let new_transform = match &self.affine_transform {
            None => transform,
            Some(old_transform) => AffineTransform::composition(*old_transform, transform),
        };
        self.affine_transform = Some(new_transform);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Sphere;
    use crate::math::matrix::Matrix;

    #[test]
    fn intersect_object() {
        let mut sphere = Object {
            affine_transform: None,
            surface: Rc::new(Sphere::new(Vec3D::ZERO, 0.5)),
            material: Rc::new(bsdf::Diffuse {
                reflectance: Vec3D::ONES,
            }),
        };

        let ray = Ray {
            origin: 10.0 * Vec3D::X,
            direction: -Vec3D::X,
        };

        let hit = sphere.hit(&ray, 0.0, 100.0).map(|(t, _)| t);
        assert_eq!(hit, Some(9.5));

        sphere.affine_transform = Some(AffineTransform {
            matrix: Matrix::identity(),
            translation: Vec3D::ZERO,
        });

        let hit = sphere.hit(&ray, 0.0, 100.0).map(|(t, _)| t);
        assert_eq!(hit, Some(9.5));

        sphere.affine_transform = Some(AffineTransform {
            matrix: Matrix::rotation(Vec3D::ONES, 1.0),
            translation: Vec3D::ZERO,
        });

        let hit = sphere.hit(&ray, 0.0, 100.0).map(|(t, _)| t);
        assert!(hit.is_some());

        sphere.affine_transform = Some(AffineTransform {
            matrix: Matrix::rotation(Vec3D::ONES, 1.0),
            translation: Vec3D::Z,
        });

        let hit = sphere.hit(&ray, 0.0, 100.0).map(|(t, _)| t);
        assert_eq!(hit, None);
    }
}
