use std::rc::Rc;

pub mod complex;
pub mod vec3d;
use vec3d::Matrix;
pub use vec3d::Vec3D;
pub mod bsdf;
pub mod bvh;
pub mod lights;
use lights::Light;
pub mod path_integrator;
pub mod primitives;
pub mod stl;
pub mod subdivision_surface;
pub mod traits;
pub mod triangle_mesh;
pub mod utils;

#[derive(Debug, Copy, Clone)]
pub struct Ray {
    pub origin: Vec3D,
    pub direction: Vec3D,
}

impl Ray {
    #[inline]
    pub fn new(origin: Vec3D, direction: Vec3D) -> Self {
        Self { origin, direction }
    }

    #[inline]
    pub fn at(&self, t: f32) -> Vec3D {
        self.origin + self.direction * t
    }
}

pub struct Object {
    surface: Rc<dyn traits::Surface>,
    material: Rc<dyn bsdf::BSDF>,
    affine_transform: Option<AffineTransform>,
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

#[derive(Debug, Clone, Copy)]
pub struct AffineTransform {
    matrix: Matrix,
    translation: Vec3D,
}

impl AffineTransform {
    #[inline]
    pub fn inverse(&self) -> Self {
        // T(v)^{-1} = M^{-1}(v - d) = M^{-1}v - M^{-1}d
        let matrix = self.matrix.inverse().unwrap();
        let translation = -matrix.mul(self.translation);
        Self {
            matrix,
            translation,
        }
    }

    #[inline]
    pub fn composition(lhs: Self, rhs: Self) -> Self {
        // B(A(v)) = B(Av + a) + b = BAv + (Ba + b)
        let matrix = Matrix::matrix_mul(lhs.matrix, rhs.matrix);
        let translation = lhs.matrix.mul(rhs.translation) + lhs.translation;
        Self {
            matrix,
            translation,
        }
    }

    #[inline]
    pub fn transform_vector(&self, vector: Vec3D) -> Vec3D {
        self.matrix.mul(vector)
    }

    #[inline]
    pub fn transform_point(&self, vector: Vec3D) -> Vec3D {
        self.matrix.mul(vector) + self.translation
    }

    #[inline]
    pub fn scale(scale: f32) -> Self {
        return Self {
            matrix: Matrix::from_diagonal(Vec3D::fill(scale)),
            translation: Vec3D::ZERO,
        };
    }

    #[inline]
    pub fn rotation(axis: Vec3D, angle: f32) -> Self {
        return Self {
            matrix: Matrix::rotation(axis, angle),
            translation: Vec3D::ZERO,
        };
    }

    #[inline]
    pub fn translate(translation: Vec3D) -> Self {
        return Self {
            matrix: Matrix::identity(),
            translation: translation,
        };
    }

    #[inline]
    pub fn rotate_translate(axis: Vec3D, angle: f32, translation: Vec3D) -> Self {
        return Self {
            matrix: Matrix::rotation(axis, angle),
            translation: translation,
        };
    }

    #[inline]
    pub fn scale_translate(scale: f32, translation: Vec3D) -> Self {
        return Self {
            matrix: Matrix::from_diagonal(Vec3D::fill(scale)),
            translation: translation,
        };
    }
}

pub struct Camera {
    origin: Vec3D,
    direction: Vec3D,
    aspect_ratio: f32,
    horizontal: Vec3D,
    vertical: Vec3D,
    lower_left_corner: Vec3D,
}

impl Camera {
    pub fn new(origin: Vec3D, direction: Vec3D, aspect_ratio: f32) -> Self {
        let mut result = Self {
            origin,
            direction,
            aspect_ratio,
            horizontal: Vec3D::default(),
            vertical: Vec3D::default(),
            lower_left_corner: Vec3D::default(),
        };
        result.set_direction(direction);
        result
    }

    pub fn get_direction(&self) -> Vec3D {
        self.direction
    }

    pub fn get_pos(&self) -> Vec3D {
        self.origin
    }

    pub fn set_pos(&mut self, point: Vec3D) {
        let translation = point - self.origin;
        self.origin = point;
        self.lower_left_corner += translation;
    }

    pub fn set_direction(&mut self, direction: Vec3D) {
        let window_height: f32 = 2.0;
        let window_width: f32 = window_height * self.aspect_ratio;
        let focal_length: f32 = 2.0; // Horizontal FOV of 46 degrees

        let direction = direction.normalise();

        self.direction = direction;
        self.horizontal = Vec3D::cross(direction, Vec3D::Z).normalise() * window_width;
        self.vertical = Vec3D::cross(self.horizontal, direction).normalise() * window_height;
        self.lower_left_corner =
            self.origin + direction * focal_length - self.horizontal / 2.0 - self.vertical / 2.0;
    }

    pub fn look_at(&mut self, point: Vec3D) {
        self.set_direction(point - self.origin);
    }

    #[inline]
    pub fn ray(&self, u: f32, v: f32) -> Ray {
        let direction =
            self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin;
        Ray::new(self.origin, direction.normalise())
    }
}

/* #[cfg(test)]
mod tests {
    use super::bsdf::BSDF;
    use super::primitives::Sphere;

    use super::*;

    #[test]
    fn intersect_object() {
        let mut sphere = Object {
            affine_transform: None,
            surface: Rc::new(Sphere::new(Vec3D::ZERO, 0.5)),
            material: Material::diffuse(Vec3D::ONES),
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
} */
