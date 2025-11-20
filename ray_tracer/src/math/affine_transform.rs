use crate::math::matrix::Matrix;
use crate::Vec3D;

#[derive(Debug, Clone, Copy)]
pub struct AffineTransform {
    pub matrix: Matrix,
    pub translation: Vec3D,
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
