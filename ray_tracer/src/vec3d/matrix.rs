use super::Vec3D;

#[derive(Debug, Clone, Copy)]
pub struct Matrix {
    data: [[f32; 3]; 3],
}

impl Matrix {
    pub fn identity() -> Self {
        Self {
            data: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    #[inline]
    pub fn from_diagonal(diag: Vec3D) -> Self {
        Self {
            data: [[diag.x, 0.0, 0.0], [0.0, diag.y, 0.0], [0.0, 0.0, diag.z]],
        }
    }

    #[inline]
    pub fn diagonal(&self) -> Vec3D {
        [0, 1, 2].map(|i| self.data[i][i]).into()
    }

    #[inline]
    pub fn mul(&self, vec: Vec3D) -> Vec3D {
        let m = self.data;
        Vec3D {
            x: m[0][0] * vec.x + m[0][1] * vec.y + m[0][2] * vec.z,
            y: m[1][0] * vec.x + m[1][1] * vec.y + m[1][2] * vec.z,
            z: m[2][0] * vec.x + m[2][1] * vec.y + m[2][2] * vec.z,
        }
    }

    /// Computes M^T * v
    #[inline]
    pub fn transpose_mul(&self, vec: Vec3D) -> Vec3D {
        let m = self.data;
        Vec3D {
            x: m[0][0] * vec.x + m[1][0] * vec.y + m[2][0] * vec.z,
            y: m[0][1] * vec.x + m[1][1] * vec.y + m[2][1] * vec.z,
            z: m[0][2] * vec.x + m[1][2] * vec.y + m[2][2] * vec.z,
        }
    }

    #[inline]
    pub fn determinant(&self) -> f32 {
        let m = self.data;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    pub fn inverse(&self) -> Option<Self> {
        let det = self.determinant();
        if det.abs() < 1e-9 {
            return None;
        }
        let inv_det = det.recip();
        let m = self.data;
        let inv_matrix = Matrix {
            data: [
                [
                    (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
                    (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
                    (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
                ],
                [
                    (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
                    (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
                    (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
                ],
                [
                    (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
                    (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
                    (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
                ],
            ],
        };

        return Some(inv_matrix);
    }

    pub fn transpose(&self) -> Self {
        let mut data = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                data[i][j] = self.data[j][i];
            }
        }

        Self { data }
    }

    pub fn matrix_mul(lhs: Self, rhs: Self) -> Self {
        let mut data = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                data[i][j] = lhs.data[i][0] * rhs.data[0][j]
                    + lhs.data[i][1] * rhs.data[1][j]
                    + lhs.data[i][2] * rhs.data[2][j];
            }
        }

        Self { data }
    }

    /// Rotation about arbitrary axis (Rodrigues' rotation formula)
    pub fn rotation(axis: Vec3D, angle: f32) -> Self {
        let axis = axis.normalise();
        let cos = angle.cos();
        let sin = angle.sin();
        let one_minus_cos = 1.0 - cos;

        let x = axis.x;
        let y = axis.y;
        let z = axis.z;

        let data = [
            [
                cos + x * x * one_minus_cos,
                x * y * one_minus_cos - z * sin,
                x * z * one_minus_cos + y * sin,
            ],
            [
                y * x * one_minus_cos + z * sin,
                cos + y * y * one_minus_cos,
                y * z * one_minus_cos - x * sin,
            ],
            [
                z * x * one_minus_cos - y * sin,
                z * y * one_minus_cos + x * sin,
                cos + z * z * one_minus_cos,
            ],
        ];

        Self { data }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse() {
        let matrix = Matrix {
            data: [[1.0, 2.0, 3.0], [0.0, 4.0, 5.0], [0.0, 0.0, 6.0]]
        };
        let inverse = matrix.inverse().unwrap();
        let identity = Matrix::matrix_mul(matrix, inverse);

        let mut test = true;
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                if (identity.data[i][j] - expected).abs() > 1e-9 {
                    test = false;
                }
            }
        }
        assert!(test);
    }
}
