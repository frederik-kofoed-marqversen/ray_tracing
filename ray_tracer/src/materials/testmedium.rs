use crate::materials::medium::{
    sample_henyey_greenstein, ConstantMajorantSegment, Medium, MediumProperties,
};
use crate::Ray;
use crate::Vec3D;
use fastrand::Rng;

/// 3D checkerboard heterogeneous medium:
/// - Space is divided into cubic cells of size `cell_size`.
/// - Even cells ("white") and odd cells ("black") have different coefficients and emission.
/// - majorant_iterator walks ray-cell intersections and yields per-cell constant majorant segments.
pub struct CheckeredMedium {
    pub cell_size: f32,

    // white cell coefficients
    pub sigma_a_white: f32,
    pub sigma_s_white: f32,
    pub emission_white: Vec3D,

    // black cell coefficients
    pub sigma_a_black: f32,
    pub sigma_s_black: f32,
    pub emission_black: Vec3D,

    // Henyey–Greenstein asymmetry parameter
    pub g: f32,
}

impl CheckeredMedium {
    pub fn new(
        cell_size: f32,
        sigma_a_white: f32,
        sigma_s_white: f32,
        emission_white: Vec3D,
        sigma_a_black: f32,
        sigma_s_black: f32,
        emission_black: Vec3D,
        g: f32,
    ) -> Self {
        Self {
            cell_size,
            sigma_a_white,
            sigma_s_white,
            emission_white,
            sigma_a_black,
            sigma_s_black,
            emission_black,
            g,
        }
    }

    #[inline]
    fn cell_index(&self, p: Vec3D) -> [i32; 3] {
        let s = self.cell_size;
        [
            (p.x / s).floor() as i32,
            (p.y / s).floor() as i32,
            (p.z / s).floor() as i32,
        ]
    }

    #[inline]
    fn is_white(&self, idx: [i32; 3]) -> bool {
        (idx[0] + idx[1] * 0 + idx[2]) & 1 == 0
    }

    #[inline]
    fn props_at_cell(&self, idx: [i32; 3]) -> (f32, f32, Vec3D) {
        if self.is_white(idx) {
            (self.sigma_a_white, self.sigma_s_white, self.emission_white)
        } else {
            (self.sigma_a_black, self.sigma_s_black, self.emission_black)
        }
    }
}

impl Medium for CheckeredMedium {
    fn evaluate(&self, p: Vec3D, _lambda: usize) -> MediumProperties {
        let idx = self.cell_index(p);
        let (sa, ss, em) = self.props_at_cell(idx);
        MediumProperties {
            absorbtion_coefficient: sa, // if you treat Vec3D spectrally, adapt this
            scattering_coefficient: ss,
            emission: em,
        }
    }

    fn majorant_iterator(
        &self,
        ray: &Ray,
        t_max: f32,
    ) -> Box<dyn Iterator<Item = ConstantMajorantSegment> + '_> {
        // Voxel (cell) traversal using a simplified 3D DDA.
        // We produce segments [t_start, t_end) each lying inside one cell, with sigma_maj
        // equal to max-component of (sigma_a + sigma_s) for that cell.
        struct Iter<'m> {
            m: &'m CheckeredMedium,
            ray: Ray,
            t: f32,
            t_max: f32,
            // current cell indices
            idx: [i32; 3],
            // stepping
            step: [i32; 3],
            // t to next boundary along each axis
            t_next: [f32; 3],
            // t delta per cell step along each axis
            t_delta: [f32; 3],
            finished: bool,
        }

        impl<'m> Iter<'m> {
            fn new(m: &'m CheckeredMedium, ray: &Ray, t_max: f32) -> Self {
                let s = m.cell_size;
                let origin = ray.origin;
                let dir = ray.direction;

                // current cell indices at t = 0
                let idx = m.cell_index(origin);

                // step direction per axis
                let step = [
                    if dir.x >= 0.0 { 1 } else { -1 },
                    if dir.y >= 0.0 { 1 } else { -1 },
                    if dir.z >= 0.0 { 1 } else { -1 },
                ];

                // compute initial t_next (distance to first boundary) and t_delta (distance per cell)
                let mut t_next = [f32::INFINITY; 3];
                let mut t_delta = [f32::INFINITY; 3];

                for a in 0..3 {
                    let (o, d) = match a {
                        0 => (origin.x, dir.x),
                        1 => (origin.y, dir.y),
                        _ => (origin.z, dir.z),
                    };
                    if d.abs() > 0.0 {
                        // boundary coordinate along axis
                        let boundary = (idx[a] + if step[a] > 0 { 1 } else { 0 }) as f32 * s;
                        let t_to_boundary = (boundary - o) / d;
                        t_next[a] = if t_to_boundary >= 0.0 {
                            t_to_boundary
                        } else {
                            0.0
                        };
                        t_delta[a] = (s / d.abs()).max(0.0);
                    }
                }

                Self {
                    m,
                    ray: *ray,
                    t: 0.0,
                    t_max,
                    idx,
                    step,
                    t_next,
                    t_delta,
                    finished: false,
                }
            }
        }

        impl<'m> Iterator for Iter<'m> {
            type Item = ConstantMajorantSegment;
            fn next(&mut self) -> Option<Self::Item> {
                if self.finished {
                    return None;
                }
                if self.t >= self.t_max {
                    self.finished = true;
                    return None;
                }

                // choose axis with smallest t_next
                let mut axis = 0;
                if self.t_next[1] < self.t_next[axis] {
                    axis = 1;
                }
                if self.t_next[2] < self.t_next[axis] {
                    axis = 2;
                }

                // next boundary along chosen axis
                let t_end = self.t_next[axis].min(self.t_max);

                // current cell props
                let (sa, ss, _em) = self.m.props_at_cell(self.idx);
                let sigma_t = sa + ss;
                let sigma_maj = sigma_t; // conservative per-cell majorant

                let seg = ConstantMajorantSegment {
                    t_start: self.t,
                    t_end,
                    sigma_majorant: sigma_maj,
                };

                // advance to next cell along chosen axis
                if self.t_next[axis].is_finite() {
                    self.idx[axis] += self.step[axis];
                    self.t_next[axis] += self.t_delta[axis];
                } else {
                    // ray parallel to this axis; push to t_max next time
                    self.t_next[axis] = self.t_max;
                }

                self.t = t_end;
                Some(seg)
            }
        }

        Box::new(Iter::new(self, ray, t_max))
    }

    fn sample_scattering_direction(&self, _: Vec3D, wo: Vec3D, rng: &mut Rng) -> (Vec3D, f32, f32) {
        sample_henyey_greenstein(self.g, wo, rng)
    }
}
