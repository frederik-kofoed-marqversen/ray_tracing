#[inline]
pub fn safe_sqrt(v: f32) -> f32 {
    if v <= 0.0 {
        0.0
    } else {
        v.sqrt()
    }
}
