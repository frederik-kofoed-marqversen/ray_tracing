mod vec3d;
mod ray;

use std::io::Write;

use ray::Ray;
use vec3d::Vec3D;

type Color = Vec3D;
type Point3D = Vec3D;

fn main() -> std::io::Result<()> {
    const HEIGHT: usize = 256;
    const WIDTH: usize = 256;

    let stdout = std::io::stdout();
    let mut lock = stdout.lock();
    let stderr = std::io::stderr();
    let mut err_lock = stderr.lock();

    write!(lock, "P3\n{WIDTH} {HEIGHT}\n255\n")?;
    for j in (0..HEIGHT).rev() {
        writeln!(err_lock, "Scanlines remaining: {j}")?;
        for i in 0..WIDTH {
            let r = i as f64 / (WIDTH - 1) as f64;
            let g = j as f64 / (HEIGHT - 1) as f64;
            let b = 0.25;

            let ir = (255.999 * r) as u32;
            let ig = (255.999 * g) as u32;
            let ib = (255.999 * b) as u32;

            writeln!(lock, "{ir} {ig} {ib}")?;
        }
    }
    writeln!(err_lock, "Done.")?;

    Ok(())
}
