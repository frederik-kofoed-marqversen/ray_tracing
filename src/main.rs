#![feature(trait_upcasting)]

extern crate ray_tracer;
use core::f32;
use ray_tracer::materials::*;
use ray_tracer::primitives::*;
use ray_tracer::subdivision_surface::loop_subdivide;
use ray_tracer::traits::BoundedSurface;
use ray_tracer::triangle_mesh::*;
use ray_tracer::{Object, Point3D, Vec3D};
use std::rc::Rc;

const ASPECT_RATIO: f32 = 16.0 / 9.0;
const IMAGE_WIDTH: usize = 400;
const SAMPLES_PER_PIXEL: u32 = 32;
const RAY_DEPTH: u32 = 2;

fn main() -> std::io::Result<()> {
    // Build a world and run the ray-tracing engine to output a bitmap by
    // running the command: `cargo +nightly run --release > renders/image.ppm`.
    // Install nightly with command: `rustup toolchain install nightly`.

    /* // World setup
    let mut world = Vec::new();

    // Small objects
    world.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(-2.0, 0.0, 0.5), 0.5)),
        material: Material::diffuse(Colour::new(1.0, 0.01, 0.01)),
    });
    world.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(-0.5, 0.0, 0.8), 0.75)),
        // material: Material::diffuse(Colour::e2())
        material: Material::metal(Colour::new(1.0, 0.5, 0.5), 0.5),
    });
    world.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(1.5, 0.0, 1.0), 1.0)),
        material: Material::diffuse(Colour::new(0.01, 0.01, 1.0)),
    });
    world.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(2.5, -2.5, 0.5), 0.5)),
        material: Material::light_source(Colour::new(0.01, 1.0, 0.01), 2.0),
    });
    world.push(Object {
        surface: Box::new(Triangle::new(
            Point3D::new(-4.5, 2.0, 0.2),
            Point3D::new(-2.0, 3.0, 0.2),
            Point3D::new(-3.0, 2.5, 3.2),
        )),
        material: Material::mirror(),
    });
    world.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(0.0, -1.6, 0.8), 0.75)),
        material: Material::dielectric(Colour::new(1.0, 1.0, 1.0), 1.5),
    });

    // Globe
    world.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(0.0, 0.0, -100.0), 100.0)),
        material: Material::diffuse(Colour::ONES * 0.5),
    });

    // Light
    world.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(0.0, 110.0, 0.0), 100.0)),
        material: Material::light_source(Colour::ONES, 3.0),
    });

    // Camera
    let camera = ray_tracer::Camera::new(
        Point3D::new(-2.5, -4.0, 3.0),
        Vec3D::new(0.5, 1.0, -0.5),
        ASPECT_RATIO,
    ); */

    /* // Use sky_colour() in engine as lightsource.
    let mesh = ray_tracer::stl::read_stl("./meshes/Baby_Yoda.stl")?;
    let mesh = Rc::new(mesh);
    let centre = mesh.vertices.iter().sum::<Vec3D>() / mesh.vertices.len() as f32;
    let ground_level = mesh
        .vertices
        .iter()
        .fold(f32::INFINITY, |min, p| min.min(p.z));

    let triangles = mesh.triangles();
    let boxed_triangles = triangles
        .into_iter()
        .map(|obj| Box::new(obj) as Box<dyn BoundedSurface>)
        .collect();
    let bvh = ray_tracer::accelerators::BVH::build(boxed_triangles);
    let baby_yoda = Object {
        surface: bvh,
        material: Material::dielectric(Colour::new(0.7, 0.7, 0.9), 1.5),
        // material: Material::diffuse(Colour::fill(0.1))
    };

    let ground = Object {
        surface: Box::new(Plane::new(Vec3D::new(0.0, 0.0, ground_level), Vec3D::Z)),
        material: Material::diffuse(Colour::ONES * 0.8),
    };

    let angle = 0.4 * std::f32::consts::PI;
    let dir = Vec3D::new_polar(1.0, std::f32::consts::PI / 2.0, angle);
    let mut camera = ray_tracer::Camera::new(
        centre + 90.0 * dir + 40.0 * Vec3D::Z,
        Vec3D::ONES,
        ASPECT_RATIO,
    );
    camera.look_at(centre);

    let up = Vec3D::Z;
    let left = Vec3D::cross(Vec3D::Z, camera.direction);
    let mut away = camera.direction.clone();
    away.z = 0.0;

    let light = Object {
        surface: Box::new(Sphere::new(
            centre + 2000.0 * up + 2000.0 * left - 1000.0 * away,
            1000.0,
        )),
        material: Material::light_source(Colour::ONES, 6.0),
    };

    let ball = Object {
        surface: Box::new(Sphere::new(Vec3D::new(-4.0, 6.0, 59.0), 5.0)),
        material: Material {
            albedo: Colour::new(0.8, 0.2, 0.2),
            roughness: 0.7,
            emission: Some((Colour::new(1.0, 0.03, 0.03), 1.0)),
            refractive_index: None,
        },
    };

    let world = vec![baby_yoda, ground, light, ball]; */

    let mesh = ray_tracer::stl::read_stl("./meshes/Tetrahedron.stl")?;
    let centre = mesh.vertices.iter().sum::<Vec3D>() / mesh.vertices.len() as f32;
    let ground_level = mesh
        .vertices
        .iter()
        .fold(f32::INFINITY, |min, p| min.min(p.z));
    let mesh = loop_subdivide(&mesh, 2);
    let mesh = Rc::new(mesh);

    let triangles = mesh.triangles();
    let boxed_triangles = triangles
        .into_iter()
        .map(|obj| Box::new(obj) as Box<dyn BoundedSurface>)
        .collect();
    let bvh = ray_tracer::accelerators::BVH::build(boxed_triangles);
    let object = Object {
        surface: bvh,
        material: Material::diffuse(Colour::fill(0.1)),
    };

    // let mut camera =
    //     ray_tracer::Camera::new(centre + 10.0 * Vec3D::e1(), Vec3D::ones(), ASPECT_RATIO);
    let mut camera = ray_tracer::Camera::new(
        centre + Vec3D::new_polar(10.0, 0.5*f32::consts::PI, 1.2 * f32::consts::PI),
        Vec3D::ONES,
        ASPECT_RATIO,
    );
    camera.look_at(centre);

    let light = Object {
        surface: Box::new(Sphere::new(
            centre + 2000.0 * Vec3D::Z + 2000.0 * Vec3D::Y,
            1000.0,
        )),
        material: Material::light_source(Colour::ONES, 6.0),
    };

    let ground = Object {
        surface: Box::new(Plane::new(
            Vec3D::new(0.0, 0.0, ground_level - 1.0),
            Vec3D::Z,
        )),
        material: Material::diffuse(Colour::fill(0.8)),
    };

    let world = vec![object, ground, light];

    // Render image
    let eng = ray_tracer::engine::Engine::new(world, camera);
    return eng.render(ASPECT_RATIO, IMAGE_WIDTH, SAMPLES_PER_PIXEL, RAY_DEPTH);
}
