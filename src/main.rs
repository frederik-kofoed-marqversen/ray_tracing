#![feature(trait_upcasting)]
#![allow(incomplete_features)]

extern crate ray_tracer;
use ray_tracer::traits::{Surface, BoundedSurface};
use ray_tracer::{Object, Point3D, Vec3D};
use ray_tracer::materials::*;
use ray_tracer::primitives::*;

const ASPECT_RATIO: f32 = 16.0 / 9.0;
const IMAGE_WIDTH: usize = 400;
const SAMPLES_PER_PIXEL: u32 = 64;
const RAY_DEPTH: u32 = 4;

fn main() -> std::io::Result<()> {
    // World setup
    let mut world = Vec::new();
    
    // Small objects
    world.push(Object{
        surface: Box::new(Sphere::new(Point3D::new(-2.0, 0.0, 0.5), 0.5)),
        material: Material::diffuse(Colour::new(1.0, 0.01, 0.01))
    });
    world.push(Object{
        surface: Box::new(Sphere::new(Point3D::new(-0.5, 0.0, 0.8), 0.75)),
        // material: Material::diffuse(Colour::e2())
        material: Material::metal(Colour::new(1.0, 0.5, 0.5), 0.5) 
    });
    world.push(Object{
        surface: Box::new(Sphere::new(Point3D::new(1.5, 0.0, 1.0), 1.0)),
        material: Material::diffuse(Colour::new(0.01, 0.01, 1.0))
    });
    world.push(Object{
        surface: Box::new(Sphere::new(Point3D::new(2.5, -2.5, 0.5), 0.5)),
        material: Material::light_source(Colour::new(0.01, 1.0, 0.01), 2.0)
    });
    world.push(Object{
        surface: Box::new(Triangle::new(
            Point3D::new(-4.5, 2.0, 0.2),
            Point3D::new(-2.0, 3.0, 0.2),
            Point3D::new(-3.0, 2.5, 3.2),
        )),
        material: Material::mirror()
    });

    // Globe
    world.push(Object{
        surface: Box::new(Sphere::new(Point3D::new(0.0, 0.0, -100.0), 100.0)),
        material: Material::diffuse(Colour::ones()*0.5)
    });

    // Light
    world.push(Object{
        surface: Box::new(Sphere::new(Point3D::new(0.0, 110.0, 0.0), 100.0)),
        material: Material::light_source(Colour::ones(), 4.0)
    });

    // Camera
    let camera = ray_tracer::Camera::new(
        Point3D::new(-2.5, -4.0, 3.0),
        Vec3D::new(0.5, 1.0, -0.5),
        ASPECT_RATIO
    );

    /* let triangles = ray_tracer::stl::read_stl("./meshes/Baby_Yoda.stl")?;
    let boxed_triangles = triangles.into_iter().map(|obj| Box::new(obj) as Box<dyn BoundedSurface>).collect();
    let bvh = ray_tracer::accelerators::BVH::build(boxed_triangles);
    let baby_yoda = Object { surface: bvh, material: Material::diffuse(Colour::new(0.3, 0.3, 0.3)) };
    let camera = ray_tracer::Camera::new(
        Point3D::new(0.0, 50.0, 70.0),
        Vec3D::new(0.0, -1.0, -0.3),
        ASPECT_RATIO
    );
    let world = vec![baby_yoda]; */

    // Render image
    let eng = ray_tracer::engine::Engine::new(world, camera);
    return eng.render(ASPECT_RATIO, IMAGE_WIDTH, SAMPLES_PER_PIXEL, RAY_DEPTH)
}
