use std::f32::consts::PI;
use std::rc::Rc;

extern crate ray_tracer;

use ray_tracer::bvh::BoundingVolumeHierarchy as BVH;
use ray_tracer::materials::*;
use ray_tracer::primitives::*;
use ray_tracer::subdivision_surface as sds;
use ray_tracer::triangle_mesh::{Triangle, TriangleMesh};
use ray_tracer::AffineTransform as Affine;
use ray_tracer::Camera;
use ray_tracer::{Object, Vec3D};

const ASPECT_RATIO: f32 = 16.0 / 9.0;
const IMAGE_WIDTH: usize = 400;
const SAMPLES_PER_PIXEL: u32 = 100;
const RAY_DEPTH: u32 = 2;

fn main() -> std::io::Result<()> {
    let (scene, camera) = scene_tetrahedra();
    // Render image
    let eng = ray_tracer::engine::Engine::new(scene, camera);
    return eng.render(ASPECT_RATIO, IMAGE_WIDTH, SAMPLES_PER_PIXEL, RAY_DEPTH);
}

#[allow(dead_code)]
fn scene_tetrahedra() -> (Vec<Object>, Camera) {
    // Compute subdivided tetrahedrons
    let tetrahedron = TriangleMesh {
        vertices: vec![
            Vec3D::new(1.0, 1.0, 1.0),
            Vec3D::new(1.0, -1.0, -1.0),
            Vec3D::new(-1.0, 1.0, -1.0),
            Vec3D::new(-1.0, -1.0, 1.0),
        ],
        triangles: vec![[0, 1, 2], [1, 3, 2], [0, 3, 1], [0, 2, 3]],
    };
    let tetrahedron_material = Material {
        albedo: Vec3D::new(1.0, 0.5, 0.5),
        roughness: 0.7,
        emission: None,
        refractive_index: None,
    };
    let subdivisions = [0, 1, 2, 6];
    let separation = 0.65 * Vec3D::Y;

    // Ground to cast shadows onto
    let ground = Object::new(
        Rc::new(Plane::new(-0.6 * Vec3D::Z, Vec3D::Z)),
        Material::diffuse(Vec3D::ONES * 0.8),
        None,
    );

    // Sun for lighting
    let light = Object::new(
        Rc::new(Sphere::new(200.0 * Vec3D::new(10.0, 2.0, 10.0), 1000.0)),
        Material::light_source(Vec3D::ONES, 6.0),
        None,
    );

    // Camera
    let loc = Vec3D::new_polar(1.8, 0.3 * PI, 0.0 * PI);
    let centre = separation * 0.5 * (subdivisions.len() - 1) as f32;
    let camera = ray_tracer::Camera::new(loc + centre, -loc, ASPECT_RATIO);

    // Prepare scene
    let mut scene = vec![ground, light];

    // Compute subdivided tetrahedron meshes
    let steps: Vec<usize> = subdivisions
        .iter()
        .zip(&subdivisions[1..])
        .map(|(a, b)| b - a)
        .collect();
    let sdsurface = sds::SubdivisionSurface::from_triangle_mesh(&tetrahedron);
    let mut tetrahedra = vec![(0..subdivisions[0]).fold(sdsurface, |res, _| res.subdivide())];
    for (i, n) in steps.into_iter().enumerate() {
        let mut sdsurface = tetrahedra[i].clone();
        for _ in 0..n {
            sdsurface = sdsurface.subdivide();
        }
        tetrahedra.push(sdsurface);
    }
    // Push to limit surface
    tetrahedra
        .iter_mut()
        .for_each(|sds| sds.push_to_limit_surface());

    // Add tetrahedron objects to scene
    for (i, sds) in tetrahedra.iter().enumerate() {
        let transform = Affine::translate(i as f32 * separation);
        let mesh = sds.to_triangle_mesh();
        let bvh = BVH::build(mesh.to_triangles().1);

        scene.push(Object::new(
            Rc::new(bvh),
            tetrahedron_material,
            Some(transform),
        ));
    }

    return (scene, camera);
}

#[allow(dead_code)]
fn scene_teapot() -> (Vec<Object>, Camera) {
    let teapot = ray_tracer::stl::read_stl("./meshes/Utah_teapot.stl").unwrap();
    let teapot = sds::loop_subdivide(&teapot, 1);
    let teapot_material = Material::diffuse(Vec3D::new(1.0, 0.5, 0.5));

    // Ground to cast shadows onto
    let ground = teapot
        .vertices
        .iter()
        .map(|v| v.z)
        .fold(f32::INFINITY, |a, b| a.min(b));
    let ground = Object::new(
        Rc::new(Plane::new(Vec3D::Z * ground, Vec3D::Z)),
        Material::diffuse(Vec3D::ONES * 0.8),
        None,
    );

    // Sun for lighting
    let light = Object::new(
        Rc::new(Sphere::new(2000.0 * Vec3D::new(1.0, 0.0, 1.0), 1000.0)),
        Material::light_source(Vec3D::ONES, 6.0),
        None,
    );

    // Camera
    let angle = 0.5 * PI;
    let loc = Vec3D::new_polar(18.0, 0.3 * PI, angle);
    let mut camera = ray_tracer::Camera::new(loc, Vec3D::ONES, ASPECT_RATIO);
    camera.look_at(teapot.volume_centroid());

    // Prepare scene for rendering
    let (_, triangles) = teapot.to_triangles();
    let bvh = BVH::build(triangles);
    let teapot = Object::new(Rc::new(bvh), teapot_material, None);
    let scene = vec![teapot, ground, light];
    return (scene, camera);
}

#[allow(dead_code)]
fn scene_baby_yoda() -> (Vec<Object>, Camera) {
    // Baby Yoda (Grogu)
    let grogu = ray_tracer::stl::read_stl("./meshes/Baby_Yoda.stl").unwrap();
    // Unfortunately the .stl file is not good and does not define a manifold
    // mesh.try_fix_mesh();
    // let grogu = loop_subdivide(&grogu, 1);
    // let grogu_material = Material::diffuse(Vec3D::fill(0.1));
    let grogu_material = Material::dielectric(Vec3D::new(0.7, 0.7, 0.9), 1.5);

    // Ground to cast shadows onto
    let ground = grogu
        .vertices
        .iter()
        .fold(f32::INFINITY, |min, v| min.min(v.z));
    let ground = Object::new(
        Rc::new(Plane::new(Vec3D::Z * ground, Vec3D::Z)),
        Material::diffuse(Vec3D::ONES * 0.8),
        None,
    );

    // Set up camera
    let centre = grogu.volume_centroid();
    let angle = 0.4 * PI;
    let dir = Vec3D::new_polar(1.0, PI / 2.0, angle);
    let mut camera = ray_tracer::Camera::new(
        centre + 90.0 * dir + 40.0 * Vec3D::Z,
        Vec3D::ONES,
        ASPECT_RATIO,
    );
    camera.look_at(centre);

    // Some directions
    let up = Vec3D::Z;
    let mut away = camera.get_direction();
    let left = Vec3D::cross(Vec3D::Z, away);
    away.z = 0.0;

    // Set up a sun to light up the scene
    let unit_sphere = Rc::new(Sphere::new(Vec3D::ZERO, 1.0));

    let light = Object::new(
        unit_sphere.clone(),
        Material::light_source(Vec3D::ONES, 6.0),
        Some(Affine::scale_translate(
            1000.0,
            centre + 2000.0 * up + 2000.0 * left - 1000.0 * away,
        )),
    );

    // An incandescent ball
    let ball = Object::new(
        unit_sphere.clone(),
        Material {
            albedo: Vec3D::new(0.8, 0.2, 0.2),
            roughness: 0.7,
            emission: Some((Vec3D::new(1.0, 0.03, 0.03), 1.0)),
            refractive_index: None,
        },
        Some(Affine::scale_translate(5.0, Vec3D::new(-4.0, 6.0, 59.0))),
    );

    // Prepare scene for rendering
    let (_, triangles) = grogu.to_triangles();
    let bvh = BVH::build(triangles);
    let grogu = Object::new(Rc::new(bvh), grogu_material, None);
    let scene = vec![grogu, ground, light, ball];

    return (scene, camera);
}

#[allow(dead_code)]
fn scene_primitives() -> (Vec<Object>, Camera) {
    // World setup
    let mut scene = Vec::new();

    // Small objects
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(-2.0, 0.0, 0.5), 0.5)),
        Material::diffuse(Vec3D::new(1.0, 0.01, 0.01)),
        None,
    ));
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(-0.5, 0.0, 0.8), 0.75)),
        // Material::diffuse(Vec3D::e2())
        Material::metal(Vec3D::new(1.0, 0.5, 0.5), 0.5),
        None,
    ));
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(1.5, 0.0, 1.0), 1.0)),
        Material::diffuse(Vec3D::new(0.01, 0.01, 1.0)),
        None,
    ));
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(2.5, -2.5, 0.5), 0.5)),
        Material::light_source(Vec3D::new(0.01, 1.0, 0.01), 2.0),
        None,
    ));
    scene.push(Object::new(
        Rc::new(Triangle::new_lonely(
            Vec3D::new(-4.5, 2.0, 0.2),
            Vec3D::new(-2.0, 3.0, 0.2),
            Vec3D::new(-3.0, 2.5, 3.2),
        )),
        Material::mirror(),
        None,
    ));
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(0.0, -1.6, 0.8), 0.75)),
        Material::dielectric(Vec3D::new(1.0, 1.0, 1.0), 1.5),
        None,
    ));

    // Globe
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(0.0, 0.0, -100.0), 100.0)),
        Material::diffuse(Vec3D::ONES * 0.5),
        None,
    ));

    // Light
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(0.0, 110.0, 0.0), 100.0)),
        Material::light_source(Vec3D::ONES, 3.0),
        None,
    ));

    // Camera
    let camera = ray_tracer::Camera::new(
        Vec3D::new(-2.5, -4.0, 3.0),
        Vec3D::new(0.5, 1.0, -0.5),
        ASPECT_RATIO,
    );

    return (scene, camera);
}
