use std::f32::consts::PI;
use std::rc::Rc;

extern crate ray_tracer;

use ray_tracer::bsdf::*;
use ray_tracer::lights::*;
use ray_tracer::primitives::*;
use ray_tracer::subdivision_surface as sds;
use ray_tracer::triangle_mesh::TriangleMesh;
use ray_tracer::AffineTransform as Affine;
use ray_tracer::Camera;
use ray_tracer::{Object, Vec3D};

const ASPECT_RATIO: f32 = 16.0 / 9.0;
const IMAGE_WIDTH: usize = 400;
const SAMPLES_PER_PIXEL: u32 = 64;
const RAY_DEPTH: u32 = 8;

fn main() -> std::io::Result<()> {
    let (scene, lights, camera) = scene_baby_yoda();
    // Render image
    let eng = ray_tracer::path_integrator::Engine::new(scene, lights, camera);
    return eng.render(ASPECT_RATIO, IMAGE_WIDTH, SAMPLES_PER_PIXEL, RAY_DEPTH);
}

#[allow(dead_code)]
fn scene_tetrahedra() -> (Vec<Object>, Vec<Rc<dyn Light>>, Camera) {
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
    let tetrahedron_material: Rc<dyn BSDF> = Rc::new(Diffuse {
        reflectance: Vec3D::new(1.0, 0.5, 0.5),
    });
    let subdivisions = [0, 1, 2, 6];
    let separation = 0.65 * Vec3D::Y;

    // Ground to cast shadows onto
    let ground = Object::new(
        Rc::new(Plane::new(-0.6 * Vec3D::Z, Vec3D::Z)),
        Rc::new(Diffuse {
            reflectance: Vec3D::ONES * 0.8,
        }),
        None,
    );

    // Sun for lighting
    let light: Rc<dyn Light> = Rc::new(PointLight {
        pos: 10.0 * Vec3D::new(10.0, 2.0, 10.0),
        intensity: Vec3D::ONES * 60000.0,
    });
    let lights = vec![light];

    // Camera
    let loc = Vec3D::new_polar(1.8, 0.3 * PI, 0.0 * PI);
    let centre = separation * 0.5 * (subdivisions.len() - 1) as f32;
    let camera = ray_tracer::Camera::new(loc + centre, -loc, ASPECT_RATIO);

    // Prepare scene
    let mut scene = vec![ground];

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
        let mesh = Rc::new(sds.to_triangle_mesh());
        let bvh = TriangleMesh::build_bvh(&mesh);
        scene.push(Object::new(
            bvh,
            tetrahedron_material.clone(),
            Some(transform),
        ));
    }

    return (scene, lights, camera);
}

#[allow(dead_code)]
fn scene_teapot() -> (Vec<Object>, Vec<Rc<dyn Light>>, Camera) {
    let teapot = ray_tracer::stl::read_stl("./meshes/utah_teapot.stl").unwrap();
    let teapot = Rc::new(sds::loop_subdivide(&teapot, 1));
    // let teapot_material = Rc::new(Dielectric {
    //     refractive_index: 1.5,
    //     roughness: None,
    // });
    let teapot_material = Rc::new(Diffuse {
        reflectance: Vec3D::new(0.9, 0.3, 0.3),
    });

    // Ground to cast shadows onto
    let ground = teapot
        .vertices
        .iter()
        .map(|v| v.z)
        .fold(f32::INFINITY, |a, b| a.min(b));
    let ground = Object::new(
        Rc::new(Plane::new(Vec3D::Z * ground, Vec3D::Z)),
        Rc::new(Diffuse {
            reflectance: Vec3D::ONES * 0.8,
        }),
        None,
    );

    // Set up camera
    let centre = teapot.volume_centroid();
    let angle = 0.5 * PI;
    let dir = Vec3D::new_polar(18.0, 0.3 * PI, angle);
    let mut camera = ray_tracer::Camera::new(centre + dir, Vec3D::ONES, ASPECT_RATIO);
    camera.look_at(centre);

    // Set up a sun to light up the scene
    let light: Rc<dyn Light> = Rc::new(PointLight {
        pos: 100.0 * Vec3D::new(1.0, 0.0, 1.0),
        intensity: Vec3D::ONES * 60000.0,
    });

    // Prepare scene for rendering
    let bvh = TriangleMesh::build_bvh(&teapot);
    let teapot = Object::new(bvh, teapot_material, None);
    let scene = vec![teapot, ground];
    let lights = vec![light];

    return (scene, lights, camera);
}

#[allow(dead_code)]
fn scene_baby_yoda() -> (Vec<Object>, Vec<Rc<dyn Light>>, Camera) {
    // Baby Yoda (Grogu)
    let grogu = Rc::new(ray_tracer::stl::read_stl("./meshes/baby_yoda.stl").unwrap());
    // Unfortunately the .stl file is not good and does not define a manifold
    // mesh.try_fix_mesh();
    // let grogu = loop_subdivide(&grogu, 1);
    // let grogu_material = Rc::new(Diffuse {
    //     reflectance: Vec3D::new(0.6, 0.6, 0.6),
    // });
    let grogu_material = Rc::new(Dielectric {
        refractive_index: 1.5,
        roughness: None,
    });

    // Ground to cast shadows onto
    let ground = grogu
        .vertices
        .iter()
        .fold(f32::INFINITY, |min, v| min.min(v.z));
    let ground = Object::new(
        Rc::new(Plane::new(Vec3D::Z * ground, Vec3D::Z)),
        Rc::new(Diffuse {
            reflectance: Vec3D::ONES * 0.8,
        }),
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

    let light: Rc<dyn Light> = Rc::new(AreaLight {
        surface: Sphere {
            center: centre + 2000.0 * up + 2000.0 * left - 1000.0 * away,
            radius: 1000.0,
        },
        emmited_radiance: Vec3D::ONES * 5.0,
    });

    // A floating ball
    let ball = Object::new(
        unit_sphere.clone(),
        Rc::new(Diffuse {
            reflectance: Vec3D::new(0.8, 0.2, 0.2) * 2.0,
        }),
        Some(Affine::scale_translate(5.0, Vec3D::new(-4.0, 6.0, 59.0))),
    );

    // Prepare scene for rendering
    let bvh = TriangleMesh::build_bvh(&grogu);
    let grogu = Object::new(bvh, grogu_material, None);
    let scene = vec![grogu, ground, ball];
    let lights = vec![light];

    return (scene, lights, camera);
}

#[allow(dead_code)]
fn scene_primitives() -> (Vec<Object>, Vec<Rc<dyn Light>>, Camera) {
    // World setup
    let mut scene = Vec::new();

    // Small objects
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(-2.0, 0.0, 0.5), 0.5)),
        Rc::new(Diffuse {
            reflectance: Vec3D::new(1.0, 0.01, 0.01),
        }),
        None,
    ));
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(-0.5, 0.0, 0.8), 0.75)),
        Rc::new(Conductor {
            refractive_index: Vec3D::new(1.0, 1.0, 1.0),
            extinction_coefficient: Vec3D::new(1.0, 1.0, 1.0),
            // roughness: Some(RoughnessModel {}),
            roughness: None,
        }),
        None,
    ));
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(1.5, 0.0, 1.0), 1.0)),
        Rc::new(Diffuse {
            reflectance: Vec3D::new(0.01, 0.01, 1.0),
        }),
        None,
    ));
    scene.push(Object::new(
        Rc::new(Triangle::new(
            Vec3D::new(-4.5, 2.0, 0.2),
            Vec3D::new(-2.0, 3.0, 0.2),
            Vec3D::new(-3.0, 2.5, 3.2),
        )),
        Rc::new(Conductor {
            refractive_index: Vec3D::new(1.0, 1.0, 1.0),
            extinction_coefficient: Vec3D::new(1.0, 1.0, 1.0),
            roughness: None,
        }),
        None,
    ));
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(0.0, -1.6, 0.8), 0.75)),
        Rc::new(Dielectric {
            refractive_index: 1.5,
            roughness: None,
        }),
        None,
    ));

    // Globe
    scene.push(Object::new(
        Rc::new(Sphere::new(Vec3D::new(0.0, 0.0, -100.0), 100.0)),
        Rc::new(Diffuse {
            reflectance: Vec3D::ONES * 0.5,
        }),
        None,
    ));

    // Lights
    let mut lights: Vec<Rc<dyn Light>> = Vec::new();
    lights.push(Rc::new(AreaLight {
        surface: Sphere {
            center: Vec3D::new(0.0, 110.0, 0.0),
            radius: 100.0,
        },
        emmited_radiance: Vec3D::ONES * 2.0,
    }));
    // lights.push(PointLight {
    //     pos: Vec3D::new(0.0, 110.0, 0.0),
    //     intensity: Vec3D::ONES * 30000.0,
    // });
    lights.push(Rc::new(AreaLight {
        surface: Sphere {
            center: Vec3D::new(2.5, -2.5, 0.5),
            radius: 0.5,
        },
        emmited_radiance: Vec3D::new(0.1, 3.0, 0.1),
    }));

    // Camera
    let camera = ray_tracer::Camera::new(
        Vec3D::new(-2.5, -4.0, 3.0),
        Vec3D::new(0.5, 1.0, -0.5),
        ASPECT_RATIO,
    );

    return (scene, lights, camera);
}
