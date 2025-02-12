use byteorder::{LittleEndian, ReadBytesExt};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Error, ErrorKind, Read, Result};

use super::triangle_mesh::TriangleMesh;
use super::Vec3D;

pub fn read_stl(file_name: &str) -> Result<TriangleMesh> {
    let (reader, num_triangles) = STLReader::new(file_name)?;
    unify_vertices(reader, num_triangles)
}

fn read_vector(reader: &mut impl Read) -> Result<Vec3D> {
    Ok(Vec3D::new(
        reader.read_f32::<LittleEndian>()?,
        reader.read_f32::<LittleEndian>()?,
        reader.read_f32::<LittleEndian>()?,
    ))
}

struct STLReader {
    remaining_triangles: usize,
    reader: BufReader<File>,
}

impl STLReader {
    fn new(file_name: &str) -> Result<(Self, usize)> {
        let file = File::open(file_name)?;
        let mut reader = BufReader::new(file);

        // Read the 80-byte file header
        let mut header = [0; 80];
        reader.read_exact(&mut header)?;
        let _header = std::str::from_utf8(&header);

        // Read the expected number of triangles
        let num_triangles = reader.read_u32::<LittleEndian>()? as usize;

        Ok((
            Self {
                remaining_triangles: num_triangles,
                reader,
            },
            num_triangles,
        ))
    }
}

impl Iterator for STLReader {
    type Item = Result<[Vec3D; 3]>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining_triangles == 0 {
            // Check that EOF has indeed been reached
            match self.reader.read(&mut [0]) {
                Err(err) => return Some(Err(err)),
                Ok(val) => {
                    if val != 0 {
                        return Some(Err(Error::new(
                            ErrorKind::InvalidData,
                            "Unexpected data after last triangle",
                        )));
                    } else {
                        return None;
                    }
                }
            }
        }

        // Read normal vector
        let normal = read_vector(&mut self.reader);
        match normal {
            Ok(_) => {}
            Err(err) => return Some(Err(err)),
        }

        // Read triangle
        let mut triangle = [Vec3D::default(); 3];
        for i in 0..3 {
            match read_vector(&mut self.reader) {
                Ok(point) => triangle[i] = point,
                Err(err) => return Some(Err(err)),
            }
        }

        self.remaining_triangles -= 1;

        // Read the 2-byte attribute data
        let attribute = self.reader.read_u16::<LittleEndian>();
        match attribute {
            Ok(_) => {}
            Err(err) => return Some(Err(err)),
        }

        return Some(Ok(triangle));
    }
}

fn to_bits(point: Vec3D) -> (u32, u32, u32) {
    (point.x.to_bits(), point.y.to_bits(), point.z.to_bits())
}

pub fn unify_vertices(
    triangle_stream: impl Iterator<Item = Result<[Vec3D; 3]>>,
    num_triangles: usize,
) -> Result<TriangleMesh> {
    // Loop through all triangles using a HashMap to unify vertices
    let mut triangles: Vec<[usize; 3]> = Vec::with_capacity(num_triangles);
    let mut vertices: Vec<Vec3D> = Vec::with_capacity(num_triangles / 2);
    let mut vertex_map: HashMap<(u32, u32, u32), usize> = HashMap::with_capacity(num_triangles / 2);

    for points in triangle_stream {
        let mut triangle = [0; 3];
        for (i, point) in points?.into_iter().enumerate() {
            let bit_repr = to_bits(point);
            if let Some(index) = vertex_map.get(&bit_repr) {
                triangle[i] = *index;
            } else {
                vertex_map.insert(bit_repr, vertices.len());
                triangle[i] = vertices.len();
                vertices.push(point);
            }
        }
        triangles.push(triangle);
    }

    return Ok(TriangleMesh {
        vertices,
        triangles,
    });
}
