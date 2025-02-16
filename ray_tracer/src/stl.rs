use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Error, ErrorKind, Read, Result};

use super::triangle_mesh::TriangleMesh;
use super::Vec3D;

pub fn read_stl(file_name: &str) -> Result<TriangleMesh> {
    let (reader, num_triangles) = STLReader::new(file_name)?;
    unify_vertices(reader, num_triangles)
}

struct STLReader {
    remaining_triangles: u32,
    reader: BufReader<File>,
}

impl STLReader {
    fn new(file_name: &str) -> Result<(Self, u32)> {
        let file = File::open(file_name)?;
        let mut reader = BufReader::new(file);

        // Read the 80-byte file header
        let mut header = [0; 80];
        reader.read_exact(&mut header)?;

        // Read the expected number of triangles
        let mut num_triangles = [0; 4];
        reader.read_exact(&mut num_triangles)?;
        let num_triangles = u32::from_le_bytes(num_triangles);

        Ok((
            Self {
                remaining_triangles: num_triangles,
                reader,
            },
            num_triangles,
        ))
    }

    fn read_vector(&mut self) -> Result<Vec3D> {
        let mut bytes = [[0; 4]; 3];

        self.reader.read_exact(&mut bytes[0])?;
        self.reader.read_exact(&mut bytes[1])?;
        self.reader.read_exact(&mut bytes[2])?;

        let floats = bytes.map(|bytes| f32::from_le_bytes(bytes));

        Ok(Vec3D::from(floats))
    }

    fn read_attribute(&mut self) -> Result<()> {
        let mut attribute = [0; 2];
        self.reader.read_exact(&mut attribute)?;
        return Ok(());
    }

    fn read_triangle(&mut self) -> Result<[Vec3D; 3]> {
        // Read normal vector
        let _normal = self.read_vector()?;

        // Read triangle
        let mut triangle = [Vec3D::default(); 3];
        for i in 0..3 {
            triangle[i] = self.read_vector()?;
        }

        // Read the 2-byte attribute data
        let _attribute = self.read_attribute()?;

        return Ok(triangle);
    }

    fn eof(&mut self) -> Result<bool> {
        // Doublecheck if EOF has indeed been reached
        if self.reader.read(&mut [0])? == 0 {
            return Ok(true);
        } else {
            return Ok(false);
        }
    }
}

impl Iterator for STLReader {
    type Item = Result<[Vec3D; 3]>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining_triangles == 0 {
            match self.eof() {
                Ok(true) => return None,
                Ok(false) => {
                    return Some(Err(Error::new(
                        ErrorKind::InvalidData,
                        "Unexpected data after last triangle",
                    )))
                }
                Err(err) => return Some(Err(err)),
            }
        }

        // Update triangle count
        self.remaining_triangles -= 1;
        // Read and return triangle
        return Some(self.read_triangle());
    }
}

#[inline]
fn to_bits(point: Vec3D) -> (u32, u32, u32) {
    (point.x.to_bits(), point.y.to_bits(), point.z.to_bits())
}

fn unify_vertices(
    triangle_stream: impl Iterator<Item = Result<[Vec3D; 3]>>,
    num_triangles: u32,
) -> Result<TriangleMesh> {
    // Loop through all triangles using a HashMap to unify vertices
    let num_triangles = num_triangles as usize;
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
