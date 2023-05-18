use std::io::Read;
use byteorder::{LittleEndian, ReadBytesExt};

use super::Vec3D;
use super::primitives::Triangle;

fn read_vector(reader: &mut impl std::io::Read) -> std::io::Result<Vec3D> {
    Ok(Vec3D::new(
        reader.read_f32::<LittleEndian>()?,
        reader.read_f32::<LittleEndian>()?,
        reader.read_f32::<LittleEndian>()?,
    ))
}

pub fn read_stl(file_name: &str) -> std::io::Result<Vec<Triangle>> {
    let file = std::fs::File::open(file_name)?;
    let mut reader = std::io::BufReader::new(file);
    
    // HEADER
    let mut header = [0; 80];
    reader.read_exact(&mut header)?;
    let _header = std::str::from_utf8(&header);

    // number of triangles
    let num_triangles = reader.read_u32::<LittleEndian>()? as usize;

    // triangles
    let mut result = Vec::with_capacity(num_triangles);
    for _ in 0..num_triangles {
        let _normal = read_vector(&mut reader)?;
        let triangle = Triangle::new(
            read_vector(&mut reader)?,
            read_vector(&mut reader)?,
            read_vector(&mut reader)?,
        );
        let _attribute = reader.read_u16::<LittleEndian>()?;

        // maybe include check that calculated normal match declared normal?

        result.push(triangle);
    }
    // check that EOF has indeed been reached
    if reader.read(&mut [0])? != 0 {
        panic!("File length does not match declared number of triangles");
    }

    return Ok(result)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read() {
        read_stl("./src/meshes/Tetrahedron.stl").unwrap();
    }
}