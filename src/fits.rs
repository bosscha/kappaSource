use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// Simple FITS image representation
#[allow(dead_code)]
pub struct FitsImage {
    pub width: usize,
    pub height: usize,
    pub data: Vec<f32>,
    pub bitpix: i32,
    pub bzero: f32,
    pub bscale: f32,
}

/// Reads a 2D float image from a standard FITS file
pub fn read_fits_image<P: AsRef<Path>>(path: P) -> std::io::Result<FitsImage> {
    let mut file = BufReader::new(File::open(path)?);

    let mut width = 0usize;
    let mut height = 0usize;
    let mut bitpix = -32i32;
    let mut bscale = 1.0f32;
    let mut bzero = 0.0f32;

    let mut header_bytes_read = 0usize;
    let mut header_done = false;
    let mut card = [0u8; 80];

    while !header_done {
        file.read_exact(&mut card)?;
        header_bytes_read += 80;

        let card_str = String::from_utf8_lossy(&card);
        let key = card_str[..8].trim();

        if key == "END" {
            header_done = true;
        } else if key == "NAXIS1" {
            if let Some(val) = parse_fits_header_val(&card_str) {
                width = val.parse::<usize>().unwrap_or(0);
            }
        } else if key == "NAXIS2" {
            if let Some(val) = parse_fits_header_val(&card_str) {
                height = val.parse::<usize>().unwrap_or(0);
            }
        } else if key == "BITPIX" {
            if let Some(val) = parse_fits_header_val(&card_str) {
                bitpix = val.parse::<i32>().unwrap_or(-32);
            }
        } else if key == "BSCALE" {
            if let Some(val) = parse_fits_header_val(&card_str) {
                bscale = val.parse::<f32>().unwrap_or(1.0);
            }
        } else if key == "BZERO" {
            if let Some(val) = parse_fits_header_val(&card_str) {
                bzero = val.parse::<f32>().unwrap_or(0.0);
            }
        }
    }

    // Seek to the beginning of the data block (multiple of 2880 bytes)
    let header_padding = (2880 - (header_bytes_read % 2880)) % 2880;
    file.seek(SeekFrom::Current(header_padding as i64))?;

    let num_pixels = width * height;
    let mut data = Vec::with_capacity(num_pixels);

    match bitpix {
        -32 => {
            let mut buf = [0u8; 4];
            for _ in 0..num_pixels {
                file.read_exact(&mut buf)?;
                let raw_val = f32::from_be_bytes(buf);
                data.push(raw_val * bscale + bzero);
            }
        }
        -64 => {
            let mut buf = [0u8; 8];
            for _ in 0..num_pixels {
                file.read_exact(&mut buf)?;
                let raw_val = f64::from_be_bytes(buf) as f32;
                data.push(raw_val * bscale + bzero);
            }
        }
        16 => {
            let mut buf = [0u8; 2];
            for _ in 0..num_pixels {
                file.read_exact(&mut buf)?;
                let raw_val = i16::from_be_bytes(buf) as f32;
                data.push(raw_val * bscale + bzero);
            }
        }
        32 => {
            let mut buf = [0u8; 4];
            for _ in 0..num_pixels {
                file.read_exact(&mut buf)?;
                let raw_val = i32::from_be_bytes(buf) as f32;
                data.push(raw_val * bscale + bzero);
            }
        }
        _ => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unsupported BITPIX: {}", bitpix),
            ));
        }
    }

    Ok(FitsImage {
        width,
        height,
        data,
        bitpix,
        bzero,
        bscale,
    })
}

/// Parses value portion of an 80-column FITS header record
fn parse_fits_header_val(card: &str) -> Option<String> {
    if let Some(eq_idx) = card.find('=') {
        let after_eq = &card[eq_idx + 1..];
        let val_part = if let Some(slash_idx) = after_eq.find('/') {
            &after_eq[..slash_idx]
        } else {
            after_eq
        };
        Some(val_part.trim().replace('\'', ""))
    } else {
        None
    }
}
