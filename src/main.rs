mod kappa;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use clap::Parser;
use rand::distributions::Distribution;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::Normal;
use rayon::prelude::*;

use kappa::{generate_n_kappa_sources, print_kappa_summary, KappaSource, Source};

/// Generate an M x M mock FITS image containing N kappa-sources with kappa <= max_kappa and max radius.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Output FITS file path
    #[arg(short, long, default_value = "mock_kappa_image.fits")]
    pub output: PathBuf,

    /// Number of kappa-sources to generate (N)
    #[arg(short = 'N', long, default_value_t = 50)]
    pub num_kappa: usize,

    /// Image dimension M for an M x M grid (NAXIS1 = NAXIS2 = M)
    #[arg(short = 'M', long, default_value_t = 4096)]
    pub size: usize,

    /// Maximum kappa multiplicity (kappa in [1, max_kappa])
    #[arg(short = 'k', long, default_value_t = 5)]
    pub max_kappa: usize,

    /// Maximum spatial clustering radius in pixels for any kappa-source
    #[arg(short = 'r', long, default_value_t = 25.0)]
    pub max_radius: f32,

    /// Detection threshold for total collective flux in units of noise RMS (e.g. 3.0 for 3xRMS)
    #[arg(short = 's', long, default_value_t = 3.0)]
    pub detection_sigma: f32,

    /// Maximum individual flux allowed for a single subcomponent in kappa >= 2 (in units of RMS)
    #[arg(long, default_value_t = 3.0)]
    pub subcomponent_max_sigma: f32,

    /// Background Gaussian noise standard deviation (RMS / sigma)
    #[arg(long, default_value_t = 1.0)]
    pub noise_sigma: f32,

    /// Background mean
    #[arg(long, default_value_t = 0.0)]
    pub noise_mean: f32,

    /// Full Width at Half Maximum (FWHM) of Gaussian point sources in pixels
    #[arg(long, default_value_t = 10.0)]
    pub fwhm: f32,

    /// Standard deviation in log(flux) for flux variations across subcomponents
    #[arg(long, default_value_t = 0.3)]
    pub flux_sigma: f32,

    /// Interpret flux as peak amplitude instead of total integrated flux
    #[arg(long, default_value_t = false)]
    pub peak_flux: bool,

    /// Maximum allowed source peak amplitude in units of noise sigma (0 to disable)
    #[arg(long, default_value_t = 100.0)]
    pub max_source_sigma: f32,

    /// Also generate a DS9 region file (.reg) for visual inspection
    #[arg(long, default_value_t = true)]
    pub save_regions: bool,

    /// Random seed for reproducibility (optional)
    #[arg(long)]
    pub seed: Option<u64>,
}

/// Formats an 80-character standard FITS header card
fn make_fits_card(key: &str, value: &str, comment: Option<&str>) -> String {
    let card = if key.is_empty() || key == "COMMENT" || key == "HISTORY" {
        if let Some(c) = comment {
            format!("{:<8} {}", key, c)
        } else {
            format!("{:<8} {}", key, value)
        }
    } else if key == "END" {
        "END".to_string()
    } else if key.starts_with("HIERARCH ") {
        if let Some(c) = comment {
            format!("{}= {} / {}", key, value, c)
        } else {
            format!("{}= {}", key, value)
        }
    } else if let Some(c) = comment {
        format!("{:<8}= {:>20} / {}", key, value, c)
    } else {
        format!("{:<8}= {:>20}", key, value)
    };

    if card.len() > 80 {
        card[..80].to_string()
    } else {
        format!("{:<80}", card)
    }
}

/// Formats a FITS string value with single quotes
fn fits_str_val(s: &str) -> String {
    format!("'{:<8}'", s)
}

/// Write a 2D float image with point sources and kappa-Sources metadata into a FITS file
pub fn write_fits_image_with_kappa(
    path: &PathBuf,
    data: &[f32],
    width: usize,
    height: usize,
    sources: &[Source],
    kappa_sources: &[KappaSource],
    noise_mean: f32,
    noise_sigma: f32,
    fwhm: f32,
    peak_flux_mode: bool,
    max_source_sigma: f32,
    max_kappa: usize,
    max_radius: f32,
    detection_sigma: f32,
    sub_max_sigma: f32,
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);

    // =========================================================================
    // 1. Primary HDU (Image + Key Metadata)
    // =========================================================================
    let mut cards = Vec::new();
    cards.push(make_fits_card("SIMPLE", "T", Some("file conforms to FITS standard")));
    cards.push(make_fits_card("BITPIX", "-32", Some("IEEE 32-bit floating point")));
    cards.push(make_fits_card("NAXIS", "2", Some("number of data axes")));
    cards.push(make_fits_card("NAXIS1", &width.to_string(), Some("length of data axis 1 (M)")));
    cards.push(make_fits_card("NAXIS2", &height.to_string(), Some("length of data axis 2 (M)")));
    cards.push(make_fits_card("EXTEND", "T", Some("FITS dataset contains extensions")));

    cards.push(make_fits_card("COMMENT", "", Some("Mock FITS image of N kappa-sources on M x M grid")));
    cards.push(make_fits_card("BG_MEAN", &format!("{:.6}", noise_mean), Some("Background noise mean")));
    cards.push(make_fits_card("BG_SIGMA", &format!("{:.6}", noise_sigma), Some("Background noise RMS (sigma)")));
    cards.push(make_fits_card("NSRC", &sources.len().to_string(), Some("Total constituent point sources")));
    cards.push(make_fits_card("SRC_FWHM", &format!("{:.4}", fwhm), Some("Point source FWHM in pixels")));
    cards.push(make_fits_card("SRC_PEAK", if peak_flux_mode { "T" } else { "F" }, Some("T if flux is peak amplitude")));
    if max_source_sigma > 0.0 {
        cards.push(make_fits_card("SRC_MAXS", &format!("{:.2}", max_source_sigma), Some("Max source peak in noise sigma units")));
    }

    // Kappa-source metadata
    cards.push(make_fits_card("NKAPPA", &kappa_sources.len().to_string(), Some("Total number of kappa-sources (N)")));
    cards.push(make_fits_card("KAP_MAX", &max_kappa.to_string(), Some("Max kappa multiplicity limit")));
    cards.push(make_fits_card("KAP_RAD", &format!("{:.2}", max_radius), Some("Max cluster radius in pixels")));
    cards.push(make_fits_card("DET_SIG", &format!("{:.2}", detection_sigma), Some("Detection sigma (total flux >= DET_SIG*RMS)")));
    cards.push(make_fits_card("SUB_MAXS", &format!("{:.2}", sub_max_sigma), Some("Max subcomponent flux in RMS units")));

    cards.push(make_fits_card("END", "", None));

    // Write primary header records
    let mut header_bytes = Vec::new();
    for card in cards {
        header_bytes.extend_from_slice(card.as_bytes());
    }
    let header_len = header_bytes.len();
    let header_padding = (2880 - (header_len % 2880)) % 2880;
    header_bytes.resize(header_len + header_padding, b' ');
    writer.write_all(&header_bytes)?;

    // Write image data block as big-endian 32-bit floats
    let mut data_bytes = Vec::with_capacity(data.len() * 4 + 2880);
    for &pixel in data {
        data_bytes.extend_from_slice(&pixel.to_be_bytes());
    }
    let data_len = data_bytes.len();
    let data_padding = (2880 - (data_len % 2880)) % 2880;
    data_bytes.resize(data_len + data_padding, 0u8);
    writer.write_all(&data_bytes)?;

    // =========================================================================
    // 2. Binary Table Extension 1: SOURCES (Point Sources Table)
    // =========================================================================
    let num_sources = sources.len();
    let src_row_size = 4 + 4 + 4 + 4 + 4 + 4 + 4 + 4 + 4; // 36 bytes

    let mut ext1_cards = Vec::new();
    ext1_cards.push(make_fits_card("XTENSION", &fits_str_val("BINTABLE"), Some("binary table extension")));
    ext1_cards.push(make_fits_card("BITPIX", "8", Some("8-bit bytes")));
    ext1_cards.push(make_fits_card("NAXIS", "2", Some("2-dimensional table")));
    ext1_cards.push(make_fits_card("NAXIS1", &src_row_size.to_string(), Some("width of table row in bytes")));
    ext1_cards.push(make_fits_card("NAXIS2", &num_sources.to_string(), Some("number of rows in table")));
    ext1_cards.push(make_fits_card("PCOUNT", "0", Some("size of special data area")));
    ext1_cards.push(make_fits_card("GCOUNT", "1", Some("one data group")));
    ext1_cards.push(make_fits_card("TFIELDS", "9", Some("number of fields per row")));
    ext1_cards.push(make_fits_card("EXTNAME", &fits_str_val("SOURCES"), Some("Constituent point sources")));

    ext1_cards.push(make_fits_card("TTYPE1", &fits_str_val("ID"), Some("Point source ID")));
    ext1_cards.push(make_fits_card("TFORM1", &fits_str_val("1J"), Some("32-bit integer")));

    ext1_cards.push(make_fits_card("TTYPE2", &fits_str_val("X"), Some("X position in pixels")));
    ext1_cards.push(make_fits_card("TFORM2", &fits_str_val("1E"), Some("32-bit float")));
    ext1_cards.push(make_fits_card("TUNIT2", &fits_str_val("pixel"), None));

    ext1_cards.push(make_fits_card("TTYPE3", &fits_str_val("Y"), Some("Y position in pixels")));
    ext1_cards.push(make_fits_card("TFORM3", &fits_str_val("1E"), Some("32-bit float")));
    ext1_cards.push(make_fits_card("TUNIT3", &fits_str_val("pixel"), None));

    ext1_cards.push(make_fits_card("TTYPE4", &fits_str_val("FLUX"), Some("Point source flux")));
    ext1_cards.push(make_fits_card("TFORM4", &fits_str_val("1E"), Some("32-bit float")));

    ext1_cards.push(make_fits_card("TTYPE5", &fits_str_val("AMPLITUDE"), Some("Gaussian peak amplitude")));
    ext1_cards.push(make_fits_card("TFORM5", &fits_str_val("1E"), Some("32-bit float")));

    ext1_cards.push(make_fits_card("TTYPE6", &fits_str_val("SIGMA"), Some("Gaussian sigma in pixels")));
    ext1_cards.push(make_fits_card("TFORM6", &fits_str_val("1E"), Some("32-bit float")));
    ext1_cards.push(make_fits_card("TUNIT6", &fits_str_val("pixel"), None));

    ext1_cards.push(make_fits_card("TTYPE7", &fits_str_val("FWHM"), Some("FWHM in pixels")));
    ext1_cards.push(make_fits_card("TFORM7", &fits_str_val("1E"), Some("32-bit float")));
    ext1_cards.push(make_fits_card("TUNIT7", &fits_str_val("pixel"), None));

    ext1_cards.push(make_fits_card("TTYPE8", &fits_str_val("KAPPA_ID"), Some("Parent kappa-Source ID")));
    ext1_cards.push(make_fits_card("TFORM8", &fits_str_val("1J"), Some("32-bit integer")));

    ext1_cards.push(make_fits_card("TTYPE9", &fits_str_val("KAPPA"), Some("Parent kappa multiplicity")));
    ext1_cards.push(make_fits_card("TFORM9", &fits_str_val("1J"), Some("32-bit integer")));

    ext1_cards.push(make_fits_card("END", "", None));

    let mut ext1_header_bytes = Vec::new();
    for card in ext1_cards {
        ext1_header_bytes.extend_from_slice(card.as_bytes());
    }
    let ext1_header_len = ext1_header_bytes.len();
    let ext1_header_padding = (2880 - (ext1_header_len % 2880)) % 2880;
    ext1_header_bytes.resize(ext1_header_len + ext1_header_padding, b' ');
    writer.write_all(&ext1_header_bytes)?;

    let mut src_table_bytes = Vec::with_capacity(num_sources * src_row_size + 2880);
    for src in sources {
        src_table_bytes.extend_from_slice(&(src.id as i32).to_be_bytes());
        src_table_bytes.extend_from_slice(&src.x.to_be_bytes());
        src_table_bytes.extend_from_slice(&src.y.to_be_bytes());
        src_table_bytes.extend_from_slice(&src.flux.to_be_bytes());
        src_table_bytes.extend_from_slice(&src.amplitude.to_be_bytes());
        src_table_bytes.extend_from_slice(&src.sigma.to_be_bytes());
        src_table_bytes.extend_from_slice(&src.fwhm.to_be_bytes());
        src_table_bytes.extend_from_slice(&(src.kappa_id as i32).to_be_bytes());
        src_table_bytes.extend_from_slice(&(src.kappa as i32).to_be_bytes());
    }
    let src_table_len = src_table_bytes.len();
    let src_table_padding = (2880 - (src_table_len % 2880)) % 2880;
    src_table_bytes.resize(src_table_len + src_table_padding, 0u8);
    writer.write_all(&src_table_bytes)?;

    // =========================================================================
    // 3. Binary Table Extension 2: KAPPA_SRCS (kappa-Sources Catalog)
    // =========================================================================
    let num_kappa = kappa_sources.len();
    let kappa_row_size = 4 + 4 + 4 + 4 + 4 + 4 + 4 + 4 + 4; // 36 bytes

    let mut ext2_cards = Vec::new();
    ext2_cards.push(make_fits_card("XTENSION", &fits_str_val("BINTABLE"), Some("binary table extension")));
    ext2_cards.push(make_fits_card("BITPIX", "8", Some("8-bit bytes")));
    ext2_cards.push(make_fits_card("NAXIS", "2", Some("2-dimensional table")));
    ext2_cards.push(make_fits_card("NAXIS1", &kappa_row_size.to_string(), Some("width of table row in bytes")));
    ext2_cards.push(make_fits_card("NAXIS2", &num_kappa.to_string(), Some("number of rows in table")));
    ext2_cards.push(make_fits_card("PCOUNT", "0", Some("size of special data area")));
    ext2_cards.push(make_fits_card("GCOUNT", "1", Some("one data group")));
    ext2_cards.push(make_fits_card("TFIELDS", "9", Some("number of fields per row")));
    ext2_cards.push(make_fits_card("EXTNAME", &fits_str_val("KAPPA_SRCS"), Some("kappa-Sources catalog")));

    ext2_cards.push(make_fits_card("TTYPE1", &fits_str_val("KAPPA_ID"), Some("kappa-Source ID")));
    ext2_cards.push(make_fits_card("TFORM1", &fits_str_val("1J"), Some("32-bit integer")));

    ext2_cards.push(make_fits_card("TTYPE2", &fits_str_val("KAPPA"), Some("kappa multiplicity (1..max_kappa)")));
    ext2_cards.push(make_fits_card("TFORM2", &fits_str_val("1J"), Some("32-bit integer")));

    ext2_cards.push(make_fits_card("TTYPE3", &fits_str_val("CEN_X"), Some("Flux-weighted centroid X (px)")));
    ext2_cards.push(make_fits_card("TFORM3", &fits_str_val("1E"), Some("32-bit float")));
    ext2_cards.push(make_fits_card("TUNIT3", &fits_str_val("pixel"), None));

    ext2_cards.push(make_fits_card("TTYPE4", &fits_str_val("CEN_Y"), Some("Flux-weighted centroid Y (px)")));
    ext2_cards.push(make_fits_card("TFORM4", &fits_str_val("1E"), Some("32-bit float")));
    ext2_cards.push(make_fits_card("TUNIT4", &fits_str_val("pixel"), None));

    ext2_cards.push(make_fits_card("TTYPE5", &fits_str_val("TOTAL_FLUX"), Some("Sum of constituent fluxes")));
    ext2_cards.push(make_fits_card("TFORM5", &fits_str_val("1E"), Some("32-bit float")));

    ext2_cards.push(make_fits_card("TTYPE6", &fits_str_val("MAX_AMP"), Some("Max constituent peak amplitude")));
    ext2_cards.push(make_fits_card("TFORM6", &fits_str_val("1E"), Some("32-bit float")));

    ext2_cards.push(make_fits_card("TTYPE7", &fits_str_val("RADIUS"), Some("Spatial extent radius (px)")));
    ext2_cards.push(make_fits_card("TFORM7", &fits_str_val("1E"), Some("32-bit float")));
    ext2_cards.push(make_fits_card("TUNIT7", &fits_str_val("pixel"), None));

    ext2_cards.push(make_fits_card("TTYPE8", &fits_str_val("SNR"), Some("Total flux / noise_sigma")));
    ext2_cards.push(make_fits_card("TFORM8", &fits_str_val("1E"), Some("32-bit float")));

    ext2_cards.push(make_fits_card("TTYPE9", &fits_str_val("N_MEMBERS"), Some("Number of member point sources")));
    ext2_cards.push(make_fits_card("TFORM9", &fits_str_val("1J"), Some("32-bit integer")));

    ext2_cards.push(make_fits_card("END", "", None));

    let mut ext2_header_bytes = Vec::new();
    for card in ext2_cards {
        ext2_header_bytes.extend_from_slice(card.as_bytes());
    }
    let ext2_header_len = ext2_header_bytes.len();
    let ext2_header_padding = (2880 - (ext2_header_len % 2880)) % 2880;
    ext2_header_bytes.resize(ext2_header_len + ext2_header_padding, b' ');
    writer.write_all(&ext2_header_bytes)?;

    let mut kappa_table_bytes = Vec::with_capacity(num_kappa * kappa_row_size + 2880);
    for ks in kappa_sources {
        kappa_table_bytes.extend_from_slice(&(ks.id as i32).to_be_bytes());
        kappa_table_bytes.extend_from_slice(&(ks.kappa as i32).to_be_bytes());
        kappa_table_bytes.extend_from_slice(&ks.centroid_x.to_be_bytes());
        kappa_table_bytes.extend_from_slice(&ks.centroid_y.to_be_bytes());
        kappa_table_bytes.extend_from_slice(&ks.total_flux.to_be_bytes());
        kappa_table_bytes.extend_from_slice(&ks.max_amplitude.to_be_bytes());
        kappa_table_bytes.extend_from_slice(&ks.radius.to_be_bytes());
        kappa_table_bytes.extend_from_slice(&ks.snr.to_be_bytes());
        kappa_table_bytes.extend_from_slice(&(ks.member_ids.len() as i32).to_be_bytes());
    }
    let kappa_table_len = kappa_table_bytes.len();
    let kappa_table_padding = (2880 - (kappa_table_len % 2880)) % 2880;
    kappa_table_bytes.resize(kappa_table_len + kappa_table_padding, 0u8);
    writer.write_all(&kappa_table_bytes)?;

    writer.flush()?;
    Ok(())
}

/// Write a DS9 region (.reg) file color-coded by kappa-Source multiplicity
pub fn write_ds9_regions_by_kappa(
    path: &PathBuf,
    sources: &[Source],
    kappa_sources: &[KappaSource],
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "# Region file format: DS9 version 4.1")?;
    writeln!(writer, "global font=\"helvetica 10 bold roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1")?;
    writeln!(writer, "image")?;

    // 1. Draw individual point sources with dashed circles
    for src in sources {
        let color = match src.kappa {
            1 => "green",
            2 => "yellow",
            3 => "orange",
            _ => "magenta",
        };
        writeln!(
            writer,
            "circle({:.3},{:.3},{:.1}) # color={} width=1 dash=1 text={{s{}: F={:.3}}}",
            src.x + 1.0,
            src.y + 1.0,
            src.fwhm / 2.0,
            color,
            src.id,
            src.flux
        )?;
    }

    // 2. Draw kappa-Source enclosing bounds and centroids
    for ks in kappa_sources {
        let color = match ks.kappa {
            1 => "green",
            2 => "yellow",
            3 => "orange",
            _ => "red",
        };

        let radius = ks.radius.max(8.0);
        writeln!(
            writer,
            "circle({:.3},{:.3},{:.1}) # color={} width=2 text={{{}-src #{}: F={:.3}, SNR={:.1}}}",
            ks.centroid_x + 1.0,
            ks.centroid_y + 1.0,
            radius,
            color,
            ks.kappa,
            ks.id,
            ks.total_flux,
            ks.snr
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();

    println!("==============================================================================");
    println!("N kappa-Sources FITS Generator (M x M Grid)");
    println!("==============================================================================");
    println!("Image Grid (M x M)       : {} x {} pixels", args.size, args.size);
    println!("Target kappa-Sources     : N = {}", args.num_kappa);
    println!("Multiplicity Range       : 1 <= kappa <= {}", args.max_kappa);
    println!("Max Cluster Radius       : <= {:.1} pixels", args.max_radius);
    println!("Gaussian Noise (RMS)     : mean = {}, sigma = {}", args.noise_mean, args.noise_sigma);
    println!("Point Source FWHM        : {} px", args.fwhm);
    println!("Detection Threshold      : Total Collective Flux >= {:.2} * RMS ({:.4})", 
        args.detection_sigma, args.detection_sigma * args.noise_sigma);
    println!("Subcomponent Limit       : Subcomponent Flux < {:.2} * RMS ({:.4}) for kappa>=2", 
        args.subcomponent_max_sigma, args.subcomponent_max_sigma * args.noise_sigma);
    if args.max_source_sigma > 0.0 {
        println!("Max Source Peak          : <= {} * noise_sigma ({:.4})", 
            args.max_source_sigma, args.max_source_sigma * args.noise_sigma);
    }
    println!("Output File              : {}", args.output.display());
    println!("==============================================================================");

    let num_pixels = args.size * args.size;
    let mut image = vec![0.0f32; num_pixels];

    // 1. Generate Gaussian background noise in parallel
    println!("Generating background Gaussian noise on {}x{} grid...", args.size, args.size);
    let chunk_size = 65536;
    let base_seed = args.seed.unwrap_or_else(|| rand::thread_rng().gen());

    image
        .par_chunks_mut(chunk_size)
        .enumerate()
        .for_each(|(chunk_idx, chunk)| {
            let chunk_seed = base_seed.wrapping_add(chunk_idx as u64 * 10007 + 1);
            let mut rng = StdRng::seed_from_u64(chunk_seed);
            let normal = Normal::new(args.noise_mean, args.noise_sigma).unwrap();

            for val in chunk.iter_mut() {
                *val = normal.sample(&mut rng);
            }
        });

    // 2. Generate N kappa-sources with collective detection >= detection_sigma and subcomponents < sub_max_sigma
    println!("Generating {} kappa-sources (1 <= kappa <= {}, radius <= {:.1} px)...", 
        args.num_kappa, args.max_kappa, args.max_radius);

    let mut source_rng = StdRng::seed_from_u64(base_seed.wrapping_add(999999));
    let (sources, kappa_sources) = generate_n_kappa_sources(
        args.num_kappa,
        args.size,
        args.size,
        args.max_kappa,
        args.max_radius,
        args.fwhm,
        args.detection_sigma,
        args.subcomponent_max_sigma,
        args.flux_sigma,
        args.peak_flux,
        args.max_source_sigma,
        args.noise_sigma,
        &mut source_rng,
    );

    // 3. Render constituent point sources onto the image
    println!("Rendering {} constituent point sources...", sources.len());
    for source in &sources {
        source.render(&mut image, args.size, args.size);
    }

    // 4. Print detailed kappa summary table
    print_kappa_summary(&kappa_sources, args.noise_sigma, args.detection_sigma);

    // 5. Write to multi-extension FITS file
    println!("Writing FITS image and metadata to {}...", args.output.display());
    write_fits_image_with_kappa(
        &args.output,
        &image,
        args.size,
        args.size,
        &sources,
        &kappa_sources,
        args.noise_mean,
        args.noise_sigma,
        args.fwhm,
        args.peak_flux,
        args.max_source_sigma,
        args.max_kappa,
        args.max_radius,
        args.detection_sigma,
        args.subcomponent_max_sigma,
    )?;

    // 6. Write DS9 region file
    if args.save_regions {
        let mut reg_path = args.output.clone();
        reg_path.set_extension("reg");
        println!("Writing DS9 regions overlay to {}...", reg_path.display());
        write_ds9_regions_by_kappa(&reg_path, &sources, &kappa_sources)?;
    }

    println!("Done! Successfully created mock FITS image with {} kappa-sources on {}x{} grid.", 
        kappa_sources.len(), args.size, args.size);
    Ok(())
}
