use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use clap::Parser;

#[path = "../fits.rs"]
mod fits;
#[path = "../kappa.rs"]
mod kappa;

use fits::read_fits_image;
use kappa::{KappaSource, Source};

/// Extract kappa-sources from a 2D FITS astronomical image
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct ExtractCli {
    /// Input FITS image file path (e.g. test.fits)
    #[arg(required = true)]
    pub input: PathBuf,

    /// Output extracted catalog FITS file (default: <input>.extracted.fits)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Maximum kappa multiplicity limit to extract (0 for no upper limit)
    #[arg(short = 'k', long, default_value_t = 0)]
    pub max_kappa: usize,

    /// Detection threshold for total collective flux in units of beam flux RMS (e.g. 3.0 for 3xRMS)
    #[arg(short = 's', long, default_value_t = 3.0)]
    pub detection_sigma: f32,

    /// Maximum spatial clustering radius in pixels from centroid
    #[arg(short = 'r', long, default_value_t = 25.0)]
    pub cluster_radius: f32,

    /// Estimated PSF FWHM in pixels
    #[arg(long, default_value_t = 10.0)]
    pub fwhm: f32,

    /// Peak detection threshold in SNR units (in matched-filtered map) to identify candidate subcomponents
    #[arg(long, default_value_t = 2.0)]
    pub peak_snr: f32,

    /// Maximum individual flux for a single subcomponent in kappa >= 2 (in units of beam RMS)
    #[arg(long, default_value_t = 3.0)]
    pub subcomponent_max_sigma: f32,

    /// Also generate a DS9 region overlay file (.reg)
    #[arg(long, default_value_t = true)]
    pub save_regions: bool,
}

/// Robust background estimation (Median and Median Absolute Deviation -> RMS)
fn estimate_background_and_rms(data: &[f32]) -> (f32, f32) {
    let mut sample: Vec<f32> = if data.len() > 500_000 {
        let step = (data.len() / 500_000).max(1);
        data.iter().step_by(step).copied().collect()
    } else {
        data.to_vec()
    };

    sample.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = sample[sample.len() / 2];

    let mut abs_diffs: Vec<f32> = sample.iter().map(|&x| (x - median).abs()).collect();
    abs_diffs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad = abs_diffs[abs_diffs.len() / 2];

    let rms = (1.4826 * mad).max(1e-6);
    (median, rms)
}

/// Apply fast separable 2D Gaussian matched filtering
fn gaussian_filter_2d(
    data: &[f32],
    width: usize,
    height: usize,
    sigma: f32,
    bg_median: f32,
) -> (Vec<f32>, f32) {
    let radius = (3.0 * sigma).ceil() as isize;
    let mut kernel = Vec::new();
    let two_sigma_sq = 2.0 * sigma * sigma;

    for i in -radius..=radius {
        kernel.push((-((i * i) as f32) / two_sigma_sq).exp());
    }
    let k_sum: f32 = kernel.iter().sum();
    for k in &mut kernel {
        *k /= k_sum;
    }

    // Kernel sum of squares for 2D kernel
    let k_sum_sq_1d: f32 = kernel.iter().map(|&v| v * v).sum();
    let k_sum_sq_2d = k_sum_sq_1d * k_sum_sq_1d;

    // Horizontal pass
    let mut temp = vec![0.0f32; width * height];
    for y in 0..height {
        let row_offset = y * width;
        for x in 0..width {
            let mut sum = 0.0f32;
            for (ki, &kval) in kernel.iter().enumerate() {
                let dx = ki as isize - radius;
                let nx = (x as isize + dx).clamp(0, width as isize - 1) as usize;
                sum += (data[row_offset + nx] - bg_median) * kval;
            }
            temp[row_offset + x] = sum;
        }
    }

    // Vertical pass
    let mut filtered = vec![0.0f32; width * height];
    for y in 0..height {
        let row_offset = y * width;
        for x in 0..width {
            let mut sum = 0.0f32;
            for (ki, &kval) in kernel.iter().enumerate() {
                let dy = ki as isize - radius;
                let ny = (y as isize + dy).clamp(0, height as isize - 1) as usize;
                sum += temp[ny * width + x] * kval;
            }
            filtered[row_offset + x] = sum;
        }
    }

    (filtered, k_sum_sq_2d)
}

/// Detects local peak subcomponents in the 2D filtered image with optimal matched filter photometry
fn detect_subcomponents(
    raw_data: &[f32],
    width: usize,
    height: usize,
    bg_median: f32,
    fwhm: f32,
    peak_snr: f32,
) -> (Vec<Source>, f32, f32) {
    let sigma_psf = fwhm / (2.0 * (2.0 * 2.0f32.ln()).sqrt());

    let (filtered, k_sum_sq_2d) = gaussian_filter_2d(raw_data, width, height, sigma_psf, bg_median);
    let (_, filtered_rms) = estimate_background_and_rms(&filtered);

    // Beam flux RMS: noise standard deviation on the integrated flux
    let flux_conv_factor = 1.0 / k_sum_sq_2d;
    let beam_flux_rms = filtered_rms * flux_conv_factor;

    let min_filtered_peak = peak_snr * filtered_rms;
    let min_peak_sep = (fwhm * 0.4).ceil().max(2.0) as isize;

    let mut peak_locs = Vec::new();

    for y in min_peak_sep as usize..(height - min_peak_sep as usize) {
        let row_offset = y * width;
        for x in min_peak_sep as usize..(width - min_peak_sep as usize) {
            let val = filtered[row_offset + x];
            if val < min_filtered_peak {
                continue;
            }

            let mut is_local_max = true;
            for dy in -min_peak_sep..=min_peak_sep {
                for dx in -min_peak_sep..=min_peak_sep {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let neighbor_val = filtered[(y as isize + dy) as usize * width + (x as isize + dx) as usize];
                    if neighbor_val > val {
                        is_local_max = false;
                        break;
                    }
                }
                if !is_local_max {
                    break;
                }
            }

            if is_local_max {
                peak_locs.push((x, y, val));
            }
        }
    }

    let mut sources = Vec::with_capacity(peak_locs.len());
    let mut sid = 1;

    for (px, py, peak_val) in peak_locs {
        // 3x3 Sub-pixel parabolic centroid refinement
        let v_c = peak_val;
        let v_l = filtered[py * width + (px - 1)];
        let v_r = filtered[py * width + (px + 1)];
        let v_u = filtered[(py - 1) * width + px];
        let v_d = filtered[(py + 1) * width + px];

        let dx = if (2.0 * v_c - v_l - v_r).abs() > 1e-5 {
            0.5 * (v_l - v_r) / (v_l - 2.0 * v_c + v_r)
        } else {
            0.0
        };
        let dy = if (2.0 * v_c - v_u - v_d).abs() > 1e-5 {
            0.5 * (v_u - v_d) / (v_u - 2.0 * v_c + v_d)
        } else {
            0.0
        };

        let sub_x = (px as f32 + dx.clamp(-0.8, 0.8)).clamp(0.0, width as f32 - 1.0);
        let sub_y = (py as f32 + dy.clamp(-0.8, 0.8)).clamp(0.0, height as f32 - 1.0);

        // Optimal unbiased flux estimate from matched filter peak
        let estimated_total_flux = peak_val * flux_conv_factor;
        let peak_amplitude = estimated_total_flux / (2.0 * std::f32::consts::PI * sigma_psf * sigma_psf);

        sources.push(Source {
            id: sid,
            x: sub_x,
            y: sub_y,
            flux: estimated_total_flux,
            amplitude: peak_amplitude,
            sigma: sigma_psf,
            fwhm,
            kappa_id: 0,
            kappa: 0,
        });
        sid += 1;
    }

    (sources, filtered_rms, beam_flux_rms)
}

/// Extract kappa-sources hierarchically with max_kappa and radius constraints
fn extract_kappa_hierarchy(
    sources: &mut [Source],
    max_radius: f32,
    max_kappa: usize,
    detection_sigma: f32,
    sub_max_sigma: f32,
    beam_flux_rms: f32,
) -> Vec<KappaSource> {
    let n = sources.len();
    if n == 0 {
        return Vec::new();
    }

    let min_detection_flux = detection_sigma * beam_flux_rms;
    let max_sub_flux = sub_max_sigma * beam_flux_rms;
    let dist_sq_thresh = (max_radius * 1.5) * (max_radius * 1.5);

    // Proximity graph
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = sources[i].x - sources[j].x;
            let dy = sources[i].y - sources[j].y;
            if dx * dx + dy * dy < dist_sq_thresh {
                adj[i].push(j);
                adj[j].push(i);
            }
        }
    }

    // Connected components
    let mut visited = vec![false; n];
    let mut raw_clusters = Vec::new();

    for start in 0..n {
        if !visited[start] {
            let mut cluster = Vec::new();
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(start);
            visited[start] = true;

            while let Some(u) = queue.pop_front() {
                cluster.push(u);
                for &v in &adj[u] {
                    if !visited[v] {
                        visited[v] = true;
                        queue.push_back(v);
                    }
                }
            }
            raw_clusters.push(cluster);
        }
    }

    let mut kappa_sources = Vec::new();

    for cluster in raw_clusters {
        let kappa = cluster.len();

        // Enforce max_kappa limit if specified (max_kappa > 0)
        if max_kappa > 0 && kappa > max_kappa {
            continue;
        }

        let mut total_flux = 0.0f32;
        let mut max_amp = 0.0f32;
        let mut weighted_x = 0.0f32;
        let mut weighted_y = 0.0f32;
        let mut member_ids = Vec::with_capacity(kappa);
        let mut all_sub_subthreshold = true;

        for &idx in &cluster {
            let s = &sources[idx];
            total_flux += s.flux;
            weighted_x += s.flux * s.x;
            weighted_y += s.flux * s.y;
            if s.amplitude > max_amp {
                max_amp = s.amplitude;
            }
            if kappa >= 2 && s.flux >= max_sub_flux {
                all_sub_subthreshold = false;
            }
            member_ids.push(s.id);
        }

        let (cen_x, cen_y) = if total_flux > 0.0 {
            (weighted_x / total_flux, weighted_y / total_flux)
        } else {
            let s = &sources[cluster[0]];
            (s.x, s.y)
        };

        let mut max_r_sq = 0.0f32;
        let fwhm = sources[cluster[0]].fwhm;
        for &idx in &cluster {
            let dx = sources[idx].x - cen_x;
            let dy = sources[idx].y - cen_y;
            let r_sq = dx * dx + dy * dy;
            if r_sq > max_r_sq {
                max_r_sq = r_sq;
            }
        }
        let radius = max_r_sq.sqrt() + fwhm / 2.0;

        if total_flux >= min_detection_flux && (kappa == 1 || (radius <= max_radius * 1.5 && all_sub_subthreshold)) {
            let snr = total_flux / beam_flux_rms;

            kappa_sources.push(KappaSource {
                id: 0,
                kappa,
                member_ids,
                centroid_x: cen_x,
                centroid_y: cen_y,
                total_flux,
                max_amplitude: max_amp,
                radius,
                snr,
            });
        }
    }

    // Hierarchical sort: 1-sources, 2-sources, 3-sources...
    kappa_sources.sort_by(|a, b| {
        a.kappa
            .cmp(&b.kappa)
            .then_with(|| b.total_flux.partial_cmp(&a.total_flux).unwrap_or(std::cmp::Ordering::Equal))
    });

    // Re-index IDs
    let mut id_to_source_idx = std::collections::HashMap::new();
    for (i, s) in sources.iter().enumerate() {
        id_to_source_idx.insert(s.id, i);
    }

    for (k_idx, ks) in kappa_sources.iter_mut().enumerate() {
        let kid = k_idx + 1;
        ks.id = kid;
        for &mid in &ks.member_ids {
            if let Some(&s_idx) = id_to_source_idx.get(&mid) {
                sources[s_idx].kappa_id = kid;
                sources[s_idx].kappa = ks.kappa;
            }
        }
    }

    kappa_sources
}

/// Print formatted breakdown summary table of extracted kappa-sources
fn print_extraction_report(
    kappa_sources: &[KappaSource],
    sources: &[Source],
    bg_median: f32,
    bg_rms: f32,
    beam_flux_rms: f32,
    detection_sigma: f32,
    max_kappa: usize,
) {
    let mut counts_by_kappa: std::collections::HashMap<usize, (usize, f32, f32, f32)> = std::collections::HashMap::new();
    let mut max_k = 0;

    for ks in kappa_sources {
        let entry = counts_by_kappa.entry(ks.kappa).or_insert((0, 0.0, 0.0, 0.0));
        entry.0 += 1;
        entry.1 += ks.total_flux;
        entry.2 += ks.radius;
        if ks.snr > entry.3 {
            entry.3 = ks.snr;
        }
        if ks.kappa > max_k {
            max_k = ks.kappa;
        }
    }

    println!("--------------------------------------------------------------------------------");
    println!("Background Level : Median = {:.4}, Pixel RMS = {:.4}, Beam Flux RMS = {:.4}", 
        bg_median, bg_rms, beam_flux_rms);
    println!("Candidate Peaks  : {} subcomponents detected", sources.len());
    println!("Detection Cutoff : Total Flux >= {:.2} * Beam RMS ({:.4})", 
        detection_sigma, detection_sigma * beam_flux_rms);
    if max_kappa > 0 {
        println!("Multiplicity Cut : kappa <= {}", max_kappa);
    }
    println!("--------------------------------------------------------------------------------");
    println!("{:<11} {:<8} {:<14} {:<12} {:<13} {:<8}", "kappa", "Count", "Sum Flux", "Mean Flux", "Mean Radius", "Max SNR");
    println!("--------------------------------------------------------------------------------");

    for k in 1..=max_k {
        if let Some(&(count, sum_flux, sum_radius, max_snr)) = counts_by_kappa.get(&k) {
            let mean_flux = sum_flux / count as f32;
            let mean_rad = sum_radius / count as f32;
            println!(
                "{:<11} {:<8} {:<14.4} {:<12.4} {:<13.1} {:<8.2}",
                format!("{}-sources", k),
                count,
                sum_flux,
                mean_flux,
                mean_rad,
                max_snr
            );
        }
    }
    println!("--------------------------------------------------------------------------------");
    println!("Total kappa-Sources Extracted: {}", kappa_sources.len());
    println!("--------------------------------------------------------------------------------");
}

/// Write DS9 region file for visual validation
fn write_ds9_regions(path: &PathBuf, kappa_sources: &[KappaSource], sources: &[Source]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "# Region file format: DS9 version 4.1")?;
    writeln!(writer, "global font=\"helvetica 10 bold roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1")?;
    writeln!(writer, "image")?;

    for src in sources {
        if src.kappa > 0 {
            writeln!(
                writer,
                "circle({:.3},{:.3},{:.1}) # color=white width=1 dash=1 text={{peak #{}: F={:.2}}}",
                src.x + 1.0,
                src.y + 1.0,
                src.fwhm / 2.0,
                src.id,
                src.flux
            )?;
        }
    }

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

fn make_fits_card(key: &str, value: &str, comment: Option<&str>) -> String {
    let card = if key.is_empty() || key == "COMMENT" || key == "HISTORY" {
        if let Some(c) = comment {
            format!("{:<8} {}", key, c)
        } else {
            format!("{:<8} {}", key, value)
        }
    } else if key == "END" {
        "END".to_string()
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

fn fits_str_val(s: &str) -> String {
    format!("'{:<8}'", s)
}

/// Write extracted kappa-sources catalog to FITS binary table
fn write_extracted_fits_catalog(
    path: &PathBuf,
    kappa_sources: &[KappaSource],
    _sources: &[Source],
    bg_median: f32,
    bg_rms: f32,
    beam_flux_rms: f32,
    detection_sigma: f32,
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);

    // Primary Header (Null Image)
    let mut p_cards = Vec::new();
    p_cards.push(make_fits_card("SIMPLE", "T", Some("file conforms to FITS standard")));
    p_cards.push(make_fits_card("BITPIX", "8", Some("Null array")));
    p_cards.push(make_fits_card("NAXIS", "0", Some("No image data in primary HDU")));
    p_cards.push(make_fits_card("EXTEND", "T", Some("Extensions present")));
    p_cards.push(make_fits_card("BG_MED", &format!("{:.6}", bg_median), Some("Background median")));
    p_cards.push(make_fits_card("BG_RMS", &format!("{:.6}", bg_rms), Some("Background pixel RMS")));
    p_cards.push(make_fits_card("FLUX_RMS", &format!("{:.6}", beam_flux_rms), Some("Beam flux RMS")));
    p_cards.push(make_fits_card("DET_SIG", &format!("{:.2}", detection_sigma), Some("Detection sigma")));
    p_cards.push(make_fits_card("NKAPPA", &kappa_sources.len().to_string(), Some("Total extracted kappa-sources")));
    p_cards.push(make_fits_card("END", "", None));

    let mut p_bytes = Vec::new();
    for c in p_cards {
        p_bytes.extend_from_slice(c.as_bytes());
    }
    let p_len = p_bytes.len();
    let p_pad = (2880 - (p_len % 2880)) % 2880;
    p_bytes.resize(p_len + p_pad, b' ');
    writer.write_all(&p_bytes)?;

    // Extension: KAPPA_SRCS
    let num_kappa = kappa_sources.len();
    let row_size = 4 + 4 + 4 + 4 + 4 + 4 + 4 + 4 + 4; // 36 bytes

    let mut ext_cards = Vec::new();
    ext_cards.push(make_fits_card("XTENSION", &fits_str_val("BINTABLE"), Some("binary table extension")));
    ext_cards.push(make_fits_card("BITPIX", "8", Some("8-bit bytes")));
    ext_cards.push(make_fits_card("NAXIS", "2", Some("2-dimensional table")));
    ext_cards.push(make_fits_card("NAXIS1", &row_size.to_string(), Some("width of table row in bytes")));
    ext_cards.push(make_fits_card("NAXIS2", &num_kappa.to_string(), Some("number of rows in table")));
    ext_cards.push(make_fits_card("PCOUNT", "0", Some("size of special data area")));
    ext_cards.push(make_fits_card("GCOUNT", "1", Some("one data group")));
    ext_cards.push(make_fits_card("TFIELDS", "9", Some("number of fields per row")));
    ext_cards.push(make_fits_card("EXTNAME", &fits_str_val("KAPPA_SRCS"), Some("Extracted kappa-sources")));

    ext_cards.push(make_fits_card("TTYPE1", &fits_str_val("KAPPA_ID"), Some("kappa-Source ID")));
    ext_cards.push(make_fits_card("TFORM1", &fits_str_val("1J"), Some("32-bit integer")));

    ext_cards.push(make_fits_card("TTYPE2", &fits_str_val("KAPPA"), Some("kappa multiplicity")));
    ext_cards.push(make_fits_card("TFORM2", &fits_str_val("1J"), Some("32-bit integer")));

    ext_cards.push(make_fits_card("TTYPE3", &fits_str_val("CEN_X"), Some("Flux-weighted centroid X (px)")));
    ext_cards.push(make_fits_card("TFORM3", &fits_str_val("1E"), Some("32-bit float")));
    ext_cards.push(make_fits_card("TUNIT3", &fits_str_val("pixel"), None));

    ext_cards.push(make_fits_card("TTYPE4", &fits_str_val("CEN_Y"), Some("Flux-weighted centroid Y (px)")));
    ext_cards.push(make_fits_card("TFORM4", &fits_str_val("1E"), Some("32-bit float")));
    ext_cards.push(make_fits_card("TUNIT4", &fits_str_val("pixel"), None));

    ext_cards.push(make_fits_card("TTYPE5", &fits_str_val("TOTAL_FLUX"), Some("Total collective flux")));
    ext_cards.push(make_fits_card("TFORM5", &fits_str_val("1E"), Some("32-bit float")));

    ext_cards.push(make_fits_card("TTYPE6", &fits_str_val("MAX_AMP"), Some("Max member peak amplitude")));
    ext_cards.push(make_fits_card("TFORM6", &fits_str_val("1E"), Some("32-bit float")));

    ext_cards.push(make_fits_card("TTYPE7", &fits_str_val("RADIUS"), Some("Spatial extent radius (px)")));
    ext_cards.push(make_fits_card("TFORM7", &fits_str_val("1E"), Some("32-bit float")));
    ext_cards.push(make_fits_card("TUNIT7", &fits_str_val("pixel"), None));

    ext_cards.push(make_fits_card("TTYPE8", &fits_str_val("SNR"), Some("Total flux / Beam RMS")));
    ext_cards.push(make_fits_card("TFORM8", &fits_str_val("1E"), Some("32-bit float")));

    ext_cards.push(make_fits_card("TTYPE9", &fits_str_val("N_MEMBERS"), Some("Number of member subcomponents")));
    ext_cards.push(make_fits_card("TFORM9", &fits_str_val("1J"), Some("32-bit integer")));

    ext_cards.push(make_fits_card("END", "", None));

    let mut ext_bytes = Vec::new();
    for c in ext_cards {
        ext_bytes.extend_from_slice(c.as_bytes());
    }
    let ext_len = ext_bytes.len();
    let ext_pad = (2880 - (ext_len % 2880)) % 2880;
    ext_bytes.resize(ext_len + ext_pad, b' ');
    writer.write_all(&ext_bytes)?;

    let mut table_bytes = Vec::with_capacity(num_kappa * row_size + 2880);
    for ks in kappa_sources {
        table_bytes.extend_from_slice(&(ks.id as i32).to_be_bytes());
        table_bytes.extend_from_slice(&(ks.kappa as i32).to_be_bytes());
        table_bytes.extend_from_slice(&ks.centroid_x.to_be_bytes());
        table_bytes.extend_from_slice(&ks.centroid_y.to_be_bytes());
        table_bytes.extend_from_slice(&ks.total_flux.to_be_bytes());
        table_bytes.extend_from_slice(&ks.max_amplitude.to_be_bytes());
        table_bytes.extend_from_slice(&ks.radius.to_be_bytes());
        table_bytes.extend_from_slice(&ks.snr.to_be_bytes());
        table_bytes.extend_from_slice(&(ks.member_ids.len() as i32).to_be_bytes());
    }
    let t_len = table_bytes.len();
    let t_pad = (2880 - (t_len % 2880)) % 2880;
    table_bytes.resize(t_len + t_pad, 0u8);
    writer.write_all(&table_bytes)?;

    writer.flush()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = ExtractCli::parse();

    println!("==============================================================================");
    println!("kappa-Source Extractor from FITS Image (kappa_extract.bin)");
    println!("==============================================================================");
    println!("Input FITS File      : {}", args.input.display());
    println!("Detection Threshold  : Total Flux >= {:.2} * Beam RMS", args.detection_sigma);
    println!("Clustering Radius    : <= {:.1} pixels", args.cluster_radius);
    if args.max_kappa > 0 {
        println!("Max Multiplicity     : kappa <= {}", args.max_kappa);
    } else {
        println!("Max Multiplicity     : Unlimited");
    }
    println!("Subcomponent Limit   : Flux < {:.2} * Beam RMS for kappa >= 2", args.subcomponent_max_sigma);
    println!("Estimated PSF FWHM   : {:.1} pixels", args.fwhm);
    println!("==============================================================================");

    // 1. Read input FITS image
    println!("Reading FITS image from {}...", args.input.display());
    let fits_img = read_fits_image(&args.input)?;
    println!("Loaded image grid: {} x {} pixels", fits_img.width, fits_img.height);

    // 2. Measure background & noise RMS
    println!("Estimating background and RMS noise level...");
    let (bg_median, bg_rms) = estimate_background_and_rms(&fits_img.data);

    // 3. Detect candidate subcomponents (local peaks) with optimal matched filter photometry
    println!("Detecting candidate point-source subcomponents with matched filtering...");
    let (mut sources, _filtered_rms, beam_flux_rms) = detect_subcomponents(
        &fits_img.data,
        fits_img.width,
        fits_img.height,
        bg_median,
        args.fwhm,
        args.peak_snr,
    );

    // 4. Extract hierarchical kappa-sources
    println!("Extracting kappa-Sources (1-sources, 2-sources, 3-sources...)...");
    let kappa_sources = extract_kappa_hierarchy(
        &mut sources,
        args.cluster_radius,
        args.max_kappa,
        args.detection_sigma,
        args.subcomponent_max_sigma,
        beam_flux_rms,
    );

    // 5. Display extraction summary report
    print_extraction_report(
        &kappa_sources,
        &sources,
        bg_median,
        bg_rms,
        beam_flux_rms,
        args.detection_sigma,
        args.max_kappa,
    );

    // 6. Save extracted catalog to FITS
    let out_fits_path = args.output.unwrap_or_else(|| {
        let mut p = args.input.clone();
        p.set_extension("extracted.fits");
        p
    });
    println!("Writing extracted catalog to {}...", out_fits_path.display());
    write_extracted_fits_catalog(&out_fits_path, &kappa_sources, &sources, bg_median, bg_rms, beam_flux_rms, args.detection_sigma)?;

    // 7. Save DS9 region overlay
    if args.save_regions {
        let mut reg_path = out_fits_path.clone();
        reg_path.set_extension("reg");
        println!("Writing DS9 region overlay to {}...", reg_path.display());
        write_ds9_regions(&reg_path, &kappa_sources, &sources)?;
    }

    println!("Extraction complete! Successfully processed {}.", args.input.display());
    Ok(())
}
