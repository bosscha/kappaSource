#[path = "psf.rs"]
pub mod psf;

use std::collections::HashMap;
use self::psf::Psf;
use rand::distributions::Distribution;
use rand::rngs::StdRng;
use rand::Rng;
use rand_distr::LogNormal;

/// Gaussian/PSF-convolved point source properties
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct Source {
    pub id: usize,
    pub x: f32,
    pub y: f32,
    pub flux: f32,
    pub amplitude: f32,
    pub sigma: f32,
    pub fwhm: f32,
    pub kappa_id: usize,
    pub kappa: usize,
}

#[allow(dead_code)]
impl Source {
    /// Render source contribution into the image buffer by convolving with the telescope PSF
    pub fn render_convolved(&self, psf: &Psf, image: &mut [f32], width: usize, height: usize) {
        psf.convolve_point_source(self.x, self.y, self.flux, image, width, height);
    }
}

/// A kappa-Source formed by the union of kappa close point sources
/// satisfying: sum(flux) >= detection_sigma * noise_sigma
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct KappaSource {
    /// Unique identifier for this kappa-source
    pub id: usize,
    /// Multiplicity kappa (number of point sources in the union)
    pub kappa: usize,
    /// Member point source IDs (1-indexed)
    pub member_ids: Vec<usize>,
    /// Flux-weighted centroid X in pixels
    pub centroid_x: f32,
    /// Flux-weighted centroid Y in pixels
    pub centroid_y: f32,
    /// Total integrated flux (sum of member fluxes)
    pub total_flux: f32,
    /// Maximum peak amplitude among member sources after PSF convolution
    pub max_amplitude: f32,
    /// Spatial extent / characteristic bounding radius in pixels
    pub radius: f32,
    /// Signal-to-noise ratio: total_flux / flux_rms
    pub snr: f32,
}

/// Generates N kappa-sources on an M x M image grid
///
/// Physical constraints:
/// - kappa multiplicity in [1, max_kappa]
/// - All constituent point sources within max_radius from centroid
/// - Total collective flux >= detection_sigma * flux_rms (default: 3.0 * RMS)
/// - For kappa >= 2: Every individual subcomponent flux < subcomponent_max_sigma * flux_rms
/// - For kappa = 1: Single source flux >= detection_sigma * flux_rms
/// - Each subcomponent is convolved by the telescope PSF
#[allow(dead_code)]
pub fn generate_n_kappa_sources(
    num_kappa_sources: usize,
    width: usize,
    height: usize,
    max_kappa: usize,
    max_radius: f32,
    psf: &Psf,
    detection_sigma: f32,
    subcomponent_max_sigma: f32,
    flux_sigma: f32,
    peak_flux_mode: bool,
    max_source_sigma: f32,
    noise_sigma: f32,
    rng: &mut StdRng,
) -> (Vec<Source>, Vec<KappaSource>) {
    let beam_flux_rms = noise_sigma * (4.0 * std::f32::consts::PI * psf.sigma * psf.sigma).sqrt();
    let min_detection_flux = if peak_flux_mode {
        detection_sigma * noise_sigma
    } else {
        detection_sigma * beam_flux_rms
    };
    let sub_max_flux_limit = if peak_flux_mode {
        subcomponent_max_sigma * noise_sigma
    } else {
        subcomponent_max_sigma * beam_flux_rms
    };

    let psf_peak_factor = psf.peak_response();

    let max_allowed_amplitude = if max_source_sigma > 0.0 {
        Some(max_source_sigma * noise_sigma)
    } else {
        None
    };

    let log_normal = if flux_sigma > 0.0 {
        Some(LogNormal::new(0.0, flux_sigma).unwrap())
    } else {
        None
    };

    let margin = (max_radius + psf.fwhm * 2.0).max(50.0);
    let effective_max_k = max_kappa.max(1);

    let mut all_sources: Vec<Source> = Vec::new();
    let mut kappa_sources: Vec<KappaSource> = Vec::with_capacity(num_kappa_sources);
    let mut point_source_id = 1;

    // Minimum separation between distinct kappa-sources to avoid accidental overlaps
    let min_separation_sq = (max_radius * 2.2) * (max_radius * 2.2);
    let mut placed_centers: Vec<(f32, f32)> = Vec::with_capacity(num_kappa_sources);

    for k_idx in 0..num_kappa_sources {
        // Multiplicity kappa in 1..=max_kappa
        let kappa = if effective_max_k == 1 {
            1
        } else {
            rng.gen_range(1..=effective_max_k)
        };

        // Find a suitable centroid position in the M x M grid
        let mut center_x = 0.0f32;
        let mut center_y = 0.0f32;
        let mut placement_attempts = 0;

        while placement_attempts < 500 {
            let cx = rng.gen_range(margin..(width as f32 - margin));
            let cy = rng.gen_range(margin..(height as f32 - margin));

            let too_close = placed_centers.iter().any(|&(px, py)| {
                let dx = cx - px;
                let dy = cy - py;
                dx * dx + dy * dy < min_separation_sq
            });

            if !too_close || placement_attempts > 400 {
                center_x = cx;
                center_y = cy;
                break;
            }
            placement_attempts += 1;
        }
        placed_centers.push((center_x, center_y));

        // Determine total target flux for this kappa-source (>= detection_sigma * beam_flux_rms)
        let total_target_flux = min_detection_flux * (1.0 + rng.gen_range(0.05..0.5));

        // Generate constituent fluxes
        let fluxes: Vec<f32> = if kappa == 1 {
            vec![total_target_flux]
        } else {
            let max_allowed_sub_flux = sub_max_flux_limit.min(total_target_flux * 0.95);

            let mut weights: Vec<f32> = (0..kappa)
                .map(|_| {
                    if let Some(dist) = &log_normal {
                        dist.sample(rng)
                    } else {
                        rng.gen_range(0.5..1.5)
                    }
                })
                .collect();

            let weight_sum: f32 = weights.iter().sum();
            for w in &mut weights {
                *w /= weight_sum;
            }

            let mut sub_fluxes: Vec<f32> = weights.iter().map(|&w| w * total_target_flux).collect();

            // Enforce that every subcomponent is strictly below sub_max_flux_limit
            for _ in 0..20 {
                let mut excess = 0.0f32;
                let mut valid_count = 0;
                for f in &mut sub_fluxes {
                    if *f > max_allowed_sub_flux {
                        excess += *f - max_allowed_sub_flux;
                        *f = max_allowed_sub_flux;
                    } else {
                        valid_count += 1;
                    }
                }
                if excess > 0.0 && valid_count > 0 {
                    let distribute = excess / valid_count as f32;
                    for f in &mut sub_fluxes {
                        if *f < max_allowed_sub_flux {
                            *f += distribute;
                        }
                    }
                } else {
                    break;
                }
            }
            sub_fluxes
        };

        // Build member point sources to be convolved by PSF
        let mut member_sources: Vec<Source> = Vec::with_capacity(kappa);
        let mut member_ids = Vec::with_capacity(kappa);
        let mut total_flux = 0.0f32;
        let mut max_amp = 0.0f32;

        for (m_idx, &flux) in fluxes.iter().enumerate() {
            let (sx, sy) = if kappa == 1 || m_idx == 0 {
                (center_x, center_y)
            } else {
                let r = max_radius * rng.gen::<f32>().sqrt();
                let theta = rng.gen_range(0.0..std::f32::consts::TAU);
                (center_x + r * theta.cos(), center_y + r * theta.sin())
            };

            let mut amplitude = if peak_flux_mode {
                flux
            } else {
                flux * psf_peak_factor
            };

            if let Some(max_amp_limit) = max_allowed_amplitude {
                if amplitude > max_amp_limit {
                    amplitude = max_amp_limit;
                }
            }

            total_flux += flux;
            if amplitude > max_amp {
                max_amp = amplitude;
            }

            let sid = point_source_id;
            point_source_id += 1;
            member_ids.push(sid);

            member_sources.push(Source {
                id: sid,
                x: sx,
                y: sy,
                flux,
                amplitude,
                sigma: psf.sigma,
                fwhm: psf.fwhm,
                kappa_id: k_idx + 1,
                kappa,
            });
        }

        // Calculate flux-weighted centroid
        let mut weighted_x = 0.0f32;
        let mut weighted_y = 0.0f32;
        for src in &member_sources {
            weighted_x += src.flux * src.x;
            weighted_y += src.flux * src.y;
        }
        let (centroid_x, centroid_y) = (weighted_x / total_flux, weighted_y / total_flux);

        // Calculate bounding radius
        let mut max_r_sq = 0.0f32;
        for src in &member_sources {
            let dx = src.x - centroid_x;
            let dy = src.y - centroid_y;
            let r_sq = dx * dx + dy * dy;
            if r_sq > max_r_sq {
                max_r_sq = r_sq;
            }
        }
        let radius = max_r_sq.sqrt() + psf.fwhm / 2.0;
        let snr = total_flux / beam_flux_rms;

        all_sources.extend(member_sources);

        kappa_sources.push(KappaSource {
            id: k_idx + 1,
            kappa,
            member_ids,
            centroid_x,
            centroid_y,
            total_flux,
            max_amplitude: max_amp,
            radius,
            snr,
        });
    }

    // Sort kappa-sources hierarchically: first 1-sources, then 2-sources, 3-sources, ...
    kappa_sources.sort_by(|a, b| {
        a.kappa
            .cmp(&b.kappa)
            .then_with(|| b.total_flux.partial_cmp(&a.total_flux).unwrap_or(std::cmp::Ordering::Equal))
    });

    // Re-assign sorted kappa_ids and update member sources
    let mut id_to_source_idx: HashMap<usize, usize> = HashMap::new();
    for (i, s) in all_sources.iter().enumerate() {
        id_to_source_idx.insert(s.id, i);
    }

    for (k_idx, k_src) in kappa_sources.iter_mut().enumerate() {
        let new_id = k_idx + 1;
        k_src.id = new_id;

        for &m_id in &k_src.member_ids {
            if let Some(&s_idx) = id_to_source_idx.get(&m_id) {
                all_sources[s_idx].kappa_id = new_id;
                all_sources[s_idx].kappa = k_src.kappa;
            }
        }
    }

    (all_sources, kappa_sources)
}

/// Print formatted breakdown summary table of extracted kappa-sources
#[allow(dead_code)]
pub fn print_kappa_summary(kappa_sources: &[KappaSource], noise_sigma: f32, detection_sigma: f32) {
    let mut counts_by_kappa: HashMap<usize, (usize, f32, f32, f32)> = HashMap::new();

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
    println!("Hierarchical kappa-Source Summary (kappa = 1, 2, ... n)");
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
    println!("Total kappa-Sources: {} (Detection Condition: Total Flux >= {:.2} * RMS (Noise σ={:.4}))", 
        kappa_sources.len(), detection_sigma, noise_sigma);
    println!("--------------------------------------------------------------------------------");
}
