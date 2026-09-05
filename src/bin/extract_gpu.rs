use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime};
use clap::Parser;
use pollster::FutureExt;
use wgpu::util::DeviceExt;

#[path = "../fits.rs"]
mod fits;
#[path = "../kappa.rs"]
mod kappa;

use fits::read_fits_image;
use kappa::{KappaSource, Source};

/// Extract kappa-sources from a 2D FITS astronomical image using GPU-accelerated two-pass multi-scale search radius
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct ExtractGpuCli {
    /// Input FITS image file path (e.g. test.fits)
    #[arg(required = true)]
    pub input: PathBuf,

    /// Output extracted catalog FITS file (default: <fitsname>_<timestamp>.extracted.fits)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Maximum kappa multiplicity limit to extract (0 for no upper limit)
    #[arg(short = 'k', long, default_value_t = 0)]
    pub max_kappa: usize,

    /// Detection threshold for total collective flux in units of beam flux RMS (e.g. 3.0 for 3xRMS)
    #[arg(short = 's', long, default_value_t = 3.0)]
    pub detection_sigma: f32,

    /// Search radius in pixels for grouping multiple subcomponents (R_search)
    #[arg(short = 'r', long = "search-radius", alias = "cluster-radius", default_value_t = 25.0)]
    pub search_radius: f32,

    /// Estimated PSF FWHM in pixels (set 0 or use --psf-auto for automatic empirical measurement)
    #[arg(long, default_value_t = 10.0)]
    pub fwhm: f32,

    /// Automatically estimate the PSF FWHM from bright isolated point sources in the image
    #[arg(long = "psf-auto", alias = "auto-psf", default_value_t = false)]
    pub psf_auto: bool,

    /// Target minimum SNR for point sources used in auto-PSF estimation (e.g. 20.0 for SNR > 20)
    #[arg(long = "min-psf-snr", default_value_t = 20.0)]
    pub min_psf_snr: f32,

    /// Maximum number of brightest sources to sample for auto-PSF calibration (default: 20)
    #[arg(long = "psf-samples", default_value_t = 20)]
    pub psf_samples: usize,

    /// Minimum SNR for individual candidate subcomponent peaks inside a search radius
    #[arg(long, default_value_t = 1.2)]
    pub min_sub_snr: f32,

    /// Candidate seed detection threshold in SNR units (alias: --peak-snr)
    #[arg(long = "seed-snr", alias = "peak-snr", default_value_t = 2.2)]
    pub seed_snr: f32,

    /// Maximum individual flux for a single subcomponent in kappa >= 2 (in units of beam RMS)
    #[arg(long, default_value_t = 3.0)]
    pub subcomponent_max_sigma: f32,

    /// Also generate a DS9 region overlay file (.reg)
    #[arg(long, default_value_t = true)]
    pub save_regions: bool,

    /// Also generate ASCII / CSV text catalogs (.cat and .csv)
    #[arg(long, default_value_t = true)]
    pub save_ascii: bool,
}

/// Generate UTC date-time timestamp string in format YYYYMMDD_HHMMSS
fn get_current_timestamp() -> String {
    let now = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    let secs_per_day = 86400;
    let days = now / secs_per_day;
    let rem_secs = now % secs_per_day;
    let hours = rem_secs / 3600;
    let mins = (rem_secs % 3600) / 60;
    let secs = rem_secs % 60;

    let mut year = 1970;
    let mut d = days;
    loop {
        let leap = if (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0) { 1 } else { 0 };
        let days_in_year = 365 + leap;
        if d < days_in_year {
            let month_days = [31, 28 + leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
            let mut month = 1;
            for &md in &month_days {
                if d < md {
                    break;
                }
                d -= md;
                month += 1;
            }
            let day = d + 1;
            return format!("{:04}{:02}{:02}_{:02}{:02}{:02}", year, month, day, hours, mins, secs);
        }
        d -= days_in_year;
        year += 1;
    }
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

/// Automatically estimate PSF FWHM by fitting 2D intensity moments to bright isolated point sources (SNR > 20 or top 20 brightest)
fn estimate_psf_fwhm_auto(
    data: &[f32],
    width: usize,
    height: usize,
    bg_median: f32,
    bg_rms: f32,
    min_psf_snr: f32,
    max_samples: usize,
) -> Option<(f32, usize, &'static str)> {
    let base_thresh = bg_median + 2.5 * bg_rms;
    let min_sep: isize = 12;
    let stamp_radius: isize = 8;

    // 1. Gather all prominent local maxima
    let mut candidate_peaks = Vec::new();
    for y in (min_sep as usize)..(height - min_sep as usize) {
        let row_offset = y * width;
        for x in (min_sep as usize)..(width - min_sep as usize) {
            let val = data[row_offset + x];
            if val < base_thresh {
                continue;
            }

            let mut is_max = true;
            for dy in -min_sep..=min_sep {
                for dx in -min_sep..=min_sep {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let ny = (y as isize + dy) as usize;
                    let nx = (x as isize + dx) as usize;
                    if data[ny * width + nx] > val {
                        is_max = false;
                        break;
                    }
                }
                if !is_max {
                    break;
                }
            }

            if is_max {
                let snr = (val - bg_median) / bg_rms;
                candidate_peaks.push((x, y, snr));
            }
        }
    }

    if candidate_peaks.is_empty() {
        return None;
    }

    // 2. Select peaks: either SNR >= min_psf_snr (e.g. > 20) or top 20 brightest
    candidate_peaks.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

    let high_snr_peaks: Vec<_> = candidate_peaks.iter().copied().filter(|p| p.2 >= min_psf_snr).collect();
    let (selected_peaks, mode_str) = if !high_snr_peaks.is_empty() {
        let limit = high_snr_peaks.len().min(max_samples);
        (&high_snr_peaks[..limit], "SNR > 20.0 (high-significance stars)")
    } else {
        let limit = candidate_peaks.len().min(max_samples);
        (&candidate_peaks[..limit], "top 20 brightest sources")
    };

    let mut fwhm_samples = Vec::new();

    for &(x, y, _) in selected_peaks {
        let mut sum_i = 0.0f32;
        let mut sum_ix = 0.0f32;
        let mut sum_iy = 0.0f32;

        for dy in -stamp_radius..=stamp_radius {
            for dx in -stamp_radius..=stamp_radius {
                let ny = (y as isize + dy) as usize;
                let nx = (x as isize + dx) as usize;
                let i_val = (data[ny * width + nx] - bg_median).max(0.0);
                sum_i += i_val;
                sum_ix += i_val * dx as f32;
                sum_iy += i_val * dy as f32;
            }
        }

        if sum_i <= 0.0 {
            continue;
        }

        let cx = sum_ix / sum_i;
        let cy = sum_iy / sum_i;

        let mut mxx = 0.0f32;
        let mut myy = 0.0f32;
        let mut mxy = 0.0f32;

        for dy in -stamp_radius..=stamp_radius {
            for dx in -stamp_radius..=stamp_radius {
                let ny = (y as isize + dy) as usize;
                let nx = (x as isize + dx) as usize;
                let i_val = (data[ny * width + nx] - bg_median).max(0.0);
                let diff_x = dx as f32 - cx;
                let diff_y = dy as f32 - cy;
                mxx += i_val * diff_x * diff_x;
                myy += i_val * diff_y * diff_y;
                mxy += i_val * diff_x * diff_y;
            }
        }

        mxx /= sum_i;
        myy /= sum_i;
        mxy /= sum_i;

        let det = mxx * myy - mxy * mxy;
        if det <= 0.0 {
            continue;
        }

        let sigma = det.powf(0.25);
        let fwhm = 2.35482 * sigma;
        let ellipticity = (mxx - myy).abs() / (mxx + myy + 1e-6);

        if fwhm >= 3.0 && fwhm <= 30.0 && ellipticity <= 0.35 {
            fwhm_samples.push(fwhm);
        }
    }

    if fwhm_samples.is_empty() {
        None
    } else {
        fwhm_samples.sort_by(|a, b| a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal));
        let median_fwhm = fwhm_samples[fwhm_samples.len() / 2];
        Some((median_fwhm, fwhm_samples.len(), mode_str))
    }
}

// =============================================================================
// GPU ACCELERATION PIPELINE (wgpu / Vulkan compute on AMD Radeon 8060S)
// =============================================================================

#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
struct GpuConvUniforms {
    width: u32,
    height: u32,
    radius: i32,
    kernel_len: u32,
    bg_median: f32,
    pad1: f32,
    pad2: f32,
    pad3: f32,
}

const SHADER_SOURCE: &str = r#"
struct Params {
    width: u32,
    height: u32,
    radius: i32,
    kernel_len: u32,
    bg_median: f32,
    pad1: f32,
    pad2: f32,
    pad3: f32,
};

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> kernel_weights: array<f32>;
@group(0) @binding(2) var<storage, read> in_data: array<f32>;
@group(0) @binding(3) var<storage, read_write> out_data: array<f32>;

@compute @workgroup_size(16, 16)
fn conv_h(@builtin(global_invocation_id) id: vec3<u32>) {
    let x = id.x;
    let y = id.y;
    if (x >= params.width || y >= params.height) {
        return;
    }

    let w = i32(params.width);
    let r = params.radius;
    let row_offset = y * params.width;

    var sum = 0.0;
    for (var k = 0u; k < params.kernel_len; k = k + 1u) {
        let dx = i32(k) - r;
        let nx = clamp(i32(x) + dx, 0, w - 1);
        let val = in_data[row_offset + u32(nx)] - params.bg_median;
        sum = sum + val * kernel_weights[k];
    }

    out_data[row_offset + x] = sum;
}

@compute @workgroup_size(16, 16)
fn conv_v(@builtin(global_invocation_id) id: vec3<u32>) {
    let x = id.x;
    let y = id.y;
    if (x >= params.width || y >= params.height) {
        return;
    }

    let h = i32(params.height);
    let r = params.radius;
    let w = params.width;

    var sum = 0.0;
    for (var k = 0u; k < params.kernel_len; k = k + 1u) {
        let dy = i32(k) - r;
        let ny = clamp(i32(y) + dy, 0, h - 1);
        let val = in_data[u32(ny) * w + x];
        sum = sum + val * kernel_weights[k];
    }

    out_data[y * w + x] = sum;
}
"#;

struct GpuContext {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline_h: wgpu::ComputePipeline,
    pipeline_v: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    pub adapter_name: String,
    pub backend_name: String,
}

impl GpuContext {
    fn new() -> Result<Self, Box<dyn std::error::Error>> {
        let instance = wgpu::Instance::default();
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
                ..Default::default()
            })
            .block_on()
            .map_err(|e| format!("Failed to find GPU adapter: {:?}", e))?;

        let info = adapter.get_info();
        let adapter_name = info.name;
        let backend_name = format!("{:?}", info.backend);

        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: Some("KappaGpuDevice"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::default(),
                memory_hints: wgpu::MemoryHints::Performance,
                ..Default::default()
            })
            .block_on()?;

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("KappaSeparableConvShader"),
            source: wgpu::ShaderSource::Wgsl(SHADER_SOURCE.into()),
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("KappaConvBGL"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("KappaPipelineLayout"),
            bind_group_layouts: &[Some(&bind_group_layout)],
            immediate_size: 0,
        });

        let pipeline_h = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("KappaConvHPipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some("conv_h"),
            compilation_options: Default::default(),
            cache: None,
        });

        let pipeline_v = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("KappaConvVPipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some("conv_v"),
            compilation_options: Default::default(),
            cache: None,
        });

        Ok(Self {
            device,
            queue,
            pipeline_h,
            pipeline_v,
            bind_group_layout,
            adapter_name,
            backend_name,
        })
    }

    /// Perform GPU-accelerated 2D separable Gaussian convolution
    fn convolve_2d_separable(
        &self,
        raw_data: &[f32],
        width: usize,
        height: usize,
        sigma: f32,
        bg_median: f32,
    ) -> Result<(Vec<f32>, f32, f64), Box<dyn std::error::Error>> {
        let t_start = Instant::now();
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

        let k_sum_sq_1d: f32 = kernel.iter().map(|&v| v * v).sum();
        let k_sum_sq_2d = k_sum_sq_1d * k_sum_sq_1d;

        let num_pixels = width * height;
        let buffer_size = (num_pixels * std::mem::size_of::<f32>()) as wgpu::BufferAddress;

        let uniforms = GpuConvUniforms {
            width: width as u32,
            height: height as u32,
            radius: radius as i32,
            kernel_len: kernel.len() as u32,
            bg_median,
            pad1: 0.0,
            pad2: 0.0,
            pad3: 0.0,
        };

        let uniform_buf = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("UniformBuf"),
            contents: bytemuck::bytes_of(&uniforms),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        let kernel_buf = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("KernelBuf"),
            contents: bytemuck::cast_slice(&kernel),
            usage: wgpu::BufferUsages::STORAGE,
        });

        let in_buf = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("InputBuf"),
            contents: bytemuck::cast_slice(raw_data),
            usage: wgpu::BufferUsages::STORAGE,
        });

        let temp_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("TempBuf"),
            size: buffer_size,
            usage: wgpu::BufferUsages::STORAGE,
            mapped_at_creation: false,
        });

        let out_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("OutBuf"),
            size: buffer_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let readback_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("ReadbackBuf"),
            size: buffer_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // Pass 1: Horizontal (in_buf -> temp_buf)
        let bind_group_h = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("BG_H"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: uniform_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: kernel_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: in_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: temp_buf.as_entire_binding() },
            ],
        });

        // Pass 2: Vertical (temp_buf -> out_buf)
        let uniforms_v = GpuConvUniforms {
            bg_median: 0.0, // already subtracted in pass 1
            ..uniforms
        };
        let uniform_buf_v = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("UniformBufV"),
            contents: bytemuck::bytes_of(&uniforms_v),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        let bind_group_v = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("BG_V"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: uniform_buf_v.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: kernel_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: temp_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: out_buf.as_entire_binding() },
            ],
        });

        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("KappaConvEncoder"),
        });

        let workgroups_x = ((width as u32) + 15) / 16;
        let workgroups_y = ((height as u32) + 15) / 16;

        {
            let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("PassHorizontal"),
                timestamp_writes: None,
            });
            cpass.set_pipeline(&self.pipeline_h);
            cpass.set_bind_group(0, &bind_group_h, &[]);
            cpass.dispatch_workgroups(workgroups_x, workgroups_y, 1);
        }

        {
            let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("PassVertical"),
                timestamp_writes: None,
            });
            cpass.set_pipeline(&self.pipeline_v);
            cpass.set_bind_group(0, &bind_group_v, &[]);
            cpass.dispatch_workgroups(workgroups_x, workgroups_y, 1);
        }

        encoder.copy_buffer_to_buffer(&out_buf, 0, &readback_buf, 0, buffer_size);
        self.queue.submit(Some(encoder.finish()));

        // Map readback buffer
        let slice = readback_buf.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        self.device.poll(wgpu::PollType::wait_indefinitely())?;
        rx.recv()??;

        let mapped = slice.get_mapped_range()?;
        let result: Vec<f32> = bytemuck::cast_slice(&mapped).to_vec();
        drop(mapped);
        readback_buf.unmap();

        let elapsed_ms = t_start.elapsed().as_secs_f64() * 1000.0;
        Ok((result, k_sum_sq_2d, elapsed_ms))
    }
}

/// Candidate cluster seed
#[derive(Debug, Clone, Copy)]
struct Seed {
    x: usize,
    y: usize,
    snr: f32,
}

/// Two-Pass Hierarchical Search Radius Extraction Engine using GPU Compute
fn extract_kappa_two_pass_gpu(
    gpu: &GpuContext,
    raw_data: &[f32],
    width: usize,
    height: usize,
    bg_median: f32,
    _bg_rms: f32,
    fwhm: f32,
    search_radius: f32,
    min_sub_snr: f32,
    seed_snr_thresh: f32,
    detection_sigma: f32,
    subcomponent_max_sigma: f32,
    max_kappa: usize,
) -> Result<(Vec<KappaSource>, Vec<Source>, f32), Box<dyn std::error::Error>> {
    let sigma_psf = fwhm / (2.0 * (2.0 * 2.0f32.ln()).sqrt());

    // 1. GPU Point-source matched filter
    let (point_filt, k_sum_sq_2d, ms_point) =
        gpu.convolve_2d_separable(raw_data, width, height, sigma_psf, bg_median)?;
    let (_, point_filt_rms) = estimate_background_and_rms(&point_filt);
    let flux_conv_factor = 1.0 / k_sum_sq_2d;
    let beam_flux_rms = point_filt_rms * flux_conv_factor;

    // 2. GPU Cluster-scale smoothed filter (matched to search radius extent)
    let sigma_cluster = (sigma_psf * sigma_psf + (search_radius / 2.0) * (search_radius / 2.0)).sqrt();
    let (cluster_filt, _, ms_cluster) =
        gpu.convolve_2d_separable(raw_data, width, height, sigma_cluster, bg_median)?;
    let (_, cluster_filt_rms) = estimate_background_and_rms(&cluster_filt);

    println!(
        "GPU Compute Time : Point Filter = {:.2} ms, Cluster Filter = {:.2} ms (Total: {:.2} ms)",
        ms_point, ms_cluster, ms_point + ms_cluster
    );

    // 3. Detect candidate subcomponent peaks down to min_sub_snr
    let min_sub_val = min_sub_snr * point_filt_rms;
    let min_peak_sep = (fwhm * 0.45).ceil().max(2.0) as isize;

    let mut sub_peaks: Vec<Source> = Vec::new();
    let mut sid = 1;

    for y in min_peak_sep as usize..(height - min_peak_sep as usize) {
        let row_offset = y * width;
        for x in min_peak_sep as usize..(width - min_peak_sep as usize) {
            let val = point_filt[row_offset + x];
            if val < min_sub_val {
                continue;
            }

            let mut is_local_max = true;
            for dy in -min_peak_sep..=min_peak_sep {
                for dx in -min_peak_sep..=min_peak_sep {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    if point_filt[(y as isize + dy) as usize * width + (x as isize + dx) as usize] > val {
                        is_local_max = false;
                        break;
                    }
                }
                if !is_local_max {
                    break;
                }
            }

            if is_local_max {
                let v_c = val;
                let v_l = point_filt[py_idx(y, 0) * width + (x - 1)];
                let v_r = point_filt[py_idx(y, 0) * width + (x + 1)];
                let v_u = point_filt[py_idx(y, -1) * width + x];
                let v_d = point_filt[py_idx(y, 1) * width + x];

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

                let sub_x = (x as f32 + dx.clamp(-0.8, 0.8)).clamp(0.0, width as f32 - 1.0);
                let sub_y = (y as f32 + dy.clamp(-0.8, 0.8)).clamp(0.0, height as f32 - 1.0);
                let flux = val * flux_conv_factor;
                let amp = flux / (2.0 * std::f32::consts::PI * sigma_psf * sigma_psf);

                sub_peaks.push(Source {
                    id: sid,
                    x: sub_x,
                    y: sub_y,
                    flux,
                    amplitude: amp,
                    sigma: sigma_psf,
                    fwhm,
                    kappa_id: 0,
                    kappa: 0,
                });
                sid += 1;
            }
        }
    }

    let mut used_sub = vec![false; sub_peaks.len()];
    let mut kappa_sources: Vec<KappaSource> = Vec::new();
    let min_detection_flux = detection_sigma * beam_flux_rms;
    let max_sub_flux = subcomponent_max_sigma * beam_flux_rms;

    // =========================================================================
    // PASS 1: Extract all genuine 1-sources (subcomponents with Flux >= detection_sigma * beam_rms)
    // =========================================================================
    for i in 0..sub_peaks.len() {
        if !used_sub[i] && sub_peaks[i].flux >= min_detection_flux {
            let s = &sub_peaks[i];
            let snr = s.flux / beam_flux_rms;

            used_sub[i] = true;
            let sx = s.x;
            let sy = s.y;

            // Lock all nearby sub-threshold peaks within 1 beam FWHM
            let beam_lock_sq = (fwhm * 0.8) * (fwhm * 0.8);
            for j in 0..sub_peaks.len() {
                if !used_sub[j] {
                    let dx = sub_peaks[j].x - sx;
                    let dy = sub_peaks[j].y - sy;
                    if dx * dx + dy * dy <= beam_lock_sq {
                        used_sub[j] = true;
                    }
                }
            }

            kappa_sources.push(KappaSource {
                id: 0,
                kappa: 1,
                member_ids: vec![s.id],
                centroid_x: s.x,
                centroid_y: s.y,
                total_flux: s.flux,
                max_amplitude: s.amplitude,
                radius: fwhm / 2.0,
                snr,
            });
        }
    }

    // =========================================================================
    // PASS 2: Extract multi-component kappa-sources (kappa >= 2) from sub-threshold peaks
    // =========================================================================
    let seed_sep = (search_radius * 0.6).round().max(10.0) as isize;
    let min_c_seed = seed_snr_thresh * cluster_filt_rms;

    let mut seeds = Vec::new();
    for y in seed_sep as usize..(height - seed_sep as usize) {
        let row_offset = y * width;
        for x in seed_sep as usize..(width - seed_sep as usize) {
            let c_val = cluster_filt[row_offset + x];
            let c_snr = c_val / cluster_filt_rms;

            if c_val >= min_c_seed && check_local_max(&cluster_filt, width, height, x, y, (search_radius * 0.6) as isize) {
                seeds.push(Seed { x, y, snr: c_snr });
            }
        }
    }

    seeds.sort_by(|a, b| b.snr.partial_cmp(&a.snr).unwrap_or(std::cmp::Ordering::Equal));
    let search_rad_sq = search_radius * search_radius;

    for seed in seeds {
        let sx = seed.x as f32;
        let sy = seed.y as f32;

        let mut member_indices = Vec::new();
        for (i, s) in sub_peaks.iter().enumerate() {
            if !used_sub[i] && s.flux < max_sub_flux {
                let dx = s.x - sx;
                let dy = s.y - sy;
                if dx * dx + dy * dy <= search_rad_sq {
                    member_indices.push(i);
                }
            }
        }

        if member_indices.len() < 2 {
            continue;
        }

        let mut total_flux = 0.0f32;
        let mut weighted_x = 0.0f32;
        let mut weighted_y = 0.0f32;
        let mut max_amp = 0.0f32;

        for &idx in &member_indices {
            let s = &sub_peaks[idx];
            total_flux += s.flux;
            weighted_x += s.flux * s.x;
            weighted_y += s.flux * s.y;
            if s.amplitude > max_amp {
                max_amp = s.amplitude;
            }
        }

        let cen_x = weighted_x / total_flux;
        let cen_y = weighted_y / total_flux;

        let mut verified_indices = Vec::new();
        let mut verified_total_flux = 0.0f32;
        let mut max_r_sq = 0.0f32;

        for &idx in &member_indices {
            let s = &sub_peaks[idx];
            let dx = s.x - cen_x;
            let dy = s.y - cen_y;
            let dist_sq = dx * dx + dy * dy;
            if dist_sq <= search_rad_sq * 1.44 {
                verified_indices.push(idx);
                verified_total_flux += s.flux;
                if dist_sq > max_r_sq {
                    max_r_sq = dist_sq;
                }
            }
        }

        let kappa = verified_indices.len();
        if kappa < 2 || (max_kappa > 0 && kappa > max_kappa) {
            continue;
        }

        let radius = max_r_sq.sqrt() + fwhm / 2.0;
        if verified_total_flux >= min_detection_flux && radius <= search_radius * 1.5 {
            for &idx in &verified_indices {
                used_sub[idx] = true;
            }

            let snr = verified_total_flux / beam_flux_rms;
            let member_ids: Vec<usize> = verified_indices.iter().map(|&i| sub_peaks[i].id).collect();

            kappa_sources.push(KappaSource {
                id: 0,
                kappa,
                member_ids,
                centroid_x: cen_x,
                centroid_y: cen_y,
                total_flux: verified_total_flux,
                max_amplitude: max_amp,
                radius,
                snr,
            });
        }
    }

    // 6. Hierarchical sort: 1-sources, 2-sources, 3-sources...
    kappa_sources.sort_by(|a, b| {
        a.kappa
            .cmp(&b.kappa)
            .then_with(|| b.total_flux.partial_cmp(&a.total_flux).unwrap_or(std::cmp::Ordering::Equal))
    });

    let mut id_to_sub_idx = std::collections::HashMap::new();
    for (i, s) in sub_peaks.iter().enumerate() {
        id_to_sub_idx.insert(s.id, i);
    }

    for (k_idx, ks) in kappa_sources.iter_mut().enumerate() {
        let kid = k_idx + 1;
        ks.id = kid;
        for &mid in &ks.member_ids {
            if let Some(&s_idx) = id_to_sub_idx.get(&mid) {
                sub_peaks[s_idx].kappa_id = kid;
                sub_peaks[s_idx].kappa = ks.kappa;
            }
        }
    }

    Ok((kappa_sources, sub_peaks, beam_flux_rms))
}

fn py_idx(y: usize, dy: isize) -> usize {
    (y as isize + dy).max(0) as usize
}

fn check_local_max(img: &[f32], width: usize, height: usize, x: usize, y: usize, radius: isize) -> bool {
    let val = img[y * width + x];
    let r = radius.max(1);
    for dy in -r..=r {
        let ny = y as isize + dy;
        if ny < 0 || ny >= height as isize {
            continue;
        }
        for dx in -r..=r {
            if dx == 0 && dy == 0 {
                continue;
            }
            let nx = x as isize + dx;
            if nx < 0 || nx >= width as isize {
                continue;
            }
            if img[ny as usize * width + nx as usize] > val {
                return false;
            }
        }
    }
    true
}

/// Print formatted breakdown summary table of extracted kappa-sources
fn print_extraction_report(
    kappa_sources: &[KappaSource],
    sources: &[Source],
    bg_median: f32,
    bg_rms: f32,
    beam_flux_rms: f32,
    detection_sigma: f32,
    search_radius: f32,
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

    let member_subcomponents_count = sources.iter().filter(|s| s.kappa > 0).count();

    println!("--------------------------------------------------------------------------------");
    println!("Background Level : Median = {:.4}, Pixel RMS = {:.4}, Beam Flux RMS = {:.4}", 
        bg_median, bg_rms, beam_flux_rms);
    println!("Constituents     : {} member peaks assigned to {} kappa-sources ({} candidates screened)", 
        member_subcomponents_count, kappa_sources.len(), sources.len());
    println!("Search Radius    : R_search <= {:.1} pixels", search_radius);
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

/// Write formatted ASCII table catalog (.cat)
fn write_ascii_catalog(
    path: &Path,
    input_filename: &str,
    timestamp: &str,
    kappa_sources: &[KappaSource],
    bg_median: f32,
    bg_rms: f32,
    beam_flux_rms: f32,
    detection_sigma: f32,
    gpu_name: &str,
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "# ==============================================================================")?;
    writeln!(writer, "# kappa-Source Extracted Catalog (GPU Accelerated)")?;
    writeln!(writer, "# Input FITS File   : {}", input_filename)?;
    writeln!(writer, "# Extraction Date   : {}", timestamp)?;
    writeln!(writer, "# GPU Acceleration  : {}", gpu_name)?;
    writeln!(writer, "# Background Median : {:.6}", bg_median)?;
    writeln!(writer, "# Background RMS    : {:.6}", bg_rms)?;
    writeln!(writer, "# Beam Flux RMS     : {:.6}", beam_flux_rms)?;
    writeln!(writer, "# Detection Sigma   : {:.2} (Total Flux >= {:.4})", detection_sigma, detection_sigma * beam_flux_rms)?;
    writeln!(writer, "# Total Extracted   : {}", kappa_sources.len())?;
    writeln!(writer, "# ==============================================================================")?;
    writeln!(writer, "#{:<7} {:<6} {:<10} {:<10} {:<14} {:<12} {:<12} {:<10} {:<10} {:<20}",
        "KAPPA_ID", "KAPPA", "CEN_X", "CEN_Y", "TOTAL_FLUX", "MAX_AMP", "RADIUS_PX", "SNR", "N_MEMBERS", "MEMBER_IDS")?;
    writeln!(writer, "# ------------------------------------------------------------------------------")?;

    for ks in kappa_sources {
        let members_str: String = ks.member_ids.iter().map(|id| id.to_string()).collect::<Vec<_>>().join(",");
        writeln!(writer, "{:<8} {:<6} {:<10.3} {:<10.3} {:<14.4} {:<12.4} {:<12.2} {:<10.2} {:<10} {:<20}",
            ks.id, ks.kappa, ks.centroid_x, ks.centroid_y, ks.total_flux, ks.max_amplitude, ks.radius, ks.snr, ks.member_ids.len(), members_str)?;
    }
    writer.flush()?;
    Ok(())
}

/// Write CSV formatted catalog (.csv)
fn write_csv_catalog(
    path: &Path,
    kappa_sources: &[KappaSource],
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "kappa_id,kappa,cen_x,cen_y,total_flux,max_amplitude,radius_px,snr,n_members,member_ids")?;

    for ks in kappa_sources {
        let members_str: String = ks.member_ids.iter().map(|id| id.to_string()).collect::<Vec<_>>().join(";");
        writeln!(writer, "{},{},{:.4},{:.4},{:.6},{:.6},{:.4},{:.4},{},\"{}\"",
            ks.id, ks.kappa, ks.centroid_x, ks.centroid_y, ks.total_flux, ks.max_amplitude, ks.radius, ks.snr, ks.member_ids.len(), members_str)?;
    }
    writer.flush()?;
    Ok(())
}

/// Write constituent subcomponents table to ASCII (.subcomponents.cat), strictly excluding unassociated peaks (kappa == 0)
fn write_subcomponents_ascii_catalog(
    path: &Path,
    input_filename: &str,
    timestamp: &str,
    sources: &[Source],
) -> std::io::Result<()> {
    let valid_sources: Vec<&Source> = sources.iter().filter(|s| s.kappa > 0).collect();

    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "# ==============================================================================")?;
    writeln!(writer, "# Constituent Subcomponents Peak Catalog (Members of Extracted kappa-Sources)")?;
    writeln!(writer, "# Input FITS File : {}", input_filename)?;
    writeln!(writer, "# Extraction Date : {}", timestamp)?;
    writeln!(writer, "# Total Members   : {}", valid_sources.len())?;
    writeln!(writer, "# ==============================================================================")?;
    writeln!(writer, "#{:<7} {:<10} {:<10} {:<12} {:<12} {:<10} {:<6}",
        "PEAK_ID", "X", "Y", "FLUX", "AMPLITUDE", "KAPPA_ID", "KAPPA")?;
    writeln!(writer, "# ------------------------------------------------------------------------------")?;

    for s in valid_sources {
        writeln!(writer, "{:<8} {:<10.3} {:<10.3} {:<12.4} {:<12.4} {:<10} {:<6}",
            s.id, s.x, s.y, s.flux, s.amplitude, s.kappa_id, s.kappa)?;
    }
    writer.flush()?;
    Ok(())
}

/// Write constituent subcomponents table to CSV (.subcomponents.csv), strictly excluding unassociated peaks (kappa == 0)
fn write_subcomponents_csv_catalog(
    path: &Path,
    sources: &[Source],
) -> std::io::Result<()> {
    let valid_sources: Vec<&Source> = sources.iter().filter(|s| s.kappa > 0).collect();

    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "peak_id,x,y,flux,amplitude,kappa_id,kappa")?;

    for s in valid_sources {
        writeln!(writer, "{},{:.4},{:.4},{:.6},{:.6},{},{}",
            s.id, s.x, s.y, s.flux, s.amplitude, s.kappa_id, s.kappa)?;
    }
    writer.flush()?;
    Ok(())
}

/// Write DS9 region file for visual validation
fn write_ds9_regions(path: &Path, kappa_sources: &[KappaSource], sources: &[Source]) -> std::io::Result<()> {
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

/// Write extracted kappa-sources catalog and constituent sources to FITS binary table
fn write_extracted_fits_catalog(
    path: &Path,
    kappa_sources: &[KappaSource],
    sources: &[Source],
    bg_median: f32,
    bg_rms: f32,
    beam_flux_rms: f32,
    detection_sigma: f32,
    gpu_name: &str,
) -> std::io::Result<()> {
    let valid_sources: Vec<&Source> = sources.iter().filter(|s| s.kappa > 0).collect();
    let mut writer = BufWriter::new(File::create(path)?);

    // Primary Header (Null Image)
    let mut p_cards = Vec::new();
    p_cards.push(make_fits_card("SIMPLE", "T", Some("file conforms to FITS standard")));
    p_cards.push(make_fits_card("BITPIX", "8", Some("Null array")));
    p_cards.push(make_fits_card("NAXIS", "0", Some("No image data in primary HDU")));
    p_cards.push(make_fits_card("EXTEND", "T", Some("Extensions present")));
    p_cards.push(make_fits_card("ACCEL", &fits_str_val("GPU"), Some("GPU acceleration enabled")));
    p_cards.push(make_fits_card("GPUNAME", &fits_str_val(&gpu_name[..gpu_name.len().min(8)]), Some("GPU Device")));
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

    // Extension 1: KAPPA_SRCS
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

    // Extension 2: SOURCES (Constituent member subcomponents only, kappa > 0)
    let num_srcs = valid_sources.len();
    let src_row_size = 4 + 4 + 4 + 4 + 4 + 4 + 4; // 28 bytes

    let mut s_cards = Vec::new();
    s_cards.push(make_fits_card("XTENSION", &fits_str_val("BINTABLE"), Some("binary table extension")));
    s_cards.push(make_fits_card("BITPIX", "8", Some("8-bit bytes")));
    s_cards.push(make_fits_card("NAXIS", "2", Some("2-dimensional table")));
    s_cards.push(make_fits_card("NAXIS1", &src_row_size.to_string(), Some("width of table row in bytes")));
    s_cards.push(make_fits_card("NAXIS2", &num_srcs.to_string(), Some("number of rows in table")));
    s_cards.push(make_fits_card("PCOUNT", "0", Some("size of special data area")));
    s_cards.push(make_fits_card("GCOUNT", "1", Some("one data group")));
    s_cards.push(make_fits_card("TFIELDS", "7", Some("number of fields per row")));
    s_cards.push(make_fits_card("EXTNAME", &fits_str_val("SOURCES"), Some("Member subcomponent peaks")));

    s_cards.push(make_fits_card("TTYPE1", &fits_str_val("ID"), Some("Peak subcomponent ID")));
    s_cards.push(make_fits_card("TFORM1", &fits_str_val("1J"), Some("32-bit integer")));

    s_cards.push(make_fits_card("TTYPE2", &fits_str_val("X"), Some("X position in pixels")));
    s_cards.push(make_fits_card("TFORM2", &fits_str_val("1E"), Some("32-bit float")));
    s_cards.push(make_fits_card("TUNIT2", &fits_str_val("pixel"), None));

    s_cards.push(make_fits_card("TTYPE3", &fits_str_val("Y"), Some("Y position in pixels")));
    s_cards.push(make_fits_card("TFORM3", &fits_str_val("1E"), Some("32-bit float")));
    s_cards.push(make_fits_card("TUNIT3", &fits_str_val("pixel"), None));

    s_cards.push(make_fits_card("TTYPE4", &fits_str_val("FLUX"), Some("Matched filter point flux")));
    s_cards.push(make_fits_card("TFORM4", &fits_str_val("1E"), Some("32-bit float")));

    s_cards.push(make_fits_card("TTYPE5", &fits_str_val("AMPLITUDE"), Some("Peak amplitude")));
    s_cards.push(make_fits_card("TFORM5", &fits_str_val("1E"), Some("32-bit float")));

    s_cards.push(make_fits_card("TTYPE6", &fits_str_val("KAPPA_ID"), Some("Parent kappa-source ID")));
    s_cards.push(make_fits_card("TFORM6", &fits_str_val("1J"), Some("32-bit integer")));

    s_cards.push(make_fits_card("TTYPE7", &fits_str_val("KAPPA"), Some("Parent kappa multiplicity")));
    s_cards.push(make_fits_card("TFORM7", &fits_str_val("1J"), Some("32-bit integer")));

    s_cards.push(make_fits_card("END", "", None));

    let mut s_bytes = Vec::new();
    for c in s_cards {
        s_bytes.extend_from_slice(c.as_bytes());
    }
    let s_len = s_bytes.len();
    let s_pad = (2880 - (s_len % 2880)) % 2880;
    s_bytes.resize(s_len + s_pad, b' ');
    writer.write_all(&s_bytes)?;

    let mut s_table_bytes = Vec::with_capacity(num_srcs * src_row_size + 2880);
    for s in valid_sources {
        s_table_bytes.extend_from_slice(&(s.id as i32).to_be_bytes());
        s_table_bytes.extend_from_slice(&s.x.to_be_bytes());
        s_table_bytes.extend_from_slice(&s.y.to_be_bytes());
        s_table_bytes.extend_from_slice(&s.flux.to_be_bytes());
        s_table_bytes.extend_from_slice(&s.amplitude.to_be_bytes());
        s_table_bytes.extend_from_slice(&(s.kappa_id as i32).to_be_bytes());
        s_table_bytes.extend_from_slice(&(s.kappa as i32).to_be_bytes());
    }
    let st_len = s_table_bytes.len();
    let st_pad = (2880 - (st_len % 2880)) % 2880;
    s_table_bytes.resize(st_len + st_pad, 0u8);
    writer.write_all(&s_table_bytes)?;

    writer.flush()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = ExtractGpuCli::parse();
    let timestamp = get_current_timestamp();

    let input_stem = args
        .input
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("image");

    let input_filename = args
        .input
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("image.fits");

    let parent_dir = args.input.parent().unwrap_or_else(|| Path::new("."));

    // 1. Initialize GPU Context
    println!("Initializing GPU device via wgpu...");
    let gpu = GpuContext::new()?;

    // 2. Read input FITS image
    println!("Reading FITS image from {}...", args.input.display());
    let fits_img = read_fits_image(&args.input)?;
    println!("Loaded image grid: {} x {} pixels", fits_img.width, fits_img.height);

    // 3. Measure background & noise RMS
    println!("Estimating background and RMS noise level...");
    let (bg_median, bg_rms) = estimate_background_and_rms(&fits_img.data);

    // 4. Automatically estimate PSF FWHM if requested or fwhm <= 0
    let effective_fwhm = if args.psf_auto || args.fwhm <= 0.0 {
        println!("Estimating PSF FWHM automatically from bright point sources (SNR > {:.1} or top {} stars)...", args.min_psf_snr, args.psf_samples);
        match estimate_psf_fwhm_auto(&fits_img.data, fits_img.width, fits_img.height, bg_median, bg_rms, args.min_psf_snr, args.psf_samples) {
            Some((measured_fwhm, count, mode)) => {
                println!("Auto-PSF: Measured median FWHM = {:.2} pixels from {} sources (using {}).", measured_fwhm, count, mode);
                measured_fwhm
            }
            None => {
                let fallback = if args.fwhm > 0.0 { args.fwhm } else { 10.0 };
                println!("Auto-PSF: No isolated stars found for calibration. Falling back to default FWHM = {:.1} px.", fallback);
                fallback
            }
        }
    } else {
        args.fwhm
    };

    println!("==============================================================================");
    println!("GPU kappa-Source Two-Pass Extractor (kappa_extract_gpu.bin)");
    println!("==============================================================================");
    println!("Active GPU Device    : {} [{}]", gpu.adapter_name, gpu.backend_name);
    println!("Input FITS File      : {}", args.input.display());
    println!("Timestamp Session    : {}", timestamp);
    println!("Search Radius (R_src): <= {:.1} pixels", args.search_radius);
    println!("Detection Threshold  : Total Flux >= {:.2} * Beam RMS", args.detection_sigma);
    println!("Subcomponent Min SNR : Peak SNR >= {:.2} in matched filter", args.min_sub_snr);
    println!("Candidate Seed SNR   : Seed SNR >= {:.2}", args.seed_snr);
    if args.max_kappa > 0 {
        println!("Max Multiplicity     : kappa <= {}", args.max_kappa);
    } else {
        println!("Max Multiplicity     : Unlimited");
    }
    println!("Subcomponent Limit   : Flux < {:.2} * Beam RMS for kappa >= 2", args.subcomponent_max_sigma);
    println!("Operating PSF FWHM   : {:.2} pixels {}", effective_fwhm, if args.psf_auto { "(Auto-Calibrated)" } else { "" });
    println!("==============================================================================");

    // 5. Run GPU-accelerated two-pass hierarchical extraction
    println!("Running GPU-Accelerated Two-Pass Hierarchical Extraction (1-sources first, then multi-sources)...");
    let (kappa_sources, sources, beam_flux_rms) = extract_kappa_two_pass_gpu(
        &gpu,
        &fits_img.data,
        fits_img.width,
        fits_img.height,
        bg_median,
        bg_rms,
        effective_fwhm,
        args.search_radius,
        args.min_sub_snr,
        args.seed_snr,
        args.detection_sigma,
        args.subcomponent_max_sigma,
        args.max_kappa,
    )?;

    // 6. Display extraction summary report
    print_extraction_report(
        &kappa_sources,
        &sources,
        bg_median,
        bg_rms,
        beam_flux_rms,
        args.detection_sigma,
        args.search_radius,
        args.max_kappa,
    );

    // 7. Determine timestamped output filenames
    let default_fits_name = format!("{}_{}.extracted.fits", input_stem, timestamp);
    let out_fits_path = args.output.unwrap_or_else(|| parent_dir.join(&default_fits_name));

    let generic_fits_path = parent_dir.join(format!("{}.extracted.fits", input_stem));

    println!("Writing timestamped FITS catalog to {}...", out_fits_path.display());
    write_extracted_fits_catalog(&out_fits_path, &kappa_sources, &sources, bg_median, bg_rms, beam_flux_rms, args.detection_sigma, &gpu.adapter_name)?;
    if out_fits_path != generic_fits_path {
        let _ = write_extracted_fits_catalog(&generic_fits_path, &kappa_sources, &sources, bg_median, bg_rms, beam_flux_rms, args.detection_sigma, &gpu.adapter_name);
    }

    // 8. Save ASCII and CSV text catalogs
    if args.save_ascii {
        let cat_path = parent_dir.join(format!("{}_{}.extracted.cat", input_stem, timestamp));
        let csv_path = parent_dir.join(format!("{}_{}.extracted.csv", input_stem, timestamp));
        let sub_cat_path = parent_dir.join(format!("{}_{}.extracted.subcomponents.cat", input_stem, timestamp));
        let sub_csv_path = parent_dir.join(format!("{}_{}.extracted.subcomponents.csv", input_stem, timestamp));

        println!("Writing ASCII kappa-sources catalog to {}...", cat_path.display());
        write_ascii_catalog(&cat_path, input_filename, &timestamp, &kappa_sources, bg_median, bg_rms, beam_flux_rms, args.detection_sigma, &gpu.adapter_name)?;

        println!("Writing CSV kappa-sources catalog to {}...", csv_path.display());
        write_csv_catalog(&csv_path, &kappa_sources)?;

        println!("Writing constituent subcomponents catalog to {}...", sub_cat_path.display());
        write_subcomponents_ascii_catalog(&sub_cat_path, input_filename, &timestamp, &sources)?;

        println!("Writing constituent subcomponents CSV to {}...", sub_csv_path.display());
        write_subcomponents_csv_catalog(&sub_csv_path, &sources)?;
    }

    // 9. Save DS9 region overlays
    if args.save_regions {
        let timestamped_reg = parent_dir.join(format!("{}_{}.extracted.reg", input_stem, timestamp));
        let generic_reg = parent_dir.join(format!("{}.extracted.reg", input_stem));
        println!("Writing DS9 region overlay to {}...", timestamped_reg.display());
        write_ds9_regions(&timestamped_reg, &kappa_sources, &sources)?;
        if timestamped_reg != generic_reg {
            let _ = write_ds9_regions(&generic_reg, &kappa_sources, &sources);
        }
    }

    println!("Extraction complete! Successfully created timestamped catalogs on {}.", gpu.adapter_name);
    Ok(())
}
