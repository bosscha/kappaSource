use std::f32::consts::PI;

/// Supported Point Spread Function (PSF) models for subcomponent convolution
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub enum PsfType {
    /// 2D Gaussian PSF: standard optical / radio clean beam
    Gaussian,
    /// 2D Moffat PSF: atmospheric seeing profile (beta ~ 4.765 or 2.5)
    Moffat,
}

/// Point Spread Function (PSF) convolution kernel and profile
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct Psf {
    pub psf_type: PsfType,
    pub fwhm: f32,
    pub sigma: f32,
    pub moffat_beta: f32,
    pub moffat_alpha: f32,
}

#[allow(dead_code)]
impl Psf {
    /// Create a new Gaussian PSF with given FWHM in pixels
    pub fn new_gaussian(fwhm: f32) -> Self {
        let sigma = fwhm / (2.0 * (2.0 * 2.0f32.ln()).sqrt());
        Self {
            psf_type: PsfType::Gaussian,
            fwhm,
            sigma,
            moffat_beta: 4.765,
            moffat_alpha: 0.0,
        }
    }

    /// Create a new Moffat PSF with given FWHM and power-law index beta
    pub fn new_moffat(fwhm: f32, beta: f32) -> Self {
        let sigma = fwhm / (2.0 * (2.0 * 2.0f32.ln()).sqrt());
        let alpha = fwhm / (2.0 * (2.0f32.powf(1.0 / beta) - 1.0).sqrt());
        Self {
            psf_type: PsfType::Moffat,
            fwhm,
            sigma,
            moffat_beta: beta,
            moffat_alpha: alpha,
        }
    }

    /// Peak amplitude of the PSF normalized to unit total integrated flux
    pub fn peak_response(&self) -> f32 {
        match self.psf_type {
            PsfType::Gaussian => 1.0 / (2.0 * PI * self.sigma * self.sigma),
            PsfType::Moffat => (self.moffat_beta - 1.0) / (PI * self.moffat_alpha * self.moffat_alpha),
        }
    }

    /// Evaluates the normalized PSF(dx, dy) at offset (dx, dy) in pixels
    #[inline(always)]
    pub fn evaluate(&self, dx: f32, dy: f32) -> f32 {
        let r_sq = dx * dx + dy * dy;
        match self.psf_type {
            PsfType::Gaussian => {
                let two_sigma_sq = 2.0 * self.sigma * self.sigma;
                let norm = 1.0 / (2.0 * PI * self.sigma * self.sigma);
                norm * (-r_sq / two_sigma_sq).exp()
            }
            PsfType::Moffat => {
                let alpha_sq = self.moffat_alpha * self.moffat_alpha;
                let norm = (self.moffat_beta - 1.0) / (PI * alpha_sq);
                norm * (1.0 + r_sq / alpha_sq).powf(-self.moffat_beta)
            }
        }
    }

    /// Convolve a delta point source of flux F at sub-pixel position (x, y) with this PSF
    /// into the 2D image raster buffer within a 5-sigma bounding box.
    pub fn convolve_point_source(
        &self,
        x: f32,
        y: f32,
        flux: f32,
        image: &mut [f32],
        width: usize,
        height: usize,
    ) {
        let cutoff_radius = match self.psf_type {
            PsfType::Gaussian => (5.0 * self.sigma).ceil() as isize,
            PsfType::Moffat => (4.0 * self.fwhm).ceil() as isize,
        };

        let x_center = x.round() as isize;
        let y_center = y.round() as isize;

        let x_min = (x_center - cutoff_radius).max(0) as usize;
        let x_max = ((x_center + cutoff_radius).min(width as isize - 1)) as usize;
        let y_min = (y_center - cutoff_radius).max(0) as usize;
        let y_max = ((y_center + cutoff_radius).min(height as isize - 1)) as usize;

        for py in y_min..=y_max {
            let dy = py as f32 - y;
            let row_offset = py * width;

            for px in x_min..=x_max {
                let dx = px as f32 - x;
                let psf_val = self.evaluate(dx, dy);
                image[row_offset + px] += flux * psf_val;
            }
        }
    }
}
