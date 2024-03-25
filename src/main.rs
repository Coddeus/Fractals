// Recovers the recordings of the cameraman who dove (with gravity) right in the Mandelbrot set towards a specific `target`.

use rug::{float::Round, ops::{AssignRound, CompleteRound, MulAssignRound, Pow}, Assign, Complex, Float};
use rayon::prelude::*;
use std::{collections::HashMap, thread::sleep, time::Duration};

extern crate ffmpeg_next as ff;
use ff::Rational;
mod video_builder;
use video_builder::{VideoBuilder, video_options::VideoOptions};

fn main() {
    let mut cameraman: MandelbrotCameraman = MandelbrotCameraman::init(
        480,
        360,
        Rational(1, 2),
        10,
        Float::with_val(128, Float::parse("-1.9997740601362903593126807559602500475710416233856384007148508574291012335").unwrap()),
        Float::with_val(128, Float::parse("-0.0000000032900403214794350534969786759266805967852946505878410088326046927").unwrap()),
        "output.mp4".to_owned(),
    );
    println!("");
    println!("Target x: {}", cameraman.r.targetx);
    println!("Target y: {}", cameraman.r.targety);
    cameraman.dive(); // Goodbye CPU
}

struct MandelbrotCameraman {
    r: Renderer,
    e: Encoder,
}

impl MandelbrotCameraman {
    fn init(width: u32, height: u32, framerate: Rational, frames: u32, targetx: Float, targety: Float, output_path: String) -> Self {
        Self {
            r: Renderer::init(width, height, targetx, targety),
            e: Encoder::init(width, height, frames, framerate, output_path),
        }
    }

    fn dive(&mut self) {
        self.e.encoder.start_encoding().unwrap();

        self.r.generate_mandelbrot();
        self.e.encoder
            .push_video_data(&self.r.pixels)
            .unwrap();
        self.e.encoder.step_encoding().unwrap();

        for frame in 1..self.e.frames {
            println!("Generated frame {}", frame);
            self.r.scale.mul_assign_round(0.7, Round::Nearest);
            // self.r.max_iter*=3;

            self.r.generate_mandelbrot();
            self.e.encoder
                .push_video_data(&self.r.pixels)
                .unwrap();
            self.e.encoder.step_encoding().unwrap();
        }
        sleep(Duration::from_secs(1));
        println!("");
        self.e.encoder.finish_encoding().unwrap();
    }
}


struct Renderer {
    // General consts
    width: u32,
    height: u32,

    // Buffer
    pixels: Vec<u8>,

    // Target-dependent
    targetx: Float,
    targety: Float,
    
    // Scale-dependent
    scale: Float,
    #[allow(unused)]
    max_iter: u32,
    #[allow(unused)]
    max_period: usize,
    precision: u32,

    // Rendering variables
    one: Float,
    rat: Float,

    cx_min: Float,
    cx_max: Float,
    cy_min: Float,
    cy_max: Float,

    pixel_width: Float,
    pixel_height: Float,
}

impl Default for Renderer {
    fn default() -> Self {
        Self {
            // General consts
            width: 360,
            height: 240,

            // Buffer
            pixels: vec![0u8; (1920 * 1080 * 3) as usize],

            // Target-dependent
            targetx: Float::new(128),
            targety: Float::new(128),

            // Scale-dependent
            scale: Float::with_val(128, 2.0),
            max_iter: 30,
            max_period: 100,
            precision: 128,

            // Rendering variables
            one: Float::new(128),
            rat: Float::new(128),

            cx_min: Float::new(128),
            cx_max: Float::new(128),
            cy_min: Float::new(128),
            cy_max: Float::new(128),

            pixel_width: Float::new(128),
            pixel_height: Float::new(128),
        }
    }
}

impl Renderer {
    fn init(width: u32, height: u32, targetx: Float, targety: Float) -> Self {
        Self {
            // General consts
            width,
            height,
        
            // Buffer
            pixels: vec![0u8; (width * height * 3) as usize],
        
            // Target-dependent
            targetx,
            targety,
        
            // Rendering variables
            one: Float::with_val(128, 1.0),
            rat: Float::with_val(128, width as f64 / height as f64),
    
            // Other non-user defaults (Scale-dependent, Rendering variables)
            ..Renderer::default()
        }
    }

    fn generate_mandelbrot(&mut self) {
        self.cx_min.assign_round(&self.targetx - &self.rat * &self.scale, Round::Up);
        self.cx_max.assign_round(&self.targetx + &self.rat * &self.scale, Round::Up);
        self.cy_min.assign_round(&self.targety - &self.one * &self.scale, Round::Up);
        self.cy_max.assign_round(&self.targety + &self.one * &self.scale, Round::Up);

        self.pixel_width .assign_round((&self.cx_max - &self.cx_min).complete(self.precision) / self.width  as f64, Round::Up);
        self.pixel_height.assign_round((&self.cy_max - &self.cy_min).complete(self.precision) / self.height as f64, Round::Up);

        self.pixels
            .par_iter_mut()
            .chunks(3)
            .enumerate()
            .for_each(|(i, mut pixel)| {
                let x = Float::with_val(self.precision, (i as u32 % self.width) as f64);
                let y = Float::with_val(self.precision, (i as u32 / self.width) as f64);
                let cx = Float::with_val(self.precision, &self.cx_min + &x * &self.pixel_width);
                let cy = Float::with_val(self.precision, &self.cy_min + &y * &self.pixel_height);
                let c = Complex::with_val(64, (cx, cy));
                let period: usize = mandelbrot_period_check(c.clone(), self.max_period);
                // *pixel[0] = (!inner) as u8 * 255;
                // *pixel[1] = (!inner) as u8 * 255;
                // *pixel[2] = (!inner) as u8 * 255;
                if period == 0 {
                    *pixel[0] = 0;
                } else {
                    let interior_dist = m_interior_distance(&Complex::with_val(128, (0, 0)), &c, period as u32);
                    *pixel[0] = (interior_dist.to_f32() * 255.0) as u8;
                }
                // let escape_time = mandelbrot_escape_time(c, self.max_iter);
                // *pixel[0] = escape_time;
                // *pixel[1] = escape_time;
                // *pixel[2] = escape_time;
        });
    }
}

#[allow(unused)]
fn mandelbrot_period_check(c: Complex, period: usize) -> usize {
    let mut z = Complex::with_val(64, (Float::with_val(128, 0), Float::with_val(128, 0)));
    let mut done: Vec<Complex> = Vec::with_capacity(period);
    let mut round: Complex;
    let mut inner: bool = false;
    for n in 0..50 {
        z = z.square() + &c;
        round = Complex::with_val(100, &z);
        // if { let norm = Float::with_val(128, z.real() * z.real() + z.imag() * z.imag()); norm<16 }  { inner = true; break; } // Another check for the distance (not optimized for really zoomed-in pics)
        for m in 0..n {
            if done[m] == round { return n }
        };
        done.push(round);
    };
    0
}

#[allow(unused)]
fn mandelbrot_escape_time(c: Complex, max_iter: u32) -> u8 {
    let mut z = Complex::with_val(128, (Float::with_val(128, 0), Float::with_val(128, 0)));
    let mut n = 0;
    while n < max_iter && { let norm = Float::with_val(128, z.real() * z.real() + z.imag() * z.imag()); norm<16} {
        z = z.square() + &c;
        n += 1;
    }
    255 - (n as f32 / max_iter as f32 * 255.0) as u8
}

fn cnorm(z: &Complex) -> Float {
    let re_squared = z.real() * z.real();
    let im_squared = z.imag() * z.imag();
    Float::with_val(128, re_squared + im_squared)
}

fn m_interior_distance(z0: &Complex, c: &Complex, p: u32) -> Float {
    let mut z = z0.clone();
    let mut dz = Complex::with_val(128, (Float::with_val(128, 1), Float::with_val(128, 0)));
    let mut dzdz = Complex::with_val(128, (Float::with_val(128, 0), Float::with_val(128, 0)));
    let mut dc = Complex::with_val(128, (Float::with_val(128, 0), Float::with_val(128, 0)));
    let mut dcdz = Complex::with_val(128, (Float::with_val(128, 0), Float::with_val(128, 0)));

    for _ in 0..p {
        dcdz.assign(Float::with_val(128, 2) * z.clone() * &dcdz + &dz * &dc);
        dc.assign((&2 * z.clone() * &dc) + Complex::with_val(128, (Float::with_val(128, 1), Float::with_val(128, 0))));
        dzdz.assign(&2 * dz.clone() * &dz + &z * &dzdz);
        dz.assign(&2 * z.clone() * &dz);
        z = z.pow(2) + c;
    }

    (Float::with_val(128, 1) - cnorm(&dz)) / cnorm(&(&dcdz + dzdz * &dc / &(Complex::with_val(128, (1, 0)) - &dz)))
}

// fn m_distance(N: u32, R: &Float, c: &Complex) -> Float {
//     let mut dc = Complex::with_val(128, (Float::with_val(128, 0), Float::with_val(128, 0)));
//     let mut z = Complex::with_val(128, (Float::with_val(128, 0), Float::with_val(128, 0)));
//     let mut m = Float::POSITIVE_INFINITY;
//     let mut p = 0;
// 
//     for n in 1..=N {
//         dc = &(&2 * &z * &dc) + Complex::with_val(128, (Float::with_val(128, 1), Float::with_val(128, 0)));
//         z = z.pow(2) + c.clone();
// 
//         if z.norm() > *R {
//             return &(&(&2 * &z.norm()) * &(&z.norm().ln()) / &dc.norm());
//         }
// 
//         if z.norm() < m {
//             m = z.norm();
//             p = n;
//             let z0 = m_interior_distance(&z, c, p);
//             let mut w = z0.clone();
//             let mut dw = Complex::with_val(128, (Float::with_val(128, 1), Float::with_val(128, 0)));
// 
//             for _ in 0..p {
//                 dw = &(&2 * &w * &dw);
//                 w = w.pow(2) + c.clone();
//             }
// 
//             if dw.norm() <= Float::with_val(128, 1) {
//                 return m_interior_distance(&z0, c, p);
//             }
//         }
//     }
// 
//     Float::with_val(128, 0)
// }


struct Encoder {
    frames: u32,
    encoder: VideoBuilder,
}
impl Encoder {
    fn init(width: u32, height: u32, frames: u32, framerate: Rational, output_path: String) -> Self {
        ffmpeg_next::init().unwrap();
        Self {
            frames,
            encoder: VideoBuilder::new(VideoOptions {
                output_path,
                metadata: HashMap::new(),
                video_time_base: framerate,
                video_codec: "libx264rgb".to_owned(),
                video_codec_params: Default::default(),
                pixel_format_in: "rgb24".to_string(),
                pixel_format_out: "rgb24".to_string(),
                resolution_in: (width, height),
                resolution_out: (width, height),
            }).unwrap(),
        }
    }
}