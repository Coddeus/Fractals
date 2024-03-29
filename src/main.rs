// Recovers the recordings of the cameraman who dove (with gravity) right in the Mandelbrot set towards a specific `target`.

use rug::{float::Round, ops::{AssignRound, CompleteRound, MulAssignRound, Pow}, Assign, Complex, Float};
use rayon::prelude::*;
use std::{collections::HashMap, time::Instant};

extern crate ffmpeg_next as ff; use ff::Rational;
mod video_builder;
use video_builder::{VideoBuilder, video_options::VideoOptions};

fn main() {
    let width: u32          = 160;
    let height: u32         = 120;
    let precision: u32      = 64;
    let max_iter: u32       = 600;
    let max_period: usize   = 20;
    let framerate = Rational(1, 24);
    let frames = 192;
    let display = false;
    let mut cameraman: MandelbrotCameraman = MandelbrotCameraman::init(
        width,
        height,
        Float::with_val(precision, Float::parse("-1.9997740601362903593126807559602500475710416233856384007148508574291012335984591928248364190215796259575718318799960175396106897245889581254834492701372949636783094955897931317174101259095891469501748126725148714587333938548443819033709904187344921523413310221887295870857771431011674873342592895504186325482220668710775749899926429101099841583206278295793058921625817004481783699245865364627140554117737774937789463895102748671351750212506004241754983473339789940659968568850689353099462034492524909310777724611601104714214019347435268544619054369865904944457792527241696528695821059623303046651934176389789308453627525109367436309636375268231073110318555064708363221007235298404379856922536028913291478442839193381367508575286692330907891402483843152933153748354825108021776358693600801782904774626935265722056455978643513448489091026679036353407968495795003386248005939867069799946547181378474054113117046900560609110812439442002663909295191705374444149326937073460052706389967886211172676612720028299452788285465688867116337489531157494508508315428488520037968118008255840569742557333862639124341116894229885253643651920014148109308402199399127712572209466874971603743536096235390414412927589954662603878558182262865151900604451937214289079939337905846647369517138325441736853526711818853134657265043099539402286244220638999824999819000131999789999857999958").unwrap()),
        Float::with_val(precision, Float::parse("-0.0000000032900403214794350534969786759266805967852946505878410088326046927853549452991056352681196631150325234171525664335353457621247922992470898021063583060218954321140472066153878996044171428801408137278072521468882260382336298800961530905692393992277070012433445706657829475924367459793505729004118759963065667029896464160298608486277109065108339157276150465318584383757554775431988245033409975361804443001325241206485033571912765723551757793318752425925728969073157628495924710926832527350298951594826689051400340011140584507852761857568007670527511272585460136585523090533629795012272916453744029579624949223464015705500594059847850617137983380334184205468184810116554041390142120676993959768153409797953194054452153167317775439590270326683890021272963306430827680201998682699627962109145863135950941097962048870017412568065614566213639455841624790306469846132055305041523313740204187090956921716703959797752042569621665723251356946610646735381744551743865516477084313729738832141633286400726001116308041460406558452004662264165125100793429491308397667995852591271957435535504083325331161340230101590756539955554407081416407239097101967362512942992702550533040602039494984081681370518238283847808934080198642728761205332894028474812918370467949299531287492728394399650466260849557177609714181271299409118059191938687461000000000000000000000000000000000000").unwrap()),
        precision,
        max_iter,
        max_period,
        Float::with_val(precision, Float::parse("4.0").unwrap()),
        Some("output.mp4".to_owned()), 
        framerate,
        frames,
        display
    );
    println!("");
    println!("Target x: {}", cameraman.r.targetx);
    println!("Target y: {}", cameraman.r.targety);
    
    let start = Instant::now();

    // for i in 0..100 {
    //     cameraman.save_image(i);
    //     println!("Saved frame {}", i);
    //     cameraman.r.scale.mul_assign_round(0.3, Round::Nearest);
    // }
    // cameraman.save_image(0);
    cameraman.dive(); // Goodbye CPU

    println!("{}", start.elapsed().as_secs_f64());
}

struct MandelbrotCameraman {
    r: Renderer,
    e: Option<Encoder>,
    // d: Option<Display>,
}

impl MandelbrotCameraman {
    fn init(width: u32, height: u32, targetx: Float, targety: Float, precision: u32, max_iter: u32, max_period: usize, base_scale: Float, output_path: Option<String>, framerate: Rational, frames: u32, display: bool) -> Self {
        let r = Renderer::init(width, height, targetx, targety, precision, max_iter, max_period, base_scale);
        let e = if output_path.is_some() { Some(Encoder::init(width, height, frames, framerate, output_path.unwrap())) } else { None };
        // let d = if display { Some(Display::init(width, height, &r.pixels)) } else { None };
        
        Self {
            r,
            e,
            // d,
        }
    }

    #[allow(unused)]
    fn dive(&mut self) {
        let mut enc = self.e.take().unwrap();

        enc.encoder.start_encoding().unwrap();
        println!("");
        for frame in 1..enc.frames+1 {
            self.r.scale.mul_assign_round(0.91, Round::Nearest);

            self.r.generate_mandelbrot();

            enc.encoder
                .push_video_data(&self.r.pixels)
                .unwrap();
            enc.encoder.step_encoding().unwrap();

            if frame%10==0 {
                println!("Generated frame {}", frame);
            }
        }
        println!("");
        enc.encoder.finish_encoding().unwrap();

        self.e = Some(enc);
    }

    #[allow(unused)]
    fn save_image(&mut self, n: u32) {
        use std::path::Path;
        use std::fs::File;
        use std::io::BufWriter;

        self.r.generate_mandelbrot();

        let filename = format!("output/{:09}.png", n);
        let path = Path::new(filename.as_str());
        let file = File::create(path).unwrap();
        let ref mut w = BufWriter::new(file);

        let mut encoder = png::Encoder::new(w, self.r.width, self.r.height);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        encoder.set_source_gamma(png::ScaledFloat::from_scaled(45455)); // 1.0 / 2.2, scaled by 100000
        encoder.set_source_gamma(png::ScaledFloat::new(1.0 / 2.2));     // 1.0 / 2.2, unscaled, but rounded
        let source_chromaticities = png::SourceChromaticities::new(     // Using unscaled instantiation here
            (0.31270, 0.32900),
            (0.64000, 0.33000),
            (0.30000, 0.60000),
            (0.15000, 0.06000)
        );
        encoder.set_source_chromaticities(source_chromaticities);
        let mut writer = encoder.write_header().unwrap();

        writer.write_image_data(&self.r.pixels).unwrap(); // Save
    }
}


struct Renderer {
    // General consts
    width: u32,
    height: u32,

    // Buffers
    escape: Vec<u32>,
    palette: Vec<u8>,
    pixels: Vec<u8>,
    // iterations: Vec<u32>,

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

// impl Default for Renderer {
//     fn default() -> Self {
//         Self {
//             // General consts
//             width: 360,
//             height: 240,
// 
//             // Buffer
//             pixels: vec![0u8; (1920 * 1080 * 3) as usize],
// 
//             // Target-dependent
//             targetx: Float::new(199),
//             targety: Float::new(199),
// 
//             // Scale-dependent
//             scale: Float::with_val(199, Float::parse("0.0000000000000000001").unwrap()),
//             max_iter: 1000,
//             max_period: 100,
//             precision: 199,
// 
//             // Rendering variables
//             one: Float::with_val(199, 1.0),
//             rat: Float::new(199),
// 
//             cx_min: Float::new(199),
//             cx_max: Float::new(199),
//             cy_min: Float::new(199),
//             cy_max: Float::new(199),
// 
//             pixel_width: Float::new(199),
//             pixel_height: Float::new(199),
//         }
//     }
// }

impl Renderer {
    fn init(width: u32, height: u32, targetx: Float, targety: Float, precision: u32, max_iter: u32, max_period: usize, base_scale: Float) -> Self {
        Self {
            // General consts
            width,
            height,
        
            // Buffer
            escape: vec![0u32; (width * height) as usize],
            palette: [255u8, 245, 235, 225, 215, 205, 195, 185, 175, 165, 155, 145, 135, 125, 115, 105, 95, 85, 75, 65, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205, 215, 225, 235, 245].into_iter().cycle().take(max_iter as usize).collect(),
            pixels: vec![0u8, 0, 0, 255].into_iter().cycle().take((width * height * 4) as usize).collect(),
            // iterations: vec![0u32; max_iter as usize + 1],
        
            // Target-dependent
            targetx,
            targety,

            // Scale-dependent
            scale: Float::with_val(precision, base_scale),
            max_iter,
            max_period,
            precision,
        
            // Rendering variables
            one: Float::with_val(precision, Float::parse("1.0").unwrap()),
            rat: Float::with_val(precision, width as f64 / height as f64),

            cx_min: Float::new(precision),
            cx_max: Float::new(precision),
            cy_min: Float::new(precision),
            cy_max: Float::new(precision),

            pixel_width: Float::new(precision),
            pixel_height: Float::new(precision),
        }
    }

    fn generate_mandelbrot(&mut self) {
        self.cx_min.assign_round(&self.targetx - &self.rat * &self.scale, Round::Up);
        self.cx_max.assign_round(&self.targetx + &self.rat * &self.scale, Round::Up);
        self.cy_min.assign_round(&self.targety - &self.one * &self.scale, Round::Up);
        self.cy_max.assign_round(&self.targety + &self.one * &self.scale, Round::Up);

        self.pixel_width .assign_round((&self.cx_max - &self.cx_min).complete(self.precision) / self.width  as f64, Round::Up);
        self.pixel_height.assign_round((&self.cy_max - &self.cy_min).complete(self.precision) / self.height as f64, Round::Up);

        // self.iterations = vec![0u32; self.max_iter as usize + 1];

        // Find escape times
        self.escape
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, esc)| {
                // if i%1000==0 { println!("Pixel {}", i) }
                let x = Float::with_val(self.precision, (i as u32 % self.width) as f64);
                let y = Float::with_val(self.precision, (i as u32 / self.width) as f64);
                let cx = Float::with_val(self.precision, &self.cx_min + &x * &self.pixel_width);
                let cy = Float::with_val(self.precision, &self.cy_min + &y * &self.pixel_height);
                let c = Complex::with_val(self.precision, (cx, cy));
                *esc = Self::mandelbrot_escape_time(c, self.max_iter, self.precision);
                // self.iterations[*esc as usize]+=1;
            }
        );


        // Colorize pixels accordingly
        // let (min, max) = (*self.escape.par_iter().min().unwrap() as f64, *self.escape.par_iter().filter(|&&u| u<self.max_iter).max().unwrap() as f64);
        // let diff = max-min;
        self.pixels
            .par_iter_mut()
            .chunks(4)
            .enumerate()
            .for_each(|(i, mut pixel)| {
                // *pixel[0] = (((self.escape[i] as f64 - min) / diff).sqrt() * 255.999f64).floor() as u8;
                let index = self.escape[i];
                if index==self.max_iter {
                    *pixel[0] = 0;
                } else {
                    *pixel[0] = self.palette[self.escape[i] as usize];// ((self.escape[i] as i32 * 200) % 511 - 255).abs() as u8; //  as f64 / total as f64 * 255.0) // self.iterations[   as usize]
                }
            });
    }
    
    #[allow(unused)]
    fn mandelbrot_escape_time(c: Complex, max_iter: u32, precision: u32) -> u32 {
        let mut z = Complex::with_val(precision, (Float::with_val(precision, 0), Float::with_val(precision, 0)));
        let mut n = 0;
        while n < max_iter && { let norm = Float::with_val(precision, z.real() * z.real() + z.imag() * z.imag()); norm<16} {
            z = z.square() + &c;
            n += 1;
        }
        // println!("{}", n); // ITERS
        n
    }

    #[allow(unused)]
    fn mandelbrot_period_check(c: Complex, max_period: usize, precision: u32) -> usize {
        let mut z = Complex::with_val(64, (Float::with_val(precision, 0), Float::with_val(precision, 0)));
        let mut done: Vec<Complex> = Vec::with_capacity(max_period);
        let mut round: Complex;
        let mut inner: bool = false;
        for n in 0..max_period {
            z = z.square() + &c;
            round = Complex::with_val(precision, &z);
            // if { let norm = Float::with_val(precision, z.real() * z.real() + z.imag() * z.imag()); norm<16 }  { inner = true; break; } // Another check for the distance (not optimized for really zoomed-in pics)
            for m in 0..n {
                if done[m] == round { return n }
            };
            done.push(round);
        };
        0
    }
    
    #[allow(unused)]
    fn cnorm(z: &Complex, precision: u32) -> Float {
        let re_squared = z.real() * z.real();
        let im_squared = z.imag() * z.imag();
        Float::with_val(precision, re_squared + im_squared)
    }
    
    #[allow(unused)]
    fn m_interior_distance(z0: &Complex, c: &Complex, p: u32, precision: u32) -> Float {
        let mut z = z0.clone();
        let mut dz = Complex::with_val(precision, (Float::with_val(precision, 1), Float::with_val(precision, 0)));
        let mut dzdz = Complex::with_val(precision, (Float::with_val(precision, 0), Float::with_val(precision, 0)));
        let mut dc = Complex::with_val(precision, (Float::with_val(precision, 0), Float::with_val(precision, 0)));
        let mut dcdz = Complex::with_val(precision, (Float::with_val(precision, 0), Float::with_val(precision, 0)));
    
        for _ in 0..p {
            dcdz.assign(Float::with_val(precision, 2) * z.clone() * &dcdz + &dz * &dc);
            dc.assign((&2 * z.clone() * &dc) + Complex::with_val(precision, (Float::with_val(precision, 1), Float::with_val(precision, 0))));
            dzdz.assign(&2 * dz.clone() * &dz + &z * &dzdz);
            dz.assign(&2 * z.clone() * &dz);
            z = z.pow(2) + c;
        }
    
        (Float::with_val(precision, 1) - Self::cnorm(&dz, precision)) / Self::cnorm(&(&dcdz + dzdz * &dc / &(Complex::with_val(precision, (1, 0)) - &dz)), precision)
    }
    
    // fn m_distance(N: u32, R: &Float, c: &Complex, precision: u32) -> Float {
    //     let mut dc = Complex::with_val(precision, (Float::with_val(precision, 0), Float::with_val(precision, 0)));
    //     let mut z = Complex::with_val(precision, (Float::with_val(precision, 0), Float::with_val(precision, 0)));
    //     let mut m = Float::POSITIVE_INFINITY;
    //     let mut p = 0;
    // 
    //     for n in 1..=N {
    //         dc = &(&2 * &z * &dc) + Complex::with_val(precision, (Float::with_val(precision, 1), Float::with_val(precision, 0)));
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
    //             let mut dw = Complex::with_val(precision, (Float::with_val(precision, 1), Float::with_val(precision, 0)));
    // 
    //             for _ in 0..p {
    //                 dw = &(&2 * &w * &dw);
    //                 w = w.pow(2) + c.clone();
    //             }
    // 
    //             if dw.norm() <= Float::with_val(precision, 1) {
    //                 return m_interior_distance(&z0, c, p);
    //             }
    //         }
    //     }
    // 
    //     Float::with_val(precision, 0)
    // }
}


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
                video_codec: "h264_nvenc".to_owned(),
                video_codec_params: Default::default(),
                pixel_format_in: "rgba".to_string(),
                pixel_format_out: "rgba".to_string(),
                resolution_in: (width, height),
                resolution_out: (width, height),
            }).unwrap(),
        }
    }
}