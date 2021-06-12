use chrono::{DateTime, NaiveDateTime, TimeZone, Utc};
use clap::{App, Arg};
use rand::Rng;
use realfft::RealFftPlanner;
use rustfft::num_complex::Complex;
use std::env;
use std::f64::consts::PI;
use std::fs::File;
use std::iter::FromIterator;
use std::path::Path;
use std::time::SystemTime;

fn main() {
    let matches = App::new("Convolve")
    .version("0.1.0")
    .author("Sam Pluta")
    .about("Convolve two sound files")
    .arg(
        Arg::with_name("file1")
        .short("f")
        .long("file1")
        .takes_value(true)
        .help("Audio File 1"),
    )
    .arg(
        Arg::with_name("file2")
        .short("g")
        .long("file2")
        .takes_value(true)
        .help("Audio File 2"),
    )
    .arg(
        Arg::with_name("file1_slice_size")
        .short("s")
        .long("slice")
        .takes_value(true)
        .help("The slice size (a power of two) used for input file 1 (optional - defaults to full file)."),
    )
    .arg(
        Arg::with_name("out_file")
        .short("o")
        .long("out_file_name")
        .takes_value(true)
        .help("The name of the output file (optional - defaults to file1_file2.wav"),
    )
    .get_matches();
    let file_name0 = matches.value_of("file1").unwrap_or("sound/tng.wav");
    println!("The input file1 is: {}", file_name0);

    let file_name1 = matches.value_of("file2").unwrap_or("sound/tng.wav");
    println!("The input file2 is: {}", file_name1);

    let s_in = matches.value_of("file1_slice_size").unwrap_or("0");
    let num_slices: usize = match s_in.parse() {
        Ok(n) => n,
        Err(_) => {
            eprintln!("error: slices argument not an integer");
            return;
        }
    };

    let path = Path::new(file_name0);
    let parent = path.parent();
    let file_stem = path.file_stem();
    let extension = path.extension();

    // let path2 = Path::new(file_name1);
    // let file_stem2 = path2.file_stem().unwrap().to_str();

    let temp_out_file = file_name0.replace(".wav", &format!("{}{}", "_convolved", ".wav"));
    let out_file = matches.value_of("out_file").unwrap_or(&temp_out_file);
    println!("The output file will be: {}", out_file);
    println!("");
    //println!("If your output file exceeds the limitations of the WAV format, the OS will think the file is too short, but all of the data will be there and it will be readable with software like Reaper (by normal import) or Audacity (via import Raw Data)");
    //let mut sound_file0 = hound::WavReader::open(file_name0).unwrap();
    //let mut sound_file1 = hound::WavReader::open(file_name1).unwrap();
    let mut channels0 = get_file_channels(file_name0);
    let mut channels1 = get_file_channels(file_name1);
    
    let npot = usize::next_power_of_two;
    let mut fft_size = std::cmp::min(npot(channels0[0].len()), npot(channels1[0].len()));

    if fft_size > usize::pow(2, 17) {
        fft_size = usize::pow(2, 17);
    }
    println!("fftSize: {}", fft_size);

    for index in 0..channels0.len() {
        make_whole (&mut channels0[index], fft_size);
    }
    for index in 0..channels1.len() {
        make_whole (&mut channels1[index], fft_size);
    } 

    let out_length = (channels0[0].len()+channels1[0].len()) as usize;
    let num_out_chans = std::cmp::max(channels0.len(), channels1.len()) as usize;
    let mut out_channels: Vec<Vec<f64>> =
        vec![vec![0.0_f64; out_length]; num_out_chans];

    println!("channels0 len: {}", channels0[0].len());
    println!("channels1 len: {}", channels1[0].len());

    let mut real_planner = RealFftPlanner::<f64>::new();
    let fft = real_planner.plan_fft_forward(fft_size*2);
    let mut spectrum0 = fft.make_output_vec();
    let mut spectrum1 = fft.make_output_vec();
    let ifft = real_planner.plan_fft_inverse(fft_size*2);
    let mut out_frame = ifft.make_output_vec();
    let mut conv = fft.make_output_vec();//vec![Complex{re:0.0, im:0.0}; fft_size*2];

    let mut aseg = vec![0.0; fft_size*2];
    let mut bseg = vec![0.0; fft_size*2];

    println!("{}", fft_size);
    for chan in 0..num_out_chans {
        println!("");
        print!("chan {}", chan);
        for i in 0..channels0[0].len()/fft_size {
            //aseg.clone_from_slice(&channels0[chan%channels0.len() as usize][fft_size*i..fft_size*(i+1)]);
            for num in 0..fft_size {
                aseg[num]=channels0[chan%channels0.len() as usize][fft_size*i + num];
                aseg[fft_size+num]=0.0;
            }
            fft.process(&mut aseg, &mut spectrum0).unwrap();
            for bi in 0..channels1[0].len()/fft_size {
                //bseg.clone_from_slice(&channels1[chan%channels0.len() as usize][fft_size*bi..fft_size*(bi+1)]);
                for num in 0..fft_size {
                    bseg[num]=channels1[chan%channels1.len() as usize][fft_size*bi + num];
                    bseg[fft_size+num]=0.0;
                }
                print!(".");
                fft.process(&mut bseg, &mut spectrum1).unwrap();
                for sp in 0..spectrum0.len() {
                  conv[sp] = spectrum0[sp]*spectrum1[sp]/(fft_size as f64*2.0);
                }
                ifft.process(&mut conv, &mut out_frame);
                for i2 in 0..out_frame.len() {
                    out_channels[chan][fft_size*(i+bi)+i2] += out_frame[i2]/(fft_size as f64*2.0).sqrt()/2.0;
                }
            }
        }
    }
    let mut max = 0.0;
    
    for i in 0..out_channels.len(){
        for i2 in 0..out_channels[0].len() {
            if out_channels[i][i2]>max {max = out_channels[i][i2]}
        }
    }
    println!("{}", max);
    for i in 0..out_channels.len(){
        for i2 in 0..out_channels[0].len() {
            out_channels[i][i2] /= max;
        }
    }

    let sound_file = hound::WavReader::open(file_name0).unwrap();
    let spec = hound::WavSpec {
        channels: out_channels.len() as u16,
        sample_rate: sound_file.spec().sample_rate,
        bits_per_sample: 32,
        sample_format: hound::SampleFormat::Float,
    };
    
    let mut writer = hound::WavWriter::create(out_file, spec).unwrap();
    for samp in 0..(out_channels[0].len()) {
        for chan in 0..spec.channels {
            writer
            .write_sample(out_channels[chan as usize][samp] as f32)
            .unwrap();
        }
    }
    writer.finalize().unwrap();
}

fn make_whole (vec: &mut Vec<f64>, fft_size: usize) {
    let missing = (vec.len() as f64 /fft_size as f64).ceil() as usize *fft_size - vec.len();
    vec.append(&mut vec![0.0; missing]);
}

fn get_file_channels(file_name: &str) -> Vec<Vec<f64>> {
    let mut sound_file = hound::WavReader::open(file_name).unwrap();
    let mut intemp = vec![0.0; 0];
    if sound_file.spec().sample_format == hound::SampleFormat::Float {
        intemp.append(
            &mut sound_file
                .samples::<f32>()
                .map(|x| x.unwrap() as f64)
                .collect::<Vec<_>>(),
        );
    } else {
        intemp.append(
            &mut sound_file
                .samples::<i32>()
                .map(|x| x.unwrap() as f64)
                .collect::<Vec<_>>(),
        );
        let bits = sound_file.spec().bits_per_sample;
        for iter in 0..intemp.len() {
            intemp[iter] = intemp[iter] / (f64::powf(2.0, bits as f64));
        }
    };
    let chunked: Vec<Vec<f64>> = intemp
        .chunks(sound_file.spec().channels as usize)
        .map(|x| x.to_vec())
        .collect();
    let mut channels = transpose(chunked);

    return channels;
}

fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>>
where
    T: Clone,
{
    assert!(!v.is_empty());
    (0..v[0].len())
        .map(|i| v.iter().map(|inner| inner[i].clone()).collect::<Vec<T>>())
        .collect()
}
