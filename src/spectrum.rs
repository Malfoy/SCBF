use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use crate::{Result, SuperCountingBloomError};

pub type ApproxSpectrum = Vec<(u64, u64)>;

pub fn write_spectrum(spectrum: &ApproxSpectrum, output: Option<&Path>) -> Result<()> {
    match output {
        Some(path) => {
            let file = File::create(path).map_err(|err| SuperCountingBloomError::Io {
                path: path.display().to_string(),
                message: err.to_string(),
            })?;
            write_spectrum_to(file, spectrum)
        }
        None => write_spectrum_to(io::stdout().lock(), spectrum),
    }
}

fn write_spectrum_to<W: Write>(writer: W, spectrum: &ApproxSpectrum) -> Result<()> {
    let mut writer = BufWriter::new(writer);
    writeln!(writer, "estimated_count\tkmers").map_err(io_error)?;
    for &(count, distinct) in spectrum {
        writeln!(writer, "{count}\t{distinct}").map_err(io_error)?;
    }
    Ok(())
}

fn io_error(err: io::Error) -> SuperCountingBloomError {
    SuperCountingBloomError::Io {
        path: "<stream>".to_string(),
        message: err.to_string(),
    }
}
