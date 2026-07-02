use std::fmt::{Display, Formatter};

#[derive(Debug)]
pub enum SuperCountingBloomError {
    InvalidConfig(String),
    FastxOpen { path: String, message: String },
    InvalidIndexFormat(String),
    ChannelClosed(String),
    WorkerPanic(&'static str),
    InternalState(String),
    Io { path: String, message: String },
}

impl Display for SuperCountingBloomError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidConfig(message) => write!(f, "invalid config: {message}"),
            Self::FastxOpen { path, message } => {
                write!(f, "failed to open FASTA/FASTQ '{path}': {message}")
            }
            Self::InvalidIndexFormat(message) => write!(f, "invalid index format: {message}"),
            Self::ChannelClosed(message) => write!(f, "worker channel closed: {message}"),
            Self::WorkerPanic(context) => write!(f, "{context} worker panicked"),
            Self::InternalState(message) => write!(f, "internal state error: {message}"),
            Self::Io { path, message } => write!(f, "I/O error for '{path}': {message}"),
        }
    }
}

impl std::error::Error for SuperCountingBloomError {}
