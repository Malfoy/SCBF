use std::io::{self, Read, Write};

pub trait Counter: Copy + Default + Ord + Send + Sync + 'static {
    fn max_value() -> Self;
    fn saturating_increment(&mut self);
    fn to_u64(self) -> u64;
}

pub trait SerializedCounter: Counter {
    fn write_slice<W: Write>(writer: &mut W, counters: &[Self]) -> io::Result<()>;
    fn read_vec<R: Read>(reader: &mut R, len: usize) -> io::Result<Vec<Self>>;
}

impl Counter for u8 {
    #[inline(always)]
    fn max_value() -> Self {
        u8::MAX
    }

    #[inline(always)]
    fn saturating_increment(&mut self) {
        *self = self.saturating_add(1);
    }

    #[inline(always)]
    fn to_u64(self) -> u64 {
        u64::from(self)
    }
}

impl SerializedCounter for u8 {
    fn write_slice<W: Write>(writer: &mut W, counters: &[Self]) -> io::Result<()> {
        writer.write_all(counters)
    }

    fn read_vec<R: Read>(reader: &mut R, len: usize) -> io::Result<Vec<Self>> {
        let mut counters = vec![0_u8; len];
        reader.read_exact(&mut counters)?;
        Ok(counters)
    }
}

impl Counter for u16 {
    #[inline(always)]
    fn max_value() -> Self {
        u16::MAX
    }

    #[inline(always)]
    fn saturating_increment(&mut self) {
        *self = self.saturating_add(1);
    }

    #[inline(always)]
    fn to_u64(self) -> u64 {
        u64::from(self)
    }
}

impl SerializedCounter for u16 {
    fn write_slice<W: Write>(writer: &mut W, counters: &[Self]) -> io::Result<()> {
        let mut buffer = Vec::with_capacity(16 * 1024);
        for chunk in counters.chunks(8 * 1024) {
            buffer.clear();
            for &counter in chunk {
                buffer.extend_from_slice(&counter.to_le_bytes());
            }
            writer.write_all(&buffer)?;
        }
        Ok(())
    }

    fn read_vec<R: Read>(reader: &mut R, len: usize) -> io::Result<Vec<Self>> {
        let mut counters = Vec::with_capacity(len);
        let mut buffer = vec![0_u8; 16 * 1024];
        let mut remaining = len;
        while remaining > 0 {
            let chunk_len = remaining.min(8 * 1024);
            let byte_len = chunk_len * 2;
            reader.read_exact(&mut buffer[..byte_len])?;
            for bytes in buffer[..byte_len].chunks_exact(2) {
                counters.push(u16::from_le_bytes([bytes[0], bytes[1]]));
            }
            remaining -= chunk_len;
        }
        Ok(counters)
    }
}

impl Counter for u32 {
    #[inline(always)]
    fn max_value() -> Self {
        u32::MAX
    }

    #[inline(always)]
    fn saturating_increment(&mut self) {
        *self = self.saturating_add(1);
    }

    #[inline(always)]
    fn to_u64(self) -> u64 {
        u64::from(self)
    }
}

impl SerializedCounter for u32 {
    fn write_slice<W: Write>(writer: &mut W, counters: &[Self]) -> io::Result<()> {
        let mut buffer = Vec::with_capacity(16 * 1024);
        for chunk in counters.chunks(4 * 1024) {
            buffer.clear();
            for &counter in chunk {
                buffer.extend_from_slice(&counter.to_le_bytes());
            }
            writer.write_all(&buffer)?;
        }
        Ok(())
    }

    fn read_vec<R: Read>(reader: &mut R, len: usize) -> io::Result<Vec<Self>> {
        let mut counters = Vec::with_capacity(len);
        let mut buffer = vec![0_u8; 16 * 1024];
        let mut remaining = len;
        while remaining > 0 {
            let chunk_len = remaining.min(4 * 1024);
            let byte_len = chunk_len * 4;
            reader.read_exact(&mut buffer[..byte_len])?;
            for bytes in buffer[..byte_len].chunks_exact(4) {
                counters.push(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]));
            }
            remaining -= chunk_len;
        }
        Ok(counters)
    }
}
