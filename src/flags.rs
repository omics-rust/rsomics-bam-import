// SAMv1 §1.4.2 FLAG bits
pub(crate) const FLAG_PAIRED: u16 = 0x1;
pub(crate) const FLAG_FUNMAP: u16 = 0x4;
pub(crate) const FLAG_MUNMAP: u16 = 0x8;
pub(crate) const FLAG_READ1: u16 = 0x40;
pub(crate) const FLAG_READ2: u16 = 0x80;

pub(crate) const SE_FLAGS: u16 = FLAG_FUNMAP;
pub(crate) const PE_R1_FLAGS: u16 = FLAG_PAIRED | FLAG_FUNMAP | FLAG_MUNMAP | FLAG_READ1;
pub(crate) const PE_R2_FLAGS: u16 = FLAG_PAIRED | FLAG_FUNMAP | FLAG_MUNMAP | FLAG_READ2;
