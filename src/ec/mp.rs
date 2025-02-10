
#[inline(always)]
fn is_digit_nonzero_ct(x: u64) -> u32 {
    const RADIX: u32 = 64; // TODO
    let result = (x | (0u64.wrapping_sub(x))) >> (RADIX - 1);
    result as u32
}

#[inline(always)]
fn is_digit_zero_ct(x: u64) -> u32 {
    1 ^ is_digit_nonzero_ct(x)
}

#[inline(always)]
fn is_digit_lessthan_ct(x: u64, y: u64) -> u32 {
    const RADIX: u32 = 64; // TODO
    let result = (x ^ ((x ^ y) | ((x.wrapping_sub(y)) ^ y))) >> (RADIX - 1);
    result as u32
}

fn subc(
    minuend: u64,
    subtrahend: u64,
    borrow_in: u32,
) -> (u64, u32) {
    let temp_reg = minuend.wrapping_sub(subtrahend);
    let borrow_reg = is_digit_lessthan_ct(minuend, subtrahend) |
                     ((borrow_in & is_digit_zero_ct(temp_reg)) as u32);
    let difference_out = temp_reg.wrapping_sub(borrow_in as u64);
    (difference_out, borrow_reg)
}

pub(crate) fn mp_sub(c: &mut [u64], a: &[u64], b: &[u64], nwords: usize) {
    let mut borrow = 0;

    for i in 0..nwords {
        let (difference, new_borrow) = subc(a[i], b[i], borrow);
        c[i] = difference;
        borrow = new_borrow;
    }
}

pub(crate) fn select_ct(a: u64, b:u64, mask: u64) -> u64 {
    ((a ^ b) & mask) ^ a
}

pub(crate) fn select_ct_arr(c: &mut [u64], a: &[u64], b: &[u64], mask: u64, nwords: usize) {
    for i in 0..nwords {
        // c[i] = ((a[i] ^ b[i]) & mask) ^ a[i];
        c[i] = select_ct(a[i], b[i], mask);
    }
}

pub(crate) fn swap_ct(a: &mut [u64], b: &mut [u64], option: u64, nwords: usize) {
    let mut temp: u64;

    for i in 0..nwords {
        temp = option & (a[i] ^ b[i]);
        a[i] = temp ^ a[i];
        b[i] = temp ^ b[i];
    }
}

fn shiftr(high_in: u64, low_in: u64, shift: u32, digit_size: u32) -> u64 {
    (low_in >> shift) ^ (high_in << (digit_size - shift))
}

pub(crate) fn mp_shiftr(x: &mut [u64], shift: u32, nwords: usize) -> u64 {
    let bit_out = x[0] & 1;
    const RADIX: u32 = 64; // TODO 

    for i in 0..nwords - 1 {
        x[i] = shiftr(x[i + 1], x[i], shift, RADIX);
    }

    x[nwords - 1] >>= shift;

    bit_out
}
