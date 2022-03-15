//! # Static Fast Fourier Transform
//!
//! Implementation developed at [my fourier notebook](https://nbviewer.org/github/unic0rn9k/fourier-notebook/blob/master/README.ipynb).
//!
//! ```rust
//! use sfft::*;
//!
//! let v = fft(&[re(0f32), re(1.), re(2.), re(3.)]);
//!
//! assert_eq!(
//!     v,
//!     [
//!         re(6f32) + im(0.),
//!         re(-1.9999999) + im(2.),
//!         re(-2.) + im(0.),
//!         re(-2.) - im(2.),
//!     ]
//! )
//! ```

pub use num::*;
use std::f32::consts::PI;

pub unsafe fn unsafe_fft<const LEN: usize>(
    x: *const Complex<f32>,
    o: *mut Complex<f32>,
    n: usize,
    ofset: usize,
) -> *const Complex<f32> {
    if n < 2 {
        return x;
    }

    let even = unsafe_fft::<LEN>(x, o.add(n / 2), n / 2, ofset * 2);
    let odd = unsafe_fft::<LEN>(x.add(ofset), o, n / 2, ofset * 2);

    for j in 0..n / 2 {
        let even = *even.add(j);
        let odd = *odd.add(j);

        *o.add(j) = im(-2. * PI * j as f32 / n as f32).exp_();
        *o.add(j + n / 2) = even - *o.add(j) * odd;
        *o.add(j) = even + *o.add(j) * odd;
    }

    o
}

pub fn fft<const LEN: usize>(x: &[Complex<f32>; LEN]) -> [Complex<f32>; LEN] {
    assert_eq!(LEN & (LEN - 1), 0);
    let mut ret = *x;
    unsafe { unsafe_fft::<LEN>(x.as_ptr(), ret.as_mut_ptr(), LEN, 1) };
    ret
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        use crate::*;

        let a = [re(0f32), re(1.), re(2.), re(3.)];
        let b = fft(&a);

        assert_eq!(
            b,
            [
                re(6f32) + im(0.),
                re(-1.9999999) + im(2.),
                re(-2.) + im(0.),
                re(-2.) - im(2.),
            ]
        )
    }
}
