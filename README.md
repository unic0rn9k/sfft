## Static Fast Fourier Transform

Implementation developed at [my fourier notebook](https://nbviewer.org/github/unic0rn9k/fourier-notebook/blob/master/README.ipynb).

```rust
use sfft::*;

let v = fft([re(0f32), re(1.), re(2.), re(3.)]);

assert_eq!(
    v,
    [
        re(6f32) + im(0.),
        re(-1.9999999) + im(2.),
        re(-2.) + im(0.),
        re(-2.) - im(2.),
    ]
)
```
