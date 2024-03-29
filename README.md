![Mandelbrot Zoom](./thumbnail.png)
*Mandelbrot set zoom to scale 0.00000000001 targetting   
(-1.99977406013629035931268, -0.00000000329004032147943)*

# Fractals

Fractals:
- Mandelbrot set

Optimization: 
- [None]

Features: 
- Encode to video
- Encode to single image

### Could be useful

```
cargo run
cargo run --release
ffmpeg -start_number 0 -framerate 24 -i output/%09d.png -vcodec h264 -frames:v 240 -pix_fmt yuv420p -crf 17 -preset veryfast -y out_concat.mp4
```