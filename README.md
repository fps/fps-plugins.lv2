# fps-plugins.lv2

A collection of LV2 plugins that probably noone besides me is interested in.

# Plugins

## relative_dynamics

A dynamics processor using the ratio of two envelope followers with different time scales.

## eq_match

An FFT-based EQ matching processor.

### How-To

1. Enable the "Analyze spectrum 1" button and input some source material which you want to have matched to a second source. Disable the button when done.

2. Enable the "Analyze spectrum 2" button. Now input some material which has the spectrum you want to match to. Disable the button when done.

3. Enable the "Apply match 1->2" button and enjoy your matched spectrum.

### Notes

- The matching filters are always recalculated when any "Analyze..." button is deactivated.

- This plugin is not meant to be automated. Due to the offloading of the response calculations to a worker thread there is a variable delay until a new match takes effect.

## stereo_decorrelation

A stereo processor that generates an exponentially decaying white noise FIR and then uses it to process the left channel with that FIR and the right channel with the negative of that FIR. This decorrelates the left and right channels to a certain degree and acts to widen the stereo field. If fed with a mono signal on both inputs it has the nice property of summing to this precise mono signal again.

### Notes

- This plugin interrupts the audio processing on parameter changes.

### Examples

https://github.com/fps/fps-plugins.lv2/raw/master/stereo_decorrelation_example.ogg

# Author

Florian Paul Schmidt (mista.tapas@gmx.net)

## Thanks

Fons Adriaensen for his invaluable tips on DSP
