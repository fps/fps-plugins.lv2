# fps-plugins.lv2

A collection of LV2 plugins that probably noone besides me is interested in.

# Plugins

## relative_dynamics

A dynamics processor using the ratio of two envelope followers with different time scales.

## eq_match

An FFT-based EQ matching processor

## stereo_decorrelation

A stereo processor that generates an exponentially decaying white noise FIR and then uses it to process the left channel with that FIR and the right channel with the negative of that FIR. This decorrelates the left and right channels to a certain degree and acts to widen the stereo field. If fed with a mono signal on both inputs it has the nice property of summing to this precise mono signal again.

### Notes:

- This plugin interrupts the audio processing on parameter changes.

# Author

Florian Paul Schmidt (mista.tapas@gmx.net)

## Thanks

Fons Adriaensen for his invaluable tips on DSP
