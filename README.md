# neuroclear-marmoset-jne

Code and data for the marmoset figures (Figures 7 and 8, Supplementary Videos 1 and 2) in
our JNE paper on the Neuro-Clear array. Everything here comes from one marmoset with a
64-contact array over MT and nearby cortex, recorded between 10 days and 9 months after
implantation.

Figure 7 is the evoked response over the first four months. Figure 8 is traveling waves
from June to August 2026.

## Running it

You need MATLAB (we used R2025a) with the Signal Processing, Statistics and Machine
Learning, and Image Processing toolboxes.

```
git clone --recursive https://github.com/bmitchell1092/neuroclear-marmoset-jne.git
```

Download the data from Zenodo (DOI to come) and unzip it so that `data/fig7` and
`data/fig8` sit at the top of the repo. Then, in MATLAB, from the repo folder:

```matlab
make_all                  % figures go to output/
make_all('Videos', true)  % also the two supplementary videos
verify                    % recompute the numbers and compare to the paper
```

`make_all` takes a few minutes, most of it detecting waves in Figure 8.

## Data

All signals are 1 kHz LFP in µV.

- `fig7/flash_<date>.mat` has single trials from full-field flash sessions: 0.5–40 Hz,
  −300 to 400 ms around the flash, baseline-subtracted. It also has the channels that
  failed QC that day.
- `fig7/darkRest_psd_<date>.mat` has resting power spectra, one per contact.
- `fig8/flash_<date>.mat` has flash trials from June 2026 on, laid out on the 8 × 10 grid
  as `cube(trial, row, col, time)`.
- `fig8/lightRest_<date>.mat` has 500 ms windows of eyes-open rest, at least 500 ms from
  any saccade, in the same layout.
- `fig8/lightRest_exemplar.mat` is the spontaneous window shown in Figure 8a, with 500 ms
  extra on each side.

A note on orientation: the figures draw the array the way the channel map supplied with
it does (row 1 = channels 57, 63, 1, 7). The grid in the data files is that map upside
down (row 1 = 53, 33, 31, 11), which is the orientation our lab pipeline uses. The figure
code flips it for display (`jneDisplayMap`), and nothing else depends on it.

Contacts that failed QC are NaN, as are the 16 corners of the grid, which have no contact.
Every file has the channel map, the electrode pitch (0.7 mm) and a note on how it was made.
`expected.mat` holds the numbers from the published figures, which `verify` checks against.

The data starts from the cut trials. Filtering, flash timing (corrected with a photodiode),
eye tracking and channel QC were done upstream in our lab pipeline. Raw recordings are
available on request.

## Code

- `code/figures` draws the panels.
- `code/data` builds what each figure needs from `data/`.
- `code/wave` is our wave detection: source-point detection after Muller lab methods, and
  speed from frequency and wavelength.

It uses the Muller lab's [wave-matlab](https://github.com/mullerlab/wave-matlab) and
[generalized-phase](https://github.com/mullerlab/generalized-phase) (submodules in
`external/`). We also copied in `smoothn` (Damien Garcia) and two CircStat functions
(Philipp Berens); their licenses are in their folders.

## Contact

Blake Mitchell, Salk Institute
