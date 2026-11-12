function Uprop = propagationField(X, gridMask, Fs, ampQuantile)
% propagationField   Per-electrode, sign-corrected traveling-wave propagation vectors.
%
% PURPOSE
% -------
% Single source of truth for the per-frame PROPAGATION direction field used by the
% wave visualizers (saccadeWaveMovie, waveGrandSnapshot). Given an analytic signal
% X (e.g. generalized_phase output) it returns the complex propagation vector at
% every electrode and time, using the SAME convention as waveGP:
%   U = conj(exp(i*phaseGradient))           (the "y-flip"/conj so angle is propagation)
%   per-frame sign from instantaneous frequency (so it is direction of travel, not a
%   180-deg-ambiguous orientation, and NOT the raw phase gradient which is ~180 off)
%
% Draw on a YDir='reverse' array axis as quiver(X,Y, real(Uprop), -imag(Uprop)) so
% arrows point where the wave actually travels (see waveFieldSnapshot). This replaces
% the older raw gradient(Ph) overlays, which were unsigned and vertically inconsistent.
%
% INPUTS
% ------
% X          nRow x nCol x T complex analytic signal (holes already filled).
% gridMask   nRow x nCol logical (true = real electrode).
% Fs         sample rate (Hz), default 1000.
% ampQuantile amplitude-mask quantile for the per-frame sign vote, default 0.20.
%
% OUTPUT
% ------
% Uprop      nRow x nCol x T complex, sign-corrected. |Uprop| = 1 on real electrodes.
%
% See also: waveGP, waveFieldSnapshot, saccadeWaveMovie, waveGrandSnapshot

if nargin < 3 || isempty(Fs), Fs = 1000; end
if nargin < 4 || isempty(ampQuantile), ampQuantile = 0.20; end
gm = logical(gridMask); T = size(X,3);

[~, pd] = phase_gradient_complex_multiplication(X, 1);
U = conj(exp(1i*pd));                                  % propagation-axis unit vectors
IF = instantaneous_frequency(X, Fs); omega = 2*pi*IF;

A = abs(X);
thr = quantile(A(repmat(gm,1,1,T)), ampQuantile);
Wt = A; Wt(A < thr) = 0; Wt = Wt .* gm;               % amplitude weight, corners zero
s_t = sign(squeeze(sum(sum(omega.*Wt,1,'omitnan'),2,'omitnan')));  % T x 1 per-frame sign
s_t(s_t==0) = 1;
Uprop = U .* reshape(s_t, 1, 1, []);
end
