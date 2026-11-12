function checkPitchDefault(fname, args, gridMask)
% checkPitchDefault   Guard against Biru's 0.7 mm pitch leaking onto another array.
%
% waveGP / waveSourcePoints / waveCriticalPoints all default 'pixelSpacing_mm' to
% 0.7 -- correct for Biru's ECoG disc, WRONG for the 0.800 mm Atlas EC32. Getting
% it wrong produces no error: every reported speed is simply ~12% low, and the
% direction/coherence results look completely normal. That is the worst kind of
% bug, so this makes it loud.
%
% Fires only when BOTH are true: the caller did not pass the pitch explicitly,
% AND the grid is not Biru's 8x10. Biru call sites are untouched.
%
% BOTH spellings count as explicit. The engine's own functions take
% 'pixelSpacing_mm'; waveMethod_explainDetection and waveEventWavelength take
% 'PixelSpacingMm'. Recognising only the first made this warning fire on calls that
% had passed 0.8 correctly -- and a guard that cries wolf is a guard that gets
% ignored on the day it is right.
%
% INPUTS
%   fname     calling function name, for the message and the warning identifier
%   args      the caller's varargin
%   gridMask  nRow x nCol logical
%
% See also: getArraySpec, accumulateSaccadeCube (returns .pitch_mm), waveGP

BIRU_GRID = [8 10];

if isequal(size(gridMask), BIRU_GRID)
    return;                                   % Biru: the 0.7 default is correct
end

if numel(args) >= 2                           % caller passed it explicitly -> trust them
    names  = args(1:2:end-1);
    isName = cellfun(@(x) ischar(x) || isstring(x), names);
    if any(ismember(lower(string(names(isName))), ["pixelspacing_mm" "pixelspacingmm"]))
        return;
    end
end

% warning()'s first argument is the message ID and may not itself contain format
% specifiers, so build it separately rather than sprintf-ing it inline.
msgID = [fname ':defaultPitch'];
warning(msgID, ...
    ['%s is using the DEFAULT pixelSpacing_mm = 0.7 mm (Biru) on a %dx%d grid, which ' ...
     'is not Biru''s 8x10. Reported SPEEDS will be wrong (directions and coherence ' ...
     'will not). Pass the array''s pitch explicitly -- accumulateSaccadeCube and ' ...
     'accumulateEventCube return it as C.pitch_mm, and getArraySpec(allData).pitch_mm ' ...
     'gives it directly.'], ...
    fname, size(gridMask,1), size(gridMask,2));
end
