function days = hms2days(varargin)
%HMS2DAYS Convert hours, minutes, and seconds to days.
%
%   DAYS = HMS2DAYS(HOUR, MINUTE, SECOND) converts the number of hours,
%   minutes, and seconds to a number of days.
%
%   The following holds (to within rounding precision):
%
%     DAYS = HOUR / 24 + MINUTE / (24 * 60) + SECOND / (24 * 60 * 60)
%          = (HOUR + (MINUTE + SECOND / 60) / 60) / 24

%   Author:      Peter J. Acklam
%   Time-stamp:  2004-09-22 08:45:33 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   error(nargchk(1, 3, nargsin));
   argv = {0 0 0};
   argv(1:nargsin) = varargin;
   [hour, minute, second] = deal(argv{:});

   days = (hour + (minute + second / 60) / 60) / 24;
