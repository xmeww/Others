function USBTTL(cmd, varargin)
% USBTTL(cmd, arg);  -- Set RTS line of a serial port low or high.
%
% USBTTL('Open' [,devName]);
% - Open a connection to a serial port, either the one specified in varagin
% or it finds one by itself
%
% USBTTL('Close');
% - Close connection to serial port.
%
% USBTTL('Set', level);
% - Set TTL level of RTS line of serial port to 'level', ie., low or high,
% for setting 'level' == 0 or 1.
% The function takes about 1 msec on a FTDI FT-232BM Serial-over-USB port
% to excecute. Operation on a native serial port may be faster.
%

% History:
% 2.10.2009  mk  Written.

persistent p;
persistent oldverbo;

if nargin < 1
    error('You must provide a "cmd" command argument.');
end

if strcmpi(cmd, 'Set')
    if isempty(p)
        error('Tried to set TTL trigger level, but connection not open!');
    end
    
    if isempty(varargin)
        error('New TTL line setting 0 or 1 missing!');
    end
    
    IOPort('ConfigureSerialPort', p, sprintf('RTS=%i', 1 - varargin{1}));
    
    return;
end

if strcmpi(cmd, 'Open')
    if ~isempty(p)
        error('Tried to open TTL trigger connection, but already open!');
    end

    if ~isempty(varargin)
        p = IOPort('OpenSerialport', FindSerialPort(char(varargin{1}), 1));
    else
        p = IOPort('OpenSerialport', FindSerialPort ([],1));
    end
    
    oldverbo = IOPort('Verbosity', 3);
    
    return;
end

if strcmpi(cmd, 'Close')
    if isempty(p)
        error('Tried to close TTL trigger connection, but none opened!');
    end
    
    IOPort('Verbosity', oldverbo);
    IOPort('Close', p);
    p = [];
    
    return;
end

error('Unknown command provided!');
