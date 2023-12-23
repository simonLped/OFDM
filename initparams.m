function [simin, nbsecs, fs] = initparams(toplay, fs)

    
    % nbsecs = 40; %number of seconds of the playback recording
    toplay = -1 + 2.*(toplay - min(toplay))./(max(toplay) - min(toplay));
    % toplay = rescale(toplay, -1, 1);    %scaling

    simin = [zeros(fs*2, 2)

             toplay', zeros(length(toplay),1);

             zeros(fs*1, 2)
             ];
    nbsecs = ceil(length(simin) / fs) + 1;

end

% 
% function [simin, nbsecs, fs] = initparams(pulse, toplay, channel, fs)
%         % Normalise the toplay and pulse signals
%         toplay = -1 + 2.*(toplay - min(toplay))./(max(toplay) - min(toplay));
% 
%         pulse = -1 + 2.*(pulse - min(pulse))./(max(pulse) - min(pulse));
% 
%         % Construct simin
%         simin = [zeros(fs*2, 2)
%                 pulse', zeros(length(pulse), 1);
%                 zeros(length(channel), 1), zeros(length(channel), 1);
%                 toplay', zeros(length(toplay), 1);
%                 zeros(fs, 2) ];
%         % Calculate nbsecs
%         nbsecs = ceil(length(simin) / fs) + 1;     
% end