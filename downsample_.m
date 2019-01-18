function [ downsampled_signal ] = downsample_( signal, factor)

%	Downsample signal by averaging neighboring values.
% 
% 	@type  signal: array_like
% 	@param signal: one-dimensional signal to be downsampled
% 
% 	@type  factor: int
% 	@param factor: this many neighboring values are averaged
% 
% 	@rtype: ndarray
% 	@return: downsampled signal
	if factor < 2
		downsampled_signal = signal;
    else
        downsampled_signal=conv(double(signal), ones(1,factor), 'same');
        downsampled_signal=downsampled_signal(1:factor:length(downsampled_signal));
    end


end

