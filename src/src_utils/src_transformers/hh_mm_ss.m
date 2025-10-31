function [ret] = hh_mm_ss(q)
	hms = floor(q);
	s = mod(hms,60);
	hm = floor(hms/60);
	m = mod(hm,60);
	h = floor(hm/60);
	ret = sprintf('%02d:%02d:%02d',h,m,s);
end