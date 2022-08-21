% -------------------------------------------------------------------------
% Copyright (c) 2013 DDE Lab, Binghamton University, NY.
% All Rights Reserved.
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The program is supplied "as is," without any
% accompanying services from DDE Lab. DDE Lab does not warrant the
% operation of the program will be uninterrupted or error-free. The
% end-user understands that the program was developed for research purposes
% and is advised not to rely exclusively on the program for any reason. In
% no event shall Binghamton University or DDE Lab be liable to any party
% for direct, indirect, special, incidental, or consequential damages,
% including lost profits, arising out of the use of this software. DDE Lab
% disclaims any warranties, and has no obligations to provide maintenance,
% support, updates, enhancements or modifications.
% -------------------------------------------------------------------------
% Author: Vojtech Holub
% -------------------------------------------------------------------------
% Contact: vojtech_holub@yahoo.com
%          fridrich@binghamton.edu
%          February 2013, http://dde.binghamton.edu
% -------------------------------------------------------------------------
% [distortion] = SI_UNIWARD(precover, cover_path, stego_path, payload, quality_factor) embeds a payload using simulated optimal
% coding for minimizing 'SI-UNIWARD' distortion - it uses default settings
%
% [distortion] = SI_UNIWARD(precover, cover_path, stego_path, payload, quality_factor, config) - same as above but using
% non-default options in 'config' structure.
%
%--------------------------------------------------------------------------
% Input format:
%   precover            - unit8                         - matrix with precover values (pixels)
%   cover_path          - string                        - string containing path to a cover image
%   stego_path          - string                        - string containing path to a stego image
%   payload             - single                        - relative payload (in bits per non zero AC coefficient)
%   quality_factor      - int32                         - Quality of JPEG compression
%
%   config.STC_h        - uint32        - default 0  	- 0 for optimal simulator, otherwise STC submatrix height (try 7-12)
%	config.seed         - int32         - default 0 	- random seed
%
% Output format:
%   distortion  - single    - value containing embedding distortion
%--------------------------------------------------------------------------
