function c = barker_code(K)
%BARKER_CODE  Binary Barker sequence of length K (values +1/-1).
%   c = barker_code(K), K in {2, 3, 4, 5, 7, 11, 13}. For K = 2 and 4 the
%   first of the two known variants is returned. The aperiodic
%   autocorrelation has peak K and sidelobes of magnitude <= 1.
switch K
    case 2,  c = [1 -1];
    case 3,  c = [1 1 -1];
    case 4,  c = [1 1 -1 1];
    case 5,  c = [1 1 1 -1 1];
    case 7,  c = [1 1 1 -1 -1 1 -1];
    case 11, c = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
    case 13, c = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
    otherwise
        error('barker_code:length', 'No Barker code of length %d exists.', K);
end
end
