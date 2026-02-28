
function pairs = Freq_Comb(T,w,lin_max)
% Freq_Comb  Generate integer frequency combinations within a bandwidth.
%
%   Constructs all integer pairs [j1, j2] such that:
%     0 <= j1*w(1) + j2*w(2) <= T   (bandwidth constraint)
%     |j2| <= lin_max               (order constraint on second frequency)
%   Rows are sorted by total integer order sum(|j|).
%
%   INPUTS:
%     T       - frequency bandwidth [rad/TU]
%     w       - [2x1] fundamental frequency vector [w0; wP]
%     lin_max - maximum allowed order of the second frequency
%
%   OUTPUT:
%     pairs   - [K x 2] integer combination matrix

w1 = w(1); w2 = w(2);
max_a = floor(T / w1);
max_b = floor(T / w2);

[A, B] = meshgrid(-max_a:1:max_a, -max_b:1:max_b);
comb   = [A(:), B(:)];
freq   = comb(:,1)*w1 + comb(:,2)*w2;

indices = (freq >= 0) & (freq <= T) & (abs(comb(:,2)) <= lin_max);
pairs   = comb(indices,:);

[~,sort_id] = sort(sum(abs(pairs), 2));
pairs = pairs(sort_id,:);
end


