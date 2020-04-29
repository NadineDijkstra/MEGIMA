function evenResult = roundEven(input)

% evenResult = roundEven(input)
%
% This function rounds the input floating point number to the nearest even
% integer. If the input number is integer, it will be converted to a
% floating point number.

if round(input) == input, input = input + 0.000000001; end

if ~mod(ceil(input),2)
    evenResult = ceil(input);
else
    evenResult = floor(input);
end

end