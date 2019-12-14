function JPEGencStream = JPEGencodeStream(img, subimg, qScale)

JPEGenc = JPEGencode(img, subimg, qScale);

JPEGencStream = [];

for count = 2:length(JPEGenc)
    JPEGencStream = [JPEGencStream , JPEGenc{count}.huffStream];
end

end