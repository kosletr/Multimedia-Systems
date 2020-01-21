load Testing/KyriaFile.mat
JPEGencKyr = JPEGenc;
load Testing/img.mat
JPEGencLet = JPEGenc;
clear JPEGenc
Tables;

JPEGencKyr{1}.ACC = JPEGencLet{1}.ACC;
JPEGencKyr{1}.DCC = JPEGencLet{1}.DCC;
JPEGencKyr{1}.ACL = JPEGencLet{1}.ACL;
JPEGencKyr{1}.DCL = JPEGencLet{1}.DCL;

JPEGenc = JPEGencKyr;