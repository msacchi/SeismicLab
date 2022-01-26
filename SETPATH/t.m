

 X = digitTestCellArrayData;

autoenc = trainAutoencoder(X,4,'MaxEpochs',400, 'DecoderTransferFunction','purelin');

Xr= predict(autoenc,X);

images([X{10},Xr{10}])
