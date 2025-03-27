“A-EOG_raw_data_example” is a short file with raw data 
(order of magnitude: couple of seconds of acquisition time, 
acquisition recorded on the computer after ).

Data is recorded by words of 32bits for each “channel”
First 8 bits (of the 32 bits word) contain the ADC channel number (0 to 13), except for channels 6 and 7 (see below)
Last 24 bits (of the 32 bits word) contain the ADC signed conversion result (negative values starting with “1” as most significant bit)
24 bits to be converted to decimal value (as these are signed, if you use a conversion instruction on 32 bits you might have to extend the MSB before conversion)
Voltage corresponding to the ADC value is given by: decimal value * 5.00 / (2^(24-1)-1)       (5.00 being the ADC reference voltage)
ADC channel “6” data is not used, and the corresponding 32 bits are replaced with time stamps in [ms] (0 at the A-EOG power-up...)
ADC channel “7” data is not used, and the corresponding 32 bits are replaced with a fixed pattern 0xAAAAAAAA (hexadecimal)
There are 4 samples per channel per millisecond (which is oversampling) –> averaging 4 samples is necessary to reduce data by 4, thereby reducing measurement noise, and getting 1 resulting sample per channel per millisecond. In case there would be short missing pieces of data this should be handled... (e.g. if only 3 samples with the same timestamp, or data missing for few milliseconds (I didn’t observe this, but who knows...))
“main.c” is the program to operate the FTDI serial to USB cable for data recoding on the PC. It asks the user a filename, and asks to press “a” to start recording and to press “z” to stop acquisition. Maybe this can be integrated in the program to be developed as well.
 
Could we imagine meeting next week one day for working on this, defining the data inputs, the data display format, ...? From Monday to Wednesday I will have to stay in the Lausanne area, but Thursday and Friday I can move either to Sion or Lausanne. Please let me know your availability.
 
Please let me know should you have further questions or needs in the meantime.
