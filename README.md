“A-EOG_raw_data_example.bin” is a 16MByte file with raw data 
(order of magnitude: 5 minutes of acquisition time on the PC).

Data is recorded by words of 32bits for each “channel”. 32bits are split into 4 bytes, and the most significant byte is the last one (reverse order) in each 32-bit word.
Most significant byte of the 32 bits word contains the ADC channel number (0 to 13), except for channels 6 and 7 (see below)
Remaining 3 bytes (24 bits) of the 32 bits word contain the ADC signed conversion result (negative values starting with “1” as most significant bit)
24 bits to be converted to decimal value (as these are signed, if you use a conversion instruction on 32 bits you might have to extend the MSB before conversion)
Voltage corresponding to the ADC value is given by: decimal value * 5.00 / (2^(24-1)-1)       (5.00 being the ADC reference voltage)
ADC channel “6” data is not used, and the corresponding 32 bits are replaced with time stamps in [ms] (0 at the A-EOG power-up...) also in bytes reverse order (most significant byte is the last one)
ADC channel “7” data is not used, and the corresponding 32 bits are replaced with a fixed pattern 0xAAAAAAAA (hexadecimal) (in reverse order too, but as this value is symmetric this is not visible)
There is 1 sample per channel per millisecond. In case there would be short missing pieces of data this should be handled... (e.g. some channels data missing, or data missing for several milliseconds)
“A-EOG_acquisition.c” is the program to operate the FTDI serial to USB cable for data recoding on the PC. It asks the user a filename, and asks to press “a” to start recording and to press “z” to stop acquisition. Maybe this can be integrated in the python program to be developed.
