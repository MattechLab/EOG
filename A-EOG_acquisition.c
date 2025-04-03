#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <conio.h>  // For _kbhit() and _getch() to capture keypresses
#include "ftd2xx.h" // Ensure this header is included for FTDI functions

#define BAUD_RATE 1000000   // 1000 kbits/s (1000000 bps)
#define BUFFER_SIZE 1024    // Buffer size for reading data

// Function declarations using explicit function pointer syntax in C
typedef FT_STATUS(WINAPI *FT_Open_func)(DWORD, FT_HANDLE*);
typedef FT_STATUS(WINAPI *FT_SetBaudRate_func)(FT_HANDLE, ULONG);
typedef FT_STATUS(WINAPI *FT_SetDataCharacteristics_func)(FT_HANDLE, UCHAR, UCHAR, UCHAR);
typedef FT_STATUS(WINAPI *FT_SetFlowControl_func)(FT_HANDLE, ULONG, UCHAR, UCHAR);
typedef FT_STATUS(WINAPI *FT_Purge_func)(FT_HANDLE, ULONG);
typedef FT_STATUS(WINAPI *FT_Read_func)(FT_HANDLE, LPVOID, DWORD, LPDWORD);
typedef FT_STATUS(WINAPI *FT_Close_func)(FT_HANDLE);

// Function to prompt the user for a file name and storage location
void getFileName(char *fileName, size_t size) {
    printf("Enter the file name to store the data: ");
    fgets(fileName, size, stdin);
    fileName[strcspn(fileName, "\n")] = '\0';  // Remove trailing newline character
}

// Function to load the FTDI DLL and get function pointers
int loadFtd2xxLibrary(HMODULE *hFtd2xx, FT_Open_func *FT_Open, FT_SetBaudRate_func *FT_SetBaudRate,
                       FT_SetDataCharacteristics_func *FT_SetDataCharacteristics, FT_SetFlowControl_func *FT_SetFlowControl,
                       FT_Purge_func *FT_Purge, FT_Read_func *FT_Read, FT_Close_func *FT_Close) {
    // Load the DLL dynamically
    *hFtd2xx = LoadLibrary("ftd2xx.dll");
    if (!(*hFtd2xx)) {
        fprintf(stderr, "Failed to load ftd2xx.dll\n");
        return 0;
    }

    // Get the function pointers
    *FT_Open = (FT_Open_func)GetProcAddress(*hFtd2xx, "FT_Open");
    *FT_SetBaudRate = (FT_SetBaudRate_func)GetProcAddress(*hFtd2xx, "FT_SetBaudRate");
    *FT_SetDataCharacteristics = (FT_SetDataCharacteristics_func)GetProcAddress(*hFtd2xx, "FT_SetDataCharacteristics");
    *FT_SetFlowControl = (FT_SetFlowControl_func)GetProcAddress(*hFtd2xx, "FT_SetFlowControl");
    *FT_Purge = (FT_Purge_func)GetProcAddress(*hFtd2xx, "FT_Purge");
    *FT_Read = (FT_Read_func)GetProcAddress(*hFtd2xx, "FT_Read");
    *FT_Close = (FT_Close_func)GetProcAddress(*hFtd2xx, "FT_Close");

    // Check if all functions were found
    if (!(*FT_Open) || !(*FT_SetBaudRate) || !(*FT_SetDataCharacteristics) || !(*FT_SetFlowControl) ||
        !(*FT_Purge) || !(*FT_Read) || !(*FT_Close)) {
        fprintf(stderr, "Failed to get function addresses\n");
        FreeLibrary(*hFtd2xx);
        return 0;
    }

    return 1;
}

// Function to setup the FTDI device
int setupFTDI(FT_HANDLE *handle, FT_Open_func FT_Open, FT_SetBaudRate_func FT_SetBaudRate,
              FT_SetDataCharacteristics_func FT_SetDataCharacteristics, FT_SetFlowControl_func FT_SetFlowControl,
              FT_Purge_func FT_Purge) {
    FT_STATUS status;

    // Open the FTDI device (0 is the first device, change if needed)
    status = FT_Open(0, handle);
    if (status != FT_OK) {
        fprintf(stderr, "Failed to open FTDI device.\n");
        return 0;
    }

    // Set the baud rate (3090 kbits/s)
    status = FT_SetBaudRate(*handle, BAUD_RATE);
    if (status != FT_OK) {
        fprintf(stderr, "Failed to set baud rate.\n");
//        FT_Close(*handle);
        return 0;
    }

    // Set data characteristics (8 data bits, 1 stop bit, no parity)
    status = FT_SetDataCharacteristics(*handle, FT_BITS_8, FT_STOP_BITS_1, FT_PARITY_NONE);
    if (status != FT_OK) {
        fprintf(stderr, "Failed to set data characteristics.\n");
//        FT_Close(*handle);
        return 0;
    }

    // Disable flow control
    status = FT_SetFlowControl(*handle, FT_FLOW_NONE, 0, 0);
    if (status != FT_OK) {
        fprintf(stderr, "Failed to set flow control.\n");
//        FT_Close(*handle);
        return 0;
    }

    // Purge any existing data in the buffer
    status = FT_Purge(*handle, FT_PURGE_RX | FT_PURGE_TX);
    if (status != FT_OK) {
        fprintf(stderr, "Failed to purge buffers.\n");
//        FT_Close(*handle);
        return 0;
    }

    return 1;
}

// Function to read and save data to a file
void readAndSaveData(FT_HANDLE handle, FT_Read_func FT_Read, const char *fileName) {
    FT_STATUS status;
    unsigned char buffer[BUFFER_SIZE];
    DWORD bytesRead = 0;

    FILE *outFile = fopen(fileName, "wb");
    if (!outFile) {
        fprintf(stderr, "Failed to open the file for writing.\n");
        return;
    }

    printf("Press 'a' to start data acquisition and 'z' to stop...\n");

    // Wait for user to press 'a' to start acquisition
    while (!_kbhit() || _getch() != 'a');  // Wait for 'a' key

    printf("Data acquisition started...\n");

    // Start acquiring data
    while (1) {
        status = FT_Read(handle, buffer, BUFFER_SIZE, &bytesRead);
        if (status != FT_OK) {
            fprintf(stderr, "Failed to read from FTDI device.\n");
            break;
        }

        if (bytesRead > 0) {
            // Write the received data to the file
            fwrite(buffer, 1, bytesRead, outFile);
            printf("Read and saved %lu bytes...\n", bytesRead);
        }

        // Check if 'z' key is pressed to stop
        if (_kbhit() && _getch() == 'z') {
            break;
        }
    }

    printf("Data acquisition stopped. File saved: %s\n", fileName);
    fclose(outFile);
}

int main() {
    HMODULE hFtd2xx;
    FT_Open_func FT_Open;
    FT_SetBaudRate_func FT_SetBaudRate;
    FT_SetDataCharacteristics_func FT_SetDataCharacteristics;
    FT_SetFlowControl_func FT_SetFlowControl;
    FT_Purge_func FT_Purge;
    FT_Read_func FT_Read;
    FT_Close_func FT_Close;

    // Load the FTDI library and get function pointers
    if (!loadFtd2xxLibrary(&hFtd2xx, &FT_Open, &FT_SetBaudRate, &FT_SetDataCharacteristics,
                           &FT_SetFlowControl, &FT_Purge, &FT_Read, &FT_Close)) {
        return 1;
    }

    // Prompt user for file name
    char fileName[256];
    getFileName(fileName, sizeof(fileName));

    // Set up the FTDI device
    FT_HANDLE handle;
    if (!setupFTDI(&handle, FT_Open, FT_SetBaudRate, FT_SetDataCharacteristics, FT_SetFlowControl, FT_Purge)) {
        FreeLibrary(hFtd2xx);
        return 1;
    }

    // Start reading and saving data
    readAndSaveData(handle, FT_Read, fileName);

    // Close the FTDI device and free the library
    FT_Close(handle);
    FreeLibrary(hFtd2xx);

    printf("Program finished. Resources freed.\n");
    return 0;
}
