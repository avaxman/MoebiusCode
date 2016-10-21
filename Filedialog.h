//
//  Filedialog.h
//  testigl
//
//  Created by Amir Vaxman on 12/08/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef testigl_Filedialog_h
#define testigl_Filedialog_h

#include <stdio.h>
#define FILE_DIALOG_MAX_BUFFER 1024

// Sets buffer to a path to an existing file
// buffer[0]=0 on cancel
//
// Usage:
//   char buffer[FILE_DIALOG_MAX_BUFFER];
//   get_open_file_path(buffer);
void get_open_file_path(char buffer[]){
#ifdef __APPLE__
    // For apple use applescript hack
    FILE * output = popen(
                          "osascript -e \""
                          "   tell application \\\"System Events\\\"\n"
                          "           activate\n"
                          "           set existing_file to choose file\n"
                          "   end tell\n"
                          "   set existing_file_path to (POSIX path of (existing_file))\n"
                          "\" 2>/dev/null | tr -d '\n' ","r");
    while ( fgets(buffer, FILE_DIALOG_MAX_BUFFER, output) != NULL ){
    }
#else
    // For every other machine type
    printf("Please enter a file path: ");
    gets(buffer);
#endif
}

// Sets buffer to a path to a new/existing file
// buffer[0]=0 on cancel
//
// Usage:
//   char buffer[FILE_DIALOG_MAX_BUFFER];
//   get_save_file_path(buffer);
void get_save_file_path(char buffer[]){
#ifdef __APPLE__
    // For apple use applescript hack
    // There is currently a bug in Applescript that strips extensions off
    // of chosen existing files in the "choose file name" dialog
    // I'm assuming that will be fixed soon :-)
    FILE * output = popen(
                          "osascript -e \""
                          "   tell application \\\"System Events\\\"\n"
                          "           activate\n"
                          "           set existing_file to choose file name\n"
                          "   end tell\n"
                          "   set existing_file_path to (POSIX path of (existing_file))\n"
                          "\" 2>/dev/null | tr -d '\n' ","r");
    while ( fgets(buffer, FILE_DIALOG_MAX_BUFFER, output) != NULL ){
    }
#else
    // For every other machine type 
    printf("Please enter a file path: ");
    gets(buffer);
#endif
}


#endif
