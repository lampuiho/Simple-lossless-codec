#include <stdio.h>
#include <stdlib.h>
#include "codec.h"

#define MAX_CHAR_NAME 150

int main()
{
    char input_name[] = "music_source.wav";
    char output_name[] = "encoded_music_source.lac";
    int choice = 0;

    do{
        printf("Decorrelation-based Lossless Audio Coder \nusing Linear Prediction and Entropy Coding\n");

        /*
        printf("\nSELECT MODE : \n\tENCODER WAV -> NEW CODEC\t1.\n\tDECODER NEW CODEC -> WAV\t2.\n");
        scanf("%d", &choice);

        printf("\nINPUT FILE (CURRENT DIRECTORY)?\n\t");
        scanf("%s", &input_name);
        printf("\nOUTPUT FILE (CURRENT DIRECTORY)?\n\t");
        scanf("%s", &output_name);
        */

        choice = 1;

        if( access( input_name, F_OK ) != -1 ) {
            switch(choice){
                case 1:
                    printf("\nEncoding file...");
                    encoder_main(input_name, output_name);
                    break;
                case 2:
                    printf("\nDecoding file...");
                    //decoder_main(input_name, output_name);
                    break;
                default:
                    printf("\nInvalid choice\n");
                    break;
            }
        } else {
            printf("INVALID INPUT FILE\n");
            int i;
            for(i = 0; i < 250000000; i++);

        }

        scanf("%d", &choice);
        // WINDOWS
        system("cls");

    }while(1);

    return 0;
}
