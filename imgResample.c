/*---------------------------------------------------------------------------
  This program reads sample PPM images and resamples them up or down by a factor
  or uses a quick 2X down sample
  
  gcc -g imgResample.c -o imgResample -lm 
  gcc -g imgResample.c -o imgResample -lm -fsanitize=address -fsanitize=undefined
  
 resample code:
  https://stackoverflow.com/questions/34622717/bicubic-interpolation-in-c
  https://pastebin.com/sQDQg7SG
----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef struct {
   unsigned char red,green,blue;
} PPMPixel;

typedef struct {
   int x, y;
   PPMPixel *data;
} PPMImage;

PPMImage *resize2(PPMImage *source_image);

#define CREATOR "FELIXKLEMM"
#define RGB_COMPONENT_COLOR 255
#define BUFFER_SIZE (256) 



// Clamps the returned value between min and max, otherwise returns the value
#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

int debug = 0;

/*---------------------------------------------------------------------------
   This function reads a PPM image and returns the binary pixel data 
   in a single 1D array.    
   
   const char *filename - File name to open
  
   Returns: PPMImage *readPPM    Pointer to a mallloced data structure
   
   Error Handling:   exits with an error code
----------------------------------------------------------------------------*/
static PPMImage *readPPM(const char *filename) {
   char buff[BUFFER_SIZE];
   PPMImage *img;
   FILE *fp;

   int c, rgb_comp_color;

   //open PPM file for reading
   fp = fopen(filename, "rb");
   if(!fp) {
      fprintf(stderr, "Unable to open file '%s'\n", filename);
      exit(1);
   }

   //read image format
   if(!fgets(buff, sizeof(buff), fp)) {
      perror(filename);
      exit(1);
   }

   //check the image format
   if(buff[0] != 'P' || buff[1] != '6') {
      fprintf(stderr, "Invalid image format (must be 'P6')\n");
      exit(1);
   }

   //alloc memory form image
   img = (PPMImage *)malloc(sizeof(PPMImage));
   if(!img) {
      fprintf(stderr, "Unable to allocate memory\n");
      exit(1);
   }

   //check for comments
   c = getc(fp);
   while (c == '#') {
      while (getc(fp) != '\n');
      c = getc(fp);
   }

   ungetc(c, fp);
   //read image size information
   if(fscanf(fp, "%d %d", &img->x, &img->y) != 2) {
      fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
      exit(1);
   }

   //read rgb component
   if(fscanf(fp, "%d", &rgb_comp_color) != 1) {
      fprintf(stderr, "Invalid rgb component (error loading '%s')\n", filename);
      exit(1);
   }

   //check rgb component depth
   if(rgb_comp_color!= RGB_COMPONENT_COLOR) {
      fprintf(stderr, "'%s' error only RGB supported\n", filename);
      exit(1);
   }

   while (fgetc(fp) != '\n');
   //memory allocation for pixel data
   img->data = (PPMPixel*)malloc(img->x * img->y * sizeof(PPMPixel));

   if(!img) {
      fprintf(stderr, "Unable to allocate memory\n");
      exit(1);
   }

   //read pixel data from file
   if(fread(img->data, 3 * img->x, img->y, fp) != img->y) {
      fprintf(stderr, "Error loading image '%s'\n", filename);
      exit(1);
   }

   fclose(fp);
   return img;
}


/*---------------------------------------------------------------------------
   This function allocates an image based on the input size and the resamples
   size.    
   
      PMImage *source   - Pointer to an open input image
      double scale      - The scale factor to use
  
   Returns: PPMImage *readPPM    Pointer to a malloced data structure
   
   Error Handling:   exits with an error code
----------------------------------------------------------------------------*/
static PPMImage *init_destination_image(PPMImage *source, double scale) {
   PPMImage *img;
   //alloc memory form image
   img = (PPMImage *)malloc(sizeof(PPMImage));
   if(!img) {

     fprintf(stderr, "Unable to allocate memory\n");
     exit(1);
   }

   //memory allocation for pixel data
    if (debug) { printf("XxY %dx%d scale %d PPM %ld dest size %ld\n", source->x, source->y, (int)scale, 
                sizeof(PPMPixel*),source->x*(int)(scale+.1)*source->y*(int)(scale+.5)*sizeof(PPMPixel));} 

   img->data = (PPMPixel*)malloc(source->x*(scale+.1)*source->y*(scale+.5)*sizeof(PPMPixel));
   if(!img) {
      fprintf(stderr, "Unable to allocate memory\n");
      exit(1);
   }
   return img;
}


/*---------------------------------------------------------------------------
   Writes a PPM format image file
      
      char *filename - The PPM file image name to write 
      PPMImage *img  - A pointer to an (PPM) image object
      
      Returns: nothing
      
      Error handling: exit with a return code
----------------------------------------------------------------------------*/
void writePPM(const char *filename, PPMImage *img) {
   FILE *fp;
   //open file for output
   fp = fopen(filename, "wb");
   if (!fp) {
       fprintf(stderr, "Unable to open file '%s'\n", filename);
       exit(1);
   }

   //write the header file as ascii data on each line
   //image format
   fprintf(fp, "P6\n");

   //comments
   fprintf(fp, "# Created by %s\n",CREATOR);

   //image size
   fprintf(fp, "%d %d\n",img->x,img->y);

   // rgb component depth
   fprintf(fp, "%d\n",RGB_COMPONENT_COLOR);

   // pixel data - binary 
   fwrite(img->data, 3 * img->x, img->y, fp);
   fclose(fp);
}


/*---------------------------------------------------------------------------
  
  
----------------------------------------------------------------------------*/
double cubic_hermite(double A, double B, double C, double D, double t) {

   double a = -A / 2.0f + (3.0f*B) / 2.0f - (3.0f*C) / 2.0f + D / 2.0f;
   double b = A - (5.0f*B) / 2.0f + 2.0f*C - D / 2.0f;
   double c = -A / 2.0f + C / 2.0f;
   double d = B;

   return a*t*t*t + b*t*t + c*t + d;
}

/*---------------------------------------------------------------------------
  This functio returns a rgb byte array for the data at the given point x,y*
  BUT will never exceed the array bounds so it handles the edge effect.
  
      PPMImage *source_image  - Pointer to an images
      int x, int y            - Image x,y coordinates
      uint8_t temp[]          - Pointer to 3 byte array to return data

   return: nothing
   
   Error handling: none
----------------------------------------------------------------------------*/
void get_pixel_clamped(PPMImage *source_image, int x, int y, uint8_t temp[])  {

   // Keep from exceeding the array index
   CLAMP(x, 0, source_image->x - 1);
   CLAMP(y, 0, source_image->y - 1);
   
   temp[0] = source_image->data[x+(source_image->x*y)].red;
   temp[1] = source_image->data[x+(source_image->x*y)].green;
   temp[2] = source_image->data[x+(source_image->x*y)].blue;
}

/*---------------------------------------------------------------------------
  

   return: nothing
   
   Error handling: none
----------------------------------------------------------------------------*/
void sample_bicubic(PPMImage *source_image, double u, double v, uint8_t sample[]) {

   double x = (u * source_image->x)-0.5;
   int xint = (int)x;
   double xfract = x-floor(x);

   double y = (v * source_image->y) - 0.5;
   int yint = (int)y;
   double yfract = y - floor(y);
   
   int i;

   uint8_t p00[3], p10[3], p20[3], p30[3];
   uint8_t p01[3], p11[3], p21[3], p31[3];
   uint8_t p02[3], p12[3], p22[3], p32[3];
   uint8_t p03[3], p13[3], p23[3], p33[3];
   
   // 1st row
   get_pixel_clamped(source_image, xint - 1, yint - 1, p00);   
   get_pixel_clamped(source_image, xint + 0, yint - 1, p10);
   get_pixel_clamped(source_image, xint + 1, yint - 1, p20);
   get_pixel_clamped(source_image, xint + 2, yint - 1, p30);
   
   // 2nd row
   get_pixel_clamped(source_image, xint - 1, yint + 0, p01);
   get_pixel_clamped(source_image, xint + 0, yint + 0, p11);
   get_pixel_clamped(source_image, xint + 1, yint + 0, p21);
   get_pixel_clamped(source_image, xint + 2, yint + 0, p31);

   // 3rd row
   get_pixel_clamped(source_image, xint - 1, yint + 1, p02);
   get_pixel_clamped(source_image, xint + 0, yint + 1, p12);
   get_pixel_clamped(source_image, xint + 1, yint + 1, p22);
   get_pixel_clamped(source_image, xint + 2, yint + 1, p32);

   // 4th row
   get_pixel_clamped(source_image, xint - 1, yint + 2, p03);
   get_pixel_clamped(source_image, xint + 0, yint + 2, p13);
   get_pixel_clamped(source_image, xint + 1, yint + 2, p23);
   get_pixel_clamped(source_image, xint + 2, yint + 2, p33);
   
   // interpolate bi-cubically!
   for (i = 0; i < 3; i++) {
      double col0 = cubic_hermite(p00[i], p10[i], p20[i], p30[i], xfract);
      double col1 = cubic_hermite(p01[i], p11[i], p21[i], p31[i], xfract);
      double col2 = cubic_hermite(p02[i], p12[i], p22[i], p32[i], xfract);
      double col3 = cubic_hermite(p03[i], p13[i], p23[i], p33[i], xfract);
  
      double value = cubic_hermite(col0, col1, col2, col3, yfract);
  
      CLAMP(value, 0.0f, 255.0f);
  
      sample[i] = (uint8_t)value;
       
   }
   if (debug) { printf("sample[]=%d %d %d\n", sample[0], sample[1], sample[2]); }
}


/*---------------------------------------------------------------------------
   This function resizes an input image to create a new destination image.
   
         PPMImage *source_image        - Input image to resize_image
         PPMImage *destination_image   - defined output images
         double scale                  - resize value
   
   returns: nothing
   
   error handling: none
----------------------------------------------------------------------------*/
void resize_image(PPMImage *source_image, PPMImage *destination_image, double scale) {
   uint8_t sample[3];
   int y, x;

   destination_image->x = (long)((double)(source_image->x)*scale);
   destination_image->y = (long)((double)(source_image->y)*scale);

   printf("Source x-width=%d | y-width=%d\n",source_image->x, source_image->y);
   printf("Dest   x-width=%d | y-width=%d\n",destination_image->x, destination_image->y);
    
   for (y = 0; y < destination_image->y; y++) {

      double v = (double)y / (double)(destination_image->y - 1);
      
      for (x = 0; x < destination_image->x; ++x) {
   
         double u = (double)x / (double)(destination_image->x - 1);
         if (debug) {printf("v=%f  u=%f\n",v, u);}
         sample_bicubic(source_image, u, v, sample);
   
         if (debug) {printf("x,y %d,%d offset %d\n", x,y, x+((destination_image->x)*y));}
          
         destination_image->data[x+((destination_image->x)*y)].red 	= sample[0];
         destination_image->data[x+((destination_image->x)*y)].green	= sample[1];  
         destination_image->data[x+((destination_image->x)*y)].blue 	= sample[2];  
      }
   }
}


/*---------------------------------------------------------------------------
   Main test program, parses command lines.  See help for documentation
  
----------------------------------------------------------------------------*/
int main(int argc, char *argv[]) {
   // Help
   if (argc != 4) {
      printf("This program resamples PPM images up or down using cubic resampling\n");
      printf("or a quick 2x down sample\n");
      printf("Syntax is  %s factor infile  outfile\n", argv[0]);
      printf("    factor - '2x' or a floating point number\n");
      printf("  eg  %s  0.5  in.ppm  out.ppm\n", argv[0]);
      printf("      %s  2x   in.ppm  out.ppm\n", argv[0]);
      return(99);
   }
   
   double scale = atof(argv[1]); 
   PPMImage *source_image;
   PPMImage *destination_image;
   printf("Starting...\n\n");
   
   if (strcmp(argv[1], "2x") && (scale <= 0.0)) { printf("error scale must be positive\n"); return(99);}
   
   if(remove(argv[3]) == 0) {	printf("Deleting old image %s...\n\n", argv[3]);}

    source_image = readPPM(argv[2]);
    if (debug) {printf("Infile x,y %dx%d\n", source_image->x, source_image->y);}
    
   // Check for quick 
   if (strcmp(argv[1], "2x") == 0) {
      printf("Using quick 2X downsample\n");
      destination_image = resize2(source_image);
   }
   else {
    destination_image = init_destination_image(source_image, scale);
    resize_image(source_image, destination_image, scale);
   }
   
    writePPM(argv[3], destination_image);
    
    // return memory
    free(source_image->data);
    free(source_image);
    source_image = NULL;
    
    free(destination_image->data);
    free(destination_image);
    destination_image = NULL;
    
    
   return(0);
}

/*---------------------------------------------------------------------------
   This is a quick function to resizes an input image down by 2
   
         PPMImage *source_image        - Input image to resize_image
   returns:  PPMImage *destination_image  
   
   error handling: none
----------------------------------------------------------------------------*/
PPMImage *resize2(PPMImage *source_image) {
   PPMImage *destination_image;
   
   destination_image = init_destination_image(source_image, 0.5);

   // fix up the size to make it always smaller
   destination_image->x = (source_image->x/2); 
   destination_image->y = (source_image->y/2); 
   
   int x,y;
   
   for (y = 0; y < destination_image->y; y++) {
      for (x = 0; x < destination_image->x; x++) {
         destination_image->data[x+((destination_image->x)*y)].red =   source_image->data[(2*x  +((source_image->x)*2*y))].red/4;
         destination_image->data[x+((destination_image->x)*y)].red +=  source_image->data[(2*x+1+((source_image->x)*2*y))].red/4;
         destination_image->data[x+((destination_image->x)*y)].red +=  source_image->data[(2*x  +(((source_image->x))*(2*y+1)))].red/4;
         destination_image->data[x+((destination_image->x)*y)].red +=  source_image->data[(2*x+1+(((source_image->x))*(2*y+1)))].red/4;
                                                                       
         destination_image->data[x+((destination_image->x)*y)].green = source_image->data[(2*x  +((source_image->x)*2*y))].green/4;
         destination_image->data[x+((destination_image->x)*y)].green +=source_image->data[(2*x+1+((source_image->x)*2*y))].green/4;
         destination_image->data[x+((destination_image->x)*y)].green +=source_image->data[(2*x  +(((source_image->x))*(2*y+1)))].green/4;
         destination_image->data[x+((destination_image->x)*y)].green +=source_image->data[(2*x+1+(((source_image->x))*(2*y+1)))].green/4;
                                                                       
         destination_image->data[x+((destination_image->x)*y)].blue =  source_image->data[(2*x  +((source_image->x)*2*y))].blue/4;
         destination_image->data[x+((destination_image->x)*y)].blue += source_image->data[(2*x+1+((source_image->x)*2*y))].blue/4;
         destination_image->data[x+((destination_image->x)*y)].blue += source_image->data[(2*x  +(((source_image->x))*(2*y+1)))].blue/4;
         destination_image->data[x+((destination_image->x)*y)].blue += source_image->data[(2*x+1+(((source_image->x))*(2*y+1)))].blue/4;
      } // End y
   } // End x
   
   return(destination_image);
}