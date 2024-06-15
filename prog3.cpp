#include "image.h"
#include <time.h>
using ImagePtr = std::shared_ptr<Image>;

/*
Some notes: 
- There shouldn't be a reason to use accumulators or re-association
*/

ImagePtr applyGamma(ImagePtr image_ptr, double gamma);
ImagePtr applyTint(ImagePtr image_ptr, const double *tints);
ImagePtr applyBlur(ImagePtr imag_ptr);
void writeImage(ImagePtr image_ptr);

void process_images(const std::vector<ImagePtr>& image_vector) {
  const double tint_array[] = {0.75, 0, 0};
  for (ImagePtr img : image_vector) {
    writeImage(img);
    img = applyGamma(img, 1.4); 
    img = applyTint(img, tint_array);
    img = applyBlur(img);
    writeImage(img);
  }
  
}
// Apply a Gamma scale to all the pixels of the input image
ImagePtr applyGamma(ImagePtr image_ptr, double gamma) {
  auto output_image_ptr = std::make_shared<Image>(image_ptr->name() + "_gamma", 
                                                  IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  //const int height = in_rows.size();
  //const int width = in_rows[1] - in_rows[0];
  for (unsigned long i = 0; i < in_rows.size(); ++i ) {
    for (int j = 0; j < in_rows[1] - in_rows[0]; ++j ) {
      //const Pixel& p = in_rows[i][j]; 
      double v = 0.3*in_rows[i][j].bgra[2] + 0.59*in_rows[i][j].bgra[1] + 0.11*in_rows[i][j].bgra[0];
      double res = pow(v, gamma);
      if(res > MAX_BGR_VALUE) res = MAX_BGR_VALUE;
      out_rows[i][j] = Pixel(res, res, res);
    }
  }
  return output_image_ptr;
}

// Apply a Gamma optimized to use code motion ()
//Complete
ImagePtr applyGammaCodeMotion(ImagePtr image_ptr, double gamma) {
  auto output_image_ptr = std::make_shared<Image>(image_ptr->name() + "_gamma", 
                                                  IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();

  //here, we can create constants for color channels
  const double red = 0.3;
  const double green = 0.59;
  const double blue = 0.11; 

  //const int height = in_rows.size();
  //const int width = in_rows[1] - in_rows[0];
  for (unsigned long i = 0; i < in_rows.size(); ++i ) {
    for (int j = 0; j < in_rows[1] - in_rows[0]; ++j ) {
      //const Pixel& p = in_rows[i][j]; 
      double v = red*in_rows[i][j].bgra[2] + green*in_rows[i][j].bgra[1] + blue*in_rows[i][j].bgra[0];
      double res = pow(v, gamma);
      if(res > MAX_BGR_VALUE) res = MAX_BGR_VALUE;
      out_rows[i][j] = Pixel(res, res, res);
    }
  }
  return output_image_ptr;
}

//apply gamma optimized to use unrolling
//the point here is to do more during a single iteration, so we would need a j+1
//THIS IS NOT COMPLETE
ImagePtr applyGammaUnrolled(ImagePtr image_ptr, double gamma) {
  auto output_image_ptr = std::make_shared<Image>(image_ptr->name() + "_gamma", 
                                                  IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  //const int height = in_rows.size();
  //const int width = in_rows[1] - in_rows[0];
  for (unsigned long i = 0; i < in_rows.size(); ++i ) {
    for (int j = 0; j < in_rows[1] - in_rows[0]; j+=2 ) {
      //const Pixel& p = in_rows[i][j]; 
      double v = 0.3*in_rows[i][j].bgra[2] + 0.59*in_rows[i][j].bgra[1] + 0.11*in_rows[i][j].bgra[0];
      double res = pow(v, gamma);
      if(res > MAX_BGR_VALUE) res = MAX_BGR_VALUE;
      out_rows[i][j] = Pixel(res, res, res);

      //here, the second would go in: 

      //processing it twice in a single iteration. 

      if (j+1 < in_rows[1] - in_rows[0]){
        double vNext = 0.3 * in_rows[i][j + 1].bgra[2] + 0.59 * in_rows[i][j + 1].bgra[1] + 0.11 * in_rows[i][j + 1].bgra[0];
        double resNext = pow(vNext, gamma);
        if (resNext > MAX_BGR_VALUE) resNext = MAX_BGR_VALUE;
        out_rows[i][j + 1] = Pixel(resNext, resNext, resNext);

      }

    }
   
  }
  return output_image_ptr;
}

// Apply a Gamma scale to all the pixels of the input image
ImagePtr applyGammaEliminateMemoryRef(ImagePtr image_ptr, double gamma) {
  auto output_image_ptr = std::make_shared<Image>(image_ptr->name() + "_gamma", 
                                                  IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  //const int height = in_rows.size();
  //const int width = in_rows[1] - in_rows[0];
  for (unsigned long i = 0; i < in_rows.size(); ++i ) {
    for (int j = 0; j < in_rows[1] - in_rows[0]; ++j ) {

      //instead of referencing the variable array each time at different indexes, initialize them

      auto pixel = in_rows[i][j].bgra[0];
      auto pixel1 = in_rows[i][j].bgra[1];
      auto pixel2 = in_rows[i][i].bgra[2];

      //const Pixel& p = in_rows[i][j]; 
      double v = 0.3*pixel2 + 0.59*pixel1 + 0.11*pixel;
      double res = pow(v, gamma);
      if(res > MAX_BGR_VALUE) res = MAX_BGR_VALUE;
      out_rows[i][j] = Pixel(res, res, res);
    }
  }
  return output_image_ptr;
}

// Apply the tint in the array tints to all the pixels of the input image
ImagePtr applyTint(ImagePtr image_ptr, const double *tints) {
  auto output_image_ptr = 
    std::make_shared<Image>(image_ptr->name() + "_tinted", 
      IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();

  for (unsigned long i = 0; i < image_ptr->rows().size(); ++i ) {
    for (int j = 0; j < image_ptr->rows()[1] - image_ptr->rows()[0]; ++j ) {
      double b = (double)in_rows[i][j].bgra[0] + (MAX_BGR_VALUE-in_rows[i][j].bgra[0])*tints[0];
      double g = (double)in_rows[i][j].bgra[1] + (MAX_BGR_VALUE-in_rows[i][j].bgra[1])*tints[1];
      double r = (double)in_rows[i][j].bgra[2] + (MAX_BGR_VALUE-in_rows[i][j].bgra[0])*tints[2];
      out_rows[i][j].bgra[0] = b > MAX_BGR_VALUE ? MAX_BGR_VALUE:b;
      out_rows[i][j].bgra[1] = g > MAX_BGR_VALUE ? MAX_BGR_VALUE:g;
      out_rows[i][j].bgra[2] = r > MAX_BGR_VALUE ? MAX_BGR_VALUE:r;
    }
  }
  return output_image_ptr;
}



// Apply the tint in the array tints to all the pixels of the input image
ImagePtr applyTintReduceMemoryRefs(ImagePtr image_ptr, const double *tints) {
  auto output_image_ptr = 
    std::make_shared<Image>(image_ptr->name() + "_tinted", 
      IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();

  for (unsigned long i = 0; i < image_ptr->rows().size(); ++i ) {
    for (int j = 0; j < image_ptr->rows()[1] - image_ptr->rows()[0]; ++j ) {
      int first[3] = {in_rows[i][j].bgra[0], in_rows[i][j].bgra[1], in_rows[i][j].bgra[2]};
      int second[3] = {out_rows[i][j].bgra[0], out_rows[i][j].bgra[1], out_rows[i][j].bgra[2]};
      double b = (double)first[0] + (MAX_BGR_VALUE - first[0]) * tints[0];
      double g = (double)first[1] + (MAX_BGR_VALUE - first[1]) * tints[1];
      double r = (double)first[2] + (MAX_BGR_VALUE - first[2]) * tints[2];
      /*
      double b = (double)in_rows[i][j].bgra[0] + (MAX_BGR_VALUE-in_rows[i][j].bgra[0])*tints[0];
      double g = (double)in_rows[i][j].bgra[1] + (MAX_BGR_VALUE-in_rows[i][j].bgra[1])*tints[1];
      double r = (double)in_rows[i][j].bgra[2] + (MAX_BGR_VALUE-in_rows[i][j].bgra[0])*tints[2];
      */
      second[0] = b > MAX_BGR_VALUE ? MAX_BGR_VALUE:b;
      second[1] = g > MAX_BGR_VALUE ? MAX_BGR_VALUE:g;
      second[2] = r > MAX_BGR_VALUE ? MAX_BGR_VALUE:r;
    }
  }
  return output_image_ptr;
}


ImagePtr applyBlur(ImagePtr image_ptr) {
  auto output_image_ptr = 
    std::make_shared<Image>(image_ptr->name() + "_blurred", 
      IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  double b, g, r;

  //int height = in_rows.size();
  //const int width = in_rows[1] - in_rows[0];
  for (unsigned long i = 0; i < in_rows.size(); ++i ) {
    for (int j = 0; j < in_rows[1] - in_rows[0]; ++j ) {
      // Average = ([i-1][j-1] + [i-1][j] + [i-1][j+1] + [i][j-1] + [i][j] + [i][j+1] + [i+1][j-1] + [i+1][j] + [i+1][j+1])/ 9
      if (i == 0) {                        /* first row */
        if (j == 0) {                     /* first row, first column */
          b = (0 + 0 + 0 + 0 + in_rows[i][j].bgra[0] + in_rows[i+1][j].bgra[0] + 0 + in_rows[i][j+1].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = (0 + 0 + 0 + 0 + in_rows[i][j].bgra[1] + in_rows[i+1][j].bgra[1] + 0 + in_rows[i][j+1].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = (0 + 0 + 0 + 0 + in_rows[i][j].bgra[2] + in_rows[i+1][j].bgra[2] + 0 + in_rows[i][j+1].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {          /* first row, last column */
          b = (0 + 0 + 0 + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + 0 + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0] + 0) / 9;
          g = (0 + 0 + 0 + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + 0 + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + 0) / 9;
          r = (0 + 0 + 0 + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + 0 + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + 0) / 9;
        } 
        else {                          /* first row, middle columns */
          b = (0 + 0 + 0 + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = (0 + 0 + 0 + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = (0 + 0 + 0 + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        }
      } 
      else if (i == in_rows.size() - 1) {        /* last row */
        if (j == 0) {             /* last row, first column */
          b = (0 + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + 0 + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + 0 + 0 + 0) / 9;
          g = (0 + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + 0 + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + 0 + 0 + 0) / 9;
          r = (0 + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + 0 + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + 0 + 0 + 0) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {      /* last row, last column */
          b = (in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j+1].bgra[0] + 0 + in_rows[i][i-1].bgra[0] + in_rows[i][j].bgra[0] + 0 + 0 + 0 + 0) / 9;
          g = (in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j+1].bgra[1] + 0 + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + 0 + 0 + 0 + 0) / 9;
          r = (in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j+1].bgra[2] + 0 + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + 0 + 0 + 0 + 0) / 9;
        } 
        else {                          /* last row, middle columns */
          b = (in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + 0 + 0 + 0) / 9;
          g = (in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + 0 + 0 + 0) / 9;
          r = (in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + 0 + 0 + 0) / 9;
        }
      } 
      else {                            /* middle rows */
        if (j == 0) {                 /* middle row, first column */
          b = ( 0 + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + 0 + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + 0 + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = ( 0 + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + 0 + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + 0 + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = ( 0 + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + 0 + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + 0 + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {      /* middle row, last column */
          b = ( in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + 0 + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + 0 + in_rows[i+1][j-1].bgra[0]+ in_rows[i+1][j].bgra[0] + 0) / 9;
          g = ( in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + 0 + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + 0 + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + 0) / 9;
          r = ( in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + 0 + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + 0 + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + 0) / 9;
        } 
        else {                          /* middle row, middle columns */
          b = ( in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = ( in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = ( in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        }
      }
      out_rows[i][j].bgra[0] = (b > MAX_BGR_VALUE)? MAX_BGR_VALUE : b;
      out_rows[i][j].bgra[1] = (g > MAX_BGR_VALUE)? MAX_BGR_VALUE : g;
      out_rows[i][j].bgra[2] = (r > MAX_BGR_VALUE)? MAX_BGR_VALUE : r;
    }
  }
  return output_image_ptr;
}

//applyBlur optimized using code motion
ImagePtr applyBlurCodeMotion(ImagePtr image_ptr) {
  auto output_image_ptr = 
    std::make_shared<Image>(image_ptr->name() + "_blurred", 
      IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  double b, g, r;
  //int height = in_rows.size();
  //const int width = in_rows[1] - in_rows[0];
  for (unsigned long i = 0; i < in_rows.size(); ++i ) {
    for (int j = 0; j < in_rows[1] - in_rows[0]; ++j ) {
      // Average = ([i-1][j-1] + [i-1][j] + [i-1][j+1] + [i][j-1] + [i][j] + [i][j+1] + [i+1][j-1] + [i+1][j] + [i+1][j+1])/ 9
      if (i == 0) {                        /* first row */
        if (j == 0) {                     /* first row, first column */
          b = (in_rows[i][j].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = (in_rows[i][j].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = (in_rows[i][j].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {          /* first row, last column */
          b = (in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0]) / 9;
          g = (in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1]) / 9;
          r = (in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2]) / 9;
        } 
        else {                          /* first row, middle columns */
          b = (in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = (in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = (in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        }
      } 
      else if (i == in_rows.size() - 1) {        /* last row */
        if (j == 0) {             /* last row, first column */
          b = (in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0]) / 9;
          g = (in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1]) / 9;
          r = (in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2]) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {      /* last row, last column */
          b = (in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][i-1].bgra[0] + in_rows[i][j].bgra[0]) / 9;
          g = (in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1]) / 9;
          r = (in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2]) / 9;
        } 
        else {                          /* last row, middle columns */
          b = (in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0]) / 9;
          g = (in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1]) / 9;
          r = (in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2]) / 9;
        }
      } 
      else {                            /* middle rows */
        if (j == 0) {                 /* middle row, first column */
          b = (in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = (in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = (in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {      /* middle row, last column */
          b = ( in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i+1][j-1].bgra[0]+ in_rows[i+1][j].bgra[0]) / 9;
          g = ( in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1]) / 9;
          r = ( in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2]) / 9;
        } 
        else {                          /* middle row, middle columns */
          b = ( in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = ( in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = ( in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        }
      }
      out_rows[i][j].bgra[0] = (b > MAX_BGR_VALUE)? MAX_BGR_VALUE : b;
      out_rows[i][j].bgra[1] = (g > MAX_BGR_VALUE)? MAX_BGR_VALUE : g;
      out_rows[i][j].bgra[2] = (r > MAX_BGR_VALUE)? MAX_BGR_VALUE : r;
    }
  }
  return output_image_ptr;
}

ImagePtr applyBlurUnRolled(ImagePtr image_ptr) {
  auto output_image_ptr = 
    std::make_shared<Image>(image_ptr->name() + "_blurred", 
      IMAGE_WIDTH, IMAGE_HEIGHT);
  auto in_rows = image_ptr->rows();
  auto out_rows = output_image_ptr->rows();
  double b, g, r;
  //int height = in_rows.size();
  //const int width = in_rows[1] - in_rows[0];
  for (unsigned long i = 0; i < in_rows.size(); ++i ) {
    for (int j = 0; j < in_rows[1] - in_rows[0]; ++j ) {
      // Average = ([i-1][j-1] + [i-1][j] + [i-1][j+1] + [i][j-1] + [i][j] + [i][j+1] + [i+1][j-1] + [i+1][j] + [i+1][j+1])/ 9
      if (i == 0) {                        /* first row */
        if (j == 0) {                     /* first row, first column */
          b = (in_rows[i][j].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = (in_rows[i][j].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = (in_rows[i][j].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {          /* first row, last column */
          b = (in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0]) / 9;
          g = (in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1]) / 9;
          r = (in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2]) / 9;
        } 
        else {                          /* first row, middle columns */
          b = (in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = (in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = (in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        }
      } 
      else if (i == in_rows.size() - 1) {        /* last row */
        if (j == 0) {             /* last row, first column */
          b = (in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0]) / 9;
          g = (in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1]) / 9;
          r = (in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2]) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {      /* last row, last column */
          b = (in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][i-1].bgra[0] + in_rows[i][j].bgra[0]) / 9;
          g = (in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1]) / 9;
          r = (in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2]) / 9;
        } 
        else {                          /* last row, middle columns */
          b = (in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0]) / 9;
          g = (in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1]) / 9;
          r = (in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2]) / 9;
        }
      } 
      else {                            /* middle rows */
        if (j == 0) {                 /* middle row, first column */
          b = (in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = (in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = (in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        } 
        else if (j == in_rows[1] - in_rows[0] - 1) {      /* middle row, last column */
          b = ( in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i+1][j-1].bgra[0]+ in_rows[i+1][j].bgra[0]) / 9;
          g = ( in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1]) / 9;
          r = ( in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2]) / 9;
        } 
        else {                          /* middle row, middle columns */
          b = ( in_rows[i-1][j-1].bgra[0] + in_rows[i-1][j].bgra[0] + in_rows[i-1][j+1].bgra[0] + in_rows[i][j-1].bgra[0] + in_rows[i][j].bgra[0] + in_rows[i][j+1].bgra[0] + in_rows[i+1][j-1].bgra[0] + in_rows[i+1][j].bgra[0] + in_rows[i+1][j+1].bgra[0]) / 9;
          g = ( in_rows[i-1][j-1].bgra[1] + in_rows[i-1][j].bgra[1] + in_rows[i-1][j+1].bgra[1] + in_rows[i][j-1].bgra[1] + in_rows[i][j].bgra[1] + in_rows[i][j+1].bgra[1] + in_rows[i+1][j-1].bgra[1] + in_rows[i+1][j].bgra[1] + in_rows[i+1][j+1].bgra[1]) / 9;
          r = ( in_rows[i-1][j-1].bgra[2] + in_rows[i-1][j].bgra[2] + in_rows[i-1][j+1].bgra[2] + in_rows[i][j-1].bgra[2] + in_rows[i][j].bgra[2] + in_rows[i][j+1].bgra[2] + in_rows[i+1][j-1].bgra[2] + in_rows[i+1][j].bgra[2] + in_rows[i+1][j+1].bgra[2]) / 9;
        }
      }
      out_rows[i][j].bgra[0] = (b > MAX_BGR_VALUE)? MAX_BGR_VALUE : b;
      out_rows[i][j].bgra[1] = (g > MAX_BGR_VALUE)? MAX_BGR_VALUE : g;
      out_rows[i][j].bgra[2] = (r > MAX_BGR_VALUE)? MAX_BGR_VALUE : r;
    }
  }
  return output_image_ptr;
}


void writeImage(ImagePtr image_ptr) {
  image_ptr->write( (image_ptr->name() + ".bmp").c_str());
}

void checkCorrectness(std::vector<ImagePtr> image_vector, std::vector<ImagePtr> image_vector_original){
  for(unsigned long k=0; k<image_vector.size(); k++){
    ImagePtr img1 = image_vector[k];
    ImagePtr img2 = image_vector_original[k];
    auto in_rows1 = img1->rows();
    auto in_rows2 = img2->rows();
    int height1 = img1->height();
    int width1 = img1->width();
    int height2 = img2->height();
    int width2 = img2->width();
    if(height1 != height2 || width1 != width2){
      printf("The two images do not have the same dimensions");
      return;
    }
    for(int i=0; i<height1; i++){
      for(int j=0; j<width1; j++){
        if(in_rows1[i][j].value != in_rows2[i][j].value){
            printf("Correctness check failed for pixels image_%ld[%d][%d] = %d and image_refrence_%ld[%d][%d] = %d\n", (k+1), i, j, in_rows1[i][j].value, (k+1), i, i, in_rows2[i][j].value);
            return;
        }
      }
    }
    printf("Correctness check passed for image %ld\n", (k+1));
  }
}
int main() {
  const double tint_array[] = {0.75, 0, 0};
  // Create two vector to hold 4 fractal images
  std::vector<ImagePtr> image_vector, image_vector_reference;
  for ( int i = 2000; i <= 2000000; i *= 10 ) {
    image_vector.push_back(makeFractalImage(i));
    image_vector_reference.push_back(makeFractalImage(i));
  }
  // Store the output of the original functions applyGamma, applyTint, applyBlur in image_vector_reference
  // The output images will be used for checking the correctness of the optimized functions

  //initialize the clock to keep count
  clock_t start, end; 

  for(ImagePtr img:image_vector_reference){
    //ALL APPLYGAMMA FUNCTIONS
    printf("Apply Gamma Functions: %lu\n");
    start = clock(); 
    img = applyGamma(img, 1.4);
    end = clock(); 
    printf("applyGamma: %lu\n", (end-start)); 

    //this is the first optimized function of six
    start = clock(); 
    img = applyGammaCodeMotion(img, 1.4);
    end = clock(); 
    printf("applyGammaCodeMotion: %lu\n", (end-start)); 

    //this is the second optimized funciton of six
    start = clock(); 
    img = applyGammaUnrolled(img, 1.4);
    end = clock(); 
    printf("applyGammaUnrolled: %lu\n", (end-start)); 

    //this is the third optimized function of six
    start = clock(); 
    img = applyGammaEliminateMemoryRef(img, 1.4);
    end = clock(); 
    printf("applyGammaEliminateMemoryRef: %lu\n", (end-start)); 


    printf("Apply Tint Functions: %lu\n");
    start = clock(); 
    img = applyTint(img, tint_array);
    end = clock(); 
    printf("applyTint: %lu\n", (end-start)); 

    //this is one of them
    start = clock(); 
    img = applyTintReduceMemoryRefs(img, tint_array);
    end = clock(); 
    printf("applyTintReduceMemoryRefs: %lu\n", (end-start)); 

    printf("Apply Blur Functions: %lu\n");
    start = clock(); 
    img = applyBlur(img);
    end = clock(); 
    printf("applyBlur: %lu\n", (end-start)); 

    start = clock(); 
    img = applyBlurCodeMotion(img);
    end = clock(); 
    printf("applyBlur: %lu\n", (end-start));

    start = clock(); 
    img = applyBlurUnRolled(img);
    end = clock(); 
    printf("applyBlur: %lu\n", (end-start));

  }
  // Process the images in the vector image_vector
  process_images(image_vector);
  // check the output images of process_images against the images processed by the original functions
  checkCorrectness(image_vector, image_vector_reference);
  return 0;
}