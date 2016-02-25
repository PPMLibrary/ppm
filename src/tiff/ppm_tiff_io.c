#ifdef __TIFF
#include "tiffio.h"

#ifdef __MPI
#include <mpi.h>
#endif

TIFF   *tif;
tmsize_t ScanlineSize; /* this is the machine addressing size type, only it's signed */
void   *bufSc;  /* image data buf ref */
void   *bufSt;  /* image data buf ref */
uint16 tdirt;  /* directory index */
char software[] = "MOSAIC Group PPM";

struct {
  uint32 width;             /* image width in pixels */
  uint32 length;            /* image height in pixels */
  uint32 depth;             /* z depth of image */
  uint16 datatype;          /* $use SampleFormat */
  uint16 bitsPerSample;     /* bits per channel (sample) */
  uint16 compression;       /* data compression technique */
  uint16 sampleFormat;      /* !data sample format */
  uint16 samplesPerPixel;   /* samples per pixel */
  uint16 photometric;       /* photometric interpretation */
  uint16 fillOrder;         /* data order within a byte */
  /* MSB2LSB  most significant -> least */
  /* LSB2MSB  least significant -> most */
  uint16 planarConfig;      /* storage organization */
  uint16 resolutionUnit;    /* units of resolutions */
  uint16 npages;            /* number of pages in a multipage file */
  float  xResolution;       /* pixels/resolution in x */
  float  yResolution;       /* pixels/resolution in y */
} ppmTIFFTAG;


int open_tiff(char *filename) {

  tif = TIFFOpen(filename,"r");

  if (!tif) {
    return -1;
  }

  ScanlineSize=TIFFScanlineSize(tif);

  bufSc = (unsigned char *)_TIFFmalloc(ScanlineSize);
  if (!bufSc){
    _TIFFfree(bufSc);
    bufSc=NULL;
    //Can't malloc strip
    return -3;
  }

  return 0;
}

int open_bigtiff(char *filename) {

  tif = TIFFOpen(filename,"r8");

  if (!tif) {
    return -1;
  }

  ScanlineSize=TIFFScanlineSize(tif);

  bufSc = (unsigned char *)_TIFFmalloc(ScanlineSize);
  if (!bufSc){
    _TIFFfree(bufSc);
    bufSc=NULL;
    //Can't malloc strip
    return -3;
  }

  return 0;
}

int open_write_tiff(char *filename) {

  tif = TIFFOpen(filename,"a+");

  if (!tif) {
    return -1;
  }
  return 0;
}

int open_write_bigtiff(char *filename) {

  tif = TIFFOpen(filename,"w8");

  if (!tif) {
    return -1;
  }
  return 0;
}

int close_tiff() {

  if (bufSc != NULL) {
    _TIFFfree(bufSc);
    bufSc=NULL;
  }

  if (bufSt != NULL) {
    _TIFFfree(bufSt);
    bufSt=NULL;
  }

  TIFFClose(tif);

  return 0;
}

#ifdef __MPI
int read_tiff_infoC(char *filename,int *ngrid,int *bitsPerSample,int *samplesPerPixel,int *rank,MPI_Fint *comm)
#else
int read_tiff_infoC(char *filename,int *ngrid,int *bitsPerSample,int *samplesPerPixel,int *rank,int *comm)
#endif
{
  if (*rank == 0) {

    tif = TIFFOpen(filename,"r");
    if (!tif) {
      return -1;
    }

    // This is where we read the header of the .TIF file and
    // obtain all required information
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,     &ppmTIFFTAG.width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH,    &ppmTIFFTAG.length);
    TIFFGetField(tif, TIFFTAG_IMAGEDEPTH,     &ppmTIFFTAG.depth);
    TIFFGetField(tif, TIFFTAG_DATATYPE,       &ppmTIFFTAG.datatype);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE,  &ppmTIFFTAG.bitsPerSample);
    TIFFGetField(tif, TIFFTAG_COMPRESSION,    &ppmTIFFTAG.compression);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL,&ppmTIFFTAG.samplesPerPixel);
    TIFFGetField(tif, TIFFTAG_PHOTOMETRIC,    &ppmTIFFTAG.photometric);
    TIFFGetField(tif, TIFFTAG_FILLORDER,      &ppmTIFFTAG.fillOrder);
    TIFFGetField(tif, TIFFTAG_PLANARCONFIG,   &ppmTIFFTAG.planarConfig);
    TIFFGetField(tif, TIFFTAG_RESOLUTIONUNIT, &ppmTIFFTAG.resolutionUnit);
    TIFFGetField(tif, TIFFTAG_XRESOLUTION,    &ppmTIFFTAG.xResolution);
    TIFFGetField(tif, TIFFTAG_YRESOLUTION,    &ppmTIFFTAG.yResolution);
    TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT,   &ppmTIFFTAG.sampleFormat);

    ppmTIFFTAG.npages = TIFFNumberOfDirectories(tif);

    // force this one for now...
    if (ppmTIFFTAG.sampleFormat == 0) ppmTIFFTAG.sampleFormat = SAMPLEFORMAT_VOID; // !untyped data
    if (ppmTIFFTAG.bitsPerSample == 0) ppmTIFFTAG.bitsPerSample = 8;
    if (ppmTIFFTAG.samplesPerPixel == 0) ppmTIFFTAG.samplesPerPixel =1;

    TIFFClose(tif);
  }

#ifdef __MPI
  MPI_Datatype ppmTIFFTAGMPI;
  MPI_Datatype type[3]={MPI_UNSIGNED,MPI_UNSIGNED_SHORT,MPI_FLOAT};
  int blocklen[3]={3,10,2};
  MPI_Aint disp[3];
  MPI_Aint base;
  MPI_Comm ccomm;
  int i;

  ccomm=MPI_Comm_f2c(*comm);

  MPI_Get_address(&ppmTIFFTAG,disp);
  MPI_Get_address(&ppmTIFFTAG.datatype,disp+1);
  MPI_Get_address(&ppmTIFFTAG.xResolution,disp+2);

  base = disp[0];
  for (i=0; i<3; i++) disp[i] -= base;

  MPI_Type_create_struct(3,blocklen,disp,type,&ppmTIFFTAGMPI);
  MPI_Type_commit(&ppmTIFFTAGMPI);
  MPI_Bcast(&ppmTIFFTAG,1,ppmTIFFTAGMPI,0,ccomm);
  MPI_Type_free(&ppmTIFFTAGMPI);
#endif

  // Read the x and y values from the image and assign those to Ngrid
  // Be careful with the x/y assignment
  ngrid[0] = ppmTIFFTAG.width;
  ngrid[1] = ppmTIFFTAG.length;
  ngrid[2] = (int) ppmTIFFTAG.npages;

  *bitsPerSample   = (int) ppmTIFFTAG.bitsPerSample;
  *samplesPerPixel = (int) ppmTIFFTAG.samplesPerPixel;

  return 0;
}

int write_tiff_header1(int *bitsPerSample,int *samplesPerPixel,int *width,int *length)
{

  if (!tif) {
    return -1;
  }

  ppmTIFFTAG.width           = (uint32) *width;
  ppmTIFFTAG.length          = (uint32) *length;
  ppmTIFFTAG.bitsPerSample   = (uint16) *bitsPerSample;
  ppmTIFFTAG.samplesPerPixel = (uint16) *samplesPerPixel;
  ppmTIFFTAG.compression     = COMPRESSION_NONE;
  ppmTIFFTAG.resolutionUnit  = RESUNIT_NONE;

  // This is where we write the header of the .TIFF file with
  // all the required information
  TIFFSetField(tif, TIFFTAG_SOFTWARE,       software);
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,     ppmTIFFTAG.width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH,    ppmTIFFTAG.length);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,  ppmTIFFTAG.bitsPerSample);
  TIFFSetField(tif, TIFFTAG_COMPRESSION,    ppmTIFFTAG.compression);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL,ppmTIFFTAG.samplesPerPixel);
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, ppmTIFFTAG.resolutionUnit);
  TIFFSetField(tif, TIFFTAG_ORIENTATION,    ORIENTATION_TOPLEFT);
  switch (ppmTIFFTAG.bitsPerSample) {
    case 32:
      ppmTIFFTAG.sampleFormat=SAMPLEFORMAT_IEEEFP;/* !IEEE floating point data */
      break;
  }

  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, ppmTIFFTAG.sampleFormat);

  switch (ppmTIFFTAG.samplesPerPixel) {
    case 1:
      ppmTIFFTAG.planarConfig = PLANARCONFIG_CONTIG; /* single image plane */
      //       ppmTIFFTAG.planarConfig = PLANARCONFIG_SEPARATE; /* separate planes of data */
      ppmTIFFTAG.photometric  = PHOTOMETRIC_MINISBLACK;
      break;

    case 3:
      ppmTIFFTAG.planarConfig = PLANARCONFIG_SEPARATE; /* separate planes of data */
      ppmTIFFTAG.photometric  = PHOTOMETRIC_RGB;
      break;

    default:
      //Unsupported image format.
      return -2;
  }

  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, ppmTIFFTAG.planarConfig);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC,  ppmTIFFTAG.photometric);

  if (bufSc != NULL) {
    _TIFFfree(bufSc);
    bufSc=NULL;
  }

  ScanlineSize=TIFFScanlineSize(tif);
  bufSc = (unsigned char *)_TIFFmalloc(ScanlineSize);
  if (!bufSc){
    _TIFFfree(bufSc);
    bufSc=NULL;
    //Can't malloc buf
    return -3;
  }

  if (bufSt != NULL) {
    _TIFFfree(bufSt);
    bufSt=NULL;
  }

  return 0;
}

int write_tiff_header2(int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  TIFFSetField(tif, TIFFTAG_SOFTWARE,       software);
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,     ppmTIFFTAG.width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH,    ppmTIFFTAG.length);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,  ppmTIFFTAG.bitsPerSample);
  TIFFSetField(tif, TIFFTAG_COMPRESSION,    ppmTIFFTAG.compression);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL,ppmTIFFTAG.samplesPerPixel);
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, ppmTIFFTAG.resolutionUnit);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,   ppmTIFFTAG.sampleFormat);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG,   ppmTIFFTAG.planarConfig);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC,    ppmTIFFTAG.photometric);

  uint16 pge = (uint16) *page;
  uint16 npg = (uint16) *npages;

  /* We are writing single page of the multipage file */
  TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE); /* one page of many */
  /* Set the page number */
  TIFFSetField(tif, TIFFTAG_PAGENUMBER, pge, npg);

  return 0;
}



int read_tiff_scanline_int1(int *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 x;

  TIFFReadScanline(tif, bufSc, row, smp);

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<wdt; x++){
        scanline[x] = *((unsigned char*)bufSc+x);
      }
      break;

    case 16:
      for (x=0; x<wdt; x++){
        scanline[x] = *((unsigned short*)bufSc+x);
      }
      break;

    case 32:
      for (x=0; x<wdt; x++){
        scanline[x] = *((unsigned int*)bufSc+x);
      }
      break;
  }

  return 0;
}

int read_tiff_scanline_long1(long *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 x;

  TIFFReadScanline(tif, bufSc, row, smp);

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<wdt; x++){
        scanline[x] = *((unsigned char*)bufSc+x);
      }
      break;

    case 16:
      for (x=0; x<wdt; x++){
        scanline[x] = *((unsigned short*)bufSc+x);
      }
      break;

    case 32:
      for (x=0; x<wdt; x++){
        scanline[x] = *((unsigned int*)bufSc+x);
      }
      break;
  }

  return 0;
}

int read_tiff_scanline_float1(float *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 x;

  TIFFReadScanline(tif, bufSc, row, smp);

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<wdt; x++){
        scanline[x] = *((unsigned char*)bufSc+x);
      }
      break;

    case 16:
      for (x=0; x<wdt; x++){
        scanline[x] = *((unsigned short*)bufSc+x);
      }
      break;

    case 32:
      for (x=0; x<wdt; x++){
        scanline[x] = *((float*)bufSc+x);
      }
      break;
  }

  return 0;
}

int read_tiff_scanline_double1(double *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 x;

  TIFFReadScanline(tif, bufSc, row, smp);

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<wdt; x++){
        scanline[x] = *((unsigned char*)bufSc+x);
      }
      break;

    case 16:
      for (x=0; x<wdt; x++){
        scanline[x] = *((unsigned short*)bufSc+x);
      }
      break;

    case 32:
      for (x=0; x<wdt; x++){
        scanline[x] = *((float*)bufSc+x);
      }
      break;
  }

  return 0;
}

int read_tiff_scanline_int2(int *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 x,y;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((unsigned char*)bufSc+x);
          k++;
        }
      }
      break;
    case 16:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((unsigned short*)bufSc+x);
          k++;
        }
      }
      break;
    case 32:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((unsigned int*)bufSc+x);
          k++;
        }
      }
      break;
  }

  return 0;
}

int read_tiff_scanline_long2(long *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 x,y;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((unsigned char*)bufSc+x);
          k++;
        }
      }
      break;
    case 16:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((unsigned short*)bufSc+x);
          k++;
        }
      }
      break;
    case 32:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((unsigned int*)bufSc+x);
          k++;
        }
      }
      break;
  }

  return 0;
}

int read_tiff_scanline_float2(float *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 x,y;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((unsigned char*)bufSc+x);
          k++;
        }
      }
      break;
    case 16:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((unsigned short*)bufSc+x);
          k++;
        }
      }
      break;
    case 32:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((float*)bufSc+x);
          k++;
        }
      }
      break;
  }

  return 0;
}

int read_tiff_scanline_double2(double *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 x,y;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((unsigned char*)bufSc+x);
          k++;
        }
      }
      break;
    case 16:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((unsigned short*)bufSc+x);
          k++;
        }
      }
      break;
    case 32:
      for (y=row;y<row+lnt;y++){
        TIFFReadScanline(tif, bufSc, y, smp);
        for (x=0; x<wdt; x++){
          scanline[k] = *((float*)bufSc+x);
          k++;
        }
      }
      break;
  }

  return 0;
}

int read_tiff_strip_int1(int *strip,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  if (bufSc != NULL) {
    _TIFFfree(bufSc);
    bufSc=NULL;
  }

  if (bufSt == NULL) {
    bufSt = (unsigned char *)_TIFFmalloc(TIFFStripSize(tif));
    if (!bufSt){
      _TIFFfree(bufSt);
      bufSt=NULL;
      //Can't malloc strip
      return -3;
    }
  }

  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 topRow;
  uint32 rowsPerStrip;
  uint32 nRowsToConvert;
  uint32 x;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsPerStrip);

  for (topRow = 0; topRow < lnt; topRow += rowsPerStrip) {
    nRowsToConvert = (topRow + rowsPerStrip > lnt ? lnt - topRow : rowsPerStrip);

    if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, topRow, smp), bufSt, nRowsToConvert*ScanlineSize) < 0) {
      _TIFFfree(bufSt);
      bufSt=NULL;
      return -4;
    }

    switch (*bitsPerSample) {
      case 8:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((unsigned char*)bufSt+x);
          k++;
        }
        break;

      case 16:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((unsigned short*)bufSt+x);
          k++;
        }
        break;

      case 32:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((unsigned int*)bufSt+x);
          k++;
        }
        break;
    }
  }

  return 0;
}

int read_tiff_strip_long1(long *strip,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  if (bufSc != NULL) {
    _TIFFfree(bufSc);
    bufSc=NULL;
  }

  if (bufSt == NULL) {
    bufSt = (unsigned char *)_TIFFmalloc(TIFFStripSize(tif));
    if (!bufSt){
      _TIFFfree(bufSt);
      bufSt=NULL;
      //Can't malloc strip
      return -3;
    }
  }

  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 topRow;
  uint32 rowsPerStrip;
  uint32 nRowsToConvert;
  uint32 x;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsPerStrip);

  for (topRow = 0; topRow < lnt; topRow += rowsPerStrip) {
    nRowsToConvert = (topRow + rowsPerStrip > lnt ? lnt - topRow : rowsPerStrip);

    if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, topRow, smp), bufSt, nRowsToConvert*ScanlineSize) < 0) {
      _TIFFfree(bufSt);
      bufSt=NULL;
      return -4;
    }

    switch (*bitsPerSample) {
      case 8:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((unsigned char*)bufSt+x);
          k++;
        }
        break;

      case 16:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((unsigned short*)bufSt+x);
          k++;
        }
        break;

      case 32:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((unsigned int*)bufSt+x);
          k++;
        }
        break;
    }
  }

  return 0;
}

int read_tiff_strip_float1(float *strip,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  if (bufSc != NULL) {
    _TIFFfree(bufSc);
    bufSc=NULL;
  }

  if (bufSt == NULL) {
    bufSt = (unsigned char *)_TIFFmalloc(TIFFStripSize(tif));
    if (!bufSt){
      _TIFFfree(bufSt);
      bufSt=NULL;
      //Can't malloc strip
      return -3;
    }
  }

  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 topRow;
  uint32 rowsPerStrip;
  uint32 nRowsToConvert;
  uint32 x;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsPerStrip);

  for (topRow = 0; topRow < lnt; topRow += rowsPerStrip) {
    nRowsToConvert = (topRow + rowsPerStrip > lnt ? lnt - topRow : rowsPerStrip);

    if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, topRow, smp), bufSt, nRowsToConvert*ScanlineSize) < 0) {
      _TIFFfree(bufSt);
      bufSt=NULL;
      return -4;
    }

    switch (*bitsPerSample) {
      case 8:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((unsigned char*)bufSt+x);
          k++;
        }
        break;

      case 16:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((unsigned short*)bufSt+x);
          k++;
        }
        break;

      case 32:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((float*)bufSt+x);
          k++;
        }
        break;
    }
  }

  return 0;
}

int read_tiff_strip_double1(double *strip,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  if (bufSc != NULL) {
    _TIFFfree(bufSc);
    bufSc=NULL;
  }

  if (bufSt == NULL) {
    bufSt = (unsigned char *)_TIFFmalloc(TIFFStripSize(tif));
    if (!bufSt){
      _TIFFfree(bufSt);
      bufSt=NULL;
      //Can't malloc strip
      return -3;
    }
  }

  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 topRow;
  uint32 rowsPerStrip;
  uint32 nRowsToConvert;
  uint32 x;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsPerStrip);

  for (topRow = 0; topRow < lnt; topRow += rowsPerStrip) {
    nRowsToConvert = (topRow + rowsPerStrip > lnt ? lnt - topRow : rowsPerStrip);

    if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, topRow, smp), bufSt, nRowsToConvert*ScanlineSize) < 0) {
      _TIFFfree(bufSt);
      bufSt=NULL;
      return -4;
    }

    switch (*bitsPerSample) {
      case 8:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((unsigned char*)bufSt+x);
          k++;
        }
        break;

      case 16:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((unsigned short*)bufSt+x);
          k++;
        }
        break;

      case 32:
        for (x=0; x<nRowsToConvert*wdt; x++){
          strip[k] = *((float*)bufSt+x);
          k++;
        }
        break;
    }
  }

  return 0;
}






int write_tiff_scanline_int1(int *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32   x;

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<wdt; x++){
        *((unsigned char*)bufSc+x) = (unsigned char) scanline[x];
      }
      break;

    case 16:
      for (x=0; x<wdt; x++){
        *((unsigned short*)bufSc+x) = (unsigned short) scanline[x];
      }
      break;

    case 32:
      for (x=0; x<wdt; x++){
        *((unsigned int*)bufSc+x) = (unsigned int) scanline[x];
      }
      break;
  }

  TIFFWriteScanline(tif, bufSc, row, smp);

  return 0;
}

int write_tiff_scanline_long1(long *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32   x;

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<wdt; x++){
        *((unsigned char*)bufSc+x) = (unsigned char) scanline[x];
      }
      break;

    case 16:
      for (x=0; x<wdt; x++){
        *((unsigned short*)bufSc+x) = (unsigned short) scanline[x];
      }
      break;

    case 32:
      for (x=0; x<wdt; x++){
        *((unsigned int*)bufSc+x) = (unsigned int) scanline[x];
      }
      break;
  }

  TIFFWriteScanline(tif, bufSc, row, smp);

  return 0;
}


int write_tiff_scanline_float1(float *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32   x;

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<wdt; x++){
        *((unsigned char*)bufSc+x) = (unsigned char) scanline[x];
      }
      break;

    case 16:
      for (x=0; x<wdt; x++){
        *((unsigned short*)bufSc+x) = (unsigned short) scanline[x];
      }
      break;

    case 32:
      for (x=0; x<wdt; x++){
        *((float*)bufSc+x) = scanline[x];
      }
      break;
  }

  TIFFWriteScanline(tif, bufSc, row, smp);

  return 0;
}

int write_tiff_scanline_double1(double *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32   x;

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<wdt; x++){
        *((unsigned char*)bufSc+x) = (unsigned char) scanline[x];
      }
      break;

    case 16:
      for (x=0; x<wdt; x++){
        *((unsigned short*)bufSc+x) = (unsigned short) scanline[x];
      }
      break;

    case 32:
      for (x=0; x<wdt; x++){
        *((float*)bufSc+x) = (float) scanline[x];
      }
      break;
  }

  TIFFWriteScanline(tif, bufSc, row, smp);

  return 0;
}

int write_tiff_scanline_int2(int *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 x,y;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((unsigned char*)bufSc+x) = (unsigned char) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
    case 16:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((unsigned short*)bufSc+x) = (unsigned short) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
    case 32:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((unsigned int*)bufSc+x) = (unsigned int) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
  }

  if (*npages > 1) {
    TIFFWriteDirectory(tif);
  }

  return 0;
}

int write_tiff_scanline_long2(long *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }


  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 x,y;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((unsigned char*)bufSc+x) = (unsigned char) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
    case 16:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((unsigned short*)bufSc+x) = (unsigned short) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
    case 32:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((unsigned int*)bufSc+x) = (unsigned int) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
  }

  if (*npages > 1) {
    TIFFWriteDirectory(tif);
  }

  return 0;
}


int write_tiff_scanline_float2(float *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 x,y;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((unsigned char*)bufSc+x) = (unsigned char) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
    case 16:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((unsigned short*)bufSc+x) = (unsigned short) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
    case 32:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((float*)bufSc+x) = scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
  }

  if (*npages > 1) {
    TIFFWriteDirectory(tif);
  }

  return 0;
}

int write_tiff_scanline_double2(double *scanline,int *rownm,int *bitsPerSample,int *Sample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  uint32 row = (uint32) *rownm;
  uint16 smp = (uint16) *Sample;
  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 x,y;
  uint32 k=0;

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((unsigned char*)bufSc+x) = (unsigned char) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
    case 16:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((unsigned short*)bufSc+x) = (unsigned short) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
    case 32:
      for (y=row;y<row+lnt;y++){
        for (x=0; x<wdt; x++){
          *((float*)bufSc+x) = (float) scanline[k];
          k++;
        }
        TIFFWriteScanline(tif, bufSc, y, smp);
      }
      break;
  }

  if (*npages > 1) {
    TIFFWriteDirectory(tif);
  }

  return 0;
}




int write_tiff_strip_int1(int *strip,int *bitsPerSample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  if (bufSc != NULL) {
    _TIFFfree(bufSc);
    bufSc=NULL;
  }

  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 pixelsPerRaster = wdt*lnt;
  uint32 pixelsPerRasterSize=pixelsPerRaster*(*bitsPerSample)/8;
  uint32 x;

  if (bufSt == NULL) {
    bufSt = (unsigned char *)_TIFFmalloc(pixelsPerRasterSize);
    if (!bufSt){
      _TIFFfree(bufSt);
      bufSt=NULL;
      //Can't malloc strip
      return -3;
    }
  }

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<pixelsPerRaster; x++){
        *((unsigned char*)bufSt+x) = (unsigned char) strip[x];
      }
      break;

    case 16:
      for (x=0; x<pixelsPerRaster; x++){
        *((unsigned short*)bufSt+x) = (unsigned short) strip[x];
      }
      break;

    case 32:
      for (x=0; x<pixelsPerRaster; x++){
        *((unsigned int*)bufSt+x) = (unsigned int) strip[x];
      }
      break;
  }

  if (TIFFWriteEncodedStrip(tif, 0, bufSt, pixelsPerRasterSize) < 0) {
    _TIFFfree(bufSt);
    bufSt=NULL;
    //-5: Can't write strip.
    return -5;
  }

  if (*npages > 1) {
    TIFFWriteDirectory(tif);
  }

  return 0;
}

int write_tiff_strip_long1(long *strip,int *bitsPerSample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  if (bufSc != NULL) {
    _TIFFfree(bufSc);
    bufSc=NULL;
  }

  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 pixelsPerRaster = wdt*lnt;
  uint32 pixelsPerRasterSize=pixelsPerRaster*(*bitsPerSample)/8;
  uint32 x;

  if (bufSt == NULL) {
    bufSt = (unsigned char *)_TIFFmalloc(pixelsPerRasterSize);
    if (!bufSt){
      _TIFFfree(bufSt);
      bufSt=NULL;
      //Can't malloc strip
      return -3;
    }
  }

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<pixelsPerRaster; x++){
        *((unsigned char*)bufSt+x) = (unsigned char) strip[x];
      }
      break;

    case 16:
      for (x=0; x<pixelsPerRaster; x++){
        *((unsigned short*)bufSt+x) = (unsigned short) strip[x];
      }
      break;

    case 32:
      for (x=0; x<pixelsPerRaster; x++){
        *((unsigned int*)bufSt+x) = (unsigned int) strip[x];
      }
      break;
  }

  if (TIFFWriteEncodedStrip(tif, 0, bufSt, pixelsPerRasterSize) < 0) {
    _TIFFfree(bufSt);
    bufSt=NULL;
    //-5: Can't write strip.
    return -5;
  }

  if (*npages > 1) {
    TIFFWriteDirectory(tif);
  }

  return 0;
}

int write_tiff_strip_float1(float *strip,int *bitsPerSample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  if (bufSc != NULL) {
    _TIFFfree(bufSc);
    bufSc=NULL;
  }

  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 pixelsPerRaster = wdt*lnt;
  uint32 pixelsPerRasterSize=pixelsPerRaster*(*bitsPerSample)/8;
  uint32 x;

  if (bufSt == NULL) {
    bufSt = (unsigned char *)_TIFFmalloc(pixelsPerRasterSize);
    if (!bufSt){
      _TIFFfree(bufSt);
      bufSt=NULL;
      //Can't malloc strip
      return -3;
    }
  }

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<pixelsPerRaster; x++){
        *((unsigned char*)bufSt+x) = (unsigned char) strip[x];
      }
      break;

    case 16:
      for (x=0; x<pixelsPerRaster; x++){
        *((unsigned short*)bufSt+x) = (unsigned short) strip[x];
      }
      break;

    case 32:
      for (x=0; x<pixelsPerRaster; x++){
        *((float*)bufSt+x) = strip[x];
      }
      break;
  }

  if (TIFFWriteEncodedStrip(tif, 0, bufSt, pixelsPerRasterSize) < 0) {
    _TIFFfree(bufSt);
    bufSt=NULL;
    //-5: Can't write strip.
    return -5;
  }

  if (*npages > 1) {
    TIFFWriteDirectory(tif);
  }

  return 0;
}

int write_tiff_strip_double1(double *strip,int *bitsPerSample,int *width,int *length,int *page,int *npages)
{

  if (!tif) {
    return -1;
  }

  if (bufSc != NULL) {
    _TIFFfree(bufSc);
    bufSc=NULL;
  }

  uint32 wdt = (uint32) *width;
  uint32 lnt = (uint32) *length;
  uint32 pixelsPerRaster = wdt*lnt;
  uint32 pixelsPerRasterSize=pixelsPerRaster*(*bitsPerSample)/8;
  uint32 x;

  if (bufSt == NULL) {
    bufSt = (unsigned char *)_TIFFmalloc(pixelsPerRasterSize);
    if (!bufSt){
      _TIFFfree(bufSt);
      bufSt=NULL;
      //Can't malloc strip
      return -3;
    }
  }

  if (*npages > 1) {
    uint16 pge = (uint16) *page;
    TIFFSetDirectory(tif, pge);
  }

  switch (*bitsPerSample) {
    case 8:
      for (x=0; x<pixelsPerRaster; x++){
        *((unsigned char*)bufSt+x) = (unsigned char) strip[x];
      }
      break;

    case 16:
      for (x=0; x<pixelsPerRaster; x++){
        *((unsigned short*)bufSt+x) = (unsigned short) strip[x];
      }
      break;

    case 32:
      for (x=0; x<pixelsPerRaster; x++){
        *((float*)bufSt+x) = (float) strip[x];
      }
      break;
  }

  if (TIFFWriteEncodedStrip(tif, 0, bufSt, pixelsPerRasterSize) < 0) {
    _TIFFfree(bufSt);
    bufSt=NULL;
    //-5: Can't write strip.
    return -5;
  }

  if (*npages > 1) {
    TIFFWriteDirectory(tif);
  }

  return 0;
}

#endif
