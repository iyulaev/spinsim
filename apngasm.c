/* APNG Assembler 2.0
 *
 * Copyright (c) 2009 Max Stepin
 * maxst at users.sourceforge.net
 *
 * GNU LGPL information
 * --------------------
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#define PNG_ZBUF_SIZE  32768

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "png.h"

#ifndef WIN32
  #define _stricmp strcasecmp
#endif

#if defined(_MSC_VER) && _MSC_VER >= 1300
inline unsigned short swap16(unsigned short data) {return(_byteswap_ushort(data));}
inline unsigned int swap32(unsigned int data) {return(_byteswap_ulong(data));}
#elif __linux__
#include <byteswap.h>
inline unsigned short swap16(unsigned short data) {return(bswap_16(data));}
inline unsigned int swap32(unsigned int data) {return(bswap_32(data));}
#else
inline unsigned short swap16(unsigned short data) {return((data >> 8) | (data << 8));}
inline unsigned int swap32(unsigned int data) {return((swap16(data) << 16) | swap16(data >> 16));}
#endif

unsigned char png_sign[8] = {137,  80,  78,  71,  13,  10,  26,  10};
unsigned char png_Software[27] = { 83, 111, 102, 116, 119, 97, 114, 101, '\0', 
                                   65,  80,  78,  71,  32, 65, 115, 115, 101, 
                                  109,  98, 108, 101, 114, 32,  50,  46,  48};

struct OP 
{
  z_stream        zstream;
  unsigned char * zbuf;
  int  x;
  int  y;
  int  w;
  int  h;
  int  valid;
} op[12];

unsigned int next_seq_num = 0;
unsigned char * row_buf;
unsigned char * sub_row;
unsigned char * up_row;
unsigned char * avg_row;
unsigned char * paeth_row;

unsigned char * LoadPNG(char * szImage, int *pWidth, int *pHeight, int *pDepth, int *pType, png_color *pPal, int *pPsize, unsigned char *pTrns, int *pTsize, int *pRes)
{
  FILE          * f;
  png_structp     png_ptr;
  png_infop       info_ptr;
  png_bytepp      row_pointers = NULL;
  unsigned char * image_data = NULL;
  png_uint_32     width, height, i, rowbytes;
  int             depth, coltype;
  int iy_returned;

  *pRes = 0;

  if ((f = fopen(szImage, "rb")) != 0)
  {
    unsigned char sig[8];

    iy_returned = fread(sig, 1, 8, f);
    if (png_sig_cmp(sig, 0, 8) == 0)
    {
      png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
      if (png_ptr != NULL)
      {
        info_ptr = png_create_info_struct(png_ptr);
        
        if (info_ptr != NULL) 
        {
          if (setjmp(png_jmpbuf(png_ptr)) == 0)
          {
            png_init_io(png_ptr, f);
            png_set_sig_bytes(png_ptr, 8);
            png_read_info(png_ptr, info_ptr);
            png_get_IHDR(png_ptr, info_ptr, &width, &height, &depth, &coltype, NULL, NULL, NULL);
            *pWidth  = width;
            *pHeight = height;
            *pDepth  = depth;
            *pType   = coltype;

            if (png_get_valid(png_ptr, info_ptr, PNG_INFO_PLTE))
            {
              memcpy(pPal, info_ptr->palette, info_ptr->num_palette * sizeof(png_color));
              *pPsize = info_ptr->num_palette;
            }
            else
              *pPsize = 0;

            if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
            {
              png_color_16p  trans_color;
              png_bytep      trans_alpha;

              png_get_tRNS(png_ptr, info_ptr, &trans_alpha, pTsize, &trans_color);

              if (coltype == PNG_COLOR_TYPE_GRAY)
              {
                pTrns[0] = trans_color->gray >> 8;
                pTrns[1] = trans_color->gray & 0xFF;
                if (depth == 16)
                {
                  pTrns[1] = pTrns[0]; pTrns[0] = 0;
                }
                *pTsize = 2;
              }
              else
              if (coltype == PNG_COLOR_TYPE_RGB)
              {
                pTrns[0] = trans_color->red >> 8;
                pTrns[1] = trans_color->red & 0xFF;
                pTrns[2] = trans_color->green >> 8;
                pTrns[3] = trans_color->green & 0xFF;
                pTrns[4] = trans_color->blue >> 8;
                pTrns[5] = trans_color->blue & 0xFF;
                if (depth == 16)
                {
                  pTrns[1] = pTrns[0]; pTrns[0] = 0;
                  pTrns[3] = pTrns[2]; pTrns[2] = 0;
                  pTrns[5] = pTrns[4]; pTrns[4] = 0;
                }
                *pTsize = 6;
              }
              else
                memcpy(pTrns, trans_alpha, *pTsize);
            }
            else
              *pTsize = 0;

            if (depth > 8)
              png_set_strip_16(png_ptr);

            if (depth < 8)
            {
              if (coltype == PNG_COLOR_TYPE_GRAY)
                png_set_expand_gray_1_2_4_to_8(png_ptr);
              else
                png_set_packing(png_ptr);
            }

            png_read_update_info(png_ptr, info_ptr);
            *pDepth  = png_get_bit_depth(png_ptr, info_ptr);

            rowbytes = png_get_rowbytes(png_ptr, info_ptr);

            if ((image_data = (unsigned char *)malloc(rowbytes*height)) != NULL) 
            {
              if ((row_pointers = (png_bytepp)malloc(height*sizeof(png_bytep))) != NULL) 
              {
                for (i=0; i<height; i++)
                  row_pointers[i] = image_data + i*rowbytes;

                png_read_image(png_ptr, row_pointers);
                free(row_pointers);
                png_read_end(png_ptr, NULL);
                png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
              }
              else
              {
                png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
                free(image_data);
                *pRes = 7;
              }
            }
            else
            {
              png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
              *pRes = 6;
            }
          }
          else
          {
            png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
            *pRes = 5;
          }
        }
        else
        {
          png_destroy_read_struct(&png_ptr, NULL, NULL);
          *pRes = 4;
        }
      }
      else
        *pRes = 3;
    }
    else
      *pRes = 2;

    fclose(f);
  }
  else
    *pRes = 1;

  return image_data;
}

void write_chunk(FILE * f, const char * name, unsigned char * data, unsigned int length)
{
  unsigned int crc = crc32(0, Z_NULL, 0);
  unsigned int len = swap32(length);

  fwrite(&len, 1, 4, f);
  fwrite(name, 1, 4, f);
  crc = crc32(crc, (const Bytef *)name, 4);

  if (memcmp(name, "fdAT", 4) == 0)
  {
    unsigned int seq = swap32(next_seq_num++);
    fwrite(&seq, 1, 4, f);
    crc = crc32(crc, (const Bytef *)(&seq), 4);
    length -= 4;
  }

  if (data != NULL && length > 0)
  {
    fwrite(data, 1, length, f);
    crc = crc32(crc, data, length);
  }

  unsigned int crc2 = swap32(crc);
  fwrite(&crc2, 1, 4, f);
}

void write_IDATs(FILE * f, int frame, unsigned char * data, unsigned int length, unsigned int idat_size)
{
  unsigned int z_cmf = data[0];
  if ((z_cmf & 0x0f) == 8 && (z_cmf & 0xf0) <= 0x70)
  {
    if (length >= 2)
    {
      unsigned int z_cinfo = z_cmf >> 4;
      unsigned int half_z_window_size = 1 << (z_cinfo + 7);
      while (idat_size <= half_z_window_size && half_z_window_size >= 256)
      {
        z_cinfo--;
        half_z_window_size >>= 1;
      }
      z_cmf = (z_cmf & 0x0f) | (z_cinfo << 4);
      if (data[0] != (unsigned char)z_cmf)
      {
        data[0] = (unsigned char)z_cmf;
        data[1] &= 0xe0;
        data[1] += (unsigned char)(0x1f - ((z_cmf << 8) + data[1]) % 0x1f);
      }
    }
  }

  while (length > 0)
  {
    unsigned int ds = length;
    if (ds > PNG_ZBUF_SIZE)
      ds = PNG_ZBUF_SIZE;

    if (frame == 0)
      write_chunk(f, "IDAT", data, ds);
    else
      write_chunk(f, "fdAT", data, ds+4);

    data += ds;
    length -= ds;
  }
}

int get_rect(int w, int h, unsigned char *pimg1, unsigned char *pimg2, unsigned char *ptemp, int *px, int *py, int *pw, int *ph, int bpp)
{
  int   i, j;
  int   x_min = w-1;
  int   y_min = h-1;
  int   x_max = 0;
  int   y_max = 0;
  int   diffnum = 0;
  int   over_is_possible = 1;

  if (bpp == 1)
  {
    unsigned char *pa = pimg1;
    unsigned char *pb = pimg2;

    for (j=0; j<h; j++)
    for (i=0; i<w; i++)
    {
      if (*pa++ != *pb++)
      {
        diffnum++;
        if (i<x_min) x_min = i;
        if (i>x_max) x_max = i;
        if (j<y_min) y_min = j;
        if (j>y_max) y_max = j;
      }
    }
    over_is_possible = 0;
  }
  else
  if (bpp == 2)
  {
    unsigned short *pa = (unsigned short *)pimg1;
    unsigned short *pb = (unsigned short *)pimg2;
    unsigned short *pc = (unsigned short *)ptemp;

    for (j=0; j<h; j++)
    for (i=0; i<w; i++)
    {
      unsigned int c1 = *pa++;
      unsigned int c2 = *pb++;
      if ((c1 != c2) && ((c1>>8) || (c2>>8)))
      {
        diffnum++;
        if ((c2 >> 8) != 0xFF) over_is_possible = 0;
        if (i<x_min) x_min = i;
        if (i>x_max) x_max = i;
        if (j<y_min) y_min = j;
        if (j>y_max) y_max = j;
      }
      else
        c2 = 0;

      *pc++ = c2;
    }
  }
  else
  if (bpp == 3)
  {
    unsigned char *pa = pimg1;
    unsigned char *pb = pimg2;

    for (j=0; j<h; j++)
    for (i=0; i<w; i++)
    {
      if ((pa[0] != pb[0]) || (pa[1] != pb[1]) || (pa[2] != pb[2]))
      {
        diffnum++;
        if (i<x_min) x_min = i;
        if (i>x_max) x_max = i;
        if (j<y_min) y_min = j;
        if (j>y_max) y_max = j;
      }
      pa += 3;
      pb += 3;
    }
    over_is_possible = 0;
  }
  else
  if (bpp == 4)
  {
    unsigned int *pa = (unsigned int *)pimg1;
    unsigned int *pb = (unsigned int *)pimg2;
    unsigned int *pc = (unsigned int *)ptemp;

    for (j=0; j<h; j++)
    for (i=0; i<w; i++)
    {
      unsigned int c1 = *pa++;
      unsigned int c2 = *pb++;
      if ((c1 != c2) && ((c1>>24) || (c2>>24)))
      {
        diffnum++;
        if ((c2 >> 24) != 0xFF) over_is_possible = 0;
        if (i<x_min) x_min = i;
        if (i>x_max) x_max = i;
        if (j<y_min) y_min = j;
        if (j>y_max) y_max = j;
      }
      else
        c2 = 0;

      *pc++ = c2;
    }
  }

  if (diffnum == 0)
  {
    *px = *py = 0;
    *pw = *ph = 1; 
  }
  else
  {
    *px = x_min;
    *py = y_min;
    *pw = x_max-x_min+1;
    *ph = y_max-y_min+1;
  }

  return over_is_possible;
}

void deflate_rect(unsigned char *pdata, int x, int y, int w, int h, int bpp, int stride, int zbuf_size, int n)
{
  int i, j, v;
  int a, b, c, pa, pb, pc, p;
  int rowbytes = w * bpp;
  unsigned char * prev = NULL;
  unsigned char * row  = pdata + y*stride + x*bpp;
  unsigned char * out;

  op[n*2].valid = 1;
  op[n*2].zstream.next_out = op[n*2].zbuf;
  op[n*2].zstream.avail_out = zbuf_size;

  op[n*2+1].valid = 1;
  op[n*2+1].zstream.next_out = op[n*2+1].zbuf;
  op[n*2+1].zstream.avail_out = zbuf_size;

  for (j=0; j<h; j++)
  {
    unsigned int    sum = 0;
    unsigned char * best_row = row_buf;
    unsigned int    mins = ((unsigned int)(-1)) >> 1;

    out = row_buf+1;
    for (i=0; i<rowbytes; i++)
    {
      v = out[i] = row[i];
      sum += (v < 128) ? v : 256 - v;
    }
    mins = sum;

    sum = 0;
    out = sub_row+1;
    for (i=0; i<bpp; i++)
    {
      v = out[i] = row[i];
      sum += (v < 128) ? v : 256 - v;
    }
    for (i=bpp; i<rowbytes; i++)
    {
      v = out[i] = row[i] - row[i-bpp];
      sum += (v < 128) ? v : 256 - v;
      if (sum > mins) break;
    }
    if (sum < mins)
    {
      mins = sum;
      best_row = sub_row;
    }

    if (prev)
    {
      sum = 0;
      out = up_row+1;
      for (i=0; i<rowbytes; i++)
      {
        v = out[i] = row[i] - prev[i];
        sum += (v < 128) ? v : 256 - v;
        if (sum > mins) break;
      }
      if (sum < mins)
      {
        mins = sum;
        best_row = up_row;
      }

      sum = 0;
      out = avg_row+1;
      for (i=0; i<bpp; i++)
      {
        v = out[i] = row[i] - prev[i]/2;
        sum += (v < 128) ? v : 256 - v;
      }
      for (i=bpp; i<rowbytes; i++)
      {
        v = out[i] = row[i] - (prev[i] + row[i-bpp])/2;
        sum += (v < 128) ? v : 256 - v;
        if (sum > mins) break;
      }
      if (sum < mins)
      { 
        mins = sum;
        best_row = avg_row;
      }

      sum = 0;
      out = paeth_row+1;
      for (i=0; i<bpp; i++)
      {
        v = out[i] = row[i] - prev[i];
        sum += (v < 128) ? v : 256 - v;
      }
      for (i=bpp; i<rowbytes; i++)
      {
        a = row[i-bpp];
        b = prev[i];
        c = prev[i-bpp];
        p = b - c;
        pc = a - c;
        pa = abs(p);
        pb = abs(pc);
        pc = abs(p + pc);
        p = (pa <= pb && pa <=pc) ? a : (pb <= pc) ? b : c;
        v = out[i] = row[i] - p;
        sum += (v < 128) ? v : 256 - v;
        if (sum > mins) break;
      }
      if (sum < mins)
      {
        best_row = paeth_row;
      }
    }

    op[n*2].zstream.next_in = row_buf;
    op[n*2].zstream.avail_in = rowbytes + 1;
    deflate(&op[n*2].zstream, Z_NO_FLUSH);

    op[n*2+1].zstream.next_in = best_row;
    op[n*2+1].zstream.avail_in = rowbytes + 1;
    deflate(&op[n*2+1].zstream, Z_NO_FLUSH);

    prev = row;
    row += stride;
  }

  deflate(&op[n*2].zstream, Z_FINISH);
  deflate(&op[n*2+1].zstream, Z_FINISH);

  op[n*2].x = op[n*2+1].x = x;
  op[n*2].y = op[n*2+1].y = y;
  op[n*2].w = op[n*2+1].w = w;
  op[n*2].h = op[n*2+1].h = h;
}

int main(int argc, char** argv)
{
  char  * szOutput;
  char  * szImage;
  char  * szOpt;
  char    szFormat[256];
  char    szNext[256];
  int i, j;
  int delay_num = -1;
  int delay_den = -1;
  int cur = 0;
  int num = 0;
  int first = 0;
  int width, height, depth, type, bpp, frames;
  int rowbytes, imagesize, idat_size, zbuf_size, zsize;
  int palsize, trnssize;
  int x0, y0, w0, h0, x1, y1, w1, h1, zero_trans, try_over;
  png_color  palette[256];
  unsigned char trns[256];
  unsigned char * zbuf;
  unsigned char * imagetemp;
  unsigned char   dop, bop, c;

  FILE * f;
    
  printf("\nAPNG Assembler 2.0\n");

  if (argc <= 2)
  {
    printf("\nUsage: apngasm.exe output.png frame001.png [1] [10] [/f]\n"
           "\n1/10 is the default delay. Use /f to skip the first frame.\n");
    return 1;
  }

  szOutput = argv[1];
  szImage  = argv[2];

  for (i=3; i<argc; i++)
  {
    szOpt = argv[i];

    if ((szOpt[0] == '/') || (szOpt[0] == '-'))
    {
      if ((szOpt[1] == 'f') || (szOpt[1] == 'F'))
        first = 1;
    }
    else
    {
      int n = atoi(szOpt);
      if ((n != 0) || (strcmp(szOpt, "0") == 0))
      {
        if (delay_num == -1) delay_num = n;
        else
        if (delay_den == -1) delay_den = n;
      }
    }
  }

  if (delay_num == -1) delay_num = 1;
  if (delay_den == -1) delay_den = 10;

  char * szExt = strrchr( szImage, '.' );
  char * szCount = szExt;

  if ((_stricmp(szExt, ".png") == 0) && (szImage<szCount))
  {
    while ((*(szCount-1) >= '0') && (*(szCount-1) <= '9'))
    {
      szCount--;
      if (szImage == szCount) break;
      if (szCount == szExt-5) break;
    }
    
    strcpy(szFormat, szImage);

    if (szCount < szExt)
    {
      cur = atoi(szCount);
      sprintf(szFormat+(szCount-szImage), "%%0%ldd%s", szExt-szCount, szExt);
    }
    else
      printf( "Error: *.png sequence not found\n" );
  }
  else
    printf( "Error: '.png' extention expected\n" );

  frames = 0;

  if ((f = fopen(szImage, "rb")) == 0)
  {
    printf("Error: can't open the file '%s'", szImage);
    return 1;
  }

  do
  {
    frames++;
    fclose(f);
    sprintf(szNext, szFormat, cur+frames);
    f = fopen(szNext, "rb");
  } 
  while (f != 0);

  unsigned char ** images = (unsigned char **)malloc(frames*sizeof(unsigned char *));
  if (images == NULL)
  {
    printf( "Error: not enough memory\n" );
    return 1;
  }

  for (i=0; i<frames; i++)
  {
    sprintf(szNext, szFormat, cur+i);
    printf("reading %s (%d of %d)\n", szNext, i-first+1, frames-first);
        
    int w, h, d, t, res;
    int ps, ts;
    png_color     pl[256];
    unsigned char tr[256];

    images[i] = LoadPNG(szNext, &w, &h, &d, &t, &pl[0], &ps, &tr[0], &ts, &res);

    if (images[i] == NULL)
    {
      printf( "Error: LoadPNG() failed\n" );
      return res;
    }

    if (i == 0)
    {
      width = w;
      height = h;
      depth = d;
      type = t;
      palsize = ps;
      trnssize = ts;
      memset(trns, 255, 256);
      if (ps) memcpy(palette, pl, ps*sizeof(png_color));
      if (ts) memcpy(trns, tr, ts);
    }
    else
    {
      if ((width != w) || (height != h) || (depth != d) || (type != t))
      {
        printf( "Error: incorrect image resolution/type\n" );
        return 1;
      }
    }
  }
    
  bpp = 1;
  if (type == 2)
    bpp = 3;
  else
  if (type == 4)
    bpp = 2;
  else
  if (type == 6)
    bpp = 4;

  rowbytes  = width * bpp;
  imagesize = rowbytes * height;
  idat_size = (rowbytes + 1) * height;
  zbuf_size = idat_size + ((idat_size + 7) >> 3) + ((idat_size + 63) >> 6) + 11;

  zero_trans = 0;

  if (type == 3)
  {
    int tcolor = 0;

    for (j=0; j<256; j++)
    if (trns[j] == 0)
    {
      tcolor = j;
      zero_trans = 1;
      break;
    }

    if (zero_trans && tcolor)
    {
      c = palette[tcolor].red;   palette[tcolor].red   = palette[0].red;   palette[0].red   = c;
      c = palette[tcolor].green; palette[tcolor].green = palette[0].green; palette[0].green = c;
      c = palette[tcolor].blue;  palette[tcolor].blue  = palette[0].blue;  palette[0].blue  = c;
      trns[tcolor] = trns[0]; trns[0] = 0;

      for (i=0; i<frames; i++)
        for (j=0; j<imagesize; j++)
        {
          if (*(images[i]+j) == tcolor) 
            *(images[i]+j) = 0;
          else
          if (*(images[i]+j) == 0) 
            *(images[i]+j) = tcolor;
        }

      if (trnssize)
        while (trnssize > 0 && trns[trnssize-1] == 255)
          trnssize--;
    }
  }
  else
  if (type == 4)
  {
    zero_trans = 1;
    for (i=0; i<frames; i++)
      for (j=0; j<width*height; j++)
        if (*(images[i]+j*2+1) == 0) 
          *(images[i]+j*2) = 0;
  }
  else
  if (type == 6)
  {
    zero_trans = 1;
    for (i=0; i<frames; i++)
      for (j=0; j<width*height; j++)
        if (*(images[i]+j*4+3) == 0) 
          *(images[i]+j*4) = *(images[i]+j*4+1) = *(images[i]+j*4+2) = 0;
  }

  for (i=0; i<12; i++)
  {
    op[i].zstream.data_type = Z_BINARY;
    op[i].zstream.zalloc = Z_NULL;
    op[i].zstream.zfree = Z_NULL;
    op[i].zstream.opaque = Z_NULL;

    if (i & 1)
      deflateInit2(&op[i].zstream, Z_BEST_COMPRESSION, 8, 15, 8, Z_FILTERED);
    else
      deflateInit2(&op[i].zstream, Z_BEST_COMPRESSION, 8, 15, 8, Z_DEFAULT_STRATEGY);

    op[i].zbuf = (unsigned char *)malloc(zbuf_size);
    if (op[i].zbuf == NULL)
    {
      printf( "Error: not enough memory\n" );
      return 1;
    }
  }

  imagetemp = (unsigned char *)malloc(imagesize);
  zbuf = (unsigned char *)malloc(zbuf_size);
  row_buf = (unsigned char *)malloc(rowbytes + 1);
  sub_row = (unsigned char *)malloc(rowbytes + 1);
  up_row = (unsigned char *)malloc(rowbytes + 1);
  avg_row = (unsigned char *)malloc(rowbytes + 1);
  paeth_row = (unsigned char *)malloc(rowbytes + 1);

  if (imagetemp && zbuf && row_buf && sub_row && up_row && avg_row && paeth_row)
  {
    row_buf[0] = 0;
    sub_row[0] = 1;
    up_row[0] = 2;
    avg_row[0] = 3;
    paeth_row[0] = 4;
  }
  else
  {
    printf( "Error: not enough memory\n" );
    return 1;
  }

  if ((f = fopen(szOutput, "wb")) != 0)
  {
    struct IHDR 
    {
      unsigned int    mWidth;
      unsigned int    mHeight;
      unsigned char   mDepth;
      unsigned char   mColorType;
      unsigned char   mCompression;
      unsigned char   mFilterMethod;
      unsigned char   mInterlaceMethod;
    } ihdr;

    struct acTL 
    {
      unsigned int    mFrameCount;
      unsigned int    mLoopCount;
    } actl;

    struct fcTL 
    {
      unsigned int    mSeq;
      unsigned int    mWidth;
      unsigned int    mHeight;
      unsigned int    mXOffset;
      unsigned int    mYOffset;
      unsigned short  mDelayNum;
      unsigned short  mDelayDen;
      unsigned char   mDisposeOp;
      unsigned char   mBlendOp;
    } fctl;

    fwrite(png_sign, 1, 8, f);

    ihdr.mWidth            = swap32(width);
    ihdr.mHeight           = swap32(height);
    ihdr.mDepth            = depth;
    ihdr.mColorType        = type;
    ihdr.mCompression      = 0;
    ihdr.mFilterMethod     = 0;
    ihdr.mInterlaceMethod  = 0;
    write_chunk(f, "IHDR", (unsigned char *)(&ihdr), 13);

    if (frames > 1)
    {
      actl.mFrameCount  = swap32(frames-first);
      actl.mLoopCount   = 0;
      write_chunk(f, "acTL", (unsigned char *)(&actl), 8);
    }
    else
      first = 0;

    if (palsize > 0)
      write_chunk(f, "PLTE", (unsigned char *)(&palette), palsize*sizeof(png_color));

    if (trnssize > 0)
      write_chunk(f, "tRNS", trns, trnssize);

    x0 = 0;
    y0 = 0;
    w0 = width;
    h0 = height;
    bop = 0;

    printf("saving frame %d of %d\n", 1-first, frames-first);
    deflate_rect(images[0], x0, y0, w0, h0, bpp, rowbytes, zbuf_size, 0);

    if (op[0].zstream.total_out <= op[1].zstream.total_out)
    {
      zsize = op[0].zstream.total_out;
      memcpy(zbuf, op[0].zbuf, zsize);
    }
    else
    {
      zsize = op[1].zstream.total_out;
      memcpy(zbuf, op[1].zbuf, zsize);
    }

    deflateReset(&op[0].zstream);
    op[0].zstream.data_type = Z_BINARY;
    deflateReset(&op[1].zstream);
    op[1].zstream.data_type = Z_BINARY;

    if (first)
    {
      write_IDATs(f, 0, zbuf, zsize, idat_size);

      printf("saving frame %d of %d\n", 1, frames-first);
      deflate_rect(images[1], x0, y0, w0, h0, bpp, rowbytes, zbuf_size, 0);

      if (op[0].zstream.total_out <= op[1].zstream.total_out)
      {
        zsize = op[0].zstream.total_out;
        memcpy(zbuf, op[0].zbuf, zsize);
      }
      else
      {
        zsize = op[1].zstream.total_out;
        memcpy(zbuf, op[1].zbuf, zsize);
      }

      deflateReset(&op[0].zstream);
      op[0].zstream.data_type = Z_BINARY;
      deflateReset(&op[1].zstream);
      op[1].zstream.data_type = Z_BINARY;
    }

    for (i=first; i<frames-1; i++)
    {
      printf("saving frame %d of %d\n", i-first+2, frames-first);
      for (j=0; j<12; j++)
        op[j].valid = 0;

      /* dispose = none */
      try_over = get_rect(width, height, images[i], images[i+1], imagetemp, &x1, &y1, &w1, &h1, bpp);
      deflate_rect(images[i+1], x1, y1, w1, h1, bpp, rowbytes, zbuf_size, 0);
      if (try_over)
        deflate_rect(imagetemp, x1, y1, w1, h1, bpp, rowbytes, zbuf_size, 1);

      /* dispose = background */
      if (zero_trans)
      {
        memcpy(imagetemp, images[i], imagesize);
        for (j=0; j<h0; j++)
          memset(imagetemp + ((j+y0)*width + x0)*bpp, 0, w0*bpp);

        try_over = get_rect(width, height, imagetemp, images[i+1], imagetemp, &x1, &y1, &w1, &h1, bpp);

        deflate_rect(images[i+1], x1, y1, w1, h1, bpp, rowbytes, zbuf_size, 2);
        if (try_over)
          deflate_rect(imagetemp, x1, y1, w1, h1, bpp, rowbytes, zbuf_size, 3);
      }

      if (i>first)
      {
        /* dispose = previous */
        try_over = get_rect(width, height, images[i-1], images[i+1], imagetemp, &x1, &y1, &w1, &h1, bpp);
        deflate_rect(images[i+1], x1, y1, w1, h1, bpp, rowbytes, zbuf_size, 4);
        if (try_over)
          deflate_rect(imagetemp, x1, y1, w1, h1, bpp, rowbytes, zbuf_size, 5);
      }

      unsigned int op_min = op[0].zstream.total_out;
      int op_best = 0;
      for (j=1; j<12; j++)
      {
        if (op[j].valid)
        {
          if (op[j].zstream.total_out < op_min)
          {
            op_min = op[j].zstream.total_out;
            op_best = j;
          }
        }
      }

      dop = op_best >> 2;

      fctl.mSeq       = swap32(next_seq_num++);
      fctl.mWidth     = swap32(w0);
      fctl.mHeight    = swap32(h0);
      fctl.mXOffset   = swap32(x0);
      fctl.mYOffset   = swap32(y0);
      fctl.mDelayNum  = swap16(delay_num);
      fctl.mDelayDen  = swap16(delay_den);
      fctl.mDisposeOp = dop;
      fctl.mBlendOp   = bop;
      write_chunk(f, "fcTL", (unsigned char *)(&fctl), 26);

      write_IDATs(f, i, zbuf, zsize, idat_size);

      if (dop == 1)
      {
        for (j=0; j<h0; j++)
          memset(images[i] + ((j+y0)*width + x0)*bpp, 0, w0*bpp);
      }
      else
      if (dop == 2)
      {
        for (j=0; j<h0; j++)
          memcpy(images[i] + ((j+y0)*width + x0)*bpp, images[i-1] + ((j+y0)*width + x0)*bpp, w0*bpp);
      }

      x0 = op[op_best].x;
      y0 = op[op_best].y;
      w0 = op[op_best].w;
      h0 = op[op_best].h;
      bop = (op_best >> 1) & 1;

      zsize = op[op_best].zstream.total_out;
      memcpy(zbuf, op[op_best].zbuf, zsize);

      for (j=0; j<12; j++)
      {
        deflateReset(&op[j].zstream);
        op[j].zstream.data_type = Z_BINARY;
      }
    }

    if (frames > 1)
    {
      fctl.mSeq       = swap32(next_seq_num++);
      fctl.mWidth     = swap32(w0);
      fctl.mHeight    = swap32(h0);
      fctl.mXOffset   = swap32(x0);
      fctl.mYOffset   = swap32(y0);
      fctl.mDelayNum  = swap16(delay_num);
      fctl.mDelayDen  = swap16(delay_den);
      fctl.mDisposeOp = 0;
      fctl.mBlendOp   = bop;
      write_chunk(f, "fcTL", (unsigned char *)(&fctl), 26);
    }

    write_IDATs(f, i, zbuf, zsize, idat_size);

    write_chunk(f, "tEXt", png_Software, 27); 
    write_chunk(f, "IEND", 0, 0);
    fclose(f);
  }
  else
  {
    printf( "Error: couldn't open file for writing\n" );
    return 1;
  }

  free(imagetemp);
  free(zbuf);
  free(row_buf);
  free(sub_row);
  free(up_row);
  free(avg_row);
  free(paeth_row);

  for (i=0; i<12; i++)
  {
    deflateEnd(&op[i].zstream);
    if (op[i].zbuf != NULL)
      free(op[i].zbuf);
  }

  for (i=0; i<frames; i++)
  {
    if (images[i] != NULL)
      free(images[i]);
  }
  free(images);

  printf("all done\n");
    
  return 0;
}
