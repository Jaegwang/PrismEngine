#pragma once

#include <GL/glew.h>
#include <stdlib.h>
#include <iostream>
#include <png.h>
#include <zlib.h>
#include "MATH_CORE.h"

class READ_WRITE_PNG
{
public:
	static void png_read_data(png_structp png_ptr, png_bytep data, png_size_t length)
	{
	   png_size_t check;

	   /* fread() returns 0 on error, so it is OK to store this in a png_size_t
		* instead of an int, which is what fread() actually returns.
		*/
	   check = (png_size_t)fread(data, (png_size_t)1, length,
		   (FILE *) png_get_io_ptr(png_ptr));

	   if (check != length)
	   {
		  png_error(png_ptr, "Read Error");
	   }
	}

	static void png_write_data(png_structp png_ptr, png_bytep data, png_size_t length)
	{
	   png_uint_32 check;

	   check = (png_uint_32)fwrite(data, 1, length, (FILE *) png_get_io_ptr(png_ptr));
	   if (check != length)
	   {
		  png_error(png_ptr, "Write Error");
	   }
	}

	static void png_flush(png_structp png_ptr)
	{
	   FILE *io_ptr;
	   io_ptr = (FILE *)png_get_io_ptr(png_ptr);
	   if (io_ptr != NULL)
		  fflush(io_ptr);
	}

	static void png_cexcept_error(png_structp png_ptr, png_const_charp msg)
	{
	   if(png_ptr)
	#ifndef PNG_NO_CONSOLE_IO
	   fprintf(stderr, "libpng error: %s\n", msg);
	#endif
	   {
		  throw msg;
	   }
	}

	// PNG image handler functions
	static bool PngLoadImage (char *pstrFileName, png_byte **ppbImageData,
		png_uint_32 *piWidth, png_uint_32 *piHeight, int *piChannels, png_color *pBkgColor)
	{
		png_structp png_ptr = NULL;
		png_infop info_ptr = NULL;

		static FILE        *pfFile;
		png_byte            pbSig[8];
		int                 iBitDepth;
		int                 iColorType;
		double              dGamma;
		png_color_16       *pBackground;
		png_uint_32         ulChannels;
		png_uint_32         ulRowBytes;
		png_byte           *pbImageData = *ppbImageData;
		static png_byte   **ppbRowPointers = NULL;
		int                 i;


		// open the PNG input file

		if (!pstrFileName)
		{
			*ppbImageData = pbImageData = NULL;
			return false;
		}

		if (!(pfFile = fopen(pstrFileName, "rb")))
		{
			*ppbImageData = pbImageData = NULL;
			return false;
		}

		// first check the eight byte PNG signature

		fread(pbSig, 1, 8, pfFile);
		if (!png_check_sig(pbSig, 8))
		{
			*ppbImageData = pbImageData = NULL;
			return false;
		}

		// create the two png(-info) structures

		png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL,
		  (png_error_ptr)png_cexcept_error, (png_error_ptr)NULL);
		if (!png_ptr)
		{
			*ppbImageData = pbImageData = NULL;
			return false;
		}

		info_ptr = png_create_info_struct(png_ptr);
		if (!info_ptr)
		{
			png_destroy_read_struct(&png_ptr, NULL, NULL);
			*ppbImageData = pbImageData = NULL;
			return false;
		}

		try
		{
        
			// initialize the png structure
        
	#if !defined(PNG_NO_STDIO)
			png_init_io(png_ptr, pfFile);
	#else
			png_set_read_fn(png_ptr, (png_voidp)pfFile, png_read_data);
	#endif
        
			png_set_sig_bytes(png_ptr, 8);
        
			// read all PNG info up to image data
        
			png_read_info(png_ptr, info_ptr);
        
			// get width, height, bit-depth and color-type
        
			png_get_IHDR(png_ptr, info_ptr, piWidth, piHeight, &iBitDepth,
				&iColorType, NULL, NULL, NULL);
        
			// expand images of all color-type and bit-dept to 3x8 bit RGB images
			// let the library process things like alpha, transparency, background
        
			if (iBitDepth == 16)
				png_set_strip_16(png_ptr);
			if (iColorType == PNG_COLOR_TYPE_PALETTE)
				png_set_expand(png_ptr);
			if (iBitDepth < 8)
				png_set_expand(png_ptr);
			if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
				png_set_expand(png_ptr);
			if (iColorType == PNG_COLOR_TYPE_GRAY ||
				iColorType == PNG_COLOR_TYPE_GRAY_ALPHA)
				png_set_gray_to_rgb(png_ptr);
        
			// set the background color to draw transparent and alpha images over.
			if (png_get_bKGD(png_ptr, info_ptr, &pBackground))
			{
				png_set_background(png_ptr, pBackground, PNG_BACKGROUND_GAMMA_FILE, 1, 1.0);
				pBkgColor->red   = (unsigned char) pBackground->red;
				pBkgColor->green = (unsigned char) pBackground->green;
				pBkgColor->blue  = (unsigned char) pBackground->blue;
			}
			else
			{
				pBkgColor = NULL;
			}
        
			// if required set gamma conversion
			if (png_get_gAMA(png_ptr, info_ptr, &dGamma))
				png_set_gamma(png_ptr, (double) 2.2, dGamma);
        
			// after the transformations have been registered update info_ptr data
        
			png_read_update_info(png_ptr, info_ptr);
        
			// get again width, height and the new bit-depth and color-type
        
			png_get_IHDR(png_ptr, info_ptr, piWidth, piHeight, &iBitDepth,
				&iColorType, NULL, NULL, NULL);
        
        
			// row_bytes is the width x number of channels
        
			ulRowBytes = (png_uint_32)png_get_rowbytes(png_ptr, info_ptr);
			ulChannels = (png_uint_32)png_get_channels(png_ptr, info_ptr);
        
			*piChannels = ulChannels;
        
			// now we can allocate memory to store the image
        
			if (pbImageData)
			{
				free (pbImageData);
				pbImageData = NULL;
			}
			if ((pbImageData = (png_byte *) malloc(ulRowBytes * (*piHeight)
								* sizeof(png_byte))) == NULL)
			{
				png_error(png_ptr, "Visual PNG: out of memory");
			}
			*ppbImageData = pbImageData;
        
			// and allocate memory for an array of row-pointers
        
			if ((ppbRowPointers = (png_bytepp) malloc((*piHeight)
								* sizeof(png_bytep))) == NULL)
			{
				png_error(png_ptr, "Visual PNG: out of memory");
			}
        
			// set the individual row-pointers to point at the correct offsets
        
			for (i = 0; i < (int)(*piHeight); i++)
				ppbRowPointers[i] = pbImageData + i * ulRowBytes;
        
			// now we can go ahead and just read the whole image
        
			png_read_image(png_ptr, ppbRowPointers);
        
			// read the additional chunks in the PNG file (not really needed)
        
			png_read_end(png_ptr, NULL);
        
			// and we're done
        
			free (ppbRowPointers);
			ppbRowPointers = NULL;
        
			// yepp, done
		}

		catch (png_const_charp msg)
		{
			png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

			*ppbImageData = pbImageData = NULL;
        
			if(ppbRowPointers)
				free (ppbRowPointers);

			fclose(pfFile);
			std::cout<< msg << std::endl;
			return false;
		}

		fclose (pfFile);

		return true;
	}

	static bool PngSaveImage (char *pstrFileName,  int width, int height)
	{
		png_structp png_ptr = NULL;
		png_infop info_ptr = NULL;

		const int           ciBitDepth = 8;
		const int           ciChannels = 3;

		 FILE        *pfFile;
		png_uint_32         ulRowBytes=0;
		static png_byte   **ppbRowPointers = NULL;
		int                 i;
		png_byte *pDiData;
		// open the PNG output file

		if (!pstrFileName)
			return false;
		 if (!(pfFile = fopen(pstrFileName, "wb")))
			return false;

		// prepare the standard PNG structures

		png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL,
		  (png_error_ptr)png_cexcept_error, (png_error_ptr)NULL);
		if (!png_ptr)
		{
			fclose(pfFile);
			return false;
		}

		info_ptr = png_create_info_struct(png_ptr);
		if (!info_ptr) {
			fclose(pfFile);
			png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
			return false;
		}

		try
		{	
			// initialize the png structure
        
	#if defined(PNG_NO_STDIO)
			png_init_io(png_ptr, pfFile);
	#else
			png_set_write_fn(png_ptr, (png_voidp)pfFile, png_write_data, png_flush);
	#endif

			// we're going to write a very simple 3x8 bit RGB image
        
			png_set_IHDR(png_ptr, info_ptr, width, height, ciBitDepth,
				PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
				PNG_FILTER_TYPE_BASE);

				png_set_compression_level(png_ptr,Z_DEFAULT_COMPRESSION);

				pDiData = (png_byte *) malloc(width * height * sizeof(GLubyte) * 3);
				glPixelStorei(GL_PACK_ALIGNMENT, 1);
				glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pDiData);


			// write the file header information
			png_write_info(png_ptr, info_ptr);
        
			// swap the BGR pixels in the DiData structure to RGB
			//png_set_bgr(png_ptr);
        
			// row_bytes is the width x number of channels
			ulRowBytes = width * ciChannels;
        
			// we can allocate memory for an array of row-pointers
        
			if ((ppbRowPointers = (png_bytepp) malloc(height * sizeof(png_bytep))) == NULL)
				throw "Visualpng: Out of memory";
        
			// set the individual row-pointers to point at the correct offsets
			// swap top and bottom
		
			int j = 0;
			for (i = height-1; i >= 0; i--)
			{	
				ppbRowPointers[j] = pDiData + i * (((ulRowBytes + 3) >> 2) << 2);
				j++;
			}
			// write out the entire image data in one call
        
			png_write_image (png_ptr, ppbRowPointers);
        
			// write the additional chunks to the PNG file (not really needed)
        
			png_write_end(png_ptr, info_ptr);
        
			// and we're done
        
			free (ppbRowPointers);
			ppbRowPointers = NULL;
        
			// clean up after the write, and free any memory allocated
        
			png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        
			// yepp, done
		}

		catch (png_const_charp msg)
		{
			png_destroy_write_struct(&png_ptr, (png_infopp) NULL);

			if(ppbRowPointers)
				free (ppbRowPointers);

			fclose(pfFile);

			std::cout<< msg << std::endl;
			return false;
		}
    
		fclose (pfFile);
		
		delete pDiData;

		return true;
	}


	static bool PngSaveImage (char *pstrFileName, const int i_res, const int j_res, const Vec3* data)
	{
		png_structp png_ptr = NULL;
		png_infop info_ptr = NULL;

		int width = i_res;
		int height = j_res;

		const int           ciBitDepth = 8;
		const int           ciChannels = 3;

		FILE        *pfFile;
		png_uint_32         ulRowBytes=0;
		static png_byte   **ppbRowPointers = NULL;
		//int                 i;
		png_byte *pDiData;
		// open the PNG output file

		if (!pstrFileName)
			return false;
		 if (!(pfFile = fopen(pstrFileName, "wb")))
			return false;

		// prepare the standard PNG structures

		png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL,
		  (png_error_ptr)png_cexcept_error, (png_error_ptr)NULL);
		if (!png_ptr)
		{
			fclose(pfFile);
			return false;
		}

		info_ptr = png_create_info_struct(png_ptr);
		if (!info_ptr) {
			fclose(pfFile);
			png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
			return false;
		}

		try
		{	
			// initialize the png structure
        
	#if defined(PNG_NO_STDIO)
			png_init_io(png_ptr, pfFile);
	#else
			png_set_write_fn(png_ptr, (png_voidp)pfFile, png_write_data, png_flush);
	#endif

			// we're going to write a very simple 3x8 bit RGB image
        
			png_set_IHDR(png_ptr, info_ptr, width, height, ciBitDepth,
				PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
				PNG_FILTER_TYPE_BASE);

				png_set_compression_level(png_ptr,Z_DEFAULT_COMPRESSION);

				pDiData = (png_byte *) malloc(width * height * sizeof(GLubyte) * 3);
			
				int count =0;
				for(int i=0;i<width;i++)
				{
					for(int j=0;j<height;j++)
					{
						pDiData[i*height*3+j*3+0]=(unsigned char)data[i*height+j].x*255;
						pDiData[i*height*3+j*3+1]=(unsigned char)data[i*height+j].y*255;
						pDiData[i*height*3+j*3+2]=(unsigned char)data[i*height+j].z*255;
					}
				}

				//glPixelStorei(GL_PACK_ALIGNMENT, 1);
				//glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pDiData);


			// write the file header information
			png_write_info(png_ptr, info_ptr);
        
			// swap the BGR pixels in the DiData structure to RGB
			//png_set_bgr(png_ptr);
        
			// row_bytes is the width x number of channels
			ulRowBytes = width * ciChannels;
        
			// we can allocate memory for an array of row-pointers
        
			if ((ppbRowPointers = (png_bytepp) malloc(height * sizeof(png_bytep))) == NULL)
				throw "Visualpng: Out of memory";
        
			// set the individual row-pointers to point at the correct offsets
			// swap top and bottom
		
			int j = 0;
			for (int i = height-1; i >= 0; i--)
			{	
				ppbRowPointers[j] = pDiData + i * (((ulRowBytes + 3) >> 2) << 2);
				j++;
			}
			// write out the entire image data in one call
        
			png_write_image (png_ptr, ppbRowPointers);
        
			// write the additional chunks to the PNG file (not really needed)
        
			png_write_end(png_ptr, info_ptr);
        
			// and we're done
        	ppbRowPointers = NULL;
        
			// clean up after the write, and free any memory allocated       
			png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
			free (ppbRowPointers);
		}

		catch (png_const_charp msg)
		{
			png_destroy_write_struct(&png_ptr, (png_infopp) NULL);

			if(ppbRowPointers)
				free (ppbRowPointers);

			fclose(pfFile);

			std::cout<< msg << std::endl;
			return false;
		}
    
		fclose (pfFile);
    
		return true;
	}


	static void savePng(int width, int height)
	{
		png_structp png_ptr;
		png_infop info_ptr;
		GLubyte *image;
		png_byte **row_pointers=NULL;
		png_ptr=NULL;info_ptr=NULL;

	
		png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL,NULL,NULL);
		if (!png_ptr){
			printf("Couldn't allocate memory for PNG file");
			return;
		}
		info_ptr= png_create_info_struct(png_ptr);
		if (!info_ptr){
			printf("Couldn't allocate image information for PNG file");
			return;
		}
		FILE *fp;
		errno_t err = fopen_s(&fp, "screenshot.png", "wb");
		if (!fp)
		{
			printf("screenshot.png file Open Failed");
		   return;
		}

		png_init_io(png_ptr, fp);

	//	png_set_write_status_fn(png_ptr, write_row_callback);

		// use default compression mode
		png_set_compression_level(png_ptr,Z_DEFAULT_COMPRESSION);
  
		png_set_IHDR(png_ptr,info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB,
					PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);


		image = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);

		row_pointers=(png_byte **)malloc(height * sizeof(png_byte*));
		int j=0;
		for (int i=height-1; i>=0; i--)
		{
			row_pointers[j++] = ((png_byte*)image) + i*width*3;
		}

		png_write_info(png_ptr, info_ptr);
		png_write_image(png_ptr, row_pointers);
		png_write_end(png_ptr, NULL);

		png_destroy_write_struct(&png_ptr,&info_ptr);
		free(row_pointers);
		free(image);
		fclose(fp);
	}
};

