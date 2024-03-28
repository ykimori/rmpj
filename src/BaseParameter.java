/*
 * RMPJ
 * ImageJ plugin for rotational morphological processing.
 * Copyright (C) 2023 Yoshitaka Kimori
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.util.Arrays;


/**
 * Common data class
 *
 * @author Yoshitaka Kimori
 */
public class BaseParameter {
    
    private static final double Coef = 2.0;
    private static final double Scl = 1.0;
    static final int MaxImgSize = 2048;
    static final int MaxArraySize = 2048 * 2 + 99 * 2; 
    
    String filter;     // Types of morphological filter
    String seShape;    // Shape of the structural element
    String seSize;     // Size of the structural element
    String rdN;	       // Number of iterations
    
    int m_bitDepth;    // Original image types (8-bit / 16-bit)
    int m_upperLimit;  // Upper limit of the pixel value (255 / 65535)
    int[][] InImg;     // Original image
    int[][] StrEL;     // Structuring element
    
    int[][] OutImg	= new int[MaxImgSize][MaxImgSize];       // Processed image
    int[][] EntRot	= new int[MaxArraySize][MaxArraySize];   // Temporary image of the rotation process
    int[][] RotImg	= new int[MaxArraySize][MaxArraySize];	 // Rotated image
    int[][] Eron	= new int[MaxArraySize][MaxArraySize];	 // Eroded image
    int[][] Diln	= new int[MaxArraySize][MaxArraySize];	 // Dilated image
    
    int ImgSizeX, ImgSizeY;	      // Size of the long side of the original image
    int realImgSizeX, realImgSizeY;   // Size of the original image
    int offsetX, offsetY;	      // Offset values
    int HfImgSizeX, HfImgSizeY;	      // Half the size of the original image
    int EntRSizeX, EntRSizeY;	      // Size of the temporary image for the rotation process
    int HfEntRSizeX, HfEntRSizeY;     // Half the size of the temporary image for the rotation process
    int DX, DY;			      // Origin coordinates of the original image to create the temporary image for rotation
    int strX, trmX, strY, trmY;	      // Range to be scanned by the structural element
    int SESizeX, SESizeY;	      // Size of the structuring element
    int HfSESizeX, HfSESizeY;	      // Half the size of the structuring element
    
    private int m_curImg;	      // Number of the current image
    private int m_totalImg;	      // Total number of images
    
    
    /**
     *	Gets the number of the current image
     *
     *	@return	Number of the current image
     */
    public int CurImg() {
	return m_curImg;
    }
    
    /**
     *	Gets the total number of images
     *
     *	@return	Total number of images
     */
    public int TotalImg() {
	return m_totalImg;
    }
    
    /**
     *	Gets the size of the structuring element
     *
     *	@return	Size of the structuring element
     */
    public int getSESize() {
	return Integer.valueOf(this.seSize);
    }
    
    /**
     *	Gets the number of iterations 
     *
     *	@return	Number of iterations 
     */
    public int getRDN() {
	return Integer.valueOf(this.rdN);
    }
    
    /**
     *	Checks if the size of the structure element is odd
     *
     *	@return	true->OK/false->NG
     */
    public boolean isSESizeOdd() {
	int size = getSESize();
	if (size % 2 == 0) {
	    return false;
	}
	return true;
    }
    
    /**
     *	Checks whether the size of the structural element is within the specified range
     *
     *	@return	true->OK/false->NG
     */
    public boolean isSESizeValid() {
	int size = getSESize();
	if (size < 3 || size > 99) {
	    return false;
	}
	return true;
    }
    
    /**
     *	Checks if the number of iterations is valid
     *
     *	@return	true->OK/false->NG
     */
    public boolean isRDNValid() {
	int dire = getRDN();
	if (dire < 1) {
	    return false;
	}
	return true;
    }
    
    
    /**
     *	Checks whether the size of the original image is within the spcified size
     *
     *	@param	w	Original image width
     *	@param	h	Original image height
     *	@return	true->OK/false->NG
     */
    public boolean isImgSizeValid(int w, int h) {
	if (w > MaxImgSize || h > MaxImgSize) {
	    return false;
	}
	return true;
    }
    
    
    /**
     *  Initializes the common data
     *
     *	@param	ip		Input image
     *	@param	curImg		Number of the current image
     *	@param	totalImg	Number of the total images
     */
    public void initialize(ImageProcessor ip, int curImg, int totalImg) {
	m_curImg = curImg;
	m_totalImg = totalImg;
	
	LoadImg(ip);
	
	HfImgSizeX = ImgSizeX / 2;
	HfImgSizeY = ImgSizeY / 2;
	
	SetSE();
	EntRSizeX = (int)(2 * (HfImgSizeX * Coef) + 2 * SESizeX);
	EntRSizeY = (int)(2 * (HfImgSizeY * Coef) + 2 * SESizeY);
	
	SetParam();
    }
    
    private void LoadImg(ImageProcessor ip) {
	realImgSizeX = ip.getWidth();
	realImgSizeY = ip.getHeight();
	
	offsetX = 0;
	offsetY = 0;

	// If the width and height of the input image are not equal in size,
	// a square image is created and the input image is set in the center
	if (realImgSizeX != realImgSizeY) {
	    if (realImgSizeX > realImgSizeY) {
		ImgSizeX = realImgSizeX;
		ImgSizeY = realImgSizeX;
		
		offsetY = (ImgSizeY - realImgSizeY) / 2;
	    } else {
		ImgSizeX = realImgSizeY;
		ImgSizeY = realImgSizeY;
		
		offsetX = (ImgSizeX - realImgSizeX) / 2;
	    }
	} else {
	    ImgSizeX = realImgSizeX;
	    ImgSizeY = realImgSizeY;
	}
	
	// Set the background pixel value to zero
	InImg = new int[ImgSizeY][ImgSizeX];
	for (int y=0; y<ImgSizeY; y++) {
	    for (int x=0; x<ImgSizeX; x++) {
		InImg[y][x] = 0;
	    }
	}
	
	if (m_bitDepth == 8) {
	    m_upperLimit = 0xff;
	    
	    byte[] pixels = (byte[])ip.getPixels();
	    
	    for (int y=0; y<realImgSizeY; y++) {
		for (int x=0; x<realImgSizeX; x++) {
		    int pos = y * realImgSizeX + x;
		    int _y = y + offsetY;
		    int _x = x + offsetX;
		    InImg[_y][_x] = pixels[pos] & 0xff;
		}
	    }
	}
	
	if (m_bitDepth == 16) {
	    m_upperLimit = 0xffff;
	    
	    short[] pixels = (short[])ip.getPixels();
	    
	    for (int y=0; y<realImgSizeY; y++) {
		for (int x=0; x<realImgSizeX; x++) {
		    int pos = y * realImgSizeX + x;
		    int _y = y + offsetY;
		    int _x = x + offsetX;
		    InImg[_y][_x] = pixels[pos] & 0xffff;
		}
	    }
	}
    }
    
    private void SetParam() {
	HfEntRSizeX = EntRSizeX / 2;
	HfEntRSizeY = EntRSizeY / 2;
	
	DX = HfEntRSizeX - HfImgSizeX;
	DY = HfEntRSizeY - HfImgSizeY;
	
	int Radi;
	if (EntRSizeY<=EntRSizeX)
	    Radi = (int)(HfEntRSizeX * Scl);
	else
	    Radi = (int)(HfEntRSizeY * Scl);
	
	strX = HfEntRSizeX - Radi;
	trmX = HfEntRSizeX + Radi;
	strY = HfEntRSizeY - Radi;
	trmY = HfEntRSizeY + Radi;
    }
    
    public void SetSE() {
	if (seShape == "Line") {
	    SESizeX = getSESize();
	    SESizeY = 1;
	}
	else { // "Square" or "Disk"
	    SESizeX = getSESize();
	    SESizeY = getSESize();
	}
	
	HfSESizeX = SESizeX/2;
	HfSESizeY = SESizeY/2;
	
	StrEL = new int[SESizeY][SESizeX];
	for (int y=0; y<SESizeY; y++) {
	    for (int x=0; x<SESizeX; x++) {
		StrEL[y][x] = 0;
	    }
	}
	
	if (this.seShape == "Disk") {
	    initDisk();
	}
    }
    
    /**
     *	Initializes the disk-shaped structuring element
     */
    private void initDisk() {
	
	// Create a disk-shaped structuring element using the Pencil Tool
	int size = getSESize();
	byte [] pixels = new byte[size * size];
	Arrays.fill(pixels, (byte)-1);
	ImageProcessor ip = new ByteProcessor(size, size, pixels);
	
	int drawPos = size / 2;
	ip.setLineWidth(size);
	ip.setColor(0);
	ip.drawDot(drawPos, drawPos);
	
	for (int y = 0; y < size; ++y) {
	    for (int x = 0; x < size; ++x) {
		int index = y * size + x;
		int pixel = pixels[index];
		StrEL[y][x] = pixel;
	    }
	}
    }
}
