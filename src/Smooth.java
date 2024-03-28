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

import java.util.Collection;
import java.util.LinkedList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ij.*;
import ij.process.*;
import ij.gui.*;


/**
 * This class implements a smoothing filter by rotational morphological processing
 *
 * @author Yoshitaka Kimori
 */
public class Smooth {
    
    private static final int MaxImgSize = 2048;
    
    private int[][] OpenImg	= new int[MaxImgSize][MaxImgSize];
    private int[][] ClosnImg	= new int[MaxImgSize][MaxImgSize];
    private int[][] OCImg	= new int[MaxImgSize][MaxImgSize];
    private int[][] COImg	= new int[MaxImgSize][MaxImgSize];
    
    private int currentIndex, finalIndex;
    
    
    /**
     *	Displays the process progress
     *
     *	@param	param	Common data
     */
    private void progress(BaseParameter param) {
	String msg = "Smoothing... " + param.CurImg() + "/" + param.TotalImg();
	IJ.showStatus(msg);
	IJ.showProgress(currentIndex, finalIndex);
	currentIndex++;
    }
    
    
    public void run(BaseParameter param) {
	currentIndex = 0;
	finalIndex = param.getRDN() * 4 + 3; 
	progress(param);
	
	RMPOpenClos(param);
	MakeOCImg(param);
	
	RMPClosOpen(param);
	MakeCOImg(param);
	
	MakeOCCOImg(param);
    }
    
    
    private void RMPOpenClos(BaseParameter param) {
	
	int rdN = param.getRDN();
	float Deg = (float)(180.0 / rdN);
	
	int [][][] IntMedImg = new int[rdN][param.ImgSizeY][param.ImgSizeX];
	
	// Opening by RMP
	for (int i=0; i<rdN; i++) {
	    MakeEnterRotImg(param);
	    Rotation(Deg * i, param);
	    
	    Opening(param);
	    
	    for (int y=0; y<param.EntRSizeY; y++) {
		for (int x=0; x<param.EntRSizeX; x++) {
		    param.EntRot[y][x] = param.Diln[y][x];
		}
	    }
	    
	    Rotation(-Deg * i, param);
	    ExtIntMedImg(IntMedImg[i], param);
	    
	    progress(param);
	}
	
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		OpenImg[y][x] = IntMedImg[0][y][x];
	    }
	}
	
	for (int i=1; i<rdN; i++) {
	    CalcMax(IntMedImg[i], param);
	}

        AFRMVO1(param);
	
	// Closing by RMP
	for (int i=0; i<rdN; i++) {
	    CopyOpenImg(param);
	    Rotation(Deg * i, param);
	    
	    Closing(param);
	    
	    for (int y=0; y<param.EntRSizeY; y++) {
		for (int x=0; x<param.EntRSizeX; x++) {
		    param.EntRot[y][x] = param.Eron[y][x];
		}
	    }
	    
	    Rotation(-Deg * i, param);
	    ExtIntMedImg(IntMedImg[i], param);
	    
	    progress(param);
	}
	
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		ClosnImg[y][x] = IntMedImg[0][y][x];
	    }
	}
	
	for (int i=1; i<rdN; i++) {
	    CalcMin(IntMedImg[i], param);
	}

	AFRMVC2(param);
    }
    
    
    private void MakeOCImg(BaseParameter param) {
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		OCImg[y][x] = ClosnImg[y][x];
	    }
	}
	
	progress(param);
    }
    
    
    private void RMPClosOpen(BaseParameter param) {
	
	int rdN = param.getRDN();
	float Deg = (float)(180.0 / rdN);
	
	int [][][] IntMedImg = new int[rdN][param.ImgSizeY][param.ImgSizeX];
	
	// Closing by RMP
	for (int i=0; i<rdN; i++) {
	    MakeEnterRotImg(param);
	    Rotation(Deg * i, param);
	    
	    Closing(param);
	    
	    for (int y=0; y<param.EntRSizeY; y++) {
		for (int x=0; x<param.EntRSizeX; x++) {
		    param.EntRot[y][x] = param.Eron[y][x];
		}
	    }
	    
	    Rotation(-Deg * i, param);
	    ExtIntMedImg(IntMedImg[i], param);
	    
	    progress(param);
	}
	
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		ClosnImg[y][x] = IntMedImg[0][y][x];
	    }
	}
	
	for (int i=1; i<rdN; i++) {
	    CalcMin(IntMedImg[i], param);
	}

	AFRMVC1(param);
	
	// Opening by RMP
	for (int i=0; i<rdN; i++) {
	    CopyClosnImg(param);
	    Rotation(Deg * i, param);
	    
	    Opening(param);
	    
	    for (int y=0; y<param.EntRSizeY; y++) {
		for (int x=0; x<param.EntRSizeX; x++) {
		    param.EntRot[y][x] = param.Diln[y][x];
		}
	    }
	    
	    Rotation(-Deg * i, param);
	    ExtIntMedImg(IntMedImg[i], param);
	    
	    progress(param);
	}
	
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		OpenImg[y][x] = IntMedImg[0][y][x];
	    }
	}
	
	for (int i=1; i<rdN; i++) {
	    CalcMax(IntMedImg[i], param);
	}

	AFRMVO2(param);
    }
    
    
    private void MakeCOImg(BaseParameter param) {
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		COImg[y][x] = OpenImg[y][x];
	    }
	}
	
	progress(param);
    }
    
    
    private void MakeOCCOImg(BaseParameter param) {
	for (int y=0; y<param.ImgSizeY; y++){
	    for (int x=0; x<param.ImgSizeX; x++) {
		int gv = OCImg[y][x] + COImg[y][x];
		int avg = gv / 2;
		param.OutImg[y][x] = avg;
	    }
	}
	
	progress(param);
    }
    
    
    private void MakeEnterRotImg(BaseParameter param) {
	for (int y=0; y<param.EntRSizeY; y++) {
	    for (int x=0; x<param.EntRSizeX; x++) {
		param.EntRot[y][x] = 0;
	    }
	}
	
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		param.EntRot[y + param.DY][x + param.DX] = param.InImg[y][x];
	    }
	}
	
	int sizeY = param.DY+param.ImgSizeY;
	int sizeX = param.DX+param.ImgSizeX;

	// Perform padding processing
	// Left direction
	for (int y=param.DY, sy=0; y<sizeY; ++y, ++sy) {
	    for (int x=0; x<param.DX; ++x) {
		param.EntRot[y][x] = param.InImg[sy][0];
	    }
	}
	// Diagonally upper left
	for (int y=0; y<param.DY; ++y) {
	    for (int x=0; x<param.DX; ++x) {
		param.EntRot[y][x] = param.InImg[0][0];
	    }
	}
	// Upward direction
	for (int y=0; y<param.DY; ++y) {
	    for (int x=param.DX, sx=0; x<sizeX; ++x, ++sx) {
		param.EntRot[y][x] = param.InImg[0][sx];
	    }
	}
	// Diagonally upper right
	for (int y=0; y<param.DY; ++y) {
	    for (int x=sizeX; x<param.EntRSizeX; ++x) {
		param.EntRot[y][x] = param.InImg[0][param.ImgSizeX-1];
	    }
	}
	// Right direction
	for (int y=param.DY, sy=0; y<sizeY; ++y, ++sy) {
	    for (int x=sizeX; x<param.EntRSizeX; ++x) {
		param.EntRot[y][x] = param.InImg[sy][param.ImgSizeX-1];
	    }
	}
	// Diagonally lower right
	for (int y=sizeY; y<param.EntRSizeY; ++y) {
	    for (int x=sizeX; x<param.EntRSizeX; ++x) {
		param.EntRot[y][x] = param.InImg[param.ImgSizeY-1][param.ImgSizeX-1];
	    }
	}
	// Downward direction
	for (int y=sizeY; y<param.EntRSizeY; ++y) {
	    for (int x=param.DX, sx=0; x<sizeX; ++x, ++sx) {
		param.EntRot[y][x] = param.InImg[param.ImgSizeY-1][sx];
	    }
	}
	// Diagonally lower left
	for (int y=sizeY; y<param.EntRSizeY; ++y) {
	    for (int x=0; x<param.DX; ++x) {
		param.EntRot[y][x] = param.InImg[param.ImgSizeY-1][0];
	    }
	}
    }
    

    // Image rotation using bilinear interpolation
    private void Rotation(float Deg, BaseParameter param) {
	
	for (int i=0; i<param.EntRSizeY; i++) {
	    for (int j=0; j<param.EntRSizeX; j++) {
		param.RotImg[i][j] = 0;
	    }
	}
	
	int Xs = param.HfEntRSizeX;
	int Ys = param.HfEntRSizeY;
	double r = Deg * 3.141592 / 180.0;
	float c = (float)Math.cos(r);
	float s = (float)Math.sin(r);
	
	for (int i=-Ys; i<=Ys; i++) {
	    for(int j=-Xs; j<=Xs; j++) {
		
		float x= j*c + i*s;
		float y= -j*s + i*c;
		
		int m, n, d;
		double cv;
		float p, q;
		
		if (y>0)
		    m = (int)y;
		else
		    m = (int)(y-1);
		if (x>0)
		    n = (int)x;
		else
		    n = (int)(x-1);
		
		q = y-m;
		p = x-n;

		if((m>=-Ys) && (m<Ys) && (n>=-Xs) && (n<Xs)){
		    cv = (double)((1.0-q)*((1.0-p)*param.EntRot[m+Ys][n+Xs]
					 +p*param.EntRot[m+Ys][n+1+Xs])
				+q*((1.0-p)*param.EntRot[m+1+Ys][n+Xs]
				    +p*param.EntRot[m+1+Ys][n+1+Xs]));
		    
		    d = (int)(cv+0.5);
		}
		else
		    d = 0;
		
		if (d < 0) d = 0;
		if (d > param.m_upperLimit) d = param.m_upperLimit;
		
		param.RotImg[i+Ys][j+Xs] = d;
	    }
	}
    }
    
    
    private void Opening(BaseParameter param) {
	// Instance for parallel-processing
	int core = Runtime.getRuntime().availableProcessors();
	ExecutorService threadPool = Executors.newFixedThreadPool(core);
	Collection<Callable<Void>> processes = new LinkedList<Callable<Void>>();
	
	// Parallel-processing parameters
	final int fstarty = param.strY+param.HfSESizeY;
	final int fendy = param.trmY-param.HfSESizeY;
	final int fstartx = param.strX+param.HfSESizeX;
	final int fendx = param.trmX-param.HfSESizeX;
	final int fhsy = param.HfSESizeY;
	final int fhsx = param.HfSESizeX;
	final int[][] strEL = param.StrEL;
	final int[][] rotImg = param.RotImg;
	final int[][] eron = param.Eron;
	final int[][] diln = param.Diln;
	
	try {
	    // Erosion
	    for (int i=0; i<param.EntRSizeY; i++) {
		for (int j=0; j<param.EntRSizeX; j++) {
		    param.Eron[i][j]=0;
		}
	    }
	    
	    for (int y=fstarty; y<fendy; y++) {
		final int fy = y;
		
		processes.add(new Callable<Void>() {
			public Void call() {
			    if (fy < 0)	return null;
			    for (int x=fstartx; x<fendx; x++) {
				if (x < 0)	continue;
				
				int min = Integer.MAX_VALUE;
				
				for (int r=-fhsy; r<=fhsy; r++) {
				    for (int c=-fhsx; c<=fhsx; c++) {
					if (strEL[r+fhsy][c+fhsx] != -1) {
					    int yn = fy + r;
					    int xn = x + c;
					    if (yn < 0 || xn < 0)	continue;
					    int gv = rotImg[yn][xn] - strEL[r+fhsy][c+fhsx];
					    
					    if (gv < min)
						min = gv;
					}
				    }
				}
				
				eron[fy][x] = min;
			    }
			    return null;
			}
		    });
	    }
	    threadPool.invokeAll(processes);
	    processes.clear();
	    
	    
	    // Dilation
	    for (int i=0; i<param.EntRSizeY; i++) {
		for (int j=0; j<param.EntRSizeX; j++) {
		    param.Diln[i][j] = 0;
		}
	    }
	    
	    for (int y=fstarty; y<fendy; y++) {
		final int fy = y;
		
		processes.add(new Callable<Void>() {
			public Void call() {
			    if (fy < 0)	return null;
			    for (int x=fstartx; x<fendx; x++) {
				if (x < 0)	continue;
				
				int max = Integer.MIN_VALUE;
				
				for (int r=-fhsy; r<=fhsy; r++) {
				    for (int c=-fhsx; c<=fhsx; c++) {
					if (strEL[r+fhsy][c+fhsx] != -1) {
					    int yn = fy - r;
					    int xn = x - c;
					    if (yn < 0 || xn < 0)	continue;
					    int gv = eron[yn][xn] + strEL[r+fhsy][c+fhsx];
					    
					    if (gv > max)
						max = gv;
					}
				    }
				}
				
				diln[fy][x] = max;
			    }
			    return null;
			}
		    });
	    }
	    threadPool.invokeAll(processes);
	    
	} catch (Exception e) {
	    throw new RuntimeException(e);
	} finally {
	    threadPool.shutdown();
	}
    }
    
    
    private void Closing(BaseParameter param) {
	// Instance for parallel-processing
	int core = Runtime.getRuntime().availableProcessors();
	ExecutorService threadPool = Executors.newFixedThreadPool(core);
	Collection<Callable<Void>> processes = new LinkedList<Callable<Void>>();
	
	// Parallel-processing parameters
	final int fstarty = param.strY+param.HfSESizeY;
	final int fendy = param.trmY-param.HfSESizeY;
	final int fstartx = param.strX+param.HfSESizeX;
	final int fendx = param.trmX-param.HfSESizeX;
	final int fhsy = param.HfSESizeY;
	final int fhsx = param.HfSESizeX;
	final int[][] strEL = param.StrEL;
	final int[][] rotImg = param.RotImg;
	final int[][] eron = param.Eron;
	final int[][] diln = param.Diln;
	
	try {
	    // Dilation
	    for (int i=0; i<param.EntRSizeY; i++) {
		for (int j=0; j<param.EntRSizeX; j++) {
		    param.Diln[i][j] = 0;
		}
	    }
	    
	    for (int y=fstarty; y<fendy; y++) {
		final int fy = y;
		
		processes.add(new Callable<Void>() {
			public Void call() {
			    if (fy < 0)	return null;
			    for (int x=fstartx; x<fendx; x++) {
				if (x < 0)	continue;
				
				int max = Integer.MIN_VALUE;
				
				for (int r=-fhsy; r<=fhsy; r++) {
				    for (int c=-fhsx; c<=fhsx; c++) {
					if (strEL[r+fhsy][c+fhsx] != -1) {
					    int yn = fy - r;
					    int xn = x - c;
					    if (yn < 0 || xn < 0)	continue;
					    int gv = rotImg[yn][xn] + strEL[r+fhsy][c+fhsx];
					    
					    if (gv > max)
						max = gv;
					}
				    }
				}
				
				diln[fy][x] = max;
			    }
			    return null;
			}
		    });
	    }
	    threadPool.invokeAll(processes);
	    processes.clear();
	    
	    
	    // Erosion
	    for (int i=0; i<param.EntRSizeY; i++) {
		for (int j=0; j<param.EntRSizeX; j++) {
		    param.Eron[i][j] = 0;
		}
	    }
	    
	    for (int y=fstarty; y<fendy; y++) {
		final int fy = y;
		
		processes.add(new Callable<Void>() {
			public Void call() {
			    if (fy < 0)	return null;
			    for (int x=fstartx; x<fendx; x++) {
				if (x < 0)	continue;
				
				int min = Integer.MAX_VALUE;
				
				for (int r=-fhsy; r<=fhsy; r++){
				    for (int c=-fhsx; c<=fhsx; c++){
					if (strEL[r+fhsy][c+fhsx] != -1) {
					    int yn = fy + r;
					    int xn = x + c;
					    if (yn < 0 || xn < 0)	continue;
					    int gv = diln[yn][xn] - strEL[r+fhsy][c+fhsx];
					    
					    if (gv < min)
						min = gv;
					}
				    }
				}
				
				eron[fy][x] = min;
			    }
			    return null;
			}
		    });
	    }
	    threadPool.invokeAll(processes);
	    
	} catch (Exception e) {
	    throw new RuntimeException(e);
	} finally {
	    threadPool.shutdown();
	}
    }
    
    
    private void ExtIntMedImg(int[][] intMedImg, BaseParameter param) {
	for (int y=0; y<param.ImgSizeY; y++) {
	    for(int x=0; x<param.ImgSizeX; x++) {
		intMedImg[y][x] = param.RotImg[y+param.DY][x+param.DX];
	    }
	}
    }
    

    private void CalcMax(int[][] intMedImg, BaseParameter param) {
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		if (OpenImg[y][x] <= intMedImg[y][x]) {
		    OpenImg[y][x] = intMedImg[y][x];
		}
	    }
	}
    }

    
     private void AFRMVO1(BaseParameter param){
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		if(OpenImg[y][x] > param.InImg[y][x])
		    OpenImg[y][x] = param.InImg[y][x];
	    }
	}
    }

    
    private void AFRMVO2(BaseParameter param){
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		if(OpenImg[y][x] > ClosnImg[y][x])
		    OpenImg[y][x] = ClosnImg[y][x];
	    }
	}
    }

    
    private void CalcMin(int[][] intMedImg, BaseParameter param) {
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		if (ClosnImg[y][x] >= intMedImg[y][x]) {
		    ClosnImg[y][x] = intMedImg[y][x];
		}
	    }
	}
    }


    private void AFRMVC1(BaseParameter param){
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		if(ClosnImg[y][x] < param.InImg[y][x])
		    ClosnImg[y][x] = param.InImg[y][x];
	    }
	}
    }
    

    private void AFRMVC2(BaseParameter param){
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		if(ClosnImg[y][x] < OpenImg[y][x])
		    ClosnImg[y][x] = OpenImg[y][x];
	    }
	}
    }
    
    
    private void CopyOpenImg(BaseParameter param) {
	for (int y=0; y<param.EntRSizeY; y++) {
	    for (int x=0; x<param.EntRSizeX; x++) {
		param.EntRot[y][x] = 0;
	    }
	}
	
	for (int y=0; y<param.ImgSizeY; y++) {
	    for (int x=0; x<param.ImgSizeX; x++) {
		param.EntRot[y+param.DY][x+param.DX] = OpenImg[y][x];
	    }
	}
	
	int sizeY = param.DY+param.ImgSizeY;
	// Mirroring (left)
	for (int y=param.DY, sy=0; y<sizeY; y++, sy++) {
	    for (int x=param.DX-1, sx=0; x>=0 && sx<param.ImgSizeX; x--, sx++) {
		param.EntRot[y][x] = OpenImg[sy][sx];
	    }
	}
	// Mirroring (right))
	for (int y=param.DY, sy=0; y<sizeY; y++, sy++) {
	    for (int x=param.DX+param.ImgSizeX, sx=param.ImgSizeX-1; x<param.EntRSizeX && sx>=0; x++, sx--) {
		param.EntRot[y][x] = OpenImg[sy][sx];
	    }
	}
	// Mirroring (top)
	for (int y=param.DY-1, sy=param.DY; y>=0; y--, sy++) {
	    for (int x=0; x<param.EntRSizeX; x++) {
		param.EntRot[y][x] = param.EntRot[sy][x];
	    }
	}
	// Mirroring (bottom)
	for (int y=param.DY+param.ImgSizeY, sy=param.DY+param.ImgSizeY-1; y<param.EntRSizeY; y++, sy--) {
	    for (int x=0; x<param.EntRSizeX; x++) {
		param.EntRot[y][x] = param.EntRot[sy][x];
	    }
	}
    }
    
    
    private void CopyClosnImg(BaseParameter param) {
	for (int y=0; y<param.EntRSizeY; y++) {
	    for (int x=0; x<param.EntRSizeX; x++) {
		param.EntRot[y][x] = 0;
	    }
	}
	
	for (int y=0; y<param.ImgSizeY; y++) {
	    for(int x=0; x<param.ImgSizeX; x++) {
		param.EntRot[y+param.DY][x+param.DX] = ClosnImg[y][x];
	    }
	}
	
	int sizeY = param.DY+param.ImgSizeY;
	// Mirroring (left)
	for (int y=param.DY, sy=0; y<sizeY; y++, sy++) {
	    for (int x=param.DX-1, sx=0; x>=0 && sx<param.ImgSizeX; x--, sx++) {
		param.EntRot[y][x] = ClosnImg[sy][sx];
	    }
	}
	// Mirroring (right)
	for (int y=param.DY, sy=0; y<sizeY; y++, sy++) {
	    for (int x=param.DX+param.ImgSizeX, sx=param.ImgSizeX-1; x<param.EntRSizeX && sx>=0; x++, sx--) {
		param.EntRot[y][x] = ClosnImg[sy][sx];
	    }
	}
	// Mirroring (top)
	for (int y=param.DY-1, sy=param.DY; y>=0; y--, sy++) {
	    for (int x=0; x<param.EntRSizeX; x++) {
		param.EntRot[y][x] = param.EntRot[sy][x];
	    }
	}
	// Mirroring (bottom)
	for (int y=param.DY+param.ImgSizeY, sy=param.DY+param.ImgSizeY-1; y<param.EntRSizeY; y++, sy--) {
	    for (int x=0; x<param.EntRSizeX; x++) {
		param.EntRot[y][x] = param.EntRot[sy][x];
	    }
	}
    }
}
