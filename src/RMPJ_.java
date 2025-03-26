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
import ij.plugin.*;
import java.awt.*;
import java.util.Arrays;

/**
* Main class
*
* @author Yoshitaka Kimori
*/

public class RMPJ_ implements PlugIn, DialogListener {
    
    private static int wndCount = 0;
    private boolean bFirstTime = true;
    private boolean stackOptions = false;
    private boolean stackApplyAll = false;
    private String stackSlices = "";
    private Integer[] stackNums;
    private boolean stackFirstTime = true;
    private Point stackPos = new Point(0, 0);
    
    /**
     *	Main processing method for the plugin
     *
     *	@param	arg	Unused
     */
    public void run(String arg) {
	
	// Avoid running the plugin twice
	if (wndCount > 0) {
	    IJ.showMessage("Warning", "RMPJ is already running!");
	    return;
	}
	wndCount++;
	
	BaseParameter param = new BaseParameter();
	
	// Set initial values
	param.filter = "Opening";
	param.seShape = "Disk";
	param.seSize = "3";
	param.rdN = "1";
	
	Point posDialog = new Point(0, 0);
	
	try {
	    while (true) {
		
		if (showRMPJDialog(param, posDialog) == false) {
		    wndCount--;
		    return;
		}
		
		ImagePlus src = WindowManager.getCurrentImage();
		if (src == null) {
		    IJ.error("No Image", "There are no images open.");
		    continue;
		}
		
		// Check original image type
		int type = src.getType();
		if (type != ImagePlus.GRAY8 &&
		    type != ImagePlus.GRAY16)
		    {
			String err = "\"RMPJ\" requires an image of type:" + System.getProperty("line.separator") + System.getProperty("line.separator");
			err += "8-bit grayscale or 16-bit grayscale";
			IJ.error(err);
			continue;
		    }
		if (type == ImagePlus.GRAY8)	param.m_bitDepth = 8;
		if (type == ImagePlus.GRAY16)	param.m_bitDepth = 16;
		
		// Check original image size
		if (param.isImgSizeValid(src.getWidth(), src.getHeight()) == false) {
		    String err = "\"RMPJ\" requires an image of size:" + System.getProperty("line.separator") + System.getProperty("line.separator");
		    err += "smaller than " + BaseParameter.MaxImgSize + " x" + BaseParameter.MaxImgSize + " pixels";
		    IJ.error(err);
		    continue;
		}
		
		ImageStack srcStack = src.getImageStack();
		ImageStack dstStack = new ImageStack(src.getWidth(), src.getHeight());
		
		int size = srcStack.getSize();
		if (size > 1 && stackOptions == false) {
		    
		    String title = "Processing comfirmation";
		    String msg = "\"" + src.getTitle() + "\"" + " has " + size + " images." + System.getProperty("line.separator") + "Are you sure you want to process?";
		    if (IJ.showMessageWithCancel(title, msg) == false) {
			continue;
		    }
		}
		
		for (int i=1; i<=size; i++) {
		    
		    if (!IsProcessImage(i)) {
			continue;
		    }
		    
		    // Initialize the basic data
		    ImageProcessor processor = srcStack.getProcessor(i);
		    param.initialize(processor, i, size);
		    
		    // Morphological operations
		    if (param.filter == "Opening") {
			Opening plugin = new Opening();
			plugin.run(param);
		    }
		    else if (param.filter == "Closing") {
			Closing plugin = new Closing();
			plugin.run(param);
		    }
		    else if (param.filter == "White top-hat") {
			Wth plugin = new Wth();
			plugin.run(param);
		    }
		    else if (param.filter == "Black top-hat") {
			Bth plugin = new Bth();
			plugin.run(param);
		    }
		    else if (param.filter == "Smoothing") {
			Smooth plugin = new Smooth();
			plugin.run(param);
		    }
		    else if (param.filter == "Enhance - type 1") {
			Enhance plugin = new Enhance();
			plugin.run(param);
		    }
		    else if (param.filter == "Enhance - type 2") {
			Enhance2 plugin = new Enhance2();
			plugin.run(param);
		    }
		    
		    // Save the processing result
		    stackResult(param, dstStack, srcStack.getShortSliceLabel(i));
		}
		
		// Display the processing result image
		String filter = param.filter;
		if (filter == "Opening") filter = "Opn";
		if (filter == "Closing") filter = "Clos";
		if (filter == "White top-hat") filter = "WTH";
		if (filter == "Black top-hat") filter = "BTH";
		if (filter == "Smoothing") filter = "MS";
		if (filter == "Enhance - type 1") filter = "MCE1";
		if (filter == "Enhance - type 2") filter = "MCE2";
		String seShape=param.seShape;
		if (seShape == "Disk") seShape = "D";
		if (seShape == "Line") seShape = "L";
		if (seShape == "Square") seShape = "S";
		String title = src.getShortTitle() + "-" + filter+ "_" + seShape + param.seSize + "N" + param.rdN;
		
		if (param.m_bitDepth == 8) {
		    ImagePlus dst = NewImage.createByteImage(title, src.getWidth(), src.getHeight(), 1, NewImage.FILL_BLACK);
		    dst.setStack(dstStack);
		    dst.show();
		}
		
		if (param.m_bitDepth == 16) {
		    ImagePlus dst = NewImage.createShortImage(title, src.getWidth(), src.getHeight(), 1, NewImage.FILL_BLACK);
		    dst.setStack(dstStack);
		    dst.show();
		}
	    }
	}
	catch (Exception e) {	
	    String err = e.getMessage() + System.getProperty("line.separator");
	    
	    StackTraceElement[] stackFrames = e.getStackTrace();
	    if (stackFrames.length >= 1) {
		String stackFrame = stackFrames[0].toString();
		err += stackFrame;
	    }
	    IJ.error("Exception", err);
	}
	
	wndCount--;
    }
    
    
    /**
     *	Displays RMPJ window
     *
     *	@param	param		Common data
     *	@param	posDialog	Window display position
     *	@return	true->OK / false->Cancel
     */
    private boolean showRMPJDialog(BaseParameter param, Point posDialog) {
	
	String[] filter	= {"Opening", "Closing", "White top-hat", "Black top-hat", "Smoothing", "Enhance - type 1", "Enhance - type 2"};
	String[] shape	= {"Disk", "Line", "Square"};
	
	stackOptions = false; 
	
	while (true) {
	    
	    NonBlockingGenericDialog gd = new NonBlockingGenericDialog("RMPJ");
	    gd.addChoice("Operation type :", filter, param.filter);
	    gd.addChoice("Shape of structuring element :", shape, param.seShape);
	    gd.addStringField("Size of structuring element :", param.seSize, 2);
	    gd.setInsets(-32, 220, 0);
	    gd.addMessage("pixels");
	    gd.addStringField("Number of image rotations :", param.rdN, 2);
	    
	    gd.addCheckbox("Stack options", stackOptions);
	    gd.addDialogListener(this);
	    
	    gd.setOKLabel("Apply");
	    
	    try {
		ImageJ ij = IJ.getInstance();
		Image img = ij.getIconImage();
		gd.setIconImage(img);
	    }
	    catch (Exception e) {
	    }
	    
	    if (bFirstTime == false) {
		gd.centerDialog(false);
		gd.setLocation(posDialog);
	    }
	    
	    gd.showDialog();
	    if (gd.wasCanceled()) {
		return false;
	    }
	    
	    posDialog.x = gd.getLocation().x;
	    posDialog.y = gd.getLocation().y;
	    if (bFirstTime == true) {
		bFirstTime = false;
	    }
	    
	    param.filter	= gd.getNextChoice();
	    param.seShape	= gd.getNextChoice();
	    param.seSize	= gd.getNextString();
	    param.rdN		= gd.getNextString();
	    
	    stackOptions = gd.getNextBoolean();
	    
	    if (param.isSESizeOdd() == false) {
		IJ.error("\"Size of the structuring element \" requires an odd value.");
		continue;
	    }
	    if (param.isSESizeValid() == false) {
		IJ.error("\"Size of the structuring element \" requires a value within 3 - 99.");
		continue;
	    }
	    if (param.isRDNValid() == false) {
		IJ.error("\"Number of iterations\" requires a value over 1 and more.");
		continue;
	    }
	    
	    return true;
	}
    }
    
    
    /**
     *	Saves the processing result
     *
     *	@param	param		Common data
     *	@param	dstStack	Where to save the processing result
     *	@param	sliceLabel	Label name of the processing result
     */
    private void stackResult(BaseParameter param, ImageStack dstStack, String sliceLabel) {
	
	int width = dstStack.getWidth();
	int height = dstStack.getHeight();
	
	if (param.m_bitDepth == 8) {
	    ImagePlus img = NewImage.createByteImage(sliceLabel, width, height, 1, NewImage.FILL_BLACK);
	    ImageProcessor ip = img.getProcessor();
	    byte[] pixels = (byte[])ip.getPixels();
	    
	    for (int y=0; y<height; ++y) {
		for (int x=0; x<width; ++x) {
		    int pos = y * width + x;
		    int _y = y + param.offsetY;
		    int _x = x + param.offsetX;
		    pixels[pos] = (byte)param.OutImg[_y][_x];
		}
	    }
	    
	    dstStack.addSlice(sliceLabel, ip);
	}
	
	if (param.m_bitDepth == 16) {
	    ImagePlus img = NewImage.createShortImage(sliceLabel, width, height, 1, NewImage.FILL_BLACK);
	    ImageProcessor ip = img.getProcessor();
	    short[] pixels = (short[])ip.getPixels();
	    
	    for (int y=0; y<height; ++y) {
		for (int x=0; x<width; ++x) {
		    int pos = y * width + x;
		    int _y = y + param.offsetY;
		    int _x = x + param.offsetX;
		    pixels[pos] = (short)param.OutImg[_y][_x];
		}
	    }
	    
	    // Reset display range
	    ip.resetMinAndMax();
	    
	    dstStack.addSlice(sliceLabel, ip);
	}
    }
    
    /**
     *	RMPJ window item-change event
     *
     *	@param	_gd	RMPJ window
     *	@param	_event	Event information
     */
    public boolean dialogItemChanged(GenericDialog _gd, AWTEvent _event) {
	
	if (_event == null) {
	    return true;
	}
	Object objStackOptions = _gd.getCheckboxes().get(0);
	if (!_event.getSource().equals(objStackOptions)) {
	    return true;
	}
	
	Checkbox cb = (Checkbox)objStackOptions;
	if (!cb.getState()) {
	    stackOptions = false;
	    return true;
	}
	
	ImagePlus image = WindowManager.getCurrentImage();
	if (image == null) {
	    IJ.error("No Image", "There are no images open.");
	    cb.setState(false);
	    return true;
	}
	
	// Display stack options window
	while (true) {
	    int stackSize = image.getImageStack().getSize();
	    GenericDialog gdStack = new GenericDialog("Stack options");
	    gdStack.setInsets(0, 80, 0);
	    gdStack.addMessage("Number of images:" + stackSize);
	    gdStack.setInsets(0, 80, 0);
	    gdStack.addCheckbox("Apply all images", stackApplyAll);
	    gdStack.setInsets(20, 50, 0);
	    gdStack.addMessage("Enter a range of images(e.g. 5-10)");
	    gdStack.setInsets(0, 130, 0);
	    gdStack.addMessage("or");
	    gdStack.setInsets(0, 70, 20);
	    gdStack.addMessage("a list of images(e.g. 1,3,5)");
	    gdStack.addStringField("Images", stackSlices, 30);
	    
	    {
		Object objSlices = gdStack.getStringFields().get(0);
		TextField tfSlices = (TextField)objSlices;
		tfSlices.setEnabled(!stackApplyAll);
	    }
	    
	    gdStack.addDialogListener(
				      new DialogListener() {
					  public boolean dialogItemChanged(GenericDialog _gdStack, AWTEvent _awtEvent)
					  {
					      if (_awtEvent == null) {
						  return true;
					      }
					      Object objApplyAllImages = _gdStack.getCheckboxes().get(0);
					      if (!_awtEvent.getSource().equals(objApplyAllImages)) {
						  return true;
					      }
					      
					      Checkbox cbApplyAllImages = (Checkbox)objApplyAllImages;
					      Object objSlices = _gdStack.getStringFields().get(0);
					      TextField tfSlices = (TextField)objSlices;
					      tfSlices.setEnabled(!cbApplyAllImages.getState());
					      
					      return true;
					  }
				      }
				      );
	    
	    if (stackFirstTime == false) {
		gdStack.centerDialog(false);
		gdStack.setLocation(stackPos);
	    }
	    gdStack.showDialog();
	    
	    stackFirstTime = false;
	    stackPos.x = gdStack.getLocation().x;
	    stackPos.y = gdStack.getLocation().y;
	    
	    if (gdStack.wasOKed()) {
		stackApplyAll = gdStack.getNextBoolean();
		stackSlices = gdStack.getNextString();
		if (!checkStackOptions(stackSize, stackApplyAll, stackSlices)) {
		    continue;
		}
	    }
	    else {
		cb.setState(false);
	    }
	    break;
	}
	return true;
    }
    
    /**
     *	Checks the configuration items in the stack options window
     *
     *	@param	imageNum	Number of total slices
     *	@param	applyAll	Process all images
     *	@param	slices		Specify images to process
     */
    private boolean checkStackOptions(int imageNum, boolean applyAll, String slices) {
	
	if (applyAll) {
	    return true;
	}
	
	try {
	    if (slices.isEmpty()) {
		throw new Exception("\"Images\" is empty.");
	    }
	    
	    int indexHyphen = slices.indexOf("-");
	    if (indexHyphen >= 1) {
		String strStart = slices.substring(0, indexHyphen);
		String strEnd = slices.substring(indexHyphen+1);
		int start = Integer.parseInt(strStart);
		int end = Integer.parseInt(strEnd);
		if (end < start) {
		    int temp = end;
		    end = start;
		    start = temp;
		}
		if (start <= 0 || end > imageNum) {
		    throw new Exception("\"Images\" is out of range.");
		}
		stackNums = new Integer[end-start+1];
		for(int n=start, i=0; n<=end; ++n, ++i) {
		    stackNums[i] = n;
		}
	    }
	    else {
		String[] arraySlice = slices.split(",");
		stackNums = new Integer[arraySlice.length];
		for(int i=0; i<arraySlice.length; ++i) {
		    int num = Integer.parseInt(arraySlice[i]);
		    if (num <= 0 || num > imageNum) {
			throw new Exception("\"Imagess\" is out of range.");
		    }
		    stackNums[i] = num;
		}
		Arrays.sort(stackNums);
	    }
	}
	catch (Exception e) {
	    IJ.error(e.getMessage());
	    return false;
	}
	return true;
    }
    
    private boolean IsProcessImage(int n) {
	
	if (stackOptions == false) {
	    return true;
	}
	if (stackApplyAll == true) {
	    return true;
	}
	
	Object target = new Integer(n);
	return Arrays.asList(stackNums).contains(target);
    }
}
