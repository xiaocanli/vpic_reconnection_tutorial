;+
; NAME:
;       TVREAD
;
; PURPOSE:
;
;       To get accurate screen dumps with the IDL command TVRD on 24-bit
;       PC and Macintosh computers, you have to be sure to set color
;       decomposition on. This program adds that capability automatically.
;       In addition, the program will optionally write BMP, GIF, JPEG,
;       PICT, PNG, and TIFF color image files of the screen dump.
;
; AUTHOR:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;
;       Graphics
;
; CALLING SEQUENCE:
;
;       image = TVREAD(xstart, ystart, ncols, nrows)
;
;       The returned image will be a 2D image on 8-bit systems and
;       a 24-bit pixel interleaved true-color image on 24-bit systems.
;       A -1 will be returned if a file output keyword is used (e.g., JPEG, TIFF, etc.).
;
; OPTIONAL INPUTS:
;
;       XSTART -- The starting column index.  By default, 0.
;
;       YSTART -- The starting row index. By default, 0.
;
;       NCOLS -- The number of columns to read. By default, !D.X_Size - XSTART
;
;       NROWS -- The number of rows to read. By default, !D.Y_Size - YSTART.
;
; KEYWORD PARAMETERS:
;
;       BMP -- Set this keyword to write the screen dump as a color BMP file.
;
;       CANCEL -- An output keyword set to 1 if the user cancels out of a
;          filename dialog. Set to 0 otherwise.
;
;       COLORS -- If a 24-bit image has to be quantized, this will set the number
;          of colors in the output image. Set to 256 by default. Applies to BMP,
;          GIF, PICT, and PNG formats written from 24-bit displays.(See the
;          COLOR_QUAN documentation for details.)
;
;       CUBE -- If this keyword is set to a value between 2 and 6 the color
;          quantization will use a cubic method of quantization. Applies to BMP,
;          GIF, PICT, and PNG formats written from 24-bit displays.(See the
;          COLOR_QUAN documentation for details.)
;
;       DITHER -- If this keyword is set the quantized image will be dithered.
;          Applies to BMP, GIF, PICT, and PNG formats written from 24-bit displays.
;          (See the COLOR_QUAN documentation for details.)
;
;       FILENAME -- The base name of the output file. (No file extensions;
;           they will be added automatically.) This name may be changed by the user.
;
;              image = TVREAD(Filename='myfile', /JPEG)
;
;           No file will be written unless a file output keyword is used
;           (e.g., JPEG, TIFF, etc.) in the call. By default the FILENAME is
;           set to "idl". The file extension will be set automatically to match
;           the type of file created.
;
;       GIF -- Set this keyword to write the screen dump as a color GIF file.
;
;       JPEG -- Set this keyword to write the screen dump as a color JPEG file.
;
;       NODIALOG -- Set this keyword if you wish to avoid the DIALOG_PICKFILE
;           dialog that asks you to name the output file. This keyword should be
;           set, for example, if you are processing screens in batch mode.
;
;       ORDER -- Set this keyword to determine the image order for reading the
;           display. Corresponds to !Order and set to such as the default.
;           
;       OVERWRITE_PROMPT -- Set this keyword if you would like to get a prompt
;           if you are overwriting a file. This applies only to operations with
;           DIALOG_PICKFILE.
;
;       PICT -- Set this keyword to write the screen dump as a color PICT file.
;
;       PNG -- Set this keyword to write the screen dump as a color PNG file.
;
;       TIFF -- Set this keyword to write the screen dump as a color TIFF file.
;
;       TRUE -- Set this keyword to the type of interleaving you want. 1 = Pixel
;           interleaved, 2 = row interleaved, 3 = band interleaved.
;
;       TYPE -- Can be set to the type of file to write. Use this instead of
;           setting BMP, GIF, JPEG, PICT, PNG, or TIFF keywords: TYPE='JPEG'. The
;           primary purpose of this is to make event handlers easier to write.
;
;       QUALITY -- This keyword sets the amount of compression for JPEG images.
;           It should be set to a value between 0 and 100. It is set to 75 by default.
;           (See the WRITE_JPEG documentation for details.)
;
;       WID -- The index number of the window to read from. The current graphics window
;           (!D.Window) is selected by default. An error is issued if no windows are
;           currently open on a device that supports windows.
;
;       _EXTRA -- Any keywords that are appropriate for the WRITE_*** routines are
;           also accepted via keyword inheritance.
;
; COMMON BLOCKS:
;
;       None
;
; RESTRICTIONS:   Requires IDL 5.2 and higher.
;
; MODIFICATION HISTORY:
;
;       Written by David W. Fanning, 9 AUG 2000.
;       Added changes to make the program more device independent. 16 SEP 2000. DWF.
;       Removed GIF file support for IDL 5.4 and above. 18 JAN 2001. DWF.
;       Added NODIALOG keyword. 28 MAR 2001. DWF.
;       Added an output CANCEL keyword. 29 AUG 2001. DWF.
;       Added ERROR_MESSAGE code to file. 17 DEC 2001. DWF.
;       Added ORDER keyword. 25 March 2002. DWF.
;       Now create 24-bit PNG files if reading from a 24-bit display. 11 May 2002. DWF.
;       Now create 24-bit BMP files if reading from a 24-bit display. 23 May 2002. DWF.
;       Removed obsolete STR_SEP and replaced with STRSPLIT. 27 Oct 2004. DWF.
;       Unleashed the GIF code for IDL 6.2 and above. 10 Nov 2005. DWF.
;       Rolled back previous change to IDL 6.1. 24 Jan 2006. DWF.
;       Fixed a problem in which 16-bit displays were treated as 24-bit displays,
;         and as a result could not produce WHITE colors. Now all 16-bit displays
;         are treated as 8-bit displays. It is an ugly trade-off. 24 Jan 2006. DWF.
;       Added TYPE keyword 27 Sep 2006. DWF.
;       Updated program to work with 24-bit Z-buffer in IDL 6.4. 11 June 2007. DWF.
;       Added OVERWRITE_PROMPT keyword. 2 Oct 2008. DWF.
;-
;
;******************************************************************************************;
;  Copyright (c) 2008, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
FUNCTION TVREAD_ERROR_MESSAGE, theMessage, Traceback=traceback, $
   NoName=noName, _Extra=extra

On_Error, 2

   ; Check for presence and type of message.

IF N_Elements(theMessage) EQ 0 THEN theMessage = !Error_State.Msg
s = Size(theMessage)
messageType = s[s[0]+1]
IF messageType NE 7 THEN BEGIN
   Message, "The message parameter must be a string.", _Extra=extra
ENDIF

   ; Get the call stack and the calling routine's name.

Help, Calls=callStack
callingRoutine = (StrSplit(StrCompress(callStack[1])," ", /Extract))[0]

   ; Are widgets supported? Doesn't matter in IDL 5.3 and higher.

widgetsSupported = ((!D.Flags AND 65536L) NE 0) OR Float(!Version.Release) GE 5.3
IF widgetsSupported THEN BEGIN
   IF Keyword_Set(noName) THEN answer = Dialog_Message(theMessage, _Extra=extra) ELSE BEGIN
      IF StrUpCase(callingRoutine) EQ "$MAIN$" THEN answer = Dialog_Message(theMessage, _Extra=extra) ELSE $
         answer = Dialog_Message(StrUpCase(callingRoutine) + ": " + theMessage, _Extra=extra)
   ENDELSE
ENDIF ELSE BEGIN
      Message, theMessage, /Continue, /NoPrint, /NoName, /NoPrefix, _Extra=extra
      Print, '%' + callingRoutine + ': ' + theMessage
      answer = 'OK'
ENDELSE

   ; Provide traceback information if requested.

IF Keyword_Set(traceback) THEN BEGIN
   Help, /Last_Message, Output=traceback
   Print,''
   Print, 'Traceback Report from ' + StrUpCase(callingRoutine) + ':'
   Print, ''
   FOR j=0,N_Elements(traceback)-1 DO Print, "     " + traceback[j]
ENDIF

RETURN, answer
END ; ----------------------------------------------------------------------------


FUNCTION TVREAD, xstart, ystart, ncols, nrows, $
   BMP=bmp, $
   Cancel=cancel, $
   Colors=colors, $
   Cube=cube, $
   Dither=dither, $
   _Extra=extra, $
   Filename=filename, $
   GIF=gif, $
   JPEG=jpeg, $
   NoDialog=nodialog, $
   Order=order, $
   Overwrite_Prompt=overwrite_prompt, $
   PICT=pict, $
   PNG=png, $
   TIFF=tiff, $
   True=true, $
   TYPE=type, $
   Quality=quality, $
   WID=wid

   ; Error handling.

Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   ok = TVRead_Error_Message(Traceback=1, /Error)
   IF N_Elements(thisWindow) EQ 0 THEN RETURN, -1
   IF thisWindow GE 0 THEN WSet, thisWindow
   RETURN, -1
ENDIF

cancel = 0

   ; Check for availability of GIF files.

thisVersion = Float(!Version.Release)
IF (thisVersion LT 5.3) OR (thisVersion GE 6.1) THEN haveGif = 1 ELSE haveGIF = 0

   ; Go to correct window.

IF N_Elements(wid) EQ 0 THEN wid =!D.Window
thisWindow = !D.Window
IF (!D.Flags AND 256) NE 0 THEN WSet, wid

   ; Check keywords and parameters. Define values if necessary.

IF N_Elements(xstart) EQ 0 THEN xstart = 0
IF N_Elements(ystart) EQ 0 THEN ystart = 0
IF N_Elements(ncols) EQ 0 THEN ncols = !D.X_Size - xstart
IF N_Elements(nrows) EQ 0 THEN nrows = !D.Y_Size - ystart
IF N_Elements(order) EQ 0 THEN order = !Order
IF N_Elements(true) EQ 0 THEN true = 1
dialog = 1 - Keyword_Set(nodialog)

   ; Do you want to write an image file instead of
   ; capturing an image?
IF N_Elements(type) NE 0 THEN BEGIN
   CASE StrUpCase(type) OF
      'BMP': bmp = 1
      'GIF': gif = 1
      'JPEG': jpeg = 1
      'JPG': jpeg = 1
      'PICT': pict = 1
      'PNG': png = 1
      'TIFF': tiff = 1
      'TIF': tif = 1
      ELSE: Message, 'Cannot write a file of type: ' + StrUpCase(type) + '.'
   ENDCASE
ENDIF
writeImage = 0
fileType = ""
extention = ""
IF Keyword_Set(bmp)THEN BEGIN
   writeImage = 1
   fileType = 'BMP'
   extension = 'bmp'
ENDIF
IF Keyword_Set(gif) THEN BEGIN
   IF havegif THEN BEGIN
      writeImage = 1
      fileType = 'GIF'
      extension = 'gif'
    ENDIF ELSE BEGIN
       ok = Dialog_Message('GIF files not supported in this IDL version. Replacing with JPEG.')
       writeImage = 1
      fileType = 'JPEG'
      extension = 'jpg'
   ENDELSE
ENDIF
IF Keyword_Set(jpeg) THEN BEGIN
   writeImage = 1
   fileType = 'JPEG'
   extension = 'jpg'
ENDIF
IF Keyword_Set(PICT) THEN BEGIN
   writeImage = 1
   fileType = 'PICT'
   extension = 'pict'
ENDIF
IF Keyword_Set(png) THEN BEGIN
   writeImage = 1
   fileType = 'PNG'
   extension = 'png'
ENDIF
IF Keyword_Set(tiff) THEN BEGIN
   writeImage = 1
   fileType = 'TIFF'
   extension = 'tif'
ENDIF

IF N_Elements(colors) EQ 0 THEN colors = 256
IF N_Elements(quality) EQ 0 THEN quality = 75
dither = Keyword_Set(dither)

   ; On 24-bit displays, make sure color decomposition is ON.

IF (!D.Flags AND 256) NE 0 THEN BEGIN
   Device, Get_Decomposed=theDecomposedState, Get_Visual_Depth=theDepth
   IF theDepth GT 8 THEN BEGIN
      Device, Decomposed=1
      IF theDepth EQ 24 THEN truecolor = true ELSE truecolor = 0
   ENDIF ELSE truecolor = 0
   IF thisWindow LT 0 THEN $
      Message, 'No currently open windows. Returning.', /NoName
ENDIF ELSE BEGIN
   truecolor = 0
   theDepth = 8
ENDELSE

; Fix for 24-bit Z-buffer.
IF (Float(!Version.Release) GE 6.4) AND (!D.NAME EQ 'Z') THEN BEGIN
   Device, Get_Decomposed=theDecomposedState, Get_Pixel_Depth=theDepth
   IF theDepth EQ 24 THEN truecolor = true ELSE truecolor = 0
ENDIF


   ; Get the screen dump. 2D image on 8-bit displays. 3D image on 24-bit displays.

image = TVRD(xstart, ystart, ncols, nrows, True=truecolor, Order=order)

   ; Need to set color decomposition back?

IF theDepth GT 8 THEN Device, Decomposed=theDecomposedState

   ; If we need to write an image, do it here.

IF writeImage THEN BEGIN

      ; Get the name of the output file.

   IF N_Elements(filename) EQ 0 THEN BEGIN
      filename = 'idl.' + StrLowCase(extension)
   ENDIF ELSE BEGIN
      filename = filename + "." + StrLowCase(extension)
   ENDELSE
   IF dialog THEN filename = Dialog_Pickfile(/Write, File=filename, OVERWRITE_PROMPT=Keyword_Set(overwrite_prompt))

   IF filename EQ "" THEN BEGIN
      cancel = 1
      RETURN, image
   ENDIF

      ; Write the file.

   CASE fileType OF

      'BMP': BEGIN
         IF truecolor THEN BEGIN
            ; BMP files assume blue, green, red planes.
            temp = image[0,*,*]
            image[0,*,*] = image[2,*,*]
            image[2,*,*] = temp
            Write_BMP, filename, image, _Extra=extra
         ENDIF ELSE BEGIN
            TVLCT, r, g, b, /Get
            Write_BMP, filename, image, r, g, b, _Extra=extra
         ENDELSE
         END

      'GIF': BEGIN
         IF truecolor THEN BEGIN
            CASE Keyword_Set(cube) OF
               0: image2D = Color_Quan(image, 1, r, g, b, Colors=colors, Dither=dither)
               1: image2D = Color_Quan(image, 1, r, g, b, Cube=2 > cube < 6)
            ENDCASE
         ENDIF ELSE BEGIN
            TVLCT, r, g, b, /Get
            image2D = image
         ENDELSE
         Write_GIF, filename, image2D, r, g, b, _Extra=extra
         END

      'JPEG': BEGIN
         IF truecolor THEN BEGIN
            image3D = image
         ENDIF ELSE BEGIN
            s = Size(image, /Dimensions)
            image3D = BytArr(3, s[0], s[1])
            TVLCT, r, g, b, /Get
            image3D[0,*,*] = r[image]
            image3D[1,*,*] = g[image]
            image3D[2,*,*] = b[image]
         ENDELSE
         Write_JPEG, filename, image3D, True=1, Quality=quality, _Extra=extra
         END

      'PICT': BEGIN
         IF truecolor THEN BEGIN
            CASE Keyword_Set(cube) OF
               0: image2D = Color_Quan(image, 1, r, g, b, Colors=colors, Dither=dither)
               1: image2D = Color_Quan(image, 1, r, g, b, Cube=2 > cube < 6)
            ENDCASE
         ENDIF ELSE BEGIN
            TVLCT, r, g, b, /Get
            image2D = image
         ENDELSE
         Write_PICT, filename, image2D, r, g, b
         END

      'PNG': BEGIN
         IF truecolor THEN BEGIN
            Write_PNG, filename, image, _Extra=extra
         ENDIF ELSE BEGIN
            TVLCT, r, g, b, /Get
            image2D = image
            Write_PNG, filename, Reverse(image2D,2), r, g, b, _Extra=extra
            Write_PNG, filename, image2D, r, g, b, _Extra=extra
         ENDELSE
         END

      'TIFF': BEGIN
         IF truecolor THEN BEGIN
            image3D = Reverse(image,3)
         ENDIF ELSE BEGIN
            s = Size(image, /Dimensions)
            image3D = BytArr(3, s[0], s[1])
            TVLCT, r, g, b, /Get
            image3D[0,*,*] = r[image]
            image3D[1,*,*] = g[image]
            image3D[2,*,*] = b[image]
            image3D = Reverse(Temporary(image3D), 3)
         ENDELSE
         Write_TIFF, filename, image3D, 1, _Extra=extra
         END
   ENDCASE
   RETURN, -1
ENDIF

   ; Return the screen dump image.

RETURN, image
END ;-------------------------------------------------------------------------------


pro read_slice,ifield,slice
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common controlplot,v
common choice, rawdata,computedquantity
common dir,directory
common list,ns,order,suffix

; Open file

nstep=0L
nstep = ns(order(v.time))
file=strcompress(directory+rawdata(ifield)+'_'+string(nstep)+'.gda',/remove_all)
if ( file_test(file) ne 1 ) then file=strcompress('success/'+rawdata(ifield)+'_'+string(nstep)+'.gda',/remove_all)
if ( file_test(file) ne 1 ) then begin
   print,"File does not exist
   exit
endif

;  Open File and read data

openr,10,file
print,"Reading Data --> ",file

; Declare array for data

slice = fltarr(nx,nz)
buff = fltarr(nx)

; Declare integer index

n = 0L

; Associate buffer space with specified data file

dfile = assoc(10,buff)

; Now loop through the 3D array in the data file and read only the
; needed y-slice of the data file - assume array is in FORTRAN ordering

j=v.ycut
for k=0,nz-1 do begin
    n = (j + ny*k)
    slice(*,k) = dfile[n]
 endfor

; Smooth data 

slice = smooth(slice,v.smoothing,/EDGE_TRUNCATE)    

; Close file

close,10

end

pro read_slice3,ifield,slice3
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common controlplot,v
common choice, rawdata,computedquantity
common dir,directory
common list,ns,order,suffix

; Open file

nstep=0L
nstep = ns(order(v.time))
file=strcompress(directory+rawdata(ifield)+'_'+string(nstep)+'.gda',/remove_all)
if ( file_test(file) ne 1 ) then file=strcompress('success/'+rawdata(ifield)+'_'+string(nstep)+'.gda',/remove_all)

;  Open File and read data

openr,10,file
print,"Reading Data --> ",file

; Declare array for data

slice3 = fltarr(nx,3,nz)
buff = fltarr(nx)

; Declare integer index

n = 0L

; Associate buffer space with specified data file

dfile = assoc(10,buff)

; Make sure we don't request data outside of array

j=v.ycut
if (v.ycut eq 0) then j = 1
if (v.ycut eq ny-1) then j = ny-2

; Now loop through the 3D array in the data file and read only the
; needed y-slice of the data file - assume array is in FORTRAN ordering

for k=0,nz-1 do begin
    n = (j - 1 + ny*k)
    slice3(*,0,k) = dfile[n]
    n = (j + ny*k)
    slice3(*,1,k) = dfile[n]
    n = (j + 1 + ny*k)
    slice3(*,2,k) = dfile[n]
 endfor

; Smooth data 

;slice = smooth(slice,v.smoothing,/EDGE_TRUNCATE)    

; Close file

close,10

end

pro beta,isp
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
if isp eq 0 then print, "Computing Beta-e"
if isp eq 1 then print, "Computing Beta-i"
if isp eq 2 then print, "Computing Beta-total"

; Declare needed variables

Pxx = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pzz = fltarr(nx,nz)
Pixx = fltarr(nx,nz)
Piyy = fltarr(nx,nz)
Pizz = fltarr(nx,nz)

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

bm = fltarr(nx,nz)

; Read needed data for contour plot

print,'--------------------------------'
print,"Time Index=",v.time

; Read magnetic field

read_slice,ibx,bx
read_slice,iby,by
read_slice,ibz,bz

; Read Pressure

if isp eq 0 then begin
   mytitle="Beta-e"
   read_slice,ipexx,Pxx
   read_slice,ipeyy,Pyy
   read_slice,ipezz,Pzz
endif
if isp eq 1 then begin
   mytitle="Beta-i"
   read_slice,ipixx,Pxx
   read_slice,ipiyy,Pyy
   read_slice,ipizz,Pzz
endif
if isp eq 2 then begin
   mytitle="Beta-Total"
   read_slice,ipixx,Pixx
   read_slice,ipiyy,Piyy
   read_slice,ipizz,Pizz
   read_slice,ipexx,Pxx
   read_slice,ipeyy,Pyy
   read_slice,ipezz,Pzz
   Pxx = Pxx + Pixx
   Pyy = Pyy + Piyy
   Pzz = Pzz + Pizz
endif

; Return beta

fulldata = 2.0*(Pxx+Pyy+Pzz)/(bx^2+by^2+bz^2)/3.0

end

pro temperature,isp
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
if isp eq 0 then print, "Computing Te"
if isp eq 1 then print, "Computing Ti"
if isp eq 2 then print, "Computing T-total"

; Declare needed variables

Pxx = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pzz = fltarr(nx,nz)
Pixx = fltarr(nx,nz)
Piyy = fltarr(nx,nz)
Pizz = fltarr(nx,nz)

nee = fltarr(nx,nz)

; Read needed data for contour plot

print,'--------------------------------'
print,"Time Index=",v.time

; Read Density

read_slice,ine,nee

; Read Pressure

if isp eq 0 then begin
   mytitle="Te"
   read_slice,ipexx,Pxx
   read_slice,ipeyy,Pyy
   read_slice,ipezz,Pzz
endif
if isp eq 1 then begin
   mytitle="Ti"
   read_slice,ipixx,Pxx
   read_slice,ipiyy,Pyy
   read_slice,ipizz,Pzz
endif
if isp eq 2 then begin
   mytitle="T-Total"
   read_slice,ipixx,Pixx
   read_slice,ipiyy,Piyy
   read_slice,ipizz,Pizz
   read_slice,ipexx,Pxx
   read_slice,ipeyy,Pyy
   read_slice,ipezz,Pzz
   Pxx = Pxx + Pixx
   Pyy = Pyy + Piyy
   Pzz = Pzz + Pizz
endif

; Return temperature

fulldata = (Pxx+Pyy+Pzz)/nee

end


pro anisotropy,isp
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
if isp eq 0 then print, "Computing Electron Anisotropy"
if isp eq 1 then print, "Computing Ion Anisotropy"

; Declare needed variables

Pxx = fltarr(nx,nz)
Pxy = fltarr(nx,nz)
Pxz = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pyz = fltarr(nx,nz)
Pzz = fltarr(nx,nz)

ppar = fltarr(nx,nz)
pper = fltarr(nx,nz)

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

bm = fltarr(nx,nz)

; Read needed data 

; Magnetic field data

read_slice,ibx,bx
read_slice,iby,by
read_slice,ibz,bz

bm = sqrt(bx^2 + by^2 + bz^2)
bx = bx/bm
by = by/bm
bz = bz/bm

; Choose either electrons or ions

if isp eq 0 then begin
   mytitle="Ane"
   read_slice,ipexx,Pxx
   read_slice,ipexy,Pxy
   read_slice,ipexz,Pxz
   read_slice,ipeyy,Pyy
   read_slice,ipeyz,Pyz
   read_slice,ipezz,Pzz
endif
if isp eq 1 then begin
   mytitle="Ani"
   read_slice,ipixx,Pxx
   read_slice,ipixy,Pxy
   read_slice,ipixz,Pxz
   read_slice,ipiyy,Pyy
   read_slice,ipiyz,Pyz
   read_slice,ipizz,Pzz
endif

; Compute parallel and perpendicular pressure

ppar = Pxx*bx^2 + Pyy*by^2 + Pzz*bz^2 + 2.0*Pxy*bx*by + 2.0*Pxz*bx*bz+ 2.0*Pyz*by*bz
pper = (Pxx + Pyy +Pzz - ppar)/2.0

; Return data

fulldata = pper/ppar

end

pro agyrotropy,isp
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
if isp eq 0 then print, "Computing Electron Agyotropy"
if isp eq 1 then print, "Computing Ion Agyotropy"

; Declare needed variables

Pxx = fltarr(nx,nz)
Pxy = fltarr(nx,nz)
Pxz = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pyz = fltarr(nx,nz)
Pzz = fltarr(nx,nz)

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

bm = fltarr(nx,nz)

Nxx = fltarr(nx,nz)
Nxy = fltarr(nx,nz)
Nxz = fltarr(nx,nz)
Nyy = fltarr(nx,nz)
Nyz = fltarr(nx,nz)
Nzz = fltarr(nx,nz)
alpha = fltarr(nx,nz)
beta = fltarr(nx,nz)

; Read needed data 

print,'--------------------------------'
print," Smoothing =",v.smoothing

; Magnetic Field data

read_slice,ibx,bx
read_slice,iby,by
read_slice,ibz,bz

bm = sqrt(bx^2 + by^2 + bz^2)
bx = bx/bm
by = by/bm
bz = bz/bm

; Choose either electrons or ions

if isp eq 0 then begin
   print, "Electron Agyrotropy"
   mytitle="A0e"
   read_slice,ipexx,Pxx
   read_slice,ipexy,Pxy
   read_slice,ipexz,Pxz
   read_slice,ipeyy,Pyy
   read_slice,ipeyz,Pyz
   read_slice,ipezz,Pzz
endif
if isp eq 1 then begin
   print, "Ion Agyrotropy"
   mytitle="A0i"
   read_slice,ipixx,Pxx
   read_slice,ipixy,Pxy
   read_slice,ipixz,Pxz
   read_slice,ipiyy,Pyy
   read_slice,ipiyz,Pyz
   read_slice,ipizz,Pzz
endif

; Compute parallel and perpendicular pressure

Nxx =  by*by*Pzz - 2.0*by*bz*Pyz + bz*bz*Pyy
Nxy = -by*bx*Pzz + by*bz*Pxz + bz*bx*Pyz - bz*bz*Pxy
Nxz =  by*bx*Pyz - by*by*Pxz - bz*bx*Pyy + bz*by*Pxy
Nyy =  bx*bx*Pzz - 2.0*bx*bz*Pxz + bz*bz*Pxx
Nyz = -bx*bx*Pyz + bx*by*Pxz + bz*bx*Pxy - bz*by*Pxx
Nzz =  bx*bx*Pyy - 2.0*bx*by*Pxy + by*by*Pxx

alpha = Nxx + Nyy + Nzz
beta = -(Nxy^2 + Nxz^2 + Nyz^2 - Nxx*Nyy - Nxx*Nzz - Nyy*Nzz)

; Return agyrotropy data

fulldata = 2.0*sqrt(alpha^2-4.0*beta)/alpha

end

pro firehose,nsp
common pdata,fulldata,mtytile
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
if nsp eq 0 then begin
    print, "Electron Firehose Condition"
    mytitle="E-Firehose"
endif
if nsp eq 1 then begin
    print, "Total Firehose Condition"
    mytitle="Firehose"
endif

; Declare needed variables

Pxx = fltarr(nx,nz)
Pxy = fltarr(nx,nz)
Pxz = fltarr(nx,nz)
Pyy = fltarr(nx,nz)
Pyz = fltarr(nx,nz)
Pzz = fltarr(nx,nz)

ppar = fltarr(nx,nz)
pper = fltarr(nx,nz)

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

bm = fltarr(nx,nz)

; Read needed data 

; Magnetic field data

read_slice,ibx,bx
read_slice,iby,by
read_slice,ibz,bz

bm = sqrt(bx^2 + by^2 + bz^2)
bx = bx/bm
by = by/bm
bz = bz/bm

; Read in Electron Pressure

read_slice,ipexx,Pxx
read_slice,ipexy,Pxy
read_slice,ipexz,Pxz
read_slice,ipeyy,Pyy
read_slice,ipeyz,Pyz
read_slice,ipezz,Pzz

; Read in Ion Pressure

if (nsp eq 1) then begin

   Pixx = fltarr(nx,nz)
   Pixy = fltarr(nx,nz)
   Pixz = fltarr(nx,nz)
   Piyy = fltarr(nx,nz)
   Piyz = fltarr(nx,nz)
   Pizz = fltarr(nx,nz)

   read_slice,ipixx,Pixx
   read_slice,ipixy,Pixy
   read_slice,ipixz,Pixz
   read_slice,ipiyy,Piyy
   read_slice,ipiyz,Piyz
   read_slice,ipizz,Pizz

; Add to compute total pressure

   Pxx = Pxx + Pixx
   Pxy = Pxy + Pixy
   Pxz = Pxz + Pixz
   Pyy = Pyy + Piyy
   Pyz = Pyz + Piyz
   Pzz = Pzz + Pizz

endif

; Compute parallel and perpendicular pressure

ppar = Pxx*bx^2 + Pyy*by^2 + Pzz*bz^2 + 2.0*Pxy*bx*by + 2.0*Pxz*bx*bz+ 2.0*Pyz*by*bz
pper = (Pxx + Pyy + Pzz - ppar)/2.0

; Compute firehose condition - this is the actual phase velocity of the
;                              Alfven wave divided by Va (isotropic)

fulldata = 1+ pper/bm^2-ppar/bm^2

end

pro vorticity,isp
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
if (isp eq 0) then begin
    print,"Computing Electron Vorticity"
    mytitle="Electron Vorticity"
endif
if (isp eq 1) then begin
    print,"Computing Total Vorticity"
    mytitle="Total Vorticity"
endif

; Declare need variables

; velocity
ux = fltarr(nx,3,nz)
uy = fltarr(nx,3,nz)
uz = fltarr(nx,3,nz)

uxe = fltarr(nx,3,nz)
uye = fltarr(nx,3,nz)
uze = fltarr(nx,3,nz)

; vorticity
wx = fltarr(nx,nz)
wy = fltarr(nx,nz)
wz = fltarr(nx,nz)

; Read data

; Electron velocity
if (isp eq 0) then begin
    read_slice3,iuex,ux
    read_slice3,iuey,uy
    read_slice3,iuez,uz
endif

; Total velocity
if (isp eq 1) then begin
    read_slice3,iuex,uxe
    read_slice3,iuey,uye
    read_slice3,iuez,uze
    read_slice3,iuix,ux
    read_slice3,iuiy,uy
    read_slice3,iuiz,uz
    ux = ux + uxe
    uy = uy + uye
    uz = uz + uze
endif

; Cell sizes

dx = 2.0*xmax/fix(nx)
dy = 2.0*ymax/fix(ny)
dz = 2.0*zmax/fix(nz)

; Now compute vorticity - edit to pick what you want

for i=1,nx-2 do begin
    for j=1,nz-2 do begin
        wx(i,j) = (uz(i,2,j)-uz(i,0,j))/dy - (uy(i,1,j+1) - uy(i,1,j-1))/dz
        wy(i,j) = (ux(i,1,j+1)-ux(i,1,j-1))/dz - (uz(i+1,1,j) - uz(i-1,1,j))/dz
        wz(i,j) = (uy(i+1,1,j)-uy(i-1,1,j))/dx - (ux(i,2,j) - ux(i,0,j))/dy
    endfor
endfor

; Now smooth if desired

wx = smooth(wx,v.smoothing,/EDGE_TRUNCATE)    
wy = smooth(wy,v.smoothing,/EDGE_TRUNCATE)    
wz = smooth(wz,v.smoothing,/EDGE_TRUNCATE)    

; Return desired component or total

fulldata = sqrt(wx^2 + wy^2 + wz^2)

end


pro curvature
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'

print,"Computing Field Line Curvature"
mytitle="Curvature"

; Declare need variables

; velocity
bx = fltarr(nx,3,nz)
by = fltarr(nx,3,nz)
bz = fltarr(nx,3,nz)

; curvature - B*grad(B)
cx = fltarr(nx,nz)
cy = fltarr(nx,nz)
cz = fltarr(nx,nz)

; Read magnetic field

read_slice3,ibx,bx
read_slice3,iby,by
read_slice3,ibz,bz

bx = smooth(bx,v.smoothing,/EDGE_TRUNCATE)    
by = smooth(by,v.smoothing,/EDGE_TRUNCATE)    
bz = smooth(bz,v.smoothing,/EDGE_TRUNCATE)    

; Cell sizes

dx = 2.0*xmax/fix(nx)
dy = 2.0*ymax/fix(ny)
dz = 2.0*zmax/fix(nz)

; Now compute curvature terms - B*grad(B)

for i=1,nx-2 do begin
    for j=1,nz-2 do begin
        cx(i,j) = bx(i,1,j)*(bx(i+1,1,j)-bx(i-1,1,j))/dx + by(i,1,j)*(bx(i,2,j)-bx(i,0,j))/dy + bz(i,1,j)*(bx(i,1,j+1)-bx(i,1,j-1))/dz 
        cy(i,j) = bx(i,1,j)*(by(i+1,1,j)-by(i-1,1,j))/dx + by(i,1,j)*(by(i,2,j)-by(i,0,j))/dy + bz(i,1,j)*(by(i,1,j+1)-by(i,1,j-1))/dz 
        cz(i,j) = bx(i,1,j)*(bz(i+1,1,j)-bz(i-1,1,j))/dx + by(i,1,j)*(bz(i,2,j)-bz(i,0,j))/dy + bz(i,1,j)*(bz(i,1,j+1)-bz(i,1,j-1))/dz 
    endfor
endfor

; Now smooth if desired

cx = smooth(cx,v.smoothing,/EDGE_TRUNCATE)    
cy = smooth(cy,v.smoothing,/EDGE_TRUNCATE)    
cz = smooth(cz,v.smoothing,/EDGE_TRUNCATE)    

; Return desired component or total

fulldata = sqrt(cx^2 + cy^2 + cz^2)/(bx^2+by^2+bz^2)

end


pro kelvin_helmholtz
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
print, "Kelvin Helmholtz Stability"
mytitle="KH Stability"

; Declare needed variables

ux = fltarr(nx,nz)
uy = fltarr(nx,nz)
uz = fltarr(nx,nz)
nee = fltarr(nx,nz)

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

; Read needed data

; Magnetic field data

read_slice,ibx,bx
read_slice,iby,by
read_slice,ibz,bz

; Ion Velocity data

read_slice,iuix,ux
read_slice,iuiy,uy
read_slice,iuiz,uz

; Electron density

read_slice,ine,nee

; Comput KH ratio
;
;  ratio = (u*u)/(2u*Va)
;
;  where u is local fluid velocity
;        Va is local Alfven velocity
;
;  Need ratio > 1 to drive KH
;
;  *** Need to specify mass ratio ***

mime=100

;  Assume "k" is in the direction of "u"
;fulldata = (ux^2 + uy^2 + uz^2)*sqrt(nee*mime)/(ux*bx + uy*by + uz*bz)/2

;  Use angle slider to pick direction

pi=3.1415927
kx = cos(v.angle*pi/180)
ky = sin(v.angle*pi/180)
kz = 0.0
print,"k=",kx,ky,kz

; Use this direction to compute KH ratio

fulldata = ((kx*ux + ky*uy + kz*uz)*sqrt(nee*mime))/(2*(kx*bx + ky*by + kz*bz))

end

pro J_magnitude
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
print, "Compute |J| "
mytitle="|J|"

; Declare needed variables

jx = fltarr(nx,nz)
jy = fltarr(nx,nz)
jz = fltarr(nx,nz)

; Current density

read_slice,ijx,jx
read_slice,ijy,jy
read_slice,ijz,jz

; Compute magntidue of total current density

fulldata = sqrt(jx^2+jy^2+jz^2)

end

pro B_magnitude
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
print, "Computing |B| "
mytitle="|B|"

; Declare needed variables

bx = fltarr(nx,nz)
by = fltarr(nx,nz)
bz = fltarr(nx,nz)

; Read needed data

; Magnetic field data

read_slice,ibx,bx
read_slice,iby,by
read_slice,ibz,bz

; Compute |B|

fulldata = sqrt(bx^2 + by^2 + bz^2)

end

pro charge_density
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
print, "Compute Net Charge Density"
mytitle="ni-ne"

; Declare needed variables

dne = fltarr(nx,nz)
dni = fltarr(nx,nz)

; Current density

read_slice,ini,dni
read_slice,ine,dne

; Compute magntidue of total current density

fulldata = dni-dne

end


pro EdotJ
common pdata,fulldata,mytitle
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct

; Screen Information

print,'--------------------------------'
print, "Computing E*J "
mytitle="E*J"

; Declare needed variables

jx = fltarr(nx,nz)
jy = fltarr(nx,nz)
jz = fltarr(nx,nz)
ex = fltarr(nx,nz)
ey = fltarr(nx,nz)
ez = fltarr(nx,nz)

; Current density

read_slice,ijx,jx
read_slice,ijy,jy
read_slice,ijz,jz

; Electric field

read_slice,iex,ex
read_slice,iey,ey
read_slice,iez,ez

; Compute magntidue of total current density

fulldata = ex*jx + ey*jy + ez*jz

end

pro doplot,outform
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct
common colortable,rgb,usecolor,red,blue,green,range1,range2,r1,r2,tmax,tmin
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common imagesize,nxpix,nypix,xoff,yoff,xpic,ypic
common pdata,fulldata,mytitle
common refresh,update
common list,ns,order,suffix
common trajectory,xp1,zp1,xp2,zp2,xp1b,zp1b,xp2b,zp2b
common presentation,csize,showlabels

; Time Slice

nstep=0L
nstep = ns(order(v.time))
print,"Simulation Time=",nstep

; First check and see if we need to read in new data.  If we are just
; re-sizing, zooming, etc - there is no need to read & recompute everthing

if ( update eq 1) then begin
   
;  If raw plot is selected - then use this.  Otherwise do the
;                            diagnostic plot

   if (v.rawplot ge 1) then begin
      read_slice,v.rawplot,fulldata
      mytitle=rawdata(v.rawplot)
   endif else begin

; Compute the desired diagnostic

; Electron Anisotropy
      if v.computed eq 0 then anisotropy,0

; Ion Anisotropy
      if v.computed eq 1 then anisotropy,1

; Electron Agyrotropy
      if v.computed eq 2 then agyrotropy,0

; Ion Agyrotropy
      if v.computed eq 3 then agyrotropy,1

; Total Firehose condition
      if v.computed eq 4 then firehose,1

; Electron Firehose condition
      if v.computed eq 5 then firehose,0

; Simple KH condition
      if v.computed eq 6 then kelvin_helmholtz

; Magnitude of current density |J|
      if v.computed eq 7 then J_magnitude

; Magnitude of current density |B|
      if v.computed eq 8 then B_magnitude

; Electron Vorticity
      if v.computed eq 9 then vorticity,0

; Ion Vorticity
      if v.computed eq 10 then vorticity,1

; E*J
      if v.computed eq 11 then EdotJ

; beta-e
      if v.computed eq 12 then beta,0

; beta-i
      if v.computed eq 13 then beta,1

; beta-total
      if v.computed eq 14 then beta,2

; Curvature
      if v.computed eq 15 then curvature

; Net charge density
      if v.computed eq 16 then charge_density

; Electron temperature
      if v.computed eq 17 then temperature,0

; Electron temperature
      if v.computed eq 18 then temperature,1

   endelse
endif

; Shift the data in y if requested

fulldata = shift(fulldata,v.shift,0)

; Select region of 2d-data to plot - This allows me to create contours just within the specified
; region and ignore the rest of the data

; Declare memory for reduced region and load region into plotting array

imin = fix(v.xmin/xmax*nx)
imax = fix(v.xmax/xmax*nx)-1
lx = imax-imin +1 
jmin = fix(v.zmin/zmax*nz)
jmax = fix(v.zmax/zmax*nz)-1
lz = jmax-jmin +1 

;print , v.zmax,zmax,imin,imax,jmin,jmax,zmax

temp = fltarr(lx,lz)
temp(*,*) = fulldata(imin:imax,jmin:jmax)

xarr = ((v.xmax-v.xmin)/lx)*(findgen(lx)+0.5) + v.xmin
zarr = v.zmin + ((v.zmax-v.zmin)/lz)*(findgen(lz)+0.5) -zmax/2

; Output info about data

print,'Maximum Value=',max(abs(temp))
print,'Minimum Value=',min(abs(temp))

; set output for either X11 or a postscript file

if ( outform eq 0 ) then begin
    set_plot, 'X'
    device,true_color=24,decomposed=0	
    !x.style=1
    !y.style=1
    !P.color =0
    !p.background=1
    !p.font = 1
endif

if (outform eq 1 ) then begin
    set_plot, 'ps'
    !p.color =1
    !p.background = 0
    !p.font = 0	
    !x.style=1
    !y.style=1
    width=12.0
    asp=0.7
    height=width*asp
    if (v.rawplot ge 1) then fname = rawdata(v.rawplot) else fname=computedquantity(v.computed)
    file=strcompress(fname+'_y'+string(v.ycut)+'_t'+string(nstep)+'.eps',/remove_all)
    print,"Writing eps file --> ",file
    device, /encapsulated, bits_per_pixel=8, filename = file, $  
            /inches, xsize=width, ysize =height, /TIMES, font_size=12, /color, set_font='Times-Roman'
endif

;  Set colors so that red is largest value

if ( outform ne 2 ) then begin
    r1=min(temp)
    r2=max(temp)
    dr = (r2-r1)*0.40
    r1=r1-dr/2
    r2=r2+dr/2
endif

; Use limited range for KH diagnostic

If (v.computed eq 6) then begin
    r1 = -1.5
    r2 = 1.5
endif

; ********************************************************************************
; *********************** Contour Plot with Colorbar ****************************
; ********************************************************************************

if ( v.plottype eq 0 ) then begin

; Set position and titles

    xoff = 0.09
    yoff = 0.11
    xpic = 0.90
    ypic = 0.92
    dx1=0.012
    dx2=0.045
    if (showlabels eq 1) then !x.title="x"
    if (showlabels eq 1) then !y.title="z"
    !p.position=[xoff,yoff,xpic,ypic] 
    shade_surf, temp, xarr, zarr, ax=90,az=0,shades=bytscl(temp,max=r2,min=r1),zstyle=4,charsize=csize,pixels=1000
    if (showlabels eq 1) then xyouts,v.xmin+(v.xmax-v.xmin)/3.1,(v.zmax)*0.51,mytitle+"     it="+strcompress(string(FORMAT='(i6)',nstep))+"    T-slice="+strcompress(v.time)+"    Y-slice="+strcompress(string(v.ycut),/remove_all),charsize=csize

;  Now add the color bar

    if (showlabels eq 1) then !x.title=""
    if (showlabels eq 1) then !y.title=""
    colorbar, Position = [xpic+dx1,yoff,xpic+dx2,ypic], /Vertical, $
              Range=[r1,r2], Format='(f6.2)', /Right, $
              Bottom=5, ncolors=251, Divisions=6, font=9, charsize=0.8

endif

; **************************************************************************
; **********************  Contour Plot with X-average **********************
; **************************************************************************

if ( v.plottype eq 1 ) then begin

    xave = fltarr(lz)
    for j=0,lz-1 do begin
        xave(j) = 0.0
        for i=0,lx-1 do begin
            xave(j) = xave(j) +  temp(i,j)
        endfor
    endfor
    xave = xave/lx

    if ( outform eq 0) then begin
        range1=min(xave)
        range2=max(xave)
        dr = (range2-range1)*0.20
        range1=range1-dr/2
        range2=range2+dr/2
    endif

; Set position and titles

    xoff = 0.09
    yoff = 0.11
    xpic = 0.70
    ypic = 0.94

    !p.position=[xoff,yoff,xpic,ypic]
    if (showlabels eq 1) then !x.title="x"
    if (showlabels eq 1) then !y.title="z"
    !p.position=[xoff,yoff,xpic,ypic]     
    shade_surf, temp, xarr, zarr, ax=90,az=0,shades=bytscl(temp,max=r2,min=r1),zstyle=4,charsize=csize,pixels=1000
    if (showlabels eq 1) then xyouts,v.xmin+(v.xmax-v.xmin)/3.1,(v.zmax)*0.51,mytitle+"    it="+strcompress(string(FORMAT='(i6)',nstep))+"    T-slice="+strcompress(v.time)+"    Y-slice="+strcompress(string(v.ycut),/remove_all),charsize=csize

    xoff2 = 0.705
    yoff2 = 0.11
    xpic2 = 0.94
    ypic2 = 0.94
    !x.title=" "
    !y.title=' '
    !p.position=[xoff2,yoff2,xpic2,ypic2] 
    plot,xave,zarr,charsize=csize,yrange=[v.zmin-zmax/2,v.zmax-zmax/2],ytickname=replicate(' ',8),xrange=[range1,range2],/NOERASE

endif

; **************************************************************************
; **********************  Contour Plot with Z-Slice ************************
; **************************************************************************

if ( v.plottype eq 2 ) then begin

; Z-slice of data (fix x and vary z)

    slice = fltarr(lx)
    islice=fix((v.zslice-v.zmin)/(v.zmax-v.zmin)*float(lz)-0.5)
    for j=0,lx-1 do begin
        slice(j) = temp(j,islice)
    endfor

    if ( outform eq 0) then begin
        range1=min(slice)
        range2=max(slice)
        dr = (range2-range1)*0.20
        range1=range1-dr/2
        range2=range2+dr/2
    endif

; Set position and titles
    
    xoff = 0.09
    yoff = 0.40
    xpic = 0.94
    ypic = 0.94

    !p.position=[xoff,yoff,xpic,ypic]
    !x.title=""
    if (showlabels eq 1) then !y.title="z"
    !p.position=[xoff,yoff,xpic,ypic]     
    shade_surf, temp, xarr, zarr, ax=90,az=0,shades=bytscl(temp,max=r2,min=r1),zstyle=4,charsize=csize,pixels=1000,xtickname=replicate(' ',8)
    if (showlabels eq 1) then xyouts,v.xmin+(v.xmax-v.xmin)/3.1,(v.zmax)*0.51,mytitle+"    it ="+strcompress(string(FORMAT='(i6)',nstep))+"    T-slice="+strcompress(v.time)+"    Y-slice="+strcompress(string(v.ycut),/remove_all),charsize=csize

; Overplot a line to show where the x-slice is located

    xs=fltarr(lx)
    for j=0,lx-1 do begin
        xs(j) = v.zslice-zmax/2
    endfor
    oplot,xarr,xs,color=0

; Add z-slice plot

    xoff2 = 0.09
    yoff2 = 0.11
    xpic2 = 0.94
    ypic2 = 0.395

    if (showlabels eq 1) then !x.title="z"
    if (showlabels eq 1) then !y.title= mytitle
    !p.position=[xoff2,yoff2,xpic2,ypic2] 
    plot,xarr,slice,charsize=csize,xrange=[v.xmin,v.xmax],yrange=[range1,range2],/NOERASE
;    oplot,xarr,slice2,color=110
;    oplot,xarr,slice3,color=220

endif

; **************************************************************************
; **********************  Contour Plot with X-Slice ************************
; **************************************************************************

if ( v.plottype eq 3 ) then begin

; X-slice of data (fix z and vary x)

    jslice=fix((v.xslice-v.xmin)/(v.xmax-v.xmin)*float(lx)-0.5)
    slice = fltarr(lz)
    slice2 = fltarr(lz)
    slice3 = fltarr(lz)
    for i=0,lz-1 do begin
        slice(i) = temp(jslice,i)
    endfor

    if ( outform eq 0) then begin
        range1=min(slice)
        range2=max(slice)
        dr = (range2-range1)*0.30
        range1=range1-dr
        range2=range2+dr*1.15
    endif

; Set position and titles

    xoff = 0.09
    yoff = 0.11
    xpic = 0.70
    ypic = 0.94

    !p.position=[xoff,yoff,xpic,ypic]
    if (showlabels eq 1) then !x.title="x"
    if (showlabels eq 1) then !y.title="z"
    !p.position=[xoff,yoff,xpic,ypic]     
    shade_surf, temp, xarr, zarr, ax=90,az=0,shades=bytscl(temp,max=r2,min=r1),zstyle=4,charsize=csize
    if (showlabels eq 1) then xyouts,v.xmin+(v.xmax-v.xmin)/3.1,(v.zmax)*0.51,mytitle+"    it="+strcompress(string(FORMAT='(i6)',nstep))+"    T-slice="+strcompress(v.time)+"    Y-slice="+strcompress(string(v.ycut),/remove_all),charsize=csize

; Overplot a line to show where the z-slice is located

    xs=fltarr(lz)
    for i=0,lz-1 do begin
        xs(i) = v.xslice
    endfor
    oplot,xs,zarr,color=0

; Add z-slice plot

    xoff2 = 0.705
    yoff2 = 0.11
    xpic2 = 0.94
    ypic2 = 0.94

    !x.title=" "
    !y.title=" "
    !p.position=[xoff2,yoff2,xpic2,ypic2] 
    plot,slice,zarr,charsize=csize,yrange=[v.zmin-zmax/2,v.zmax-zmax/2],ytickname=replicate(' ',8),xrange=[range1,range2],/NOERASE
    oplot,slice2,zarr,color=110
    oplot,slice3,zarr,color=50

endif

; **************************************************************************
; ***************  Contour Plot with Arbitrary Trajectory ******************
; **************************************************************************
;
;  Use right mouse button to select points

if ( v.plottype eq 4 ) then begin

; Place trajectory plot to the left
;   xoff = 0.09
;   yoff = 0.11
;   xpic = 0.70
;   ypic = 0.94
; Place trajectory plot below

   xoff = 0.09
   yoff = 0.44
   xpic = 0.92
   ypic = 0.94
   dx1=0.012
   dx2=0.04

   if (showlabels eq 1) then !x.title="x"
   if (showlabels eq 1) then !y.title="z"
   !p.position=[xoff,yoff,xpic,ypic]
   shade_surf, temp, xarr, zarr, ax=90,az=0,shades=bytscl(temp,max=r2,min=r1),zstyle=4,charsize=csize
   if (showlabels eq 1) then xyouts,v.xmin+(v.xmax-v.xmin)/3.1,(v.zmax)*0.51,mytitle+"    it="+strcompress(string(FORMAT='(i6)',nstep))+"    T-slice="+strcompress(v.time)+"    Y-slice="+strcompress(string(v.ycut),/remove_all),charsize=csize

; Overplot a line to show trajectory

   num = lx
   x=fltarr(num)
   z=fltarr(num)
   s=fltarr(num)
   smag = sqrt((xp1-xp2)^2+(zp1-zp2)^2)
   s = smag*findgen(num)/num 
   x = xp1 + (xp2-xp1)*findgen(num)/num
   z = zp1 + (zp2-zp1)*findgen(num)/num 
   oplot,x,z,color=0

;  Now add the color bar

    if (showlabels eq 1) then !x.title=""
    if (showlabels eq 1) then !y.title=""
    colorbar, Position = [xpic+dx1,yoff,xpic+dx2,ypic], /Vertical, $
              Range=[r1,r2], Format='(f6.2)', /Right, $
              Bottom=5, ncolors=251, Divisions=6, font=9, charsize=0.8

; Now interpolate data along this line
; Use bilinear interpolation

   ix=0L
   ix1=0L
   iz=0L
   iz1=0L
   trajectory = fltarr(num)
   zero = fltarr(num)
   hx = double(xmax/nx)
   hz = double(zmax/nz)
   for j=0,num-1 do begin
      rx = (x(j)-v.xmin)/hx
      rz = (z(j)-v.zmin + zmax/2)/hz
      ix = long(rx)
      ix1 = ix + 1
      iz = long(rz)
      iz1 = iz + 1
      if (ix gt 0 and ix1 le lx-1 and iz gt 0 and iz1 le lz-1) then begin
         wx1 = rx  - double(ix)
         wx = 1.0 - wx1
         wz1 = rz - double(iz)
         wz = 1.0 - wz1
         trajectory(j) = temp(ix,iz)*wx*wz + temp(ix1,iz)*wx1*wz + temp(ix,iz1)*wx*wz1 + temp(ix1,iz1)*wx1*wz1
         zero(j)=0.0
      endif 
   endfor
   range1=min(trajectory)
   range2=max(trajectory)

; Plot data along trajectory

   if (showlabels eq 1) then !x.title="s"
   if (showlabels eq 1) then !y.title=mytitle
; Place trajectory plot to the left
;   !p.position=[0.735,0.11,0.995,0.94]
; Place trajectory plot below
    !p.position=[0.09,0.11,0.92,0.38]

    plot,s,trajectory,charsize=csize,/NOERASE
    oplot,s,zero,linestyle=1
endif

; **************************************************************************
; ***************   Contour Plot with 2 Trajectories   *********************
; **************************************************************************
;
;  Use right mouse button to select points

if ( v.plottype eq 5 ) then begin

; Place trajectory plot to the left
;   xoff = 0.09
;   yoff = 0.11
;   xpic = 0.70
;   ypic = 0.94
; Place trajectory plot below
   xoff = 0.09
   yoff = 0.44
   xpic = 0.92
   ypic = 0.94
   dx1=0.012
   dx2=0.04

   if (showlabels eq 1) then !x.title="x"
   if (showlabels eq 1) then !y.title="z"
   !p.position=[xoff,yoff,xpic,ypic]
   shade_surf, temp, xarr, zarr, ax=90,az=0,shades=bytscl(temp,max=r2,min=r1),zstyle=4,charsize=csize
   if (showlabels eq 1) then xyouts,v.xmin+(v.xmax-v.xmin)/3.1,(v.zmax)*0.51,mytitle+"    it="+strcompress(string(FORMAT='(i6)',nstep))+"    T-slice="+strcompress(v.time)+"    Y-slice="+strcompress(string(v.ycut),/remove_all),charsize=csize

; Overplot a line to show trajectory

   num = lx
   x=fltarr(num)
   z=fltarr(num)
   s=fltarr(num)
   smag = sqrt((xp1-xp2)^2+(zp1-zp2)^2)
   s = smag*findgen(num)/num 
   x = xp1 + (xp2-xp1)*findgen(num)/num
   z = zp1 + (zp2-zp1)*findgen(num)/num 
   oplot,x,z,color=0

   xb=fltarr(num)
   zb=fltarr(num)
   sb=fltarr(num)
   smag = sqrt((xp1b-xp2b)^2+(zp1b-zp2b)^2)
   sb = smag*findgen(num)/num 
   xb = xp1b + (xp2b-xp1b)*findgen(num)/num
   zb = zp1b + (zp2b-zp1b)*findgen(num)/num 
   oplot,xb,zb,color=0,linestyle=2

; Now interpolate data along this line
; Use bilinear interpolation

   ix=0L
   ix1=0L
   iz=0L
   iz1=0
   zero = fltarr(num)
   trajectory = fltarr(num)
   hx = double(xmax/nx)
   hz = double(zmax/nz)
   for j=0,num-1 do begin
      rx = (x(j)-v.xmin)/hx
      rz = (z(j)-v.zmin + zmax/2)/hz
      ix = long(rx)
      ix1 = ix + 1
      iz = long(rz)
      iz1 = iz + 1
      if (ix gt 0 and ix1 le lx-1 and iz gt 0 and iz1 le lz-1) then begin
         wx1 = rx  - double(ix)
         wx = 1.0 - wx1
         wz1 = rz - double(iz)
         wz = 1.0 - wz1
         trajectory(j) = temp(ix,iz)*wx*wz + temp(ix1,iz)*wx1*wz + temp(ix,iz1)*wx*wz1 + temp(ix1,iz1)*wx1*wz1
         zero(j) = 0.0
      endif 
   endfor
   range1=min(trajectory)
   range2=max(trajectory)

   trajectoryb = fltarr(num)
   for j=0,num-1 do begin
      rx = (xb(j)-v.xmin)/hx
      rz = (zb(j)-v.zmin + zmax/2)/hz
      ix = long(rx)
      ix1 = ix + 1
      iz = long(rz)
      iz1 = iz + 1
      if (ix gt 0 and ix1 le lx-1 and iz gt 0 and iz1 le lz-1) then begin
         wx1 = rx  - double(ix)
         wx = 1.0 - wx1
         wz1 = rz - double(iz)
         wz = 1.0 - wz1
         trajectoryb(j) = temp(ix,iz)*wx*wz + temp(ix1,iz)*wx1*wz + temp(ix,iz1)*wx*wz1 + temp(ix1,iz1)*wx1*wz1
      endif 
   endfor

;  Now add the color bar

    if (showlabels eq 1) then !x.title=""
    if (showlabels eq 1) then !y.title=""
    colorbar, Position = [xpic+dx1,yoff,xpic+dx2,ypic], /Vertical, $
              Range=[r1,r2], Format='(f6.2)', /Right, $
              Bottom=5, ncolors=251, Divisions=6, font=9, charsize=0.8

; Plot data along trajectory

   if (showlabels eq 1) then !x.title="s"
   if (showlabels eq 1) then !y.title=mytitle
; Place trajectory plot to the left
;   !p.position=[0.735,0.11,0.995,0.94]
; Place trajectory plot below
    !p.position=[0.09,0.11,0.92,0.38]

    plot,s,trajectory,charsize=csize,/NOERASE
    oplot,sb,trajectoryb,linestyle=2
    oplot,sb,zero,linestyle=1
endif


; Close postscript device

if (outform eq 1 ) then device, /close

; Also dump as TIFF
;if (outform eq 1 ) then begin
;   if (v.rawplot ge 1) then fname = rawdata(v.rawplot) else fname=computedquantity(v.computed)
;    file=strcompress(fname+'_y'+string(v.ycut)+'_t'+string(nstep)+'.eps',/remove_all)
;    print,"Writing eps file --> ",file
;    image = TVRead(Filename=file,/NODIALOG,/TIFF,quality=100)
;endif    

; Set refresh to false (i.e. won't read back in data unless needed)

update = 0

end


pro handle_event, ev
common controlplot,v
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq, anames
common imagesize,nxpix,nypix,xoff,yoff,xpic,ypic
common oldposition,xmin0,zmax0
common refresh,update
common trajectory,xp1,zp1,xp2,zp2,xp1b,zp1b,xp2b,zp2b

widget_control, ev.id, get_uval = whichbutton
case whichbutton of
    'done' : begin
        close,/All
        free_lun, 1
        widget_control, ev.top, /destroy
    end
    'plot' : begin
        widget_control,ev.top,get_uvalue=v

        v.zmin=0.
        v.zmax=zmax
        v.xmin=0.
        v.xmax=xmax
        widget_control,ev.top,set_uvalue=v    
        doplot,0
    end
    'render' : begin
        widget_control,ev.top,get_uvalue=v
        doplot,1
    end
    'animate' : begin
        widget_control,ev.top,get_uvalue=v
        create_movie
    end
    'region' : begin
        widget
    end

    'mouse' : begin

       widget_control,ev.top,get_uvalue=v

; left mouse button
       if (ev.press eq 1) then begin
          xmin0=v.xmin
          zmax0=v.zmax
          v.xmin=(float(ev.x)/float(nxpix)-xoff)*(v.xmax-v.xmin)/(xpic-xoff)+v.xmin
          v.zmax=(float(ev.y)/float(nypix)-yoff)*(v.zmax-v.zmin)/(ypic-yoff)+v.zmin
          if (v.zmax gt zmax) then v.zmax=zmax
          if (v.xmin lt 0.0) then v.xmin=0.
          widget_control,ev.top,set_uvalue=v    
          print,"xmin=",v.xmin," zmax=",v.zmax-zmax/2
       endif
       if (ev.release eq 1) then begin
          v.xmax=(float(ev.x)/float(nxpix)-xoff)*(v.xmax-xmin0)/(xpic-xoff)+xmin0
          v.zmin=(float(ev.y)/float(nypix)-yoff)*(zmax0-v.zmin)/(ypic-yoff)+v.zmin
          if (v.zmin lt 0 ) then v.zmin=0
          if (v.xmax gt xmax) then v.xmax=xmax
          widget_control,ev.top,set_uvalue=v  
          print,"xmax=",v.xmax," zmin=",v.zmin-zmax/2  
          doplot,0
       endif

; Right button
       if (ev.press eq 2) then begin
          xp1b =(float(ev.x)/float(nxpix)-xoff)*(v.xmax-v.xmin)/(xpic-xoff) + v.xmin
          zp1b =(float(ev.y)/float(nypix)-yoff)*(v.zmax-v.zmin)/(ypic-yoff) + v.zmin - zmax/2
          if (xp1b lt v.xmin) then xp1b = v.xmin
          if (xp1b gt v.xmax) then xp1b = v.xmax
          if (zp1b + zmax/2 gt v.zmax) then zp1b = v.zmax - zmax/2
          if (zp1b + zmax/2 lt v.zmin) then zp1b = v.zmin - zmax/2
          print,"Press middle Button --> x=",xp1b,"  z=",zp1b
       endif
       
       if (ev.release eq 2) then begin
          xp2b=(float(ev.x)/float(nxpix)-xoff)*(v.xmax-v.xmin)/(xpic-xoff) + v.xmin
          zp2b=(float(ev.y)/float(nypix)-yoff)*(v.zmax-v.zmin)/(ypic-yoff) + v.zmin - zmax/2

          if (xp2b lt v.xmin) then xp2b = v.xmin
          if (xp2b gt v.xmax) then xp2b = v.xmax
          if (zp2b + zmax/2 gt v.zmax) then zp2b = v.zmax - zmax/2
          if (zp2b + zmax/2 lt v.zmin) then zp2b = v.zmin - zmax/2

          print,"Release middle button --> x=",xp2b,"  z=",zp2b
          doplot,0
       endif
       
; Right button
       if (ev.press eq 4) then begin
          xp1=(float(ev.x)/float(nxpix)-xoff)*(v.xmax-v.xmin)/(xpic-xoff) + v.xmin
          zp1=(float(ev.y)/float(nypix)-yoff)*(v.zmax-v.zmin)/(ypic-yoff) + v.zmin - zmax/2
          if (xp1 lt v.xmin) then xp1 = v.xmin
          if (xp1 gt v.xmax) then xp1 = v.xmax
          if (zp1 + zmax/2 gt v.zmax) then zp1 = v.zmax - zmax/2
          if (zp1 + zmax/2 lt v.zmin) then zp1 = v.zmin - zmax/2
          print,"Press right button --> x=",xp1,"  z=",zp1
       endif
       
       if (ev.release eq 4) then begin
          xp2=(float(ev.x)/float(nxpix)-xoff)*(v.xmax-v.xmin)/(xpic-xoff) + v.xmin
          zp2=(float(ev.y)/float(nypix)-yoff)*(v.zmax-v.zmin)/(ypic-yoff) + v.zmin - zmax/2

          if (xp2 lt v.xmin) then xp2 = v.xmin
          if (xp2 gt v.xmax) then xp2 = v.xmax
          if (zp2 + zmax/2 gt v.zmax) then zp2 = v.zmax - zmax/2
          if (zp2 + zmax/2 lt v.zmin) then zp2 = v.zmin - zmax/2

          print,"Release right button --> x=",xp2,"  z=",zp2
          doplot,0
       endif
       
    
    end

    'computed' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.computed=ev.index
        widget_control,ev.top,set_uvalue=v    
        update=1
    end

    'angle' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.angle=ev.value
        widget_control,ev.top,set_uvalue=v    
        update=1
    end

    'rawplot' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.rawplot=ev.index
        widget_control,ev.top,set_uvalue=v    
        update=1
    end

    'options' : begin
        widget_control,ev.top,get_uvalue=v                	
        print," Selected Option = ",ev.value       
        if ( ev.value ge 2 ) and ( ev.value le 9 ) then v.smoothing=ev.value-1
        if ( ev.value ge 11 ) and ( ev.value le 30) then v.contours=(ev.value*2)-10
        if ( ev.value eq 31 ) then v.data=0 
        if ( ev.value eq 32 ) then v.data=1
        if ( ev.value ge 35 ) then v.map=35-ev.value
        if ( ev.value eq 36 ) then xloadct,bottom=5
        widget_control,ev.top,set_uvalue=v    
        update=1        
    end

    'ptype' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.plottype=ev.index
        widget_control,ev.top,set_uvalue=v
        doplot,0
     end

    'time' : begin
        widget_control,ev.top,get_uvalue=v                	 	
        v.time=ev.value
        widget_control,ev.top,set_uvalue=v 
        update=1
        doplot,0
    end

    'yplane' : begin
        widget_control,ev.top,get_uvalue=v                	 	
        v.ycut=ev.value
        print,"v.ycut=",v.ycut
        widget_control,ev.top,set_uvalue=v 
        update=1
        doplot,0
    end

    'shift' : begin
        widget_control,ev.top,get_uvalue=v                	 	
        v.shift=ev.value
        widget_control,ev.top,set_uvalue=v 
        doplot,0
    end

    'zmax' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.zmax=ev.value + zmax/2
        v.zmin=zmax/2  - ev.value
        widget_control,ev.top,set_uvalue=v    
    end

    'xslice' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.xslice=ev.value
        widget_control,ev.top,set_uvalue=v    
    end

    'zslice' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.zslice=ev.value+zmax/2
        widget_control,ev.top,set_uvalue=v    
    end

    'xmax' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.xmax=ev.value
        widget_control,ev.top,set_uvalue=v    
    end

    'xmin' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.xmin=ev.value
        widget_control,ev.top,set_uvalue=v    
    end

    'ymax' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.ymax=ev.value
        widget_control,ev.top,set_uvalue=v    
    end

    'ymin' : begin
        widget_control,ev.top,get_uvalue=v                	
        v.ymin=ev.value
        widget_control,ev.top,set_uvalue=v    
    end


endcase

end

; NAME: COLORBAR
; PURPOSE:
;       The purpose of this routine is to add a color bar to the current
;       graphics window.
; AUTHOR:
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   2642 Bradbury Court
;   Fort Collins, CO 80521 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
; CATEGORY:
;       Graphics, Widgets.
; CALLING SEQUENCE:
;       COLORBAR
; INPUTS:
;       None.
; KEYWORD PARAMETERS:
;       BOTTOM:   The lowest color index of the colors to be loaded in
;                 the bar.
;       CHARSIZE: The character size of the color bar annotations. Default is 1.0.
;
;       COLOR:    The color index of the bar outline and characters. Default
;                 is !P.Color..
;       DIVISIONS: The number of divisions to divide the bar into. There will
;                 be (divisions + 1) annotations. The default is 6.
;       FONT:     Sets the font of the annotation. Hershey: -1, Hardware:0, 
;                 True-Type: 1.
;       FORMAT:   The format of the bar annotations. Default is '(I5)'.
;       MAXRANGE: The maximum data value for the bar annotation. Default is
;                 NCOLORS.
;       MINRANGE: The minimum data value for the bar annotation. Default is 0.
;       MINOR:    The number of minor tick divisions. Default is 2.
;       NCOLORS:  This is the number of colors in the color bar.
;
;       POSITION: A four-element array of normalized coordinates in the same
;                 form as the POSITION keyword on a plot. Default is
;                 [0.88, 0.15, 0.95, 0.95] for a vertical bar and
;                 [0.15, 0.88, 0.95, 0.95] for a horizontal bar.
;       RANGE:    A two-element vector of the form [min, max]. Provides an
;                 alternative way of setting the MINRANGE and MAXRANGE keywords.
;       RIGHT:    This puts the labels on the right-hand side of a vertical
;                 color bar. It applies only to vertical color bars.
;       TITLE:    This is title for the color bar. The default is to have
;                 no title.
;       TOP:      This puts the labels on top of the bar rather than under it.
;                 The keyword only applies if a horizontal color bar is rendered.
;       VERTICAL: Setting this keyword give a vertical color bar. The default
;                 is a horizontal color bar.
; COMMON BLOCKS:
;       None.
; SIDE EFFECTS:
;       Color bar is drawn in the current graphics window.
; RESTRICTIONS:
;       The number of colors available on the display device (not the
;       PostScript device) is used unless the NCOLORS keyword is used.
; EXAMPLE:
;       To display a horizontal color bar above a contour plot, type:
;       LOADCT, 5, NCOLORS=100
;       CONTOUR, DIST(31,41), POSITION=[0.15, 0.15, 0.95, 0.75], $
;          C_COLORS=INDGEN(25)*4, NLEVELS=25
;       COLORBAR, NCOLORS=100, POSITION=[0.15, 0.85, 0.95, 0.90]
; MODIFICATION HISTORY:
;       Written by: David Fanning, 10 JUNE 96.
;       10/27/96: Added the ability to send output to PostScript. DWF
;       11/4/96: Substantially rewritten to go to screen or PostScript
;           file without having to know much about the PostScript device
;           or even what the current graphics device is. DWF
;       1/27/97: Added the RIGHT and TOP keywords. Also modified the
;            way the TITLE keyword works. DWF
;       7/15/97: Fixed a problem some machines have with plots that have
;            no valid data range in them. DWF
;       12/5/98: Fixed a problem in how the colorbar image is created that
;            seemed to tickle a bug in some versions of IDL. DWF.
;       1/12/99: Fixed a problem caused by RSI fixing a bug in IDL 5.2. 
;                Sigh... DWF.
;       3/30/99: Modified a few of the defaults. DWF.
;       3/30/99: Used NORMAL rather than DEVICE coords for positioning bar. DWF.
;       3/30/99: Added the RANGE keyword. DWF.
;       3/30/99: Added FONT keyword. DWF
;       5/6/99: Many modifications to defaults. DWF.
;       5/6/99: Removed PSCOLOR keyword. DWF.
;       5/6/99: Improved error handling on position coordinates. DWF.
;       5/6/99. Added MINOR keyword. DWF.
;       5/6/99: Set Device, Decomposed=0 if necessary. DWF.
;       2/9/99: Fixed a problem caused by setting BOTTOM keyword, but not 
;               NCOLORS. DWF.
;       8/17/99. Fixed a problem with ambiguous MIN and MINOR keywords. DWF
;       8/25/99. I think I *finally* got the BOTTOM/NCOLORS thing sorted out. 
;                :-( DWF.
;       10/10/99. Modified the program so that current plot and map 
;                 coordinates are
;            saved and restored after the colorbar is drawn. DWF.
;       3/18/00. Moved a block of code to prevent a problem with color 
;                decomposition. DWF.
;       4/28/00. Made !P.Font default value for FONT keyword. DWF.
;###########################################################################
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright Å© 2000 Fanning Software Consulting.
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
; 3. This notice may not be removed or altered from any source distribution.
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;###########################################################################
PRO COLORBAR, BOTTOM=bottom, CHARSIZE=charsize, COLOR=color, $
   DIVISIONS=divisions, $
   FORMAT=format, POSITION=position, MAXRANGE=maxrange, MINRANGE=minrange, $
   NCOLORS=ncolors, $
   TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right, MINOR=minor, $
   RANGE=range, FONT=font, TICKLEN=ticklen, _EXTRA=extra

   ; Return to caller on error.
On_Error, 2

   ; Save the current plot state.

bang_p = !P
bang_x = !X
bang_Y = !Y
bang_Z = !Z
bang_Map = !Map

   ; Is the PostScript device selected?

postScriptDevice = (!D.NAME EQ 'PS' OR !D.NAME EQ 'PRINTER')

   ; Which release of IDL is this?

thisRelease = Float(!Version.Release)

    ; Check and define keywords.

IF N_ELEMENTS(ncolors) EQ 0 THEN BEGIN

   ; Most display devices to not use the 256 colors available to
   ; the PostScript device. This presents a problem when writing
   ; general-purpose programs that can be output to the display or
   ; to the PostScript device. This problem is especially bothersome
   ; if you don't specify the number of colors you are using in the
   ; program. One way to work around this problem is to make the
   ; default number of colors the same for the display device and for
   ; the PostScript device. Then, the colors you see in PostScript are
   ; identical to the colors you see on your display. Here is one way to
   ; do it.

   IF postScriptDevice THEN BEGIN
      oldDevice = !D.NAME

         ; What kind of computer are we using? SET_PLOT to appropriate
         ; display device.

      thisOS = !VERSION.OS_FAMILY
      thisOS = STRMID(thisOS, 0, 3)
      thisOS = STRUPCASE(thisOS)
      CASE thisOS of
         'MAC': SET_PLOT, thisOS
         'WIN': SET_PLOT, thisOS
         ELSE: SET_PLOT, 'X'
      ENDCASE

         ; Here is how many colors we should use.

      ncolors = !D.TABLE_SIZE
      SET_PLOT, oldDevice
    ENDIF ELSE ncolors = !D.TABLE_SIZE
ENDIF
IF N_ELEMENTS(bottom) EQ 0 THEN bottom = 0B
IF N_ELEMENTS(charsize) EQ 0 THEN charsize = 1.0
IF N_ELEMENTS(format) EQ 0 THEN format = '(I5)'
IF N_ELEMENTS(color) EQ 0 THEN color = !P.Color
IF N_ELEMENTS(minrange) EQ 0 THEN minrange = 0
IF N_ELEMENTS(maxrange) EQ 0 THEN maxrange = ncolors
IF N_ELEMENTS(ticklen) EQ 0 THEN ticklen = 0.2
IF N_ELEMENTS(minor) EQ 0 THEN minor = 2
IF N_ELEMENTS(range) NE 0 THEN BEGIN
   minrange = range[0]
   maxrange = range[1]
ENDIF
IF N_ELEMENTS(divisions) EQ 0 THEN divisions = 6
IF N_ELEMENTS(font) EQ 0 THEN font = !P.Font
IF N_ELEMENTS(title) EQ 0 THEN title = ''

IF KEYWORD_SET(vertical) THEN BEGIN
   bar = REPLICATE(1B,20) # BINDGEN(ncolors)
   IF N_ELEMENTS(position) EQ 0 THEN BEGIN
      position = [0.88, 0.1, 0.95, 0.9]
   ENDIF ELSE BEGIN
      IF position[2]-position[0] GT position[3]-position[1] THEN BEGIN
         position = [position[1], position[0], position[3], position[2]]
      ENDIF
      IF position[0] GE position[2] THEN Message, "Position coordinates can't be reconciled."
      IF position[1] GE position[3] THEN Message, "Position coordinates can't be reconciled."
   ENDELSE
ENDIF ELSE BEGIN
   bar = BINDGEN(ncolors) # REPLICATE(1B, 20)
   IF N_ELEMENTS(position) EQ 0 THEN BEGIN
      position = [0.1, 0.88, 0.9, 0.95]
   ENDIF ELSE BEGIN
      IF position[3]-position[1] GT position[2]-position[0] THEN BEGIN
         position = [position[1], position[0], position[3], position[2]]
      ENDIF
      IF position[0] GE position[2] THEN Message, "Position coordinates can't be reconciled."
      IF position[1] GE position[3] THEN Message, "Position coordinates can't be reconciled."
   ENDELSE
ENDELSE

   ; Scale the color bar.

 bar = BYTSCL(bar, TOP=(ncolors-1 < (255-bottom))) + bottom

   ; Get starting locations in NORMAL coordinates.

xstart = position(0)
ystart = position(1)

   ; Get the size of the bar in NORMAL coordinates.

xsize = (position(2) - position(0))
ysize = (position(3) - position(1))

   ; Display the color bar in the window. Sizing is
   ; different for PostScript and regular display.

IF postScriptDevice THEN BEGIN

   TV, bar, xstart, ystart, XSIZE=xsize, YSIZE=ysize, /Normal

ENDIF ELSE BEGIN

   bar = CONGRID(bar, CEIL(xsize*!D.X_VSize), CEIL(ysize*!D.Y_VSize), /INTERP)

        ; Decomposed color off if device supports it.

   CASE  StrUpCase(!D.NAME) OF
        'X': BEGIN
            IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
            Device, Decomposed=0
            ENDCASE
        'WIN': BEGIN
            IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
            Device, Decomposed=0
            ENDCASE
        'MAC': BEGIN
            IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
            Device, Decomposed=0
            ENDCASE
        ELSE:
   ENDCASE

   TV, bar, xstart, ystart, /Normal

      ; Restore Decomposed state if necessary.

   CASE StrUpCase(!D.NAME) OF
      'X': BEGIN
         IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
         ENDCASE
      'WIN': BEGIN
         IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
         ENDCASE
      'MAC': BEGIN
         IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
         ENDCASE
      ELSE:
   ENDCASE

ENDELSE

   ; Annotate the color bar.

IF KEYWORD_SET(vertical) THEN BEGIN

   IF KEYWORD_SET(right) THEN BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=1, $
         YTICKS=divisions, XSTYLE=1, YSTYLE=9, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', YTICKLEN=ticklen , $
         YRANGE=[minrange, maxrange], FONT=font, _EXTRA=extra, YMINOR=minor

      AXIS, YAXIS=1, YRANGE=[minrange, maxrange], YTICKFORMAT=format, YTICKS=divisions, $
         YTICKLEN=ticklen, YSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         FONT=font, YTITLE=title, _EXTRA=extra, YMINOR=minor

   ENDIF ELSE BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=1, $
         YTICKS=divisions, XSTYLE=1, YSTYLE=9, YMINOR=minor, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT=format, XTICKFORMAT='(A1)', YTICKLEN=ticklen , $
         YRANGE=[minrange, maxrange], FONT=font, YTITLE=title, _EXTRA=extra

      AXIS, YAXIS=1, YRANGE=[minrange, maxrange], YTICKFORMAT='(A1)', YTICKS=divisions, $
         YTICKLEN=ticklen, YSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         FONT=font, _EXTRA=extra, YMINOR=minor

   ENDELSE

ENDIF ELSE BEGIN

   IF KEYWORD_SET(top) THEN BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=divisions, $
         YTICKS=1, XSTYLE=9, YSTYLE=1, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', XTICKLEN=ticklen, $
         XRANGE=[minrange, maxrange], FONT=font, _EXTRA=extra, XMINOR=minor

      AXIS, XTICKS=divisions, XSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         XTICKFORMAT=format, XTICKLEN=ticklen, XRANGE=[minrange, maxrange], XAXIS=1, $
         FONT=font, XTITLE=title, _EXTRA=extra, XCHARSIZE=charsize, XMINOR=minor

   ENDIF ELSE BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=divisions, $
         YTICKS=1, XSTYLE=1, YSTYLE=1, TITLE=title, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT=format, XTICKLEN=ticklen, $
         XRANGE=[minrange, maxrange], FONT=font, XMinor=minor, _EXTRA=extra

    ENDELSE

ENDELSE

   ; Restore the previous plot and map system variables.

!P = bang_p
!X = bang_x
!Y = bang_y
!Z = bang_z
!Map = bang_map

END


pro create_movie
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,zcenter,numq,anames
common picdata, field, struct
common colortable,rgb,usecolor,red,blue,green,range1,range2,r1,r2,tmax,tmin
common controlplot,v
common imagesize,nxpix,nypix,xoff,yoff,xpic,ypic
common refresh,update

;  Loop over all time slices

; for itime = 0,3,1 do begin
for itime = 0,tslices do begin
    v.time = itime
    update =1

; Create Plot

    doplot,2

; Save JPEG frame

    if (v.rawplot ge 1) then fname = rawdata(v.rawplot) else fname=computedquantity(v.computed)
    file=strcompress('temp/'+fname+'-y'+string(v.ycut)+'.'+string(FORMAT='(i4.4)',itime),/remove_all)
    print,"Saving Time slice=",itime," File name = ",file
    image = TVRead(Filename=file,/NODIALOG,/JPEG,quality=100)

; end of time loop

endfor

;  end of subroutine
end

pro diagnostic3D
common choice, rawdata,computedquantity
common parameters,nx,ny,nz,tslices,xmax,ymax,zmax,xcenter,numq, anames
common pdata,fulldata,mytitle
common picdata, field, struct
common colortable,rgb,usecolor,red,blue,green,range1,range2,r1,r2,tmax,tmin
common controlplot,v
common labels,ipexx,ipexy,ipexz,ipeyy,ipeyz,ipezz,ipixx,ipixy,ipixz,ipiyy,ipiyz,ipizz,ibx,iby,ibz,ine,ini,iuix,iuiy,iuiz,iuex,iuey,iuez,ijx,ijy,ijz,iex,iey,iez
common imagesize,nxpix,nypix,xoff,yoff,xpic,ypic
common refresh,update
common dir,directory
common list,ns,order,suffix
common trajectory,xp1,zp1,xp2,zp2,xp1b,zp1b,xp2b,zp2b
common presentation,csize,showlabels

; Alter presentation for X vs eps

csize = 1.5
showlabels = 0

; First time starting

update = 1

; First determine the screen size

dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
print,"Dimensions of screen=",dimensions

; Pick size of viewer based on fraction of screen

nxpix = dimensions(0)*0.80
nypix = dimensions(1)*0.60

; Hard code in your choice (better for movies)

;nxpix=1100
;nypix=800
;nxpix=900
;nypix=540

; Declare structure for controlling the plots

v={quantity:0,time:0,xmin:0.0,xmax:0.0,ymin:0.0,ymax:0.0,xslice:0.0,zmin:0.0,zmax:0.0,zslice:0.0,smoothing:1,contours:24,plottype:0,rawplot:0,shift:0,data:0,map:0,xcut:1,ycut:0,zcut:1,computed:0,angle:0.0}

; Declare strings for menu

ptype = strarr(6)  
porient = strarr(3)
anames = strarr(3,2)

; Read in my color table

nc=256
rgb=fltarr(3,nc)
usecolor=fltarr(3,nc)
red=fltarr(nc)
blue=fltarr(nc)
green=fltarr(nc)
openr, unit, '~/Color_Tables/c5.tbl', /get_lun
readf, unit, rgb
close, unit
free_lun, unit
red=rgb(0,*)
green=rgb(1,*)
blue=rgb(2,*)
tvlct,red,green,blue

; Declare long integers 

nx=0L
ny=0L
nz=0L
numq=0L

; Determine if there are multiple data directories and allow user to select

datafiles = file_search('dat*',count=numq) 
if (numq gt 1) then directory = dialog_pickfile(/directory,filter='data*',TITLE='Choose directory with data') else  directory = 'data/'
print,"Reading data in",directory
;directory='data'    


; open binary file for problem description 

if ( file_test(directory+'info') eq 1 ) then begin
    print," *** Found Data Information File *** "
endif else begin
    print," *** Error - File info is missing ***"
endelse

;  First Try little endian then switch to big endian if this fails

little = 0
on_ioerror, switch_endian
;openr,unit,directory+'info',/f77_unformatted, /get_lun,/swap_if_big_endian
openr,unit,directory+'/info', /get_lun,/swap_if_big_endian
readu,unit,nx,ny,nz
little=1
switch_endian: if (not little) then begin
    print, " ** Little endian failed --> Switch to big endian"
    close,unit
    free_lun,unit
    openr,unit,directory+'/info',/f77_unformatted, /get_lun,/swap_if_little_endian
    readu,unit,nx,ny,nz
endif
;nx = 1024
;ny = 1024
;nz = 512
;xmax = 700.
;ymax = 700.
;zmax = 350.
; Read the problem desciption

on_ioerror, halt_error1
readu,unit,xmax,ymax,zmax
close,unit
free_lun, unit

; Find the names of data files in the data directory

;datafiles = file_search(directory+'/*.gda',count=numq) 

; Default points for directory

zp1=0.0
xp1=0.0
zp2=0.0
xp2=xmax
zp1b=0.0
xp1b=0.0
zp2b=0.0
xp2b=xmax

; Do some magic and find time slices

suffix = file_search(directory+'bx*.gda',count=tslices) 
if (tslices eq 0) then suffix = file_search(directory+'bx*.gda',count=tslices) 
print, suffix
ns = lonarr(tslices)
for i=0,tslices-1 do begin
;    print,i," ",suffix(i)
    suffix(i) = strmid(suffix(i),3+strlen(directory))
;    print,"suffix=",suffix(i)
    split = strsplit(suffix(i),'.',/EXTRACT)
    ns(i) = string(split(0))
;    print,i," suffix=",suffix(i),ns(i)
endfor
order = sort(ns)
;for i=0,tslices-1 do begin
;    print,i," suffix=",suffix(order(i))
;endfor

; Check and see if there is a filelist file that specifies names for
; the data files - we will use this if available

if ( file_test('filelist') eq 1 ) then begin
    openr,5,'filelist'
    numq=0
    name ='name for the data file'
    while not eof(5) do begin
        readf,5,name
        numq = numq+1
    endwhile
    close,5
    print,"Found the filelist file -->",numq," filenames"
    openr,5,'filelist'
    rawdata = strarr(numq+1)  
    rawdata(0) = "None"
    for i=1,numq do begin
        readf,5,name
        rawdata(i) = name
        print,i," ",rawdata(i)
    endfor
endif else begin

; Manutally specify base file names

    numq=32
    rawdata = strarr(numq+1)  
    rawdata(0) = "None"
    rawdata(1) = "bx"
    rawdata(2) = "by"
    rawdata(3) = "bz"
    rawdata(4) = "ex"
    rawdata(5) = "ey"
    rawdata(6) = "ez"
    rawdata(7) = "ni"
    rawdata(8) = "ex"
    rawdata(9) = "ey"
    rawdata(10) = "ez"
    rawdata(11) = "jx"
    rawdata(12) = "jy"
    rawdata(13) = "jz"
    rawdata(14) = "uix"
    rawdata(15) = "uiy"
    rawdata(16) = "uiz"
    rawdata(17) = "uex"
    rawdata(18) = "uey"
    rawdata(19) = "uez"
    rawdata(20) = "pe-xx"
    rawdata(21) = "pe-xy"
    rawdata(22) = "pe-xz"
    rawdata(23) = "pe-yy"
    rawdata(24) = "pe-yz"
    rawdata(25) = "pe-zz"
    rawdata(26) = "pi-xx"
    rawdata(27) = "pi-xy"
    rawdata(28) = "pi-xz"
    rawdata(29) = "pi-yy"
    rawdata(30) = "pi-yz"
    rawdata(31) = "pi-zz"
endelse

; Echo file names & match to specified qnantities for analysis routines
print," Number of files=",numq
instring='     '
for i=1,numq do begin
    print,"i=",i," --> ",rawdata(i)
    if ( rawdata(i) eq 'Pe-xx' or rawdata(i) eq 'pe-xx') then ipexx =i
    if ( rawdata(i) eq 'Pe-xy' or rawdata(i) eq 'pe-xy') then ipexy =i
    if ( rawdata(i) eq 'Pe-xz' or rawdata(i) eq 'pe-xz') then ipexz =i
    if ( rawdata(i) eq 'Pe-yy' or rawdata(i) eq 'pe-yy') then ipeyy =i
    if ( rawdata(i) eq 'Pe-yz' or rawdata(i) eq 'pe-yz') then ipeyz =i
    if ( rawdata(i) eq 'Pe-zz' or rawdata(i) eq 'pe-zz') then ipezz =i
    if ( rawdata(i) eq 'Pi-xx' or rawdata(i) eq 'pi-xx') then ipixx =i
    if ( rawdata(i) eq 'Pi-xy' or rawdata(i) eq 'pi-xy') then ipixy =i
    if ( rawdata(i) eq 'Pi-xz' or rawdata(i) eq 'pi-xz') then ipixz =i
    if ( rawdata(i) eq 'Pi-yy' or rawdata(i) eq 'pi-yy') then ipiyy =i
    if ( rawdata(i) eq 'Pi-yz' or rawdata(i) eq 'pi-yz') then ipiyz =i
    if ( rawdata(i) eq 'Pi-zz' or rawdata(i) eq 'pi-zz') then ipizz =i
    if ( rawdata(i) eq 'Bx' or rawdata(i) eq 'bx' ) then ibx =i
    if ( rawdata(i) eq 'By' or rawdata(i) eq 'by' ) then iby =i
    if ( rawdata(i) eq 'Bz' or rawdata(i) eq 'bz') then ibz =i
    if ( rawdata(i) eq 'Ex' or rawdata(i) eq 'ex') then iex =i
    if ( rawdata(i) eq 'Ey' or rawdata(i) eq 'ey') then iey =i
    if ( rawdata(i) eq 'Ez' or rawdata(i) eq 'ez' ) then iez =i
    if ( rawdata(i) eq 'Jx' or rawdata(i) eq 'jx' ) then ijx =i
    if ( rawdata(i) eq 'Jy' or rawdata(i) eq 'jy' ) then ijy =i
    if ( rawdata(i) eq 'Jz' or rawdata(i) eq 'jz' ) then ijz =i
    if ( rawdata(i) eq 'ni' ) then ini =i
    if ( rawdata(i) eq 'ne' ) then ine =i
    if ( rawdata(i) eq 'Uix' or rawdata(i) eq 'uix' ) then iuix =i
    if ( rawdata(i) eq 'Uiy' or rawdata(i) eq 'uiy') then iuiy =i
    if ( rawdata(i) eq 'Uiz' or rawdata(i) eq 'uiz') then iuiz =i
    if ( rawdata(i) eq 'Uex' or rawdata(i) eq 'uex') then iuex =i
    if ( rawdata(i) eq 'Uey' or rawdata(i) eq 'uey') then iuey =i
    if ( rawdata(i) eq 'Uez' or rawdata(i) eq 'uez') then iuez =i
 endfor

; Check that we have found all the files needed

  if (n_elements(ibx) eq 0) then print," ****** WARNING - Missing Bx data"
  if (n_elements(iby) eq 0) then print," ****** WARNING - Missing By data"
  if (n_elements(ibz) eq 0) then print," ****** WARNING - Missing Bz data"
  if (n_elements(ine) eq 0) then print," ****** WARNING - Missing ne data"
  if (n_elements(ipexx) eq 0) then print," ****** WARNING - Missing Pe-xx data"
  if (n_elements(ipexy) eq 0) then print," ****** WARNING - Missing Pe-xy data"
  if (n_elements(ipexz) eq 0) then print," ****** WARNING - Missing Pe-xz data"
  if (n_elements(ipeyy) eq 0) then print," ****** WARNING - Missing Pe-yy data"
  if (n_elements(ipeyz) eq 0) then print," ****** WARNING - Missing Pe-yz data"
  if (n_elements(ipezz) eq 0) then print," ****** WARNING - Missing Pe-zz data"
  if (n_elements(ipixx) eq 0) then print," ****** WARNING - Missing Pi-xx data"
  if (n_elements(ipixy) eq 0) then print," ****** WARNING - Missing Pi-xy data"
  if (n_elements(ipixz) eq 0) then print," ****** WARNING - Missing Pi-xz data"
  if (n_elements(ipiyy) eq 0) then print," ****** WARNING - Missing Pi-yy data"
  if (n_elements(ipiyz) eq 0) then print," ****** WARNING - Missing Pi-yz data"
  if (n_elements(ipizz) eq 0) then print," ****** WARNING - Missing Pi-zz data"
  if (n_elements(iuix) eq 0) then print," ****** WARNING - Missing Uix data"
  if (n_elements(iuiy) eq 0) then print," ****** WARNING - Missing Uiy data"
  if (n_elements(iuiz) eq 0) then print," ****** WARNING - Missing Uiz data"
  if (n_elements(iuex) eq 0) then print," ****** WARNING - Missing Uex data"
  if (n_elements(iuey) eq 0) then print," ****** WARNING - Missing Uey data"
  if (n_elements(iuez) eq 0) then print," ****** WARNING - Missing Uez data"
  if (n_elements(ijx) eq 0) then print," ****** WARNING - Missing Jx data"
  if (n_elements(ijy) eq 0) then print," ****** WARNING - Missing Jy data"
  if (n_elements(ijz) eq 0) then print," ****** WARNING - Missing Jz data"

; Define different types of reduced diagnostic

computedquantity = strarr(19)  
computedquantity(0) = "Electron Anisotropy"
computedquantity(1) = "Ion Anisotropy"
computedquantity(2) = "Electron Agyrotopy"
computedquantity(3) = "Ion Agyrotopy"
computedquantity(4) = "Firehose Condition"
computedquantity(5) = "Electron Firehose"
computedquantity(6) = "KH Condition"
computedquantity(7) = "|J|"
computedquantity(8) = "|B|"
computedquantity(9) = "Electron Vorticity"
computedquantity(10) = "Total Vorticity"
computedquantity(11) = "E*J"
computedquantity(12) = "beta-e"
computedquantity(13) = "beta-i"
computedquantity(14) = "beta-total"
computedquantity(15) = "Curvature"
computedquantity(16) = "Charge Density"
computedquantity(17) = "Te"
computedquantity(18) = "Ti"

; Close the input file

close,unit
free_lun, unit

; Define twci

; twci = toutput

; Define zcenter 

zcenter=zmax/2.0

; Echo information

print,'nx=',nx,'  ny=',ny,'  nz=',nz
print,'xmax=',xmax,'  ymax=',ymax,'  zmax=',zmax

; Plotting array

fulldata = fltarr(nx,nz)

; Define structure and Association data to use direct access

;struct = {data:fltarr(nx,ny,nz),time:0.0,it:500000}
;struct = {data:fltarr(nx,ny,nz)}
;field = assoc(1,struct)

; Determine number of time slices

;information=file_info(datafiles(0))
;record_length = 4L*(nx*ny*nz+2L)
;tslices=information.size/record_length-1
;print,"File Size=",information.size
;print,"Record Length=",record_length
;if (tslices lt 1) then tslices=1
;print,"Time Slices",tslices

; Plot type

ptype(0)='Contour'
ptype(1)='Contour+X-Average'
ptype(2)='Contour+Z-Slice'
ptype(3)='Contour+X-Slice'
ptype(4)='Contour+1-Trajectory'
ptype(5)='Contour+2-Trajectories'

; plane orientation

porient(0)='XZ'
porient(1)='XY'
porient(2)='ZY'

; axis names
anames(0,0) = 'Z'
anames(0,1) = 'X'

anames(1,0) = 'Y'
anames(1,1) = 'X'

anames(2,0) = 'Z'
anames(2,1) = 'Y'


; Options menu

desc =[ '1\Options' , $
        '1\Smoothing' , $
        '0\1' , $
        '0\2' , $
        '0\3' , $
        '0\4' , $
        '0\5' , $
        '0\6' , $
        '0\7' , $
        '2\8' , $
        '1\Contours' , $
        '0\12' , $
        '0\14' , $
        '0\16' , $
        '0\18' , $
        '0\20' , $
        '0\22' , $
        '0\24' , $
        '0\26' , $
        '0\28' , $
        '0\30' , $
        '0\32' , $
        '0\34' , $
        '0\36' , $
        '0\38' , $
        '0\40' , $
        '0\42' , $
        '0\44' , $
        '0\46' , $
        '2\48' , $
        '1\Transform Data' , $
        '0\Full Data' , $
        '0\Perturbed', $
        '2\',  $
        '1\Color Map', $
        '0\Default' , $
        '2\Load Table' ]

; Setup widgets

base = widget_base(/column)
row1 = widget_base(base,/row,scr_ysize=nypix/7)
row2 = widget_base(base,/row)
row3 = widget_base(base,/row)
button1 = widget_button(row1, value = '  Done  ', uvalue = 'done',/sensitive)
button2 = widget_button(row1, value = '  Plot  ', uvalue = 'plot',/sensitive)
button3 = widget_button(row1, value = 'Render PS', uvalue = 'render')
button4 = widget_button(row1, value = 'Animate', uvalue = 'animate')
list5=widget_list(row1,value=ptype,uvalue='ptype')
list4=widget_droplist(row1,title=' Raw Plot',value=rawdata,uvalue='rawplot')
list1=widget_droplist(row1,title=' Diagnostic Plot',value=computedquantity,uvalue='computed')
;list6=widget_droplist(row1,title=' Orientation',value=porient,uvalue='orientation')
opt = cw_pdmenu(row1,desc,/return_index,uvalue='options')
slider10=cw_fslider(row1,title='KH Angle',value=0.0,format='(f4.0)',uvalue='angle',minimum=0,maximum=180.0,/drag,scroll=1.0)
slider2=cw_fslider(row2,title='X-min',value=0.0,format='(f7.1)',uvalue='xmin',minimum=0,maximum=xmax,/drag)
slider3=cw_fslider(row2,title='X-max',value=xmax,format='(f7.1)',uvalue='xmax',minimum=0,maximum=xmax,/drag)
slider2b=cw_fslider(row2,title='Y-min',value=0.0,format='(f7.1)',uvalue='ymin',minimum=0,maximum=ymax,/drag)
slider3b=cw_fslider(row2,title='Y-max',value=ymax,format='(f7.1)',uvalue='ymax',minimum=0,maximum=ymax,/drag)
slider6=cw_fslider(row2,title='Z-max',value=zcenter,format='(f7.1)',uvalue='zmax',minimum=0.0,maximum=zcenter,/drag)
slider4=cw_fslider(row2,title='X-Slice',value=xmax/2,format='(f7.1)',uvalue='xslice',minimum=0,maximum=xmax,/drag,scroll=1.0)
slider5=cw_fslider(row2,title='Z-Slice',value=0,format='(f7.1)',uvalue='zslice',minimum=-zmax/2,maximum=zmax/2,/drag,scroll=1.0)
slider7=widget_slider(row2,title=' X-Shift',scrol=1,value=0,uvalue='shift',minimum=-nx/2,maximum=nx/2)

;slider8  =widget_slider(row3,title='x-plane',scrol=1,value=0,uvalue='xplane',minimum=0,maximum=nx-1,scr_xsize=nxpix/3.6)

; Set this slider to work with 2D or 3D
nymax = ny-1
if (ny eq 1) then nymax = 1 
slider9  =widget_slider(row3,title='y-plane',scrol=1,value=0,uvalue='yplane',minimum=0,maximum=nymax,scr_xsize=nxpix/3.6)
tsmax=tslices-1
if (tslices eq 1) then tsmax=1
slider1=widget_slider(row3,title='Time Slice',scrol=1,value=0,uvalue='time',minimum=0,maximum=tsmax,scr_xsize=nxpix/3.6)
;slider10 =widget_slider(row3,title='z-plane',scrol=1,value=0,uvalue='zplane',minimum=0,maximum=nz-1,scr_xsize=nxpix/3.6)

draw = widget_draw(base,retain=2, xsize = nxpix, ysize = nypix,/button_events,uvalue='mouse')
widget_control, base, /realize
widget_control,list5,set_list_select=0
widget_control,list1,set_droplist=0
widget_control, base, set_uvalue={quantity:1,time:0,xmin:0.0,xmax:xmax,ymin:0.0,ymax:ymax,xslice:xmax/2,zmin:0.0,zmax:zmax,zslice:zcenter,smoothing:2,contours:24,plottype:0,rawplot:0,shift:0,data:0,map:0,orientation:0,xcut:1,ycut:0,zcut:1,computed:0,angle:0.0}
widget_control, draw, get_value = index
xmanager, 'handle', base
halt_error1: print, ' *** Halting Program ***'
end


