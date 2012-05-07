APNG Assembler 2.0

by Max Stepin 
maxst@users.sourceforge.net
---------------------------
Creates APNG animation from PNG image sequence.

Usage: 
apngasm.exe output.png frame001.png [1] [10] [/f]

1/10 sec is the default delay. Use /f to skip the first frame.

----------------------------------------------------------------
Example 1:

Let's say you have following frame sequence:
frame01.png
frame02.png
frame03.png

And you want to have 3/4 seconds delay between frames.
The correct command will be

apngasm output.png frame01.png 3 4

----------------------------------------------------------------
Example 2: 

The same as above, but you added "invisible" frame00.png :
frame00.png - invisible
frame01.png
frame02.png
frame03.png

The correct command will be

apngasm output.png frame00.png 3 4 /f

That way APNG supported browsers and image viewers 
will show frame01-frame02-frame03 animation, 
while IE and Chrome will display static frame00.png image.

Some people like to put "your browser sucks" message in that 
frame00.png but keep in mind that you will see that message 
in thumbnails when you are browsing image folders, 
or using google image search. So be careful.

----------------------------------------------------------------
Example 3: 

apngasm output.png frame01.png

That way you'll get 1/10 sec delay.

----------------------------------------------------------------
Linux binaries compiled under Puppy, like that:

gcc -O2 -o apngasm apngasm.c -lz -lpng 

