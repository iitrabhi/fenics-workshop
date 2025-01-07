## Install Salome
Install salome from <https://www.salome-platform.org/>
Go to downloads and select your OS.

## Open SALOME
After you install and open salome, the first interface should look something like this:
![Salome](../salome_tutorial/image.png)

* Click on the shaper
![shaper](../salome_tutorial/image-1.png)

* You should see a window like this:
![alt text](../salome_tutorial/image-2.png)

## Now we will try to draw a rectangle with a hole in it

* Click on the "Sketch" tab or selection sketch from the "Sketch" in the upper toolbar
![alt text](../salome_tutorial/image-3.png)

* After clicking it you will be asked to choose a plane. Select ZX plane, which should be green in colour. There are 3 colours Red Green and Blue. Lets select Green now. The orientation might look different.

![alt text](../salome_tutorial/image-4.png)

* After selecting the Green rectangle, also click "Set plane view"

![alt text](../salome_tutorial/image-5.png)

* Now you will see the ZX plane as shown:

![alt text](../salome_tutorial/image-6.png)

* Here we will draw a rectangle using "Rectangle" tool as shown:

![alt text](../salome_tutorial/image-7.png)

* After clicking on the rectangle, take your mouse closer to the center, a blue dot will be highlighted, i.e. your mouse will now snap to the center. Click here.
![alt text](../salome_tutorial/image-8.png)

* After clicking move your mouse any where in the first quadrant and click for the second time. You will see something like this:

![alt text](../salome_tutorial/image-9.png)

* Now we will constraint the dimensions parametrically. Then click length in the top panel as shown

![alt text](../salome_tutorial/image-11.png)

* Now click on any horizontal line

![alt text](../salome_tutorial/image-12.png)

* Now drag a little up and click. After you click the text box will be shown. in the text box type "L=50"

![alt text](../salome_tutorial/image-13.png)

* Press enter. Now notice two things one the dimension has been set, and next a new parameter "L" has been created whose value can be changed parametrically and whose change will be reflected.

![alt text](../salome_tutorial/image-14.png)

* Now click on the green tick arrow for confirmation

![alt text](../salome_tutorial/image-15.png)

* Repeat the process for vertical lines, in the text box that appears put "H=40". Press "Enter" and click the green arrow for confirmation

![alt text](../salome_tutorial/image-16.png)

* We have our rectangle ready. Now lets draw a circle. Click on the circle icon as shown:

![alt text](../salome_tutorial/image-17.png)

* After that click on any point inside the rectangle and then mouse your mouse around for specifying for radius. Right now just keep it whatever you want. Get something line this:

![alt text](../salome_tutorial/image-18.png)

You can see that the circle is red. Red in salome is unconstrained. Green is constrained. Now we will constrain our circle. First we will specify the radius parametrically.

* Click on the radius as shown:

![alt text](../salome_tutorial/image-19.png)

* After clicking on circle you will see something like:
![alt text](../salome_tutorial/image-20.png)

* Now click again so that the text box appears: In the box put r = 10

![alt text](../salome_tutorial/image-21.png)

Press Enter and tick the green arrow like before.

* Great work. Now we have circle and a rectangle whose dimensions can be parametrically changed.

* The final thing is to contrain the circle in the centre of rectangle to do that draw a diagonal as follows:

Click the line:
![alt text](../salome_tutorial/image-22.png)

Then draw a diagonal by clicking two ends of the rectangle as:

![alt text](../salome_tutorial/image-23.png)

* Now we will use mid point contraint to fix the center of circle to the center of diagonal.

* To do that first click on the "Middle point" constraint:

![alt text](../salome_tutorial/image-24.png)

* You will see a side panel like:
![alt text](../salome_tutorial/image-25.png)

* For the first object, now click the center of circle. After clicking you will see:

![alt text](../salome_tutorial/image-26.png)

* Then select the diagonal

![alt text](../salome_tutorial/image-27.png)

* Doing this will immediately move the circle to the center and everything will be green as above. Then click the tick mark.

![alt text](../salome_tutorial/image-28.png)

* Now we have drawn the sketch completely. A sketch consists of lines only. We make use of this lines to draw further shapes.

To mark completion click green tick as shown:

![alt text](../salome_tutorial/image-29.png)

* After click tick everything will be greyed out and look like this:

![alt text](../salome_tutorial/image-30.png)

* Now try chaning the parameters and see if there is any effect:

* Double click on L as shown:
![alt text](../salome_tutorial/image-31.png)

* Change the 50 to 60 as shown:
![alt text](../salome_tutorial/image-32.png)

* Click the green tick:
![alt text](../salome_tutorial/image-33.png)

* You will see that the figure stretches:

![alt text](../salome_tutorial/image-34.png)

* Now lets add faces. To add faces click "Face":
![alt text](../salome_tutorial/image-35.png)

* Just hover around the upper half to see:
![alt text](../salome_tutorial/image-36.png)

* Click on the upper right corner area and then "SHIFT+ CLICK" on the bottom left to get:

![alt text](../salome_tutorial/image-37.png)

* Then click the green tick

![alt text](../salome_tutorial/image-38.png)

* You will see something like:

![alt text](../salome_tutorial/image-39.png)

* Finally we will merge those two halves into one for meshing. Click on "Fuse"

![alt text](../salome_tutorial/image-40.png)

* Click on Face_1_1 first, then shift click on Face_1_2

![alt text](../salome_tutorial/image-41.png)

* Then Click on the green tick

![alt text](../salome_tutorial/image-42.png)

* You will se a fuse item here:

![alt text](../salome_tutorial/image-43.png)
