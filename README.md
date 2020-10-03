# interpol

General interpolation inside a polygon, and color maps with >1 input value. Generalisation of linear colormap, which has (at most) 2 extremes, to triangular colormap, which has 3 extremes.

## Interpolation inside a shape: `triangles` and `polygon`
There are 2 functions to interpolate inside a shape.

* `triangles` tesselates the plane with triangles based on the provided anchorpoints, and does 'standard' barycentric interpolation within each triangle. Disadvantage: a maximum of 3 anchors is used for any point, which does not always look good. Also, extrapolation, i.e. to points that do not lay within the convex hull around the anchorpoints, is not possible. If wanted, the function tries nonetheless, but results are often poor.
* `polygon` is a more general function, that gives better results. The only drawback is that the anchorpoints must be given in order, as a chain describing the boundary of the polygon - there may be no internal points. 

Here a comparison between the two. The first column shows how the first function interpolates; note how the tesselation triangles are quite visible. The second column shows how the second function interpolates, which is much better. (The dots indicate the anchor points.)
![image with squares](in_shape_square.png)
![image with circles](in_shape_circle.png)

## Background: colormaps
The interpolation functions above can be used for colormaps with 3 axes. For a better understanding, we go into colormaps with 1 and 2 axes first, even though these don't use the interpolation functions.

### With 1 axis
A common colormap turns any value in the interval [0, 1] into its corresponding color, for example the 'terrain' colormap in matplotlib: 

![cmap](cmap_notdiverging.png)

Here the values for which a color has been explicitly specified (the anchor points), are marked with dots. The color for any other value is interpolated from the nearest on either side).

If the values that we want to turn into colors don't live on the interval [0, 1], some mapping of the input domain onto [0, 1] is necessary, but this is left out of consideration here.

### With 2 axes, showing relative size: `cmap2`
If we have a 'diverging' colormap, like the red-green one below, we can see this as a special case of the 1-axis colormap we just saw. But we can also imagine that there are in fact 2 axes (A on the left and B on the right), that meet in the middle. If values for *a* and *b*, both in the interval [0, 1], are provided, we are able to pick a single color on the bar. We do this by taking difference d = b - a. If b is the larger value, we have d > 0, and the color that belongs to this point (a, b) is the color that we find on axis B at the value d. If a is the larger value, we have d < 0, and the corresponding color is found on axis A at the value -d. 

An example would be when a is my expenditure in a given time period, and b is my income in the same time period. The resulting diverging colormap shows the 'income balance':

![cmap2 for income and expenditure](cmap2_income.png)

(Anchor points are again indicated; here, 2 colors were used on axis A, 3 on axis B, and also the center color was specified.)

It might seem a bit silly to turn what is basically a single variable (*d*) into two variables (*a* and *b*), especially since our color map is really only 1-dimensional. However, this will help us generalise onto three axes, further below. 

It's important that we see that this color map actually displays the *relative size* of *a* and *b*:
* It shows us, *which is the largest* of the 2 values, from the principle color (red or green) - if *a* is larger, the color comes from axis A. 
* It shows us, *how much larger* the largest of the two values is - the further the color is towards the end of the axis, the larger their difference.

Note that this is only possible if the values can be meaningfully subtracted from each other. 

Another example might be the time gained by switching between driving and biking to a certain location. Either biking (left) or driving (right) is faster, and the larger the time difference between the two transport options, the further away from the center the color is picked:

![cmap2 for transport modes](cmap2_transport.png)

Using this final example, we see that one of the values (*a*, *b*) is always 0. That's because we are showing the 'time gained'. The slowest of the 2 options becomes the reference for the other one, i.e., (*a*, *b*) = (something, 0) if the car is the slowest option, which translates to a color from the left side of the color map. Vice versa, when biking is slower, its 'time gained' is 0, so (*a*, *b*) = (0, something). 

Note that this means that the color that corresponds to (bike takes 4 minutes, car takes 7 minutes) is the same as the one that corresponds to (bike takes 5 minutes, car takes 8 minutes), which underlines that this colormap can be used to show *relative size*: the gain by going by bike is 3 minutes in both cases.

The color mapping function to find the color of a point (*a*, *b*), is returned by `cmap2`. The arguments to this function are the colors on either axis, as well as the one in the middle. It's also possible to include arguments that handle the mapping of the input domain, in (the likely) case the extremes on the axes do *not* correspond to a difference between *a* and *b* of exactly 1.


## Colormap with 3 axes: `cmap3`
Now, what if we have another transport mode in that last example? What if we also have the option of taking public transport, and want to be able to pick a color, based on the differences in time needed to get somewhere by bike, car, or public transport?

Having 3 values means 3 axes: A, B, C. These still meet in a point, but at 120 degree angles:

![3 axes, partial](cmap3_transport_partial.png)

This image, however, only gives us the color when of the 3 values (*a*, *b*, *c*), *two* are zero. For example, when biking has an equal time gain over car and bus, we're in the orange 'arm' of the image using the point (*a*, *b*, *c*) = (something, 0, 0).

In order to also find the colors when only *one* of the coordinates is 0 - remember that at least one of the coordinates will always be 0 - we need to fill in more colors, which is where the interpolation function shown above comes in:

![3 axes, full](cmap3_transport_full.png)

(Here, 2 colors are specified for each of the axes, and none is specified for the center.)

The location on the triangle that a color is selected from, indicates the relationship between the values *a*, *b* and *c*. In this image, we focus on *a*:

![3 axes, explainer](cmap3_transport_special.png)

The color mapping function to find the color of a point (*a*, *b*, *c*), is returned by `cmap3`. The arguments to this function are the colors on each axis, as well as the one in the middle. Here, too, it's possible to include arguments that handle the mapping of the input domain.
