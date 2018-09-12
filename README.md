# Sturm sequences

A way to find real roots of polynomial equations.

## Example of usage:

```
Please enter order of polynomial: 6

Please enter coefficient number 6: -1
Please enter coefficient number 5: 1
Please enter coefficient number 4: 2
Please enter coefficient number 3: -2
Please enter coefficient number 2: 1
Please enter coefficient number 1: 2
Please enter coefficient number 0: -1

Sturm sequence for:
-1.000000 1.000000 2.000000 -2.000000 1.000000 2.000000 -1.000000 

-1.000000 1.000000 2.000000 -2.000000 1.000000 2.000000 -1.000000 
-1.000000 0.833333 1.333333 -1.000000 0.333333 0.333333 
-1.000000 0.965517 -0.620690 -2.137931 1.172414 
-1.000000 -0.667969 0.304687 -0.097656 
1.000000 0.764977 -0.502304 
1.000000 0.180000 
0.607600 

4 distinct real roots for x: -1.246980 -1.000000 0.445042 1.801938 
```

## References
* [C code](https://webdocs.cs.ualberta.ca/~graphics/books/GraphicsGems/gems/Sturm/) from Graphic Gems by D.G. Hook and P. R. McAree.
* [Roots of Polynomials](polyroots.pdf) - Yan-Bin Jia 