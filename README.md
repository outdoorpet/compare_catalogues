# compare_catalogues
GUI (written with PyQt) for comparing event catalogues from different sources (via obspy FDSN client).

![Imgur](http://i.imgur.com/JrqtrW2.png)
![Imgur](http://i.imgur.com/jAo0ZNL.png)

### Dependencies

**The GUI has the following dependencies:**

* `Python 2.7`
* `PyQt4`
* `qdarkstyle`
* `pyqtgraph`
* `obspy`
* `pandas`
* `numpy`

Currently the User must edit the code in lines 209-213 to specify the 
desired earthquake catalogues to compare and the start and end date/times for the comparison.

This will be built into the GUI in the future.
