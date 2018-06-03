from time import sleep
from random import random
import numpy as np
import matplotlib.pyplot as plt

class DynamicPlot(object):
    def __init__(self, fid_name):
        self.fig = plt.figure(fid_name)
        plt.clf()
        self.ax = self.fig.add_subplot(111)
        self.line1, = self.ax.plot([1], [0], 'ko')
        self.ax.axis((0, 2, 0, 1))
        self.fig.canvas.draw()
        plt.show()

        # def _update_callback(sender, **kwargs):
        #     self.line1.set_ydata(kwargs['y'])
        #     self.fig.canvas.draw()

    def update_callback(self, y):
        self.line1.set_ydata(y)
        self.line1.set_xdata(2*y)
        self.fig.canvas.draw()


        # plot_update.connect(_update_callback, weak=False)
        # plot_update.connect(tmp(self))


if __name__ == "__main__":

    d1 = DynamicPlot("Not it Works!")

    for i in range(10):
        d1.update_callback(i*np.arange(10))
