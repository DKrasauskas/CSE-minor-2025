import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np

class Renderer:

    def __init__(self):
        self.fig = plt.figure(figsize=(8, 6))
        #self.gs = self.fig.add_gridspec(1, 2, height_ratios=[2, 1], hspace=0.3, wspace=0.3)
        self.gs = self.fig.add_gridspec(1, 2, hspace=0.3, wspace=0.3)
        self.ax1 = self.fig.add_subplot(self.gs[0, 0], title ="u(x, y)")
        self.ax2 = self.fig.add_subplot(self.gs[0, 1], title ="v(x, y)")
        #self.ax3 = self.fig.add_subplot(self.gs[1, :], title="u(x, y, t)")
        self.dynamic_data1 = None
        self.dynamic_data2 = None
        self.img1 = None
        self.img2 = None
        self.T = 0
        self.dT = 0

    def render(self):
        if self.dynamic_data1 is not None:
            if self.img1 is None:
                self.img1 = self.ax1.imshow(self.dynamic_data1)
                plt.colorbar(self.img1, orientation='horizontal')
                self.txt = self.ax1.text(0, -20, f"time = {0.245:.2f}", color='white',
                                         fontsize=12, ha='left', va='center',
                                         bbox=dict(facecolor='black', alpha=0.5)
                 )
            else:
                # Flip the data along y-axis
                flipped_data = np.flipud(self.dynamic_data1)

                # Display the data without rescaling, but set axes 0→10
                self.img1.set_data(flipped_data)
                self.img1.set_clim(flipped_data.min(), flipped_data.max())

                # Set the extent so that axes go from 0 to 10
                # Format: [x_min, x_max, y_min, y_max]
                self.img1.set_extent([0, 4, 0, 4])

                # Update text
                self.txt.set_text(f"T {self.T:2f} s")

                # Update your text

        if self.dynamic_data2 is not None:
            if self.img2 is None:
                self.img2 = self.ax2.imshow(self.dynamic_data2)
                plt.colorbar(self.img2, orientation='horizontal')
                self.txt2 = self.ax2.text(0, -20, f"Δt = {0.245:.2f}", color='white',
                                         fontsize=12, ha='left', va='center',
                                         bbox=dict(facecolor='black', alpha=0.5)
                 )
            else:
                # Flip the data along y-axis
                flipped_data = np.flipud(self.dynamic_data2)

                # Display the data without rescaling, but set axes 0→10
                self.img2.set_data(flipped_data)
                self.img2.set_clim(flipped_data.min(), flipped_data.max())

                # Set the extent so that axes go from 0 to 10
                # Format: [x_min, x_max, y_min, y_max]
                self.img2.set_extent([0, 4, 0, 4])

                # Update text
                self.txt2.set_text(f"Δt {self.dT:2f} s")


class DebugRenderer:

    def __init__(self):
        self.fig = plt.figure(figsize=(8, 6))
        self.gs = self.fig.add_gridspec(2, 2, height_ratios=[2, 1], hspace=0.3, wspace=0.3)
        self.ax1 = self.fig.add_subplot(self.gs[0, 0], title ="u(x, y)")
        self.ax2 = self.fig.add_subplot(self.gs[0, 1], title ="v(x, y)")
        self.ax3 = self.fig.add_subplot(self.gs[1, :], title="NR iteration count / update")
        self.dynamic_data1 = None
        self.dynamic_data2 = None
        self.dynamic_data3 = None
        self.img1 = None
        self.img2 = None
        self.img3 = None
        self.it = None
        self.T = 0
        self.dT = 0

    def render(self):
        if self.dynamic_data1 is not None:
            if self.img1 is None:
                self.img1 = self.ax1.imshow(self.dynamic_data1)
                plt.colorbar(self.img1, orientation='horizontal')
                self.txt = self.ax1.text(0, -20, f"time = {0.245:.2f}", color='white',
                                         fontsize=12, ha='left', va='center',
                                         bbox=dict(facecolor='black', alpha=0.5)
                 )
            else:
                flipped_data = np.flipud(self.dynamic_data2)

                # Display the data without rescaling, but set axes 0→10
                self.img1.set_data(flipped_data)
                self.img1.set_clim(flipped_data.min(), flipped_data.max())

                # Set the extent so that axes go from 0 to 10
                # Format: [x_min, x_max, y_min, y_max]
                self.img1.set_extent([0, 4, 0, 4])

                # Update text
                self.txt.set_text(f"T {self.T:2f} s")

        if self.dynamic_data2 is not None:
            if self.img2 is None:
                self.img2 = self.ax2.imshow(self.dynamic_data2)
                plt.colorbar(self.img2, orientation='horizontal')
                self.txt2 = self.ax2.text(0, -20, f"time = {0.245:.2f}", color='white',
                                         fontsize=12, ha='left', va='center',
                                         bbox=dict(facecolor='black', alpha=0.5)
                 )
            else:
                flipped_data = np.flipud(self.dynamic_data2)

                # Display the data without rescaling, but set axes 0→10
                self.img2.set_data(flipped_data)
                self.img2.set_clim(flipped_data.min(), flipped_data.max())

                # Set the extent so that axes go from 0 to 10
                # Format: [x_min, x_max, y_min, y_max]
                self.img2.set_extent([0, 4, 0, 4])

                # Update text
                self.txt2.set_text(f"Δt {self.dT:2f} s")
        if self.dynamic_data3 is not None:
            self.img3 = self.ax3.plot(self.it, self.dynamic_data3)
