import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

class Renderer:

    def __init__(self, u, f_plot, k_plot, params):
        self.fig = plt.figure(figsize=(8, 6))
        self.gs = self.fig.add_gridspec(2, 2, height_ratios=[1, 2], hspace=0.3, wspace=0.3)
        self.ax1 = self.fig.add_subplot(self.gs[0, 0], title ="f(x, y, t)")
        self.ax2 = self.fig.add_subplot(self.gs[0, 1], title ="k(x, y, t)")
        self.ax3 = self.fig.add_subplot(self.gs[1, :])
        self.forcing_plot = self.ax1.imshow(f_plot)
        self.k_plot = self.ax2.imshow(k_plot)
        plt.colorbar(self.forcing_plot, orientation='vertical')
        plt.colorbar(self.k_plot, orientation='vertical')
        self.forcing_plot.set_clim(-1, 1)
        self.img = self.ax3.imshow(
            u,
            extent=[*params],
            interpolation='None'
        )
        self.txt = self.ax3.text(11, 0, f"time = {0.245:.2f}", color='white',
               fontsize=12, ha='left', va='center',
               bbox=dict(facecolor='black', alpha=0.5)
        )
       # plt.gca().invert_xaxis()
        plt.colorbar(self.img, orientation='horizontal')

    def new_frame(self, u, T, F = None):
        self.img.set_data(u)
        self.img.set_clim(u.min(), u.max())
        self.txt.set_text(f"T {T:2f} s")
        if F is not None:
            self.forcing_plot.set_data(F)
        plt.draw()
        plt.pause(0.01)
