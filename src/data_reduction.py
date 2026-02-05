import numpy as np 
from src.options import * 
from src.frames import * 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm
from matplotlib.widgets import Button
import time
from astropy.coordinates import SkyCoord, EarthLocation
from scipy.optimize import curve_fit
from scipy.constants import speed_of_light

grey_cmap = copy.copy(matplotlib.cm.get_cmap('Greys_r'))
grey_cmap.set_bad((0,0,0))	

def gaussian(x, a, mean, std_dev, h):
    return a / (std_dev * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mean) / std_dev) ** 2) + h

def markov_gaussian(x, amp, mean, std):
    return amp * np.exp(-(x - mean) ** 2 / (2 * std ** 2))

def gaussian_with_tilt(x, a, mean, std_dev, h, b):
    return a / (std_dev * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mean) / std_dev) ** 2) + h + x * b


def lo(x, gamma, x_0):
    return 1 / np.pi * (gamma / 2) / ((x - x_0) ** 2 + (gamma / 2) ** 2)


log_two = np.log(2)


def ga(x, gamma, x_0):
    sigma = gamma / (2 * np.sqrt(2 * log_two))
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp((-(x - x_0) ** 2) / (2 * sigma ** 2))


def pseudo_voigt(x, scaling, gamma, shift, slope, height, eta):
    g = ga(x, gamma, shift)
    l = lo(x, gamma, shift)
    return -scaling * (eta * g + (1 - eta) * l) + slope * x + height


def line(x, w, v, u, m, n):
    return w * x ** 4 + v * x ** 3 + u * x ** 2 + x * m + n

def lowpoly(x, u, m, n):
    return u * x ** 2 + x * m + n


def wlshift(wl, vel_corr):
    # wl_shift = vel_corr/speed_of_light * wl
    # return wl+wl_shift
    return wl / (1 - (vel_corr / (speed_of_light / 1000)))


def get_central_wavelength(gra_angle, cam_angle, d_grating):
    return (np.sin(gra_angle * 2 * np.pi / 360) + np.sin((cam_angle - gra_angle) * 2 * np.pi / 360)) / (
            d_grating * 1.e-7)


def get_flux(image, x_ind, y_ind, width, boxcut=True):

    col = int(round(x_ind))
    ysize = image.shape[0]

    if boxcut:
        y_start = int(np.ceil(y_ind - width))
        y_end = int(np.floor(y_ind + width))
        fluxsum = np.sum(image[y_start:y_end+1, col])
        
        # fractional pixel contribution from edges
        if y_start > 0:
            upperfraction = image[y_start - 1, col] * (
                np.abs(y_start - (y_ind - width)))
        else:
            upperfraction = 0
        
        if y_end + 1 < ysize:
            lowerfraction = image[y_end + 1, col] * (
                np.abs((y_ind + width) - y_end))
        else:
            lowerfraction = 0
        
        fluxsum += upperfraction + lowerfraction
        return fluxsum

    else:
        sigma = width / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to sigma
        y_min = int(np.floor(y_ind - width))
        y_max = int(np.ceil(y_ind + width))

        # Clip bounds to image limits
        y_min = max(y_min, 0)
        y_max = min(y_max, ysize - 1)

        y = np.arange(y_min, y_max + 1)
        weights = np.exp(-0.5 * ((y - y_ind) / sigma)**2)

        fluxsum = np.sum(image[y, col] * weights)
        return fluxsum


def fluxstatistics(wl, flux):
    med = median_filter(flux, 5)
    flux_norm = flux / med - 1
    std = pd.Series(flux_norm).rolling(min_periods=1, window=20, center=True).std().to_numpy()

    flux = flux[flux_norm < 3 * std]
    wl = wl[flux_norm < 3 * std]

    med = median_filter(flux, 5)
    flux_norm = flux / med - 1
    std = pd.Series(flux_norm).rolling(min_periods=1, window=20, center=True).std().to_numpy()

    flx_std = flux * std

    return wl, flux, flx_std


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import median_filter
import matplotlib.colors as colors

class InteractiveTraceSelector:
    def __init__(self, image):
        self.image = image
        self.clicks = []
        self.lines = []
        self.bounds = None
        self.reset_pressed = False
        self.finished = False
        
    def setup_plot(self):
        self.fig, self.ax = plt.subplots(1, 1, figsize=(4.8 * 16 / 9, 4.8))
        
        # Display the image
        norm_l = np.percentile(self.image, 1)
        norm_h = np.percentile(self.image, 99.9)
        if norm_l == 0:
            norm_l = 1
            
        self.ax.imshow(self.image, cmap=grey_cmap, norm=colors.LogNorm(norm_l, norm_h), aspect="auto")
        self.ax.set_title("Click twice to select Y-bounds for trace fitting\nPress 'R' to reset, Enter or close to proceed")
        
        # Connect event handlers
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.fig.canvas.mpl_connect('close_event', self.on_close)
        
        plt.tight_layout()
        
    def on_click(self, event):
        if event.inaxes != self.ax or len(self.clicks) >= 2:
            return
            
        y_pos = event.ydata
        if y_pos is None:
            return
            
        self.clicks.append(y_pos)
        line = self.ax.axhline(y=y_pos, color='red', linestyle='--', linewidth=2)
        self.lines.append(line)
        
        if len(self.clicks) == 2:
            self.ax.set_title("Bounds selected. Press Enter or close to proceed, 'R' to reset")
        else:
            self.ax.set_title(f"Click #{len(self.clicks)+1} to select Y-bounds for trace fitting\nPress 'R' to reset, Enter or close to proceed")
            
        self.fig.canvas.draw()
        
    def on_key(self, event):
        if event.key == 'r' or event.key == 'R':
            # Reset clicks and lines
            self.clicks = []
            for line in self.lines:
                line.remove()
            self.lines = []
            self.ax.set_title("Click twice to select Y-bounds for trace fitting\nPress 'R' to reset, Enter or close to proceed")
            self.fig.canvas.draw()
        elif event.key == 'enter':
            self.finish_selection()
            
    def on_close(self, event):
        self.finish_selection()
        
    def finish_selection(self):
        if len(self.clicks) >= 2:
            self.bounds = (min(self.clicks), max(self.clicks))
        else:
            # Use default bounds if insufficient clicks
            self.bounds = (self.image.shape[0] * 1/4, self.image.shape[0] * 3/4)
        self.finished = True
        plt.close(self.fig)
        
    def get_bounds(self):
        self.setup_plot()
        plt.show()
        
        # Wait for user interaction to complete
        while not self.finished:
            plt.pause(0.1)
            
        return self.bounds

class TraceValidator:
    def __init__(self, image, xcenters, ycenters, params, width, reduction_options):
        self.image = image
        self.xcenters = xcenters
        self.ycenters = ycenters
        self.params = params
        self.width = width
        self.reduction_options = reduction_options
        self.accepted = False
        self.finished = False
        
    def setup_plot(self):
        self.fig, self.axs = plt.subplots(2, 1, figsize=(4.8 * 16 / 9, 4.8))
        
        # Plot trace fit
        xspace = np.linspace(0, self.image.shape[1], 1000)
        self.axs[0].plot(xspace, lowpoly(xspace, *self.params), zorder=1)
        self.axs[0].scatter(self.xcenters, self.ycenters, color="red", marker="x", zorder=5)
        self.axs[0].set_title("Fitted trace - Accept or Reject?")
        
        # Plot image with trace
        norm_l = np.percentile(self.image, 1)
        norm_h = np.percentile(self.image, 99.9)
        if norm_l == 0:
            norm_l = 1
            
        self.axs[1].imshow(self.image, cmap=grey_cmap, norm=colors.LogNorm(norm_l, norm_h), aspect="auto")
        self.axs[1].plot(xspace, lowpoly(xspace, *self.params), color="lime", linewidth=0.5)
        self.axs[1].plot(xspace, lowpoly(xspace, *self.params) - self.reduction_options.skyfluxsep, color="red", linestyle="--", linewidth=0.5)
        self.axs[1].plot(xspace, lowpoly(xspace, *self.params) + self.reduction_options.skyfluxsep, color="red", linestyle="--", linewidth=0.5)
        self.axs[1].plot(xspace, lowpoly(xspace, *self.params) + self.width, color="lime", linestyle="--", linewidth=0.5)
        self.axs[1].plot(xspace, lowpoly(xspace, *self.params) - self.width, color="lime", linestyle="--", linewidth=0.5)
        
        # Add title with instructions
        self.axs[0].set_title("Fitted trace - Press 'Y' to Accept, 'N' to Reject, or close window")
        
        # Connect event handlers
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.fig.canvas.mpl_connect('close_event', self.on_close)
        
        plt.tight_layout()
        
    def on_key(self, event):
        if event.key in ['y', 'Y']:
            self.accepted = True
            self.finished = True
            plt.close(self.fig)
        elif event.key in ['n', 'N']:
            self.accepted = False
            self.finished = True
            plt.close(self.fig)
            
    def on_close(self, event):
        # Default to rejection if window is closed without decision
        if not self.finished:
            self.accepted = True
            self.finished = True
        
    def validate(self):
        self.setup_plot()
        plt.show()
        
        # Wait for user decision
        while not self.finished:
            plt.pause(0.1)
            
        return self.accepted

class ManualTraceSelector:
    """Interactive manual point selection for low quality data trace fitting"""
    
    def __init__(self, image, reduction_options):
        # Disable default matplotlib keybindings that conflict
        plt.rcParams['keymap.home'].remove('r') if 'r' in plt.rcParams['keymap.home'] else None
        plt.rcParams['keymap.fullscreen'].remove('f') if 'f' in plt.rcParams['keymap.fullscreen'] else None

        self.image = image
        self.reduction_options = reduction_options
        self.points = []  # List of (x, y) tuples
        self.widths = []  # Corresponding widths at each point
        self.scatter_points = None
        self.trace_line = None
        self.upper_width_line = None
        self.lower_width_line = None
        self.current_params = None
        self.current_width = None
        self.finished = False
        self.accepted = False
        self.cursor_x = None
        self.cursor_y = None
        
        # Set up the figure
        self.fig, self.ax = plt.subplots(figsize=(14, 8))
        
        # Display image
        norm_l = np.percentile(image, 1)
        norm_h = np.percentile(image, 99.9)
        if norm_l <= 0:
            norm_l = 0.1
        
        self.ax.imshow(image, cmap='gray', norm=colors.LogNorm(norm_l, norm_h), 
                       aspect='auto', origin='lower')
        self.ax.set_xlabel('X (pixels)')
        self.ax.set_ylabel('Y (pixels)')
        self.ax.set_title('A: Add point | R: Remove nearest | F: Force add | Enter: Accept | Esc: Cancel')
        
        # Initialize empty scatter plot for points
        self.scatter_points = self.ax.scatter([], [], c='red', s=100, marker='x', 
                                               zorder=10, linewidths=2)
        
        # Connect events
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        
        # Add instruction text
        self.instruction_text = self.ax.text(
            0.02, 0.98, f'Points: 0 | Need at least 3 points for trace fit',
            transform=self.ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        )
        
        plt.tight_layout()
        plt.show(block=True)
    
    def on_click(self, event):
        """Handle mouse click events"""
        if event.inaxes != self.ax:
            return
        
        if event.button == 1:  # Left click - add point
            self.add_point(event.xdata, event.ydata)
        elif event.button == 3:  # Right click - remove nearest point
            self.remove_nearest_point(event.xdata, event.ydata)

    def on_motion(self, event):
        """Track cursor position"""
        if event.inaxes == self.ax:
            self.cursor_x = event.xdata
            self.cursor_y = event.ydata
    
    def on_key(self, event):
        """Handle keyboard events"""
        if event.key == 'enter':
            if len(self.points) >= 3 and self.current_params is not None:
                self.accepted = True
                self.finished = True
                plt.close(self.fig)
            else:
                print("Need at least 3 points with a valid trace fit to accept!")
        elif event.key == 'escape':
            self.accepted = False
            self.finished = True
            plt.close(self.fig)
        elif event.key == 'a':
            if self.cursor_x is not None and self.cursor_y is not None:
                self.add_point(self.cursor_x, self.cursor_y, force=False)
        elif event.key == 'f':
            if self.cursor_x is not None and self.cursor_y is not None:
                self.add_point(self.cursor_x, self.cursor_y, force=True)
        elif event.key == 'r':
            if self.cursor_x is not None and self.cursor_y is not None:
                self.remove_nearest_point(self.cursor_x, self.cursor_y)
    
    def add_point(self, x, y, force=False):
        """Add a new point and fit width at that location"""
        x = int(round(x))
        y = float(y)
        
        # Fit gaussian to get width at this x position
        width = self.fit_width_at_point(x, y)
        
        if width is not None:
            self.points.append((x, y))
            self.widths.append(width)
            print(f"Added point at ({x}, {y:.1f}) with width {width:.2f}")
            self.update_display()
        elif force:
            # Use median of existing widths, or default if none exist
            if self.widths:
                fallback_width = np.nanmedian(self.widths)
            else:
                fallback_width = np.nan  # Default fallback: NaN
            self.points.append((x, y))
            self.widths.append(fallback_width)
            print(f"Force-added point at ({x}, {y:.1f}) with fallback width {fallback_width:.2f}")
            self.update_display()
        else:
            print(f"Could not fit width at ({x}, {y:.1f}) - use 'F' to force add")
    
    def fit_width_at_point(self, x, y, half_window=15):
        """Fit a gaussian to vertical slice at x position with center constrained near y"""
        # Extract vertical slice
        x_lo = max(0, x - half_window)
        x_hi = min(self.image.shape[1], x + half_window)
        
        if x_hi - x_lo < 5:
            return None
        
        data = np.median(self.image[:, x_lo:x_hi], axis=1)
        data = median_filter(data, 5)
        xarr = np.arange(len(data))
        
        if np.all(data < 5):
            return None
        
        try:
            # Constrain center to be within ±1 pixel of clicked y position
            params, *_ = curve_fit(
                gaussian,
                xarr,
                data,
                p0=[np.max(data) - np.median(data), y, 10, np.median(data)],
                bounds=[
                    [0, y - 1, 1, -np.inf],      # Lower bounds: center at y-1
                    [np.inf, y + 1, 50, np.inf]  # Upper bounds: center at y+1
                ],
                maxfev=10000
            )
            return params[2]  # Return width (sigma)
        except (ValueError, RuntimeError) as e:
            print(f"Width fit failed at x={x}: {e}")
            return None
    
    def remove_nearest_point(self, x, y):
        """Remove the point nearest to the click location"""
        if not self.points:
            return
        
        # Find nearest point
        distances = [np.sqrt((px - x)**2 + (py - y)**2) for px, py in self.points]
        nearest_idx = np.argmin(distances)
        
        # Only remove if reasonably close (within 50 pixels)
        if distances[nearest_idx] < 50:
            removed = self.points.pop(nearest_idx)
            self.widths.pop(nearest_idx)
            print(f"Removed point at ({removed[0]}, {removed[1]:.1f})")
            self.update_display()
        else:
            print("No point close enough to remove")
    
    def update_display(self):
        """Update the display with current points and fitted trace"""
        # Update scatter points
        if self.points:
            xs, ys = zip(*self.points)
            self.scatter_points.set_offsets(np.column_stack([xs, ys]))
        else:
            self.scatter_points.set_offsets(np.empty((0, 2)))
        
        # Remove old trace lines
        if self.trace_line is not None:
            self.trace_line.remove()
            self.trace_line = None
        if self.upper_width_line is not None:
            self.upper_width_line.remove()
            self.upper_width_line = None
        if self.lower_width_line is not None:
            self.lower_width_line.remove()
            self.lower_width_line = None
        
        # Fit and display trace if we have enough points
        if len(self.points) >= 3:
            self.fit_trace()
            
            if self.current_params is not None:
                xspace = np.linspace(0, self.image.shape[1], 500)
                trace_y = lowpoly(xspace, *self.current_params)
                
                self.trace_line, = self.ax.plot(xspace, trace_y, 'lime', 
                                                 linewidth=1.5, label='Trace')
                self.upper_width_line, = self.ax.plot(xspace, trace_y + self.current_width, 
                                                       'lime', linewidth=1, linestyle='--',
                                                       label=f'Width (±{self.current_width:.1f})')
                self.lower_width_line, = self.ax.plot(xspace, trace_y - self.current_width, 
                                                       'lime', linewidth=1, linestyle='--')
        
        # Update instruction text
        status = f'Points: {len(self.points)}'
        if len(self.points) < 3:
            status += ' | Need at least 3 points for trace fit'
        elif self.current_params is not None:
            status += f' | Width: {self.current_width:.1f}px | Press Enter to accept'
        else:
            status += ' | Trace fit failed - add more points'
        
        self.instruction_text.set_text(status)
        
        self.fig.canvas.draw_idle()
    
    def fit_trace(self):
        """Fit polynomial trace to current points"""
        if len(self.points) < 3:
            self.current_params = None
            self.current_width = None
            return
        
        xs, ys = zip(*self.points)
        xs = np.array(xs)
        ys = np.array(ys)
        
        try:
            # Sort by x for proper fitting
            sort_idx = np.argsort(xs)
            xs_sorted = xs[sort_idx]
            ys_sorted = ys[sort_idx]
            
            params, *_ = curve_fit(
                lowpoly,
                xs_sorted,
                ys_sorted,
                p0=[0, np.mean(np.diff(ys_sorted) / np.diff(xs_sorted)) if len(xs_sorted) > 1 else 0, 
                    np.mean(ys_sorted)]
            )
            self.current_params = params
            self.current_width = 2 * np.nanmedian(self.widths)
            
        except (ValueError, RuntimeError, TypeError) as e:
            print(f"Trace fit failed: {e}")
            self.current_params = None
            self.current_width = None
    
    def get_results(self):
        """Return the fitting results"""
        if not self.accepted or self.current_params is None:
            return None, None, None, None
        
        xs, ys = zip(*self.points)

        if np.isnan(self.current_width):
            self.current_width = 10.0

        return self.current_params, self.current_width, np.array(xs), np.array(ys)


# Main trace fitting function (enhanced version of your original code)
def fit_spectrum_trace(image, reduction_options, progress_window):
    """
    Enhanced trace fitting function with multi-trace and low quality data support
    """
    
    if reduction_options.low_quality_mode:
        # Manual point selection for low quality data
        return fit_lowquality_manual(image, reduction_options, progress_window)
    elif reduction_options.multitrace:
        # Multi-trace interactive behavior
        return fit_multitrace_interactive(image, reduction_options, progress_window)
    else:
        # Original single-trace behavior
        return fit_single_trace(image, reduction_options, progress_window)

def fit_single_trace(image, reduction_options, progress_window):
    """Your original trace fitting code (unchanged)"""
    ycenters = []
    xcenters = []
    width = []
    progress_window.update_current("Fitting Trace...", 0.24)
    
    for i in np.linspace(25, image.shape[1] - 10, 31):
        if 1013 < i < 1017:  # Ignore bad column
            continue
        data = np.median(image[:, int(i - 15):int(i + 15)], axis=1)
        xarr = np.arange(len(data))
        data = median_filter(data, 7)
        if np.all(data < 5):
            continue
        try:
            params, *_ = curve_fit(gaussian,
                                  xarr,
                                  data,
                                  p0=[
                                      5 * np.max(data), xarr[np.argmax(data)], 10, np.median(data)
                                  ],
                                  bounds=[
                                      [0, len(data) * 1 / 4, 1, -np.inf],
                                      [np.inf, len(data) * 3 / 4, len(xarr)/4, np.inf]
                                  ],
                                  maxfev=10000)
        except (ValueError, RuntimeError):
            continue
            
        width.append(params[2])
        xcenters.append(int(i))
        ycenters.append(params[1])
        
    width = 2 * np.median(width)
    
    try:
        params, *_ = curve_fit(lowpoly,
                              xcenters,
                              ycenters,
                              p0=[0, np.mean(np.diff(ycenters) / np.diff(xcenters)), np.mean(ycenters)])
    except (ValueError, RuntimeError, TypeError) as e:
        print(f"Error fitting polynomial to trace: {e}")
        return None, None, None, None
    
    xcenters = np.array(xcenters)
    ycenters = np.array(ycenters)
    
    # Check if params is valid before using it
    if params is None:
        print("Error: Polynomial fit returned None")
        return None, None, None, None
        
    resids = ((ycenters - lowpoly(xcenters, *params))**2)/lowpoly(xcenters, *params)
    outsidestd = resids > 1 * np.std(resids)+np.mean(resids)
    
    if np.sum(outsidestd.astype(int)) > 0 and not np.sum(outsidestd.astype(int)) > 0.5 * len(xcenters):
        try:
            params, *_ = curve_fit(lowpoly,
                                  xcenters[~outsidestd],
                                  ycenters[~outsidestd],
                                  p0=[0, np.mean(np.diff(ycenters) / np.diff(xcenters)), np.mean(ycenters)])
        except (ValueError, RuntimeError, TypeError) as e:
            print(f"Error in outlier-filtered polynomial fit: {e}")
            # Keep the original params if the filtered fit fails
    
    # Original debug plotting code
    if reduction_options.debugimages:
        xspace = np.linspace(0, image.shape[1], 1000)
        fig, axs = plt.subplots(2, 1, figsize=(4.8 * 16 / 9, 4.8))
        axs[0].plot(xspace, lowpoly(xspace, *params), zorder=1)
        axs[0].scatter(xcenters, ycenters, color="red", marker="x", zorder=5)
        norm_l = np.percentile(image, 1)
        norm_h = np.percentile(image, 99.9)
       
        if norm_l == 0:
            norm_l = 1
               
        axs[1].imshow(image, cmap=grey_cmap, norm=colors.LogNorm(norm_l, norm_h), aspect="auto")
        axs[1].plot(xspace, lowpoly(xspace, *params), color="lime", linewidth=0.5)
        axs[1].plot(xspace, lowpoly(xspace, *params) - reduction_options.skyfluxsep, color="red", linestyle="--", linewidth=0.5)
        axs[1].plot(xspace, lowpoly(xspace, *params) + reduction_options.skyfluxsep, color="red", linestyle="--", linewidth=0.5)
        axs[1].plot(xspace, lowpoly(xspace, *params) + width, color="lime", linestyle="--", linewidth=0.5)
        axs[1].plot(xspace, lowpoly(xspace, *params) - width, color="lime", linestyle="--", linewidth=0.5)
        plt.tight_layout()
        plt.show()
    
    return params, width, xcenters, ycenters

def fit_multitrace_interactive(image, reduction_options, progress_window):
    """New multi-trace interactive fitting function with retry capability"""
    
    while True:  # Loop until user accepts a trace or exits
        # Interactive bound selection
        print("Starting interactive trace selection...")
        selector = InteractiveTraceSelector(image)
        y_min, y_max = selector.get_bounds()
        print(f"Selected Y bounds: {y_min:.1f} to {y_max:.1f}")
        
        # Fit trace with constrained bounds
        ycenters = []
        xcenters = []
        width = []
        progress_window.update_current("Fitting Trace...", 0.24)
        
        for i in np.linspace(25, image.shape[1] - 10, 31):
            if 1013 < i < 1017:  # Ignore bad column
                continue
            data = np.min(image[:, int(i - 5):int(i + 5)], axis=1)
            data = median_filter(data, 5)
            if np.all(data < 5):
                continue
            xarr = np.arange(len(data))
            try:
                params, *_ = curve_fit(gaussian,
                                      xarr,
                                      data,
                                      p0=[
                                          5 * np.max(data), (y_min + y_max) / 2, 10, np.median(data)  # Fixed initial guess
                                      ],
                                      bounds=[
                                          [0, y_min, 1, -np.inf],  # Use interactive bounds
                                          [np.inf, y_max, len(xarr)/4, np.inf]  # Use interactive bounds
                                      ],
                                      maxfev=10000)
            except (ValueError, RuntimeError):
                continue
                
            width.append(params[2])
            xcenters.append(int(i))
            ycenters.append(params[1])
        
        if not ycenters:
            print("Warning: No valid trace points found with selected bounds!")
            print("Please try selecting different bounds.")
            continue  # Go back to bound selection
            
        width = 2 * np.median(width)
        
        try:
            params, *_ = curve_fit(lowpoly,
                                  xcenters,
                                  ycenters,
                                  p0=[0, np.mean(np.diff(ycenters) / np.diff(xcenters)), np.mean(ycenters)])
        except (ValueError, RuntimeError, TypeError) as e:
            print(f"Error fitting polynomial to trace: {e}")
            print("Please try selecting different bounds.")
            continue  # Go back to bound selection
        
        xcenters = np.array(xcenters)
        ycenters = np.array(ycenters)
        
        # Check if params is valid before using it
        if params is None:
            print("Error: Polynomial fit returned None")
            print("Please try selecting different bounds.")
            continue  # Go back to bound selection
            
        resids = ((ycenters - lowpoly(xcenters, *params))**2)/lowpoly(xcenters, *params)
        outsidestd = resids > 1 * np.std(resids)+np.mean(resids)
        
        if np.sum(outsidestd.astype(int)) > 0 and not np.sum(outsidestd.astype(int)) > 0.5 * len(xcenters):
            try:
                params, *_ = curve_fit(lowpoly,
                                      xcenters[~outsidestd],
                                      ycenters[~outsidestd],
                                      p0=[0, np.mean(np.diff(ycenters) / np.diff(xcenters)), np.mean(ycenters)])
            except (ValueError, RuntimeError, TypeError) as e:
                print(f"Error in outlier-filtered polynomial fit: {e}")
                # Keep the original params if the filtered fit fails
        
        # Interactive validation
        print("Validating fitted trace...")
        validator = TraceValidator(image, xcenters, ycenters, params, width, reduction_options)
        accepted = validator.validate()
        
        if accepted:
            print("Trace accepted!")
            return params, width, xcenters, ycenters
        else:
            print("Trace rejected! Returning to bound selection...")
            # Continue the loop to restart the process


def fit_lowquality_manual(image, reduction_options, progress_window):
    """Manual trace fitting for low quality data"""
    
    while True:
        progress_window.update_current("Manual trace selection...", 0.24)
        
        print("Starting manual trace selection for low quality data...")
        print("Left-click to add points along the spectrum trace")
        print("Right-click to remove points")
        print("Press Enter to accept, Escape to cancel")
        
        selector = ManualTraceSelector(image, reduction_options)
        params, width, xcenters, ycenters = selector.get_results()
        
        if params is None:
            # User cancelled or no valid fit
            response = input("No valid trace. Retry? (y/n): ").strip().lower()
            if response == 'y':
                continue
            else:
                return None, None, None, None
        
        # Show validation
        print("Validating fitted trace...")
        validator = TraceValidator(image, xcenters, ycenters, params, width, reduction_options)
        accepted = validator.validate()
        
        if accepted:
            print("Trace accepted!")
            return params, width, xcenters, ycenters
        else:
            print("Trace rejected! Restarting manual selection...")
            continue

def get_montecarlo_results(reduction_options):
    i = 0
    data_list = []
    while os.path.isfile(f".temp/mcmkc_output{i}.txt"):
        data_list.append(np.loadtxt(f".temp/mcmkc_output{i}.txt", delimiter=",", dtype=float))
        i += 1

    data = np.concatenate(data_list)

    threshold = 0
    pct = 0.1
    d= []
    if len(data) > 100000:
        while np.sum(data[:, -1] < threshold) < 100000:
            threshold = np.percentile(data[:, -1], pct)
            d = data[data[:, -1] < threshold]
            pct += 0.05
    else:
        threshold = np.percentile(data[:, -1], pct)
        d = data[data[:, -1] < threshold]

    params = []
    nbins = int(np.ceil(2 * (len(data[:, -1]) ** (1 / 3))))

    for i in range(4):
        hist, bin_edges = np.histogram(data[:, i], weights=1 / data[:, -1], bins=nbins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Fit the Gaussian to the histogram data
        popt, pcov = curve_fit(markov_gaussian, bin_centers, hist,
                               p0=[np.max(hist), bin_centers[np.argmax(hist)], np.std(data[:, i])], maxfev=1000000)

        # Extract the fitting parameters and their errors
        amp, mean, std = popt
        amp_err, mean_err, std_err = np.sqrt(np.diag(pcov))

        params.append(mean)
        if reduction_options.debugimages:
            plt.hist(data[:, i], weights=1 / data[:, -1], bins=nbins, alpha=0.6, label='Data')
            x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)
            y_fit = markov_gaussian(x_fit, *popt)
            plt.plot(x_fit, y_fit, color='red', label='Gaussian fit')
            plt.title(f"MCMC result histogram: {['Offset', 'Linear', 'Quadratic', 'Cubic'][i]} Parameter")
            plt.xlabel('Data')
            plt.ylabel('Frequency')
            plt.legend()
            plt.tight_layout()
            plt.show()
    return params


def extract_spectrum(frame, master_flat, complamplist, frame_config, reduction_options, progress_window):
    comp_header = frame_config.comparisonframes.master_frame.header

    if reduction_options.offset_zero == "Auto":
        if "GRATING" in frame_config.comparisonframes.master_frame.header:
            if "930" in comp_header["GRATING"]:
                d_grating = 930.
            elif "2100" in comp_header["GRATING"]:
                d_grating = 2100.
            else:
                d_grating = int(re.search(r'\d+', comp_header["GRATING"]).group())
        else:
            raise ValueError("No GRATING Keyword found in Comparison Lamp header but 'Auto' set as offset_zero! Specify a reasonable offset_zero manually!")

        try:
           central_wl = get_central_wavelength(comp_header["GRT_ANG"], comp_header["CAM_ANG"], d_grating)
        except:
            raise ValueError("Could not determine central Wavelength! Is the header compatible? Try to manually specify offset_zero!")
    else:
        reduction_options.offset_zero = float(reduction_options.offset_zero)
        central_wl = None

    progress_window.update_current("Correcting Spectral Image...", 0.21)
    frame.apply_corrections(master_flat, frame_config)
    image = frame.data
     
    traceparams, width, xcenters, ycenters = fit_spectrum_trace(image, reduction_options, progress_window)

    resids = ((ycenters - lowpoly(xcenters, *traceparams))**2)/lowpoly(xcenters, *traceparams)
    outsidestd = resids > 1 * np.std(resids)+np.mean(resids)
    if np.sum(outsidestd.astype(int)) > 0 and not np.sum(outsidestd.astype(int)) > 0.5 * len(xcenters):
        traceparams, _ = curve_fit(lowpoly,
                              xcenters[~outsidestd],
                              ycenters[~outsidestd],
                              p0=[0, np.mean(np.diff(ycenters) / np.diff(xcenters)), np.mean(ycenters)])

    #if reduction_options.debugimages:
    #    xspace = np.linspace(0, image.shape[1], 1000)
    #    fig, axs = plt.subplots(2, 1, figsize=(4.8 * 16 / 9, 4.8))
    #    axs[0].plot(xspace, lowpoly(xspace, *traceparams), zorder=1)
    #    axs[0].scatter(xcenters, ycenters, color="red", marker="x", zorder=5)
    #    norm_l = np.percentile(image, 1)
    #    norm_h = np.percentile(image, 99.9)
    #   
    #    if norm_l == 0:
    #        norm_l = 1
               
    #    axs[1].imshow(image, cmap=grey_cmap, norm=colors.LogNorm(norm_l, norm_h), aspect="auto")
    #    axs[1].plot(xspace, lowpoly(xspace, *traceparams), color="lime", linewidth=0.5)
    #    axs[1].plot(xspace, lowpoly(xspace, *traceparams) - reduction_options.skyfluxsep, color="red", linestyle="--", linewidth=0.5)
    #    axs[1].plot(xspace, lowpoly(xspace, *traceparams) + reduction_options.skyfluxsep, color="red", linestyle="--", linewidth=0.5)
    #    axs[1].plot(xspace, lowpoly(xspace, *traceparams) + width, color="lime", linestyle="--", linewidth=0.5)
    #    axs[1].plot(xspace, lowpoly(xspace, *traceparams) - width, color="lime", linestyle="--", linewidth=0.5)
    #    plt.tight_layout()
    #    plt.show()

    progress_window.update_current("Reading Trace...", 0.33)

    image64 = image.astype(np.float64)
    master_comp64 = complamplist.master_frame.data.astype(np.float64)

    pixel = np.arange(image.shape[1]).astype(np.float64)
    flux = np.array([get_flux(image64, p, lowpoly(p, *traceparams), width, reduction_options.use_boxcut) for p in pixel])
    compflux = np.array([get_flux(master_comp64, p, lowpoly(p, *traceparams), width, reduction_options.use_boxcut) for p in pixel])
    uskyflx = np.array([get_flux(image64, p, lowpoly(p, *traceparams) + reduction_options.skyfluxsep, width, reduction_options.use_boxcut) for p in pixel])
    lskyflx = np.array([get_flux(image64, p, lowpoly(p, *traceparams) - reduction_options.skyfluxsep, width, reduction_options.use_boxcut) for p in pixel])

    skyflx = np.minimum(uskyflx, lskyflx)
    flux -= skyflx

    compflux_cont = minimum_filter(compflux, 10)
    compflux -= compflux_cont

    if reduction_options.cosmicrejection:
        flux[flux > 3 * np.median(flux)] = np.nan

    lines = np.genfromtxt(reduction_options.linelist)[:, 0]

    progress_window.update_current("Solving Wavelength Array...", 0.45)

    complamplist.get_compparams(traceparams)

    if complamplist.compparams is None:
        def call_fitlines_markov(compspec_x, compspec_y):
            print("Finding wavelength solution, this may take some time...")
            compspec_x = np.array(compspec_x, dtype=np.double)
            compspec_y = np.array(compspec_y, dtype=np.double)

            compspec_y /= maximum_filter(compspec_y, reduction_options.lampfilterwindow)

            compspec_y = gaussian_filter(compspec_y, 1)

            if not os.path.isdir("./.temp"):
                os.mkdir(".temp")
            else:
                shutil.rmtree("./.temp")
                os.mkdir(".temp")

            offset_zero = 0
            if central_wl is None:
                offset_zero = reduction_options.offset_zero
            else:
                offset_zero = central_wl - reduction_options.extent_zero / 2

            np.savetxt("./.temp/compspec_x.txt", compspec_x, fmt="%.9e")
            np.savetxt("./.temp/compspec_y.txt", compspec_y, fmt="%.9e")
            np.savetxt("./.temp/lines.txt", lines, fmt="%.9e")
            np.savetxt("./.temp/arguments.txt", np.array([len(lines), len(compspec_x), reduction_options.sampleamt, offset_zero,
                                                         reduction_options.extent_zero / len(compspec_x), reduction_options.quad_zero, reduction_options.cube_zero,
                                                         reduction_options.offset_stepsize, reduction_options.linear_stepsize, reduction_options.quad_stepsize, reduction_options.cube_stepsize, 
                                                         reduction_options.c_cov, reduction_options.s_cov, reduction_options.q_cov, reduction_options.cub_cov,
                                                         reduction_options.accept_param]), fmt="%.9e")

            

            process = subprocess.Popen(
                "./linefit .temp/compspec_x.txt .temp/compspec_y.txt .temp/lines.txt .temp/arguments.txt 1",
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True)
            for line in process.stdout:
                print(line, end='')  # end='' prevents adding extra newlines

            # Wait for the process to finish and capture any remaining output or errors
            stdout, stderr = process.communicate()

            progress_window.update_current("Reading Solution...", 0.75)

            result = get_montecarlo_results(reduction_options)

            return result

        result = call_fitlines_markov(pixel, compflux)

        compparams = [result[0], result[1], result[2], result[3]]

        complamplist.compparams = compparams
        complamplist.save_compparams(traceparams)
    final_wl_arr = line(pixel, 0.0, *complamplist.compparams[::-1])

    if reduction_options.debugimages:
        plt.plot(final_wl_arr, compflux.min() - np.nanmax(flux) + flux, linewidth=1, color="gray")
        plt.plot(final_wl_arr, gaussian_filter(compflux, 1), color="darkred")
        for b in lines:
            plt.axvline(b, linestyle="--", color="darkgreen", zorder=-5)
        plt.show()

    progress_window.update_current("Applying Barycentric Correction...", 0.85)

    sc = SkyCoord(ra=frame.star_info["ra"] * u.deg, dec=frame.star_info["dec"] * u.deg)
    barycorr = sc.radial_velocity_correction(obstime=frame.obstime, location=reduction_options.earth_location)
    barycorr = barycorr.to(u.km / u.s)
    barycorr = barycorr.value

    final_wl_arr = wlshift(final_wl_arr, barycorr)

    final_wl_arr, flux, flx_std = fluxstatistics(final_wl_arr, flux)

    return final_wl_arr, flux, flx_std


def get_matching_complamps(frame, frame_config, reduction_options):
    mode = reduction_options.comparison_match_mode
    pos = reduction_options.comparison_match_position
    value = reduction_options.comparison_match_value

    compframes = frame_config.comparisonframes
    filtered = []

    try:
        if mode == "count":
            # Filter first, then sort
            if pos == "before":
                filtered = [f for f in compframes if f.obstime < frame.obstime]
                filtered.sort(key=lambda f: f.obstime, reverse=True)
                filtered = filtered[:value]
                filtered.reverse()  # maintain chronological order
            elif pos == "after":
                filtered = [f for f in compframes if f.obstime > frame.obstime]
                filtered.sort(key=lambda f: f.obstime)
                filtered = filtered[:value]
            elif pos == "around":
                # Get all, sort by proximity to the science frame
                filtered = sorted(compframes, key=lambda f: abs((f.obstime - frame.obstime).sec))
                filtered = filtered[:value]

        elif mode == "time":
            dt = TimeDelta(value * 60, format='sec')  # convert minutes to seconds
            if pos == "before":
                filtered = [f for f in compframes if frame.obstime - dt <= f.obstime < frame.obstime]
            elif pos == "after":
                filtered = [f for f in compframes if frame.obstime < f.obstime <= frame.obstime + dt]
            elif pos == "around":
                filtered = [f for f in compframes if abs((f.obstime - frame.obstime).sec) <= dt.sec]

            filtered.sort(key=lambda f: f.obstime)

        else:
            raise ValueError("Invalid comparison_match_mode. Expected 'count' or 'time'.")
    except TypeError:
        raise ValueError("Could not find DATE-OBS field in header. Did you select the right Frames? Did you select the right header HDU index?")

    return filtered


def save_to_ascii(frame, wl, flx, flx_std, reduction_options, frame_config):
    ordered_cols = ["name", "source_id", "ra", "dec",
                        "file", "SPEC_CLASS", "bp_rp", "gmag", "nspec",
                        "pmra", "pmra_error", "pmdec", "pmdec_error", "parallax", "parallax_error"]
    trow_new = dict([(a, b) for a, b in frame.star_info.items() if a in ordered_cols])
    trow = dict(sorted(trow_new.items(), key=lambda pair: ordered_cols.index(pair[0])))
    
    outtable = None
    if os.path.isfile(reduction_options.outputfile_path):
        outtable = pd.read_csv(reduction_options.outputfile_path)
        if not trow["file"] in outtable["file"].to_list():
            outtable = pd.concat([outtable, pd.DataFrame([trow])])
    else:
        outtable = pd.DataFrame([trow])

    if not os.path.isdir(reduction_options.output_dict):
        os.mkdir(reduction_options.output_dict)
    fname = trow["file"].replace(".fits", "_01.txt")
    fname = reduction_options.output_dict + "/" + fname

    with open(fname.replace("_01.", "_mjd."), "w") as datefile:
        datefile.write(str(frame.obstime.mjd))

    outtable.to_csv(reduction_options.outputfile_path, index=False)
    outdata = np.stack((wl, flx, flx_std), axis=-1)
    np.savetxt(fname, outdata, fmt='%1.4f')


def reduce_data(reduction_options, frame_config, progress_window):
    print("Starting data reduction...")

    progress_window.update_overall("Rotating Frames...", 0.0)
    progress_window.update_current("Applying Rotation to Biases...", 0.0)
    for frame in frame_config.biases:
        frame.rotate()

    progress_window.update_current("Applying Rotation to Flats...", 0.2)
    for frame in frame_config.flats:
        frame.rotate()

    progress_window.update_current("Applying Rotation to Shifted Flats...", 0.4)
    for frame in frame_config.shiftedflats:
        frame.rotate()

    progress_window.update_current("Applying Rotation to Science Frames...", 0.6)
    for frame in frame_config.scienceframes:
        frame.rotate()

    progress_window.update_current("Applying Rotation to Comparison Frames...", 0.8)
    for frame in frame_config.comparisonframes:
        frame.rotate()
    
    progress_window.update_overall("Cropping Frames...", 0.02)
    progress_window.update_current("Creating Master Flat...", 0.0)	
    # Determine Image Crop and apply to all frames	

    if len(frame_config.shiftedflats.frames) != 0:
        master_flat = Frame(frame_config.flats.create_master_flat(second_image_list = frame_config.shiftedflats), reduction_options)
    else:
        master_flat = Frame(frame_config.flats.create_master_flat(), reduction_options)

    progress_window.update_current("Detecting Spectral Area...", 0.33)
    if not reduction_options.manual_crop:
        crop = master_flat.detect_spectral_area()
    else:
        crop = [[reduction_options.x_lo, master_flat.data.shape[1]-reduction_options.x_hi], [reduction_options.y_hi, master_flat.data.shape[0]-reduction_options.y_lo]]

    progress_window.update_current("Determined Cropping!", 0.33)

    if reduction_options.debugimages:
        plt.imshow(master_flat.data, cmap="Greys_r", zorder=1)
        if crop[0][0] >= 0:
            plt.axvline(crop[0][0], color="lime", zorder=5)
        if crop[0][1] <= master_flat.data.shape[1]:
            plt.axvline(crop[0][1], color="lime", zorder=5)
        if crop[1][0] >= 0:
            plt.axhline(crop[1][0], color="lime", zorder=5)
        if crop[1][1] <= master_flat.data.shape[0]:
            plt.axhline(crop[1][1], color="lime", zorder=5)
        plt.axis("off")
        plt.title("Detected Spectral Area of Frame")
        plt.tight_layout()
        plt.show()

    progress_window.update_current("Applying Cropping to Biases...", 0.66)

    for frame in frame_config.biases:
        frame.crop(crop)

    progress_window.update_current("Applying Cropping to Flats...", 0.66+0.066)
    for frame in frame_config.flats:
        frame.crop(crop)

    progress_window.update_current("Applying Cropping to Shifted Flats...", 0.66+2*0.066)
    for frame in frame_config.shiftedflats:
        frame.crop(crop)

    progress_window.update_current("Applying Cropping to Science Frames...", 0.66+3*0.066)
    for frame in frame_config.scienceframes:
        frame.crop(crop)

    progress_window.update_current("Applying Cropping to Comparison Frames...", 0.66+4*0.066)
    for frame in frame_config.comparisonframes:
        frame.crop(crop)

    progress_window.update_current("Cropping Applied!", 1.0)

    progress_window.update_overall("Creating Master Bias...", 0.05)
    progress_window.update_current("Adding Bias Frames...", 0.0)

    master_bias = Frame(frame_config.biases.create_master_image()[0], reduction_options)

    progress_window.update_overall("Master Bias Created!", 0.1)
    progress_window.update_current("Added Bias Frames!", 1.0)

    if reduction_options.debugimages:

        norm = colors.LogNorm(vmin=np.nanmin(master_bias.data[master_bias.data != 0]), vmax=np.nanmax(master_bias.data))

        plt.imshow(master_bias.data, cmap=grey_cmap, interpolation='none', norm=norm, zorder=1)
        plt.title("Master Bias Image")
        plt.axis("off")
        plt.tight_layout()
        plt.show()

    progress_window.update_overall("Creating Master Flat...", 0.1)
    progress_window.update_current("Adding Flat Frames...", 0.0)

    if len(frame_config.shiftedflats.frames) != 0:
        master_flat, master_continuum = frame_config.flats.create_master_flat(second_image_list=frame_config.shiftedflats, master_bias=master_bias, bounds=crop)
    else:
        master_flat, master_continuum = frame_config.flats.create_master_flat(second_image_list=None, master_bias=master_bias, bounds=crop)

    master_flat = Frame(master_flat, reduction_options)
    master_continuum = Frame(master_continuum, reduction_options)

    progress_window.update_overall("Master Flat Created!", 0.15)
    progress_window.update_current("Added Flat Frames!", 1.0)

    if reduction_options.debugimages:
        plt.imshow(master_flat.data, cmap="Greys_r", zorder=1)
        plt.title("Master Flat Image")
        plt.axis("off")
        plt.tight_layout()
        plt.show()

        plt.imshow(master_continuum.data, cmap="Greys_r", zorder=1)
        plt.title("Master Continuum Image")
        plt.axis("off")
        plt.tight_layout()
        plt.show()


    location = reduction_options.populate_earth_location()
    n_frames = len(frame_config.scienceframes.frames)
    for i, frame in enumerate(frame_config.scienceframes):
        progress_window.update_overall(f"Processing spectrum {i+1}/{n_frames}...", 0.15+0.85*(i / n_frames))
        progress_window.update_current("Processing Arc Frames...", 0.0)

        filtered_complamps = ComplampList(get_matching_complamps(frame, frame_config, reduction_options), reduction_options)
        if len(filtered_complamps.frames) == 0:
            print(f"WARNING: DID NOT FIND MATCHING COMPLAMP FRAMES WITHIN CONSTRAINTS FOR {frame.filepath}, skipping!")
            continue

        master_comp, master_comp_header = filtered_complamps.create_master_image(master_bias=master_bias, master_continuum=master_continuum, return_header=True)
        frame_config.comparisonframes.master_frame = Frame(master_comp, reduction_options)
        frame_config.comparisonframes.master_frame.header = master_comp_header


        if reduction_options.debugimages:
            plt.imshow(master_comp, cmap=grey_cmap, norm=colors.LogNorm(vmin=np.nanmin(master_comp[master_comp != 0]), vmax=np.nanmax(master_comp)), zorder=1)
            plt.title("Master Arc Lamp Image")
            plt.axis("off")
            plt.tight_layout()
            plt.show()

        progress_window.update_current("Getting Star Info...", 0.15)
        frame.get_star_info()
        progress_window.update_current("Extracting Spectrum...", 0.20)

        wl, flx, flx_std = extract_spectrum(frame, master_flat, filtered_complamps, frame_config, reduction_options, progress_window)

        progress_window.update_current("Saving Solution...", 0.9)	
        save_to_ascii(frame, wl, flx, flx_std, reduction_options, frame_config)

        progress_window.update_current("", 0.0)

    progress_window.update_overall("Reduction complete!", 1.0)
    time.sleep(0.5)
    progress_window.destroy()