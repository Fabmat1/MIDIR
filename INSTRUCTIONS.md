# MIDIR Data Reduction Instructions

Welcome to **MIDIR**, the Multi-Instrument Data Input Reduce-inator! This is a modular data reduction pipeline built specifically for **long slit spectrographs**. MIDIR is designed to be flexible and extensible, supporting multiple instruments and providing a clean graphical interface for configuring and executing the full data reduction workflow.

---

## Data Reduction Step-by-Step

### 1. Load a Preset (If Available)

If you're working with one of the supported instruments — **SOAR**, **EFOSC**, or **ALFOSC** — you can make your life easier by loading a preset configuration:

* Navigate to the **General Settings** tab.
* Choose the appropriate instrument from the **Preset** dropdown.

> ⚠️ **Important**: Load the preset **before** selecting any files. The preset determines critical settings like the HDU indices. These must be set correctly before data files are loaded.

### 2. Select Input Files

Next, go to the **File Selection** tab. Here you will load your calibration and science data:

* Use the buttons to select files for each category:

  * **Bias Frames**
  * **Flat Frames**
  * **Shifted Flats** *(optional; only relevant for SOAR to eliminate the Littrow ghost)*
  * **Arc Lamps**
  * **Science Frames**
* The number next to each button shows how many files are currently selected in that category.
* Files can be in `.fits, .fits.Z or .fits.gz` formats

### 3. Set Comparison Lamp Matching Mode

After selecting your files, choose how MIDIR should assign **comparison lamp frames** to each science frame. This is configured with two dropdowns and a numeric field:

* **Match Mode** (Dropdown 1):

  * `Time cutoff`: selects frames within a time window.
  * `Number cutoff`: selects a fixed number of comparison frames.
* **Relative Position** (Dropdown 2):

  * `Before`: selects frames that occur before the science frame.
  * `After`: selects frames that occur after the science frame.
  * `Around`: selects frames on both sides — this means:

    * For **time cutoff**: a ± time window.
    * For **number cutoff**: the *n* closest frames.
* **N/ Minutes**: enter the number of frames or the time window size (depending on the mode).

### 4. Start the Reduction

Once everything is configured:

* Click the **Start Reduction** button.
* MIDIR will begin processing all selected science frames with the chosen calibration data and settings.

### 5. Check the Output

MIDIR is not a perfect tool and can make mistakes - and so can you! To ensure proper data reduction try to enable **Show Plots** in the General Settings tab and check if the reduction steps work properly.

> ⚠️ When **Show Plots** is enabled the results of the MCMC wavelength calibration will be shown as histograms. If any of these histograms do not look approximately Gaussian, your wavelength solution might be unreliable!

---

## Advanced Configuration

If needed, you can modify existing presets or create your own to support additional long slit spectrographs. Here's how to configure MIDIR in more detail:

### General Settings

This tab contains essential options for handling input data:

* **Rotate Images**: If your input data has the spectral trace oriented vertically or the wavelength does not increase along the x-axis, select the appropriate rotation angle to correct this.

* **Sky Flux Separation**: MIDIR subtracts sky flux by sampling regions above and below the spectral trace. Use this setting to define the pixel offset between these sky traces and the spectral trace.

* **Telescope Location**: Select your telescope’s location. If it is not available in the dropdown, open an issue or submit a pull request.

* **Crop Frames Manually**: When disabled, MIDIR will attempt automatic cropping. This may be unreliable for instruments other than SOAR. Enable this to manually set image margins (in pixels) used for cropping around the spectral region.

* **Attempt Cosmic Rejection**: MIDIR rejects flux values that are >3σ from the spectrum’s standard deviation by default. Disable this if you're working with emission line spectra, where this behavior might remove valid data.

* **Use Boxcut Extraction**: By default, MIDIR uses a Gaussian-weighted extraction. Enable this option to instead use a simpler boxcut method.

### Wavelength Reduction Settings

These settings determine how MIDIR calibrates wavelength solutions along the spectral trace:

* **Linelists**: MIDIR includes default linelists for supported instruments. To use a custom list, add a whitespace-separated two-column ASCII file (wavelengths in angstroms) to the `linelists/` directory. Your file will appear in the linelist dropdown automatically.

* **Amount of MC Samples**: Sets the number of samples per CPU core used in MCMC sampling. The total number of samples is scaled by the number of CPU cores available.

* **Acceptance Rate Modifier**: MIDIR uses a modified Metropolis-Hastings algorithm for fitting. This value controls the maximum allowed χ² ratio between a new proposal and the current solution. Ratios above this threshold are always rejected.

* **MCMC Polynomial Parameter Configuration**: This table defines the initial guesses, step sizes, and bounds for each polynomial coefficient used to transform pixels into wavelengths. Ensure the bounds are broad enough to explore valid solutions, but not too wide to avoid degeneracies.

> ⚠️ Ensure the printed acceptance rate (visible in the terminal during MCMC) is between **7% and 30%** for good convergence.

> ⚠️ Due to a legacy design choice, the zero point of the linear wavelength term defines the **total wavelength range** (in Å), not its literal intercept. This will be revised in a future version.

* **Arc Lamp Maximum Filter Size**: Specifies the width (in pixels) of a maximum filter applied to the arc lamp trace, normalizing it for MCMC use. The filter size should be such that usable emission lines scale up to 1.0.

---

Enjoy using MIDIR!
