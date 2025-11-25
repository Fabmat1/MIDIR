import numpy as np

def fit_pixel_wavelength(pixels, wavelengths):
    """
    Fit polynomial: wavelength = a + b*pixel + c*pixel^2 + d*pixel^3
    
    Returns coefficients and covariance matrix for uncertainty estimation
    """
    pixels = np.array(pixels, dtype=float)
    wavelengths = np.array(wavelengths, dtype=float)
    
    # polyfit with covariance for uncertainty estimation
    coeffs, cov = np.polyfit(pixels, wavelengths, 3, cov=True)
    
    d, c, b, a = coeffs
    # Covariance diagonal gives variances (in descending order: d, c, b, a)
    uncertainties = np.sqrt(np.diag(cov))
    d_err, c_err, b_err, a_err = uncertainties
    
    return {
        'a': a, 'b': b, 'c': c, 'd': d,
        'a_err': a_err, 'b_err': b_err, 'c_err': c_err, 'd_err': d_err,
        'coeffs_descending': coeffs,
        'covariance': cov
    }

def format_value(val):
    """Format value: scientific notation only for |val| < 1e-4"""
    if val == 0:
        return "0"
    elif abs(val) < 1e-4:
        return f"{val:.6e}"
    else:
        # Determine appropriate decimal places based on magnitude
        if abs(val) >= 1000:
            return f"{val:.4f}"
        elif abs(val) >= 1:
            return f"{val:.6f}"
        else:
            return f"{val:.8f}"

def calculate_mcmc_parameters(result, pixels, wavelengths, n_pixels_total):
    """
    Calculate MCMC stepsizes and bounds based on fit results
    """
    a, b, c, d = result['a'], result['b'], result['c'], result['d']
    a_err, b_err, c_err, d_err = result['a_err'], result['b_err'], result['c_err'], result['d_err']
    
    pixel_range = max(pixels) - min(pixels)
    wavelength_range_data = max(wavelengths) - min(wavelengths)
    
    # Wavelength range for legacy code (b * n_pixels)
    wavelength_range_legacy = b * n_pixels_total
    
    # Stepsizes: based on uncertainties, scaled for efficient MCMC exploration
    step_a = (a_err / 2)
    step_b = (b_err / 2)
    step_c = (c_err / 2)
    step_d = (d_err / 2)
    
    # Symmetric bounds margins (distance from initial guess to bound)
    margin_a = (10 * a_err)
    margin_b = (10 * b_err)
    margin_c = (10 * c_err)
    margin_d = (10 * d_err)
    
    return {
        'wavelength_range_legacy': wavelength_range_legacy,
        'n_pixels_total': n_pixels_total,
        'steps': {'a': step_a, 'b': step_b, 'c': step_c, 'd': step_d},
        'margins': {'a': margin_a, 'b': margin_b, 'c': margin_c, 'd': margin_d}
    }

def print_results(result, mcmc_params, pixels, wavelengths):
    """Print the results, residuals, and MCMC parameters"""
    a, b, c, d = result['a'], result['b'], result['c'], result['d']
    a_err, b_err, c_err, d_err = result['a_err'], result['b_err'], result['c_err'], result['d_err']
    
    print("\n" + "=" * 70)
    print("POLYNOMIAL FIT RESULTS")
    print("wavelength = a + b*pixel + c*pixel² + d*pixel³")
    print("=" * 70)
    
    print(f"\n{'Parameter':<12} {'Value':>20} {'Uncertainty':>20}")
    print("-" * 54)
    print(f"{'a':<12} {format_value(a):>20} {format_value(a_err):>20}")
    print(f"{'b':<12} {format_value(b):>20} {format_value(b_err):>20}")
    print(f"{'c':<12} {format_value(c):>20} {format_value(c_err):>20}")
    print(f"{'d':<12} {format_value(d):>20} {format_value(d_err):>20}")
    
    print("\n" + "-" * 70)
    print("VERIFICATION (Pixel → Calculated vs Actual Wavelength)")
    print("-" * 70)
    print(f"{'Pixel':>10} {'Actual (Å)':>14} {'Calculated (Å)':>16} {'Residual (Å)':>14}")
    
    total_residual_sq = 0
    for px, wl_actual in zip(pixels, wavelengths):
        wl_calc = a + b*px + c*px**2 + d*px**3
        residual = wl_actual - wl_calc
        total_residual_sq += residual**2
        print(f"{px:>10.1f} {wl_actual:>14.2f} {wl_calc:>16.2f} {residual:>14.4f}")
    
    rms = np.sqrt(total_residual_sq / len(pixels))
    print(f"\nRMS Error: {rms:.4f} Angstrom")
    
    # MCMC Parameters
    print("\n" + "=" * 70)
    print("MCMC PARAMETERS")
    print("=" * 70)
    
    print(f"\n--- Legacy Wavelength Range (b × n_pixels) ---")
    print(f"n_pixels_total: {mcmc_params['n_pixels_total']}")
    print(f"wavelength_range = {format_value(mcmc_params['wavelength_range_legacy'])} Angstrom")
    
    print(f"\n--- Initial Guesses ---")
    print(f"a = {format_value(a)}")
    print(f"b = {format_value(b)}  (or wavelength_range = {format_value(mcmc_params['wavelength_range_legacy'])})")
    print(f"c = {format_value(c)}")
    print(f"d = {format_value(d)}")
    
    print(f"\n--- Suggested Stepsizes ---")
    steps = mcmc_params['steps']
    print(f"step_a = {format_value(steps['a'])}")
    print(f"step_b = {format_value(steps['b'])}")
    print(f"step_c = {format_value(steps['c'])}")
    print(f"step_d = {format_value(steps['d'])}")
    
    print(f"\n--- Suggested Bounds (± margin from initial guess) ---")
    margins = mcmc_params['margins']
    print(f"a: ± {format_value(margins['a'])}")
    print(f"b: ± {format_value(margins['b'])}")
    print(f"c: ± {format_value(margins['c'])}")
    print(f"d: ± {format_value(margins['d'])}")

def display_entries(pixels, wavelengths):
    """Display current entries in a table"""
    if not pixels:
        print("\n  (No entries yet)")
        return
    
    print(f"\n{'Index':<8} {'Pixel':<15} {'Wavelength (Å)':<15}")
    print("-" * 40)
    for i, (px, wl) in enumerate(zip(pixels, wavelengths)):
        print(f"{i:<8} {px:<15.2f} {wl:<15.2f}")
    print(f"\nTotal entries: {len(pixels)}")

def interactive_mode():
    """Interactive mode for entering pixel-wavelength pairs"""
    pixels = []
    wavelengths = []
    n_pixels_total = None
    
    print("\n" + "=" * 70)
    print("INTERACTIVE PIXEL-WAVELENGTH CALIBRATION")
    print("=" * 70)
    print("\nCommands:")
    print("  [Enter pixel,wavelength] - Add a new pair (e.g., '1802,5875.62')")
    print("  'list' or 'l'            - Show all current entries")
    print("  'delete N' or 'd N'      - Delete entry at index N")
    print("  'edit N' or 'e N'        - Edit entry at index N")
    print("  'clear'                  - Clear all entries")
    print("  'npix N'                 - Set total number of pixels (for legacy output)")
    print("  'fit' or 'f'             - Perform fit and show results")
    print("  'done' or 'q'            - Finish and show final results")
    print("  'help' or 'h'            - Show this help message")
    print("-" * 70)
    
    while True:
        try:
            user_input = input("\n> ").strip()
            
            if not user_input:
                continue
            
            # Help
            if user_input.lower() in ['help', 'h']:
                print("\nCommands:")
                print("  [pixel,wavelength] - Add pair (e.g., '1802,5875.62' or '1802 5875.62')")
                print("  'list' / 'l'       - Show all entries")
                print("  'delete N' / 'd N' - Delete entry at index N")
                print("  'edit N' / 'e N'   - Edit entry at index N")
                print("  'clear'            - Clear all entries")
                print("  'npix N'           - Set total pixel count")
                print("  'fit' / 'f'        - Perform fit")
                print("  'done' / 'q'       - Finish")
                continue
            
            # List entries
            if user_input.lower() in ['list', 'l']:
                display_entries(pixels, wavelengths)
                if n_pixels_total:
                    print(f"Total pixels (for legacy): {n_pixels_total}")
                continue
            
            # Clear all
            if user_input.lower() == 'clear':
                confirm = input("Clear all entries? (y/n): ").strip().lower()
                if confirm == 'y':
                    pixels = []
                    wavelengths = []
                    print("All entries cleared.")
                continue
            
            # Set n_pixels_total
            if user_input.lower().startswith('npix'):
                parts = user_input.split()
                if len(parts) == 2:
                    try:
                        n_pixels_total = int(parts[1])
                        print(f"Total pixels set to: {n_pixels_total}")
                    except ValueError:
                        print("Invalid number. Usage: npix 2048")
                else:
                    print("Usage: npix <number> (e.g., 'npix 2048')")
                continue
            
            # Delete entry
            if user_input.lower().startswith(('delete ', 'd ')):
                parts = user_input.split()
                if len(parts) == 2:
                    try:
                        idx = int(parts[1])
                        if 0 <= idx < len(pixels):
                            removed_px = pixels.pop(idx)
                            removed_wl = wavelengths.pop(idx)
                            print(f"Deleted entry {idx}: pixel={removed_px}, wavelength={removed_wl}")
                        else:
                            print(f"Invalid index. Valid range: 0-{len(pixels)-1}")
                    except ValueError:
                        print("Invalid index. Usage: delete 3")
                else:
                    print("Usage: delete <index> (e.g., 'delete 3' or 'd 3')")
                continue
            
            # Edit entry
            if user_input.lower().startswith(('edit ', 'e ')):
                parts = user_input.split()
                if len(parts) == 2:
                    try:
                        idx = int(parts[1])
                        if 0 <= idx < len(pixels):
                            print(f"Current: pixel={pixels[idx]}, wavelength={wavelengths[idx]}")
                            new_input = input("Enter new pixel,wavelength (or press Enter to cancel): ").strip()
                            if new_input:
                                # Parse new values
                                for sep in [',', ' ', '\t']:
                                    if sep in new_input:
                                        new_parts = [p.strip() for p in new_input.split(sep) if p.strip()]
                                        if len(new_parts) == 2:
                                            new_px = float(new_parts[0])
                                            new_wl = float(new_parts[1])
                                            pixels[idx] = new_px
                                            wavelengths[idx] = new_wl
                                            print(f"Updated entry {idx}: pixel={new_px}, wavelength={new_wl}")
                                            break
                                else:
                                    print("Could not parse input. Use format: pixel,wavelength")
                            else:
                                print("Edit cancelled.")
                        else:
                            print(f"Invalid index. Valid range: 0-{len(pixels)-1}")
                    except ValueError:
                        print("Invalid input.")
                else:
                    print("Usage: edit <index> (e.g., 'edit 3' or 'e 3')")
                continue
            
            # Fit
            if user_input.lower() in ['fit', 'f', 'done', 'q']:
                if len(pixels) < 4:
                    print(f"Need at least 4 data points for cubic fit. Currently have {len(pixels)}.")
                    if user_input.lower() in ['done', 'q']:
                        confirm = input("Exit anyway? (y/n): ").strip().lower()
                        if confirm == 'y':
                            print("Exiting without fit.")
                            return None, None
                    continue
                
                if n_pixels_total is None:
                    while True:
                        npix_input = input("Enter total number of pixels (for legacy wavelength_range): ").strip()
                        try:
                            n_pixels_total = int(npix_input)
                            break
                        except ValueError:
                            print("Please enter a valid integer.")
                
                # Perform fit
                result = fit_pixel_wavelength(pixels, wavelengths)
                mcmc_params = calculate_mcmc_parameters(result, pixels, wavelengths, n_pixels_total)
                print_results(result, mcmc_params, pixels, wavelengths)
                
                if user_input.lower() in ['done', 'q']:
                    return result, mcmc_params
                continue
            
            # Try to parse as pixel,wavelength pair
            parsed = False
            for sep in [',', ' ', '\t']:
                if sep in user_input:
                    parts = [p.strip() for p in user_input.split(sep) if p.strip()]
                    if len(parts) == 2:
                        try:
                            px = float(parts[0])
                            wl = float(parts[1])
                            pixels.append(px)
                            wavelengths.append(wl)
                            print(f"Added: pixel={px}, wavelength={wl} Å  [Entry {len(pixels)-1}]")
                            parsed = True
                            break
                        except ValueError:
                            pass
            
            if not parsed:
                print("Could not parse input. Enter 'pixel,wavelength' or type 'help'.")
                
        except KeyboardInterrupt:
            print("\n\nInterrupted. Type 'done' to finish or 'q' to quit.")
        except Exception as e:
            print(f"Error: {e}")

def hardcoded_mode():
    """Run with hardcoded values"""
    # ============================================================
    # INPUT YOUR DATA HERE
    # ============================================================
    
    pixels = [
        1802,
        1604,
        1424,
        1287,
        1228,
        1070,
        912,
        702
    ]
    
    wavelengths = [
        5875.62,
        5460.74,
        5085.82,
        4799.91,
        4678.15,
        4358.33,
        4046.56,
        3650.15
    ]
    
    n_pixels_total = 2048
    
    # ============================================================
    
    result = fit_pixel_wavelength(pixels, wavelengths)
    mcmc_params = calculate_mcmc_parameters(result, pixels, wavelengths, n_pixels_total)
    print_results(result, mcmc_params, pixels, wavelengths)
    
    return result, mcmc_params

def main():
    print("\n" + "=" * 70)
    print("WAVELENGTH CALIBRATION - POLYNOMIAL FITTER FOR MCMC")
    print("=" * 70)
    print("\nSelect mode:")
    print("  1. Interactive mode (enter data manually)")
    print("  2. Hardcoded mode (use values in script)")
    
    while True:
        choice = input("\nEnter choice (1 or 2): ").strip()
        if choice == '1':
            interactive_mode()
            break
        elif choice == '2':
            hardcoded_mode()
            break
        else:
            print("Please enter 1 or 2")

if __name__ == "__main__":
    main()