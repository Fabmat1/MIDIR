# src/font_loader.py
"""
Cross-platform TTF font loader for tkinter.
Registers bundled Inter static TTFs so every OS renders the same UI.
"""

import os
import platform
import subprocess
import tkinter as tk
import tkinter.font as tkfont
import tkinter.ttk as ttk
from pathlib import Path

FONT_DIR = Path(__file__).resolve().parent.parent / "fonts"

# Only the weights we actually use — all from the 18pt optical-size family
FONT_FILES = {
    "regular":  "Inter_18pt-Regular.ttf",
    "medium":   "Inter_18pt-Medium.ttf",
    "semibold": "Inter_18pt-SemiBold.ttf",
    "bold":     "Inter_18pt-Bold.ttf",
}

# The family name embedded in the 18pt static TTFs
# (You can verify with: fc-query Inter_18pt-Regular.ttf | grep family)
FONT_FAMILY = "Inter 18pt"


# ──────────────────────────────────────────────
# Platform-specific registration backends
# ──────────────────────────────────────────────

def _register_linux():
    """Install into user fontconfig dir and rebuild cache."""
    user_font_dir = Path.home() / ".local" / "share" / "fonts" / "midir"
    user_font_dir.mkdir(parents=True, exist_ok=True)

    changed = False
    for filename in FONT_FILES.values():
        src = FONT_DIR / filename
        dst = user_font_dir / filename
        if src.exists() and not dst.exists():
            try:
                dst.symlink_to(src.resolve())
            except OSError:
                import shutil
                shutil.copy2(src, dst)
            changed = True

    if changed:
        try:
            subprocess.run(
                ["fc-cache", "-f", str(user_font_dir)],
                capture_output=True, timeout=15,
            )
            print(f"[Font] Installed fonts → {user_font_dir}")
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            print(f"[Font] fc-cache failed: {e}")


def _register_windows():
    """Add fonts to the current GDI session (process-scoped)."""
    import ctypes
    FR_PRIVATE = 0x10
    for filename in FONT_FILES.values():
        path = str(FONT_DIR / filename)
        if os.path.exists(path):
            added = ctypes.windll.gdi32.AddFontResourceExW(path, FR_PRIVATE, 0)
            if added:
                print(f"[Font] Loaded {filename} (GDI)")


def _register_macos():
    """Register with Core Text for the current process."""
    try:
        import ctypes, ctypes.util
        ct = ctypes.cdll.LoadLibrary(ctypes.util.find_library("CoreText"))
        cf = ctypes.cdll.LoadLibrary(ctypes.util.find_library("CoreFoundation"))

        for filename in FONT_FILES.values():
            path = str((FONT_DIR / filename).resolve())
            if not os.path.exists(path):
                continue
            url_bytes = f"file://{path}".encode("utf-8")

            cf.CFStringCreateWithCString.restype = ctypes.c_void_p
            cf.CFStringCreateWithCString.argtypes = [
                ctypes.c_void_p, ctypes.c_char_p, ctypes.c_uint32,
            ]
            cf_str = cf.CFStringCreateWithCString(None, url_bytes, 0x08000100)

            cf.CFURLCreateWithString.restype = ctypes.c_void_p
            cf.CFURLCreateWithString.argtypes = [
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
            ]
            cf_url = cf.CFURLCreateWithString(None, cf_str, None)

            ct.CTFontManagerRegisterFontsForURL.argtypes = [
                ctypes.c_void_p, ctypes.c_uint32, ctypes.c_void_p,
            ]
            ct.CTFontManagerRegisterFontsForURL(cf_url, 1, None)
            cf.CFRelease(cf_str)
            cf.CFRelease(cf_url)

        print("[Font] Registered via CoreText")
    except Exception as e:
        print(f"[Font] macOS registration failed: {e}")


def _try_tkextrafont() -> bool:
    """
    Best option: tkextrafont injects directly into Tk's runtime.
    pip install tkextrafont
    """
    try:
        from tkextrafont import Font as ExtraFont
        for filename in FONT_FILES.values():
            path = str(FONT_DIR / filename)
            if os.path.exists(path):
                ExtraFont(file=path)
                print(f"[Font] Loaded {filename} (tkextrafont)")
        return True
    except ImportError:
        return False


# ──────────────────────────────────────────────
# Resolve the actual Tk family name
# ──────────────────────────────────────────────

def _resolve_family(root: tk.Tk) -> str:
    """
    Different Tk versions / platforms may register the family under
    slightly different names.  Try the most likely ones.
    """
    available = set(tkfont.families(root=root))

    # Candidates in order of likelihood
    for candidate in (
        "Inter 18pt",       # what fc-query reports for the static 18pt files
        "Inter",            # some Tk versions drop the optical-size suffix
        "Inter 18",
    ):
        if candidate in available:
            return candidate

    # Fuzzy match: anything starting with "Inter"
    for fam in sorted(available):
        if fam.lower().startswith("inter"):
            print(f"[Font] Fuzzy-matched family: '{fam}'")
            return fam

    return ""


# ──────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────

_resolved_family: str = ""


def register_fonts(root: tk.Tk):
    """
    Register bundled TTFs with the OS/toolkit.
    Call ONCE right after creating the root Tk() window.
    """
    global _resolved_family

    if not FONT_DIR.is_dir():
        print(f"[Font] Font directory missing: {FONT_DIR}")
        return

    # Try tkextrafont first (most reliable), then platform fallback
    if not _try_tkextrafont():
        system = platform.system()
        if system == "Linux":
            _register_linux()
        elif system == "Windows":
            _register_windows()
        elif system == "Darwin":
            _register_macos()

    # On Linux with fc-cache, Tk may need a nudge to reload
    # Force Tk to refresh its internal font list
    try:
        root.tk.call("font", "families")
    except tk.TclError:
        pass

    _resolved_family = _resolve_family(root)
    if _resolved_family:
        print(f"[Font] Resolved family: '{_resolved_family}'")
    else:
        print("[Font] Could not find Inter in Tk font families")


def apply_font_globally(root: tk.Tk, size: int = 10):
    """
    Override all default Tk/ttk fonts to use the bundled family.
    Gracefully falls back to the best available system font.
    """
    global _resolved_family

    family = _resolved_family

    if not family:
        # Fallback chain: try common good-looking fonts
        available = set(tkfont.families(root=root))
        for fallback in (
            "Segoe UI",
            "Helvetica Neue",
            "Noto Sans",
            "DejaVu Sans",
            "Liberation Sans",
            "Cantarell",
            "Ubuntu",
            "Helvetica",
        ):
            if fallback in available:
                family = fallback
                print(f"[Font] Falling back to: {family}")
                break
        else:
            print("[Font] No preferred font found, using Tk defaults")
            return

    _resolved_family = family  # cache for get_font()

    # ── Override Tk named fonts (controls ALL default widget fonts) ──
    font_overrides = {
        "TkDefaultFont":      {"family": family, "size": size},
        "TkTextFont":         {"family": family, "size": size},
        "TkMenuFont":         {"family": family, "size": size},
        "TkHeadingFont":      {"family": family, "size": size, "weight": "bold"},
        "TkCaptionFont":      {"family": family, "size": size},
        "TkSmallCaptionFont": {"family": family, "size": size - 1},
        "TkIconFont":         {"family": family, "size": size},
        "TkTooltipFont":      {"family": family, "size": size - 1},
    }

    for name, config in font_overrides.items():
        try:
            tkfont.nametofont(name).configure(**config)
        except tk.TclError:
            pass

    # ── Override ttk styles ──
    style = ttk.Style(root)
    default_font = (family, size)
    bold_font = (family, size, "bold")

    for widget_style in (
        ".", "TButton", "TLabel", "TCheckbutton",
        "TRadiobutton", "TNotebook.Tab", "TCombobox",
        "TEntry", "TLabelframe.Label", "TMenubutton",
    ):
        style.configure(widget_style, font=default_font)

    style.configure("Treeview", font=default_font)
    style.configure("Treeview.Heading", font=bold_font)


def get_font(size: int = 10, weight: str = "normal") -> tuple:
    """
    Convenience helper for explicit font tuples.

    Usage:
        label = ttk.Label(root, text="Title", font=get_font(14, "bold"))
    """
    family = _resolved_family or FONT_FAMILY
    return (family, size, weight)