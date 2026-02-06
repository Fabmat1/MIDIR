# src/updater.py
import json
import os
import subprocess
import threading
import tkinter as tk
from tkinter import ttk, messagebox
from urllib.request import urlopen, Request
from urllib.error import URLError
from packaging import version  # pip install packaging

from version import __version__

GITHUB_API_URL = f"https://api.github.com/repos/Fabmat1/MIDIR/releases/latest"
CONFIG_FILE = os.path.join(os.path.dirname(__file__), "..", "update_config.json")


def load_update_config():
    """Load user preferences for update reminders."""
    defaults = {
        "skip_version": None,       # version string the user chose to skip
        "dont_remind": False,       # global "don't remind me again"
        "last_check": None,         # ISO timestamp of last check
    }
    if os.path.exists(CONFIG_FILE):
        try:
            with open(CONFIG_FILE, "r") as f:
                saved = json.load(f)
            defaults.update(saved)
        except (json.JSONDecodeError, IOError):
            pass
    return defaults


def save_update_config(config: dict):
    """Persist user preferences."""
    with open(CONFIG_FILE, "w") as f:
        json.dump(config, f, indent=2)


def fetch_latest_release() -> dict | None:
    """
    Query the GitHub API for the latest release.
    Returns dict with keys: tag_name, html_url, body (release notes)
    Returns None on failure.
    """
    try:
        req = Request(GITHUB_API_URL, headers={"Accept": "application/vnd.github.v3+json"})
        with urlopen(req, timeout=5) as resp:
            data = json.loads(resp.read().decode())
        return {
            "tag_name": data["tag_name"].lstrip("v"),
            "html_url": data["html_url"],
            "body": data.get("body", ""),
        }
    except (URLError, KeyError, json.JSONDecodeError) as e:
        print(f"[Updater] Failed to check for updates: {e}")
        return None


def is_newer(remote_version: str, local_version: str) -> bool:
    """Compare semantic versions."""
    try:
        return version.parse(remote_version) > version.parse(local_version)
    except version.InvalidVersion:
        return False


def pull_and_build(progress_callback=None) -> tuple[bool, str]:
    """
    Run `git pull` then `make` in the project root.
    Returns (success: bool, message: str).
    """
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    try:
        # --- git pull ---
        if progress_callback:
            progress_callback("Pulling latest changes from GitHub...")

        result = subprocess.run(
            ["git", "pull", "--ff-only"],
            cwd=project_root,
            capture_output=True,
            text=True,
            timeout=60,
        )
        if result.returncode != 0:
            return False, f"git pull failed:\n{result.stderr}"

        # --- make ---
        if progress_callback:
            progress_callback("Rebuilding C++ components (make)...")

        result = subprocess.run(
            ["make", "-C", project_root],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0:
            return False, f"make failed:\n{result.stderr}"

        if progress_callback:
            progress_callback("Update complete!")

        return True, result.stdout

    except FileNotFoundError as e:
        return False, f"Command not found: {e}"
    except subprocess.TimeoutExpired:
        return False, "Command timed out."


class UpdateDialog(tk.Toplevel):
    """
    Modal dialog that informs the user about a new version
    and optionally pulls + rebuilds.
    """

    def __init__(self, parent, remote_version: str, release_url: str,
                 release_notes: str):
        super().__init__(parent)
        self.title("Update Available")
        self.geometry("520x370")
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()

        self.result = None          # "update" | "skip" | "later" | "dismiss"
        self.remote_version = remote_version
        self._parent = parent

        # ---- Header ----
        header = ttk.Label(
            self,
            text=f"A new version is available!",
            font=("Segoe UI", 13, "bold"),
        )
        header.pack(pady=(18, 4))

        version_label = ttk.Label(
            self,
            text=f"Current: v{__version__}   →   Latest: v{remote_version}",
            font=("Segoe UI", 10),
        )
        version_label.pack(pady=(0, 10))

        # ---- Release notes (scrollable) ----
        notes_frame = ttk.LabelFrame(self, text="Release Notes")
        notes_frame.pack(fill="both", expand=True, padx=20, pady=(0, 5))

        notes_text = tk.Text(notes_frame, wrap="word", height=8, relief="flat")
        notes_text.insert("1.0", release_notes or "No release notes provided.")
        notes_text.configure(state="disabled")
        notes_text.pack(fill="both", expand=True, padx=5, pady=5)

        # ---- "Don't remind me" checkbox ----
        self.dont_remind_var = tk.BooleanVar(value=False)
        self.skip_version_var = tk.BooleanVar(value=False)

        check_frame = ttk.Frame(self)
        check_frame.pack(fill="x", padx=20, pady=(5, 0))

        ttk.Checkbutton(
            check_frame,
            text="Don't remind me again",
            variable=self.dont_remind_var,
        ).pack(side="left")

        ttk.Checkbutton(
            check_frame,
            text=f"Skip v{remote_version}",
            variable=self.skip_version_var,
        ).pack(side="left", padx=(15, 0))

        # ---- Buttons ----
        btn_frame = ttk.Frame(self)
        btn_frame.pack(fill="x", padx=20, pady=(10, 15))

        ttk.Button(btn_frame, text="Update Now",
                    command=self._on_update).pack(side="right", padx=(5, 0))
        ttk.Button(btn_frame, text="Remind Me Later",
                    command=self._on_later).pack(side="right", padx=(5, 0))
        ttk.Button(btn_frame, text="View on GitHub",
                    command=lambda: __import__("webbrowser").open(release_url)
                    ).pack(side="left")

        # ---- Progress widgets (hidden initially) ----
        self.progress_label = ttk.Label(self, text="")
        self.progress_bar = ttk.Progressbar(self, mode="indeterminate",
                                             length=460)

        self.protocol("WM_DELETE_WINDOW", self._on_later)

    # ------------------------------------------------------------------
    def _save_prefs(self):
        config = load_update_config()
        if self.dont_remind_var.get():
            config["dont_remind"] = True
        if self.skip_version_var.get():
            config["skip_version"] = self.remote_version
        save_update_config(config)

    def _on_later(self):
        self._save_prefs()
        self.result = "later"
        self.destroy()

    def _on_update(self):
        self._save_prefs()
        # Show progress widgets
        self.progress_label.pack(padx=20, pady=(0, 2))
        self.progress_bar.pack(padx=20, pady=(0, 10))
        self.progress_bar.start(15)

        def _do_update():
            success, msg = pull_and_build(
                progress_callback=lambda t: self.after(
                    0, lambda: self.progress_label.config(text=t))
            )
            self.after(0, lambda: self._update_finished(success, msg))

        threading.Thread(target=_do_update, daemon=True).start()

    def _update_finished(self, success, msg):
        self.progress_bar.stop()

        if success:
            self.progress_label.config(text="Update successful!")
            messagebox.showinfo(
                "Update Complete",
                f"Updated to v{self.remote_version}.\n"
                "Please restart the application for changes to take effect.",
                parent=self,
            )
            # Reset skip/dismiss for the now-current version
            config = load_update_config()
            config["skip_version"] = None
            config["dont_remind"] = False
            save_update_config(config)
        else:
            self.progress_label.config(text="Update failed.")
            messagebox.showerror("Update Failed", msg, parent=self)

        self.result = "update" if success else "failed"
        self.destroy()


# ------------------------------------------------------------------
# Public API – call this once at startup
# ------------------------------------------------------------------
def check_for_updates(parent: tk.Tk, force: bool = False):
    """
    Non-blocking update check.  Runs the network call in a thread
    and pops up the dialog on the main thread only when needed.
    """
    config = load_update_config()

    if not force and config.get("dont_remind"):
        print("[Updater] User opted out of reminders.")
        return

    def _background():
        release = fetch_latest_release()
        if release is None:
            return

        remote_ver = release["tag_name"]

        if not is_newer(remote_ver, __version__):
            print(f"[Updater] Up to date (v{__version__}).")
            return

        if not force and config.get("skip_version") == remote_ver:
            print(f"[Updater] User skipped v{remote_ver}.")
            return

        # Schedule dialog on the main thread
        parent.after(0, lambda: UpdateDialog(
            parent,
            remote_version=remote_ver,
            release_url=release["html_url"],
            release_notes=release["body"],
        ))

    threading.Thread(target=_background, daemon=True).start()